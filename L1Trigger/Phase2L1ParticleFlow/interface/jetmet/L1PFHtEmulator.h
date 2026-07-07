#ifndef L1Trigger_Phase2L1ParticleFlow_HTMHT_h
#define L1Trigger_Phase2L1ParticleFlow_HTMHT_h

#include "DataFormats/L1TParticleFlow/interface/jets.h"
#include "DataFormats/L1TParticleFlow/interface/sums.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/dbgPrintf.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/jetmet/L1SeedConePFJetEmulator.h"

#ifndef CMSSW_GIT_HASH
#include "hls_math.h"
#endif

#include <vector>
#include <numeric>
#include <algorithm>
#include "ap_int.h"
#include "ap_fixed.h"

namespace P2L1HTMHTEmu {
  typedef l1ct::pt_t pt_t;          // Type for pt/ht 1 unit = 0.25 GeV; max = 16 TeV
  typedef l1ct::glbeta_t etaphi_t;  // Type for eta & phi

  typedef ap_fixed<12, 3> radians_t;
  typedef ap_fixed<9, 2> cossin_t;
  typedef ap_fixed<16, 13> pxy_t;
  static constexpr int Fin = 6;  // Number of decimal bits/precision of met squared i.e. 2*16 - 2*13
  static constexpr int Fout = pt_t::width - pt_t::iwidth;  // Number of decimal bits/precision of output met

  static constexpr int N_TABLE = 2048;

  // Constants for CORDIC atan2 implementation
  // Based on those used by hls implementation
  static constexpr int CORDIC_INPUT_W = pxy_t::width;
  static constexpr int CORDIC_INPUT_I = pxy_t::iwidth;
  static constexpr int CORDIC_GUARD_BITS = 7;
  static constexpr int CORDIC_WORKING_W =
      CORDIC_INPUT_W + CORDIC_GUARD_BITS;                // Larger internal working bitwidth with "enough" guard bits
  static constexpr int CORDIC_WORKING_I = 3;             // Number of integer bits in CORDIC working representation
  static constexpr int CORDIC_NITER = CORDIC_WORKING_W;  // Number of iterations to perform

  // Fixed constants, in precision used by the HLS implementation and their exact values
  static const ap_fixed<CORDIC_INPUT_W + 1, CORDIC_WORKING_I> pi_ap("0x3.243F6A8885A308D3");
  static const ap_fixed<CORDIC_INPUT_W + 2, CORDIC_WORKING_I> pi2_ap("0x1.921FB54442D1846");    // pi/2
  static const ap_fixed<CORDIC_INPUT_W + 1, CORDIC_WORKING_I> pi4_ap("0x0.C90FDAA22168C23");    // pi/4
  static const ap_fixed<CORDIC_INPUT_W + 1, CORDIC_WORKING_I> pi3n_ap("-0x2.5B2F8FE6643A469");  // -3*pi/4

  // Useful typedefs for CORDIC implementation
  typedef ap_fixed<CORDIC_INPUT_W + 1, CORDIC_INPUT_I + 1> cordic_abs_t;
  typedef ap_fixed<CORDIC_INPUT_W + 1, 2> cordic_scaled_input_t;
  typedef ap_fixed<CORDIC_WORKING_W, CORDIC_WORKING_I> cordic_working_t;

  // The exact LUT values used by the HLS implementation
  static const std::array<ap_ufixed<128, 2>, 23> atan_lut = {
      "0x0.C90FDAA22168C234C4C6628B80DC1CD0", "0x0.76B19C1586ED3DA2B7F222F65E1D4680",
      "0x0.3EB6EBF25901BAC55B71E7BD7DE885F8", "0x0.1FD5BA9AAC2F6DC65912F313E7D111DC",
      "0x0.0FFAADDB967EF4E36CB2792DC0E2E0D4", "0x0.07FF556EEA5D892A13BCEBBB6ED46310",
      "0x0.03FFEAAB776E5356EF9E31590057DD80", "0x0.01FFFD555BBBA972D00C46A3F77CC15C",
      "0x0.00FFFFAAAADDDDB94BB12AFB6B6D4F7C", "0x0.007FFFF55556EEEEA5CA6ADEAB02251C",
      "0x0.003FFFFEAAAAB77776E52E5A019FBCE8", "0x0.001FFFFFD55555BBBBBA972976256248",
      "0x0.000FFFFFFAAAAAADDDDDDB94B94D5BD4", "0x0.0007FFFFFF5555556EEEEEEA5CA5CB40",
      "0x0.0003FFFFFFEAAAAAAB7777776E52E52C", "0x0.0001FFFFFFFD5555555BBBBBBBA97294",
      "0x0.0000FFFFFFFFAAAAAAAADDDDDDDDB948", "0x0.00007FFFFFFFF555555556EEEEEEEEA4",
      "0x0.00003FFFFFFFFEAAAAAAAAB777777774", "0x0.00001FFFFFFFFFD555555555BBBBBBB8",
      "0x0.00000FFFFFFFFFFAAAAAAAAAADDDDDDC", "0x0.000007FFFFFFFFFF55555555556EEEEC",
      "0x0.000003FFFFFFFFFFEAAAAAAAAAAB7774"};

  // Class for intermediate variables
  class PtPxPy {
  public:
    pt_t pt = 0.;
    pxy_t px = 0.;
    pxy_t py = 0.;

    PtPxPy operator+(const PtPxPy &b) const {
      PtPxPy c;
      c.pt = this->pt + b.pt;
      c.px = this->px + b.px;
      c.py = this->py + b.py;
      return c;
    }
  };

  namespace Scales {
    const ap_fixed<12, -4> scale_degToRad = M_PI / 180.;
  };  // namespace Scales

  template <class data_T, class table_T, int N>
  void init_sinphi_table(table_T table_out[N]) {
    for (int i = 0; i < N; i++) {
      double x = i * (M_PI / 180.) / 2.;
      table_T sin_x = std::sin(x);
      table_out[i] = sin_x;
    }
  }
  template <class in_t, class table_t, int N>
  table_t sine_with_conversion(etaphi_t hwPhi) {
    table_t sin_table[N];
    init_sinphi_table<in_t, table_t, N>(sin_table);
    table_t out = sin_table[hwPhi];
    return out;
  }
  
  inline void init_atan_table( std::array<ap_ufixed<128, 2>, 23> &atan_lut) {
    for (int i = 0; i < CORDIC_NITER; ++i) {
      atan_lut[i] = cordic_working_t(ap_ufixed<128, 2>(std::atan(std::ldexp(1.0, -i))));
    }
  }

  // Software emulation of hls::atan2(pxy_t, pxy_t) = generic_atan2<W=16,I=13>.
  // Replicates the fixed-point CORDIC arithmetic bit-exactly so that the CMSSW
  // emulator matches the HLS firmware/csim output.
  inline ap_fixed<12, 3> atan2_cordic(pxy_t in1, pxy_t in2) {
    // Encode the sign of the inputs (0=negative, 1=zero, 2=positive)
    const ap_uint<2> signin1 = (in1 > 0) ? 2 : (in1 == 0) ? 1 : 0;
    const ap_uint<2> signin2 = (in2 > 0) ? 2 : (in2 == 0) ? 1 : 0;

    // Special cases (match generic_atan2)
    // If any inputs are zero, no need to run CORDIC
    if (signin1 == 1 && signin2 == 2)
      return 0;
    if (signin1 == 1 && signin2 == 0)
      return pi_ap;
    if (signin1 == 2 && signin2 == 1)
      return pi2_ap;
    if (signin1 == 0 && signin2 == 1)
      return -pi2_ap;
    // If inputs are equal, return +/- pi/4 or -3pi/4 depending on the signs
    if (in1 == in2) {
      if (signin1 == 2)
        return pi4_ap;
      if (signin1 == 1)
        return 0;
      return pi3n_ap;
    }

    // Absolute values of inputs
    // Widen by one bit to ensure -in1/2 is representable
    cordic_abs_t in1abs = (signin1 == 0) ? cordic_abs_t(-in1) : cordic_abs_t(in1);
    cordic_abs_t in2abs = (signin2 == 0) ? cordic_abs_t(-in2) : cordic_abs_t(in2);

    // Bit reinterpretation
    // CORDIC prefers working with ~2 integer bits and many fractional bits
    cordic_scaled_input_t in1abs_sft, in2abs_sft;
    in1abs_sft.range() = in1abs.range();
    in2abs_sft.range() = in2abs.range();

    // Ensure cx >= cy for CORDIC, swap in2 and in1 if necessary
    // CORDIC then operates in 0-pi/4 range
    const bool swap = (in1abs < in2abs);
    cordic_working_t cx = swap ? in2abs_sft : in1abs_sft;
    cordic_working_t cy = swap ? in1abs_sft : in2abs_sft;
    cordic_working_t cz = 0;  // Initial angle accumulator

    // CORDIC iterations
    // Each iteration rotates the vector (cx, cy) by atan(2^-i) towards the x-axis, accumulating the angle in cz.
    // Sign check of cy determines the direction of rotation for the current iteration. i.e. if cy<0, the previous iteration overshot the x-axis and the next iteration rotates back towards the x-axis.
    // After all iterations, cy~0 and cz contains the angle of the original vector (in1, in2) in radians.
    for (int i = 0; i < CORDIC_NITER; ++i) {
      cordic_working_t cx_new, cy_new, cz_new;

      static std::array<ap_ufixed<128, 2>, 23> atan_lut;
      init_atan_table(atan_lut);
      cordic_working_t angle = cordic_working_t(atan_lut[i]);  // convert once

      if (cy[CORDIC_WORKING_W - 1] == 0) {
        cx_new = cx + (cy >> i);
        cy_new = cy - (cx >> i);
        cz_new = cz + angle;
      } else {
        cx_new = cx - (cy >> i);
        cy_new = cy + (cx >> i);
        cz_new = cz - angle;
      }

      cx = cx_new;
      cy = cy_new;
      cz = cz_new;
    }

    // Map cz back to the original quadrant based on the signs of the inputs and whether they were swapped.
    if (!swap)
      cz = pi2_ap - cz;

    if (signin2 == 0 && signin1 == 2)
      return pi_ap - cz;
    else if (signin2 == 0 && signin1 == 0)
      return cz - pi_ap;
    else if (signin2 == 2 && signin1 == 0)
      return -cz;
    else
      return cz;
  }

  inline etaphi_t phi_cordic(pxy_t y, pxy_t x) {
#ifdef CMSSW_GIT_HASH
    ap_fixed<12, 3> phi = atan2_cordic(y, x);
#else
    ap_fixed<12, 3> phi = hls::atan2(y, x);
#endif
    ap_fixed<16, 9> etaphiscale = (float)l1ct::Scales::INTPHI_PI / M_PI;  // radians to hwPhi
    return phi * etaphiscale;
  }

  inline PtPxPy mht_compute(l1ct::Jet jet) {
    // Add an extra bit to px/py for the sign, and one additional bit to improve precision (pt_t is ap_ufixed<14, 12>)
    PtPxPy v_pxpy;

    //Initialize table once
    static cossin_t sin_table[N_TABLE];
    init_sinphi_table<etaphi_t, cossin_t, N_TABLE>(sin_table);

    cossin_t sinphi;
    cossin_t cosphi;
    bool sign = jet.hwPhi.sign();

    etaphi_t hwphi = jet.hwPhi;

    // Reduce precision of hwPhi
    ap_int<10> phi;
    phi.V = hwphi(11, 1);
    phi = (phi > 0) ? phi : (ap_int<10>)-phi;  //Only store values for positive phi, pick up sign later

    sinphi = sin_table[phi];

    sinphi = (sign > 0) ? (cossin_t)(-sign * sinphi) : sinphi;  // Change sign bit if hwPt is negative, sin(-x)=-sin(x)
    cosphi = sin_table[phi + 90 * 2];  //cos(x)=sin(x+90). Do nothing with sign, cos(-θ) = cos θ,

    v_pxpy.pt = jet.hwPt;
    v_pxpy.py = jet.hwPt * sinphi;
    v_pxpy.px = jet.hwPt * cosphi;

    return v_pxpy;
  }
}  // namespace P2L1HTMHTEmu

//TODO replace with l1ct::Jet
inline l1ct::Sum htmht(std::vector<l1ct::Jet> jets) {
  // compute jet px, py
  std::vector<P2L1HTMHTEmu::PtPxPy> ptpxpy;
  ptpxpy.resize(jets.size());
  std::transform(
      jets.begin(), jets.end(), ptpxpy.begin(), [](const l1ct::Jet &jet) { return P2L1HTMHTEmu::mht_compute(jet); });

  // Sum pt, px, py over jets
  P2L1HTMHTEmu::PtPxPy hthxhy = std::accumulate(ptpxpy.begin(), ptpxpy.end(), P2L1HTMHTEmu::PtPxPy());

  // Compute the MHT magnitude and direction
  l1ct::Sum ht;
  ht.hwSumPt = hthxhy.pt;

  // this emulates the following firmware function, since hls_math.h is not available in CMSSW
  // ht.hwPt = hls::sqrt(((hthxhy.px * hthxhy.px) + (hthxhy.py * hthxhy.py)));
  double d = std::sqrt(((hthxhy.px * hthxhy.px) + (hthxhy.py * hthxhy.py)).to_double());
  // emulate hls::sqrt internal rounding
  double rounded = std::round(d * (1 << P2L1HTMHTEmu::Fin)) / (1 << P2L1HTMHTEmu::Fin);
  // emulate AP_TRN conversion to output type
  double truncated = std::floor(rounded * (1 << P2L1HTMHTEmu::Fout)) / (1 << P2L1HTMHTEmu::Fout);
  P2L1HTMHTEmu::pt_t hwPt_hls = truncated;
  ht.hwPt = hwPt_hls;

  ht.hwPhi = P2L1HTMHTEmu::phi_cordic(hthxhy.py, hthxhy.px);
  return ht;
}

#endif
