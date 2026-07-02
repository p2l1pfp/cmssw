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

  // Class for intermediate variables
  class PtPxPy {
  public:
    pt_t pt = 0.;
    pxy_t px = 0.;
    pxy_t py = 0.;

    PtPxPy operator+(const PtPxPy& b) const {
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

  // Software emulation of hls::atan2(pxy_t, pxy_t) = generic_atan2<W=16,I=13>.
  // Replicates the fixed-point CORDIC arithmetic bit-exactly so that the CMSSW
  // emulator matches the HLS firmware/csim output.
  inline ap_fixed<12, 3> atan2_cordic(pxy_t in1, pxy_t in2) {
    static constexpr int W = pxy_t::width;
    static constexpr int I = pxy_t::iwidth;
    static constexpr int WC = W + 7;
    static constexpr int NITER = WC - 3;

    static const ap_fixed<W + 1, 3> pi_ap = M_PI;
    static const ap_fixed<WC, 3> pi2_ap = M_PI / 2.0;
    static const ap_fixed<W + 1, 3> pi4_ap = M_PI / 4.0;
    static const ap_fixed<W + 1, 3> pi3n_ap = -3.0 * M_PI / 4.0;

    // Construct once. Values are truncated into ap_fixed<WC,3>, matching HLS.
    static const auto atan_lut = [] {
      std::array<ap_fixed<WC, 3>, NITER> lut{};
      for (int i = 0; i < NITER; ++i)
        lut[i] = std::atan(std::ldexp(1.0, -i));  // atan(2^-i)
      return lut;
    }();

    const ap_uint<2> signin1 = (in1 > 0) ? 2 : (in1 == 0) ? 1 : 0;
    const ap_uint<2> signin2 = (in2 > 0) ? 2 : (in2 == 0) ? 1 : 0;

    ap_fixed<W, 3> out;

    // Special cases (match generic_atan2)
    if (signin1 == 1 && signin2 == 2)
      return 0;
    if (signin1 == 1 && signin2 == 0)
      return pi_ap;
    if (signin1 == 2 && signin2 == 1)
      return pi2_ap;
    if (signin1 == 0 && signin2 == 1)
      return -pi2_ap;

    if (in1 == in2) {
      if (signin1 == 2)
        return pi4_ap;
      if (signin1 == 1)
        return 0;
      return pi3n_ap;
    }

    // Absolute values
    ap_fixed<W + 1, I + 1> in1abs = (signin1 == 0) ? (ap_fixed<W + 1, I + 1>)(-in1) : (ap_fixed<W + 1, I + 1>)(in1);
    ap_fixed<W + 1, I + 1> in2abs = (signin2 == 0) ? (ap_fixed<W + 1, I + 1>)(-in2) : (ap_fixed<W + 1, I + 1>)(in2);

    // Bit reinterpretation (not a numeric conversion!)
    ap_fixed<W + 1, 2> in1abs_sft, in2abs_sft;
    in1abs_sft.range() = in1abs.range();
    in2abs_sft.range() = in2abs.range();

    const bool swap = (in1abs <= in2abs);

    ap_fixed<WC, 3> cx = swap ? in2abs_sft : in1abs_sft;
    ap_fixed<WC, 3> cy = swap ? in1abs_sft : in2abs_sft;
    ap_fixed<WC, 3> cz = 0;

    for (int i = 0; i < NITER; ++i) {
      ap_fixed<WC, 3> cx_new, cy_new, cz_new;

      if (cy >= 0) {
        cx_new = cx + (cy >> i);
        cy_new = cy - (cx >> i);
        cz_new = cz + atan_lut[i];
      } else {
        cx_new = cx - (cy >> i);
        cy_new = cy + (cx >> i);
        cz_new = cz - atan_lut[i];
      }

      cx = cx_new;
      cy = cy_new;
      cz = cz_new;
    }

    if (!swap)
      cz = pi2_ap - cz;

    if (signin2 == 0 && signin1 == 2)
      out = pi_ap - cz;
    else if (signin2 == 0 && signin1 == 0)
      out = cz - pi_ap;
    else if (signin2 == 2 && signin1 == 0)
      out = -cz;
    else
      out = cz;

    return out;
  }

  inline etaphi_t phi_cordic(pxy_t y, pxy_t x) {
#ifdef CMSSW_GIT_HASH
    // ap_fixed<12, 3> phi = atan2(y.to_double(), x.to_double());  // hls_math.h not available yet in CMSSW
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
    cossin_t sin_table[N_TABLE];
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
      jets.begin(), jets.end(), ptpxpy.begin(), [](const l1ct::Jet& jet) { return P2L1HTMHTEmu::mht_compute(jet); });

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
