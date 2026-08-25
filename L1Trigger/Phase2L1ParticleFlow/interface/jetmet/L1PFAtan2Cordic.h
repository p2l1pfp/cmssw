#ifndef L1Trigger_Phase2L1ParticleFlow_Atan2Cordic_h
#define L1Trigger_Phase2L1ParticleFlow_Atan2Cordic_h

#include "DataFormats/L1TParticleFlow/interface/jets.h"
#include "DataFormats/L1TParticleFlow/interface/sums.h"

#include <array>
#include <cmath>

#include "ap_fixed.h"
#include "ap_int.h"

#ifndef CMSSW_GIT_HASH
#include "hls_math.h"
#endif

namespace P2L1ATanCordicEmu {
  template <typename IN_T>
  struct CordicConfig {
    // Constants for CORDIC atan2 implementation
    // Based on those used by hls implementation
    static constexpr int CORDIC_INPUT_W = IN_T::width;
    static constexpr int CORDIC_INPUT_I = IN_T::iwidth;
    static constexpr int CORDIC_GUARD_BITS = 7;
    static constexpr int CORDIC_WORKING_W = CORDIC_INPUT_W + CORDIC_GUARD_BITS; // Larger internal working bitwidth with "enough" guard bits
    static constexpr int CORDIC_WORKING_I = 3; // Number of integer bits in CORDIC working representation
    static constexpr int CORDIC_NITER = CORDIC_WORKING_W; // Number of iterations to perform

    // Useful typedefs for CORDIC implementation
    typedef ap_fixed<CORDIC_INPUT_W + 1, CORDIC_INPUT_I + 1> cordic_abs_t;
    typedef ap_fixed<CORDIC_INPUT_W + 1, 2> cordic_scaled_input_t;
    typedef ap_fixed<CORDIC_WORKING_W, CORDIC_WORKING_I> cordic_working_t;
    typedef std::array<ap_ufixed<128, 2>, CORDIC_NITER> cordic_lut_t;
  };

  template <typename IN_T>
  inline void init_atan_table(typename CordicConfig<IN_T>::cordic_lut_t& atan_lut) {
    for (int i = 0; i < CordicConfig<IN_T>::CORDIC_NITER; ++i) {
      atan_lut[i] = typename CordicConfig<IN_T>::cordic_working_t(
          ap_ufixed<128, 2>(std::atan(std::ldexp(1.0, -i))));
    }
  }

  // Software emulation of hls::atan2(pxy_t, pxy_t) = generic_atan2<W=16,I=13>.
  // Replicates the fixed-point CORDIC arithmetic bit-exactly so that the CMSSW
  // emulator matches the HLS firmware/csim output.
  template <typename OUT_T, typename IN_T>
  inline OUT_T atan2_cordic(IN_T in1, IN_T in2) {
    typedef CordicConfig<IN_T> cfg_t;
    typedef typename cfg_t::cordic_abs_t cordic_abs_t;
    typedef typename cfg_t::cordic_scaled_input_t cordic_scaled_input_t;
    typedef typename cfg_t::cordic_working_t cordic_working_t;

    // Fixed constants, in precision used by the HLS implementation
    const cordic_working_t pi_ap(M_PI);
    const cordic_working_t pi2_ap(M_PI_2);
    const cordic_working_t pi4_ap(M_PI_4);
    const cordic_working_t pi3n_ap(-3 * M_PI_4);

    // Encode the sign of the inputs (0=negative, 1=zero, 2=positive)
    const ap_uint<2> signin1 = (in1 > 0) ? 2 : (in1 == 0) ? 1 : 0;
    const ap_uint<2> signin2 = (in2 > 0) ? 2 : (in2 == 0) ? 1 : 0;

    // Special cases (match generic_atan2)
    // If any inputs are zero, no need to run CORDIC
    if (signin1 == 1 && signin2 == 2)
      return OUT_T(0);
    if (signin1 == 1 && signin2 == 0)
      return OUT_T(pi_ap);
    if (signin1 == 2 && signin2 == 1)
      return OUT_T(pi2_ap);
    if (signin1 == 0 && signin2 == 1)
      return OUT_T(-pi2_ap);
    // If inputs are equal, return +/- pi/4 or -3pi/4 depending on the signs
    if (in1 == in2) {
      if (signin1 == 2)
        return OUT_T(pi4_ap);
      if (signin1 == 1)
        return OUT_T(0);
      return OUT_T(pi3n_ap);
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
    // Also initialize atan LUT once
    static typename cfg_t::cordic_lut_t atan_lut;
    static const bool atan_lut_init = []() {
      init_atan_table<IN_T>(atan_lut);
      return true;
    }();
    (void)atan_lut_init;  // Does nothing, just to avoid unused variable warning

    for (int i = 0; i < cfg_t::CORDIC_NITER; ++i) {
      cordic_working_t cx_new, cy_new, cz_new;
      const cordic_working_t angle = cordic_working_t(atan_lut[i]);

      if (cy[cfg_t::CORDIC_WORKING_W - 1] == 0) {
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
      return OUT_T(pi_ap - cz);
    else if (signin2 == 0 && signin1 == 0)
      return OUT_T(cz - pi_ap);
    else if (signin2 == 2 && signin1 == 0)
      return OUT_T(-cz);
    else
      return OUT_T(cz);
  }
}  // namespace P2L1ATanCordicEmu

#endif
