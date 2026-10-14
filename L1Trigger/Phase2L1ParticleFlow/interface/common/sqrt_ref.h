#ifndef L1Trigger_Phase2L1ParticleFlow_sqrt_ref_h
#define L1Trigger_Phase2L1ParticleFlow_sqrt_ref_h

#include "ap_fixed.h"
#include <cmath>

namespace l1ct {
  // Emulate fixed-point hls::sqrt for nonnegative MET/MHT squared magnitudes.
  // Round at the input fractional precision before applying the output type's quantization and saturation.
  template <class out_t, class in_t>
  inline out_t sqrt_ref(const in_t& input) {
    const ap_ufixed<in_t::width, in_t::iwidth, AP_RND_ZERO> rounded = std::sqrt(input.to_double());
    return out_t(rounded);
  }
}  // namespace l1ct

#endif
