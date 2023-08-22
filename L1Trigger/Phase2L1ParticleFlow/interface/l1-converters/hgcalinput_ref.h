#ifndef L1Trigger_Phase2L1ParticleFlow_newfirmware_hgcalinput_ref_h
#define L1Trigger_Phase2L1ParticleFlow_newfirmware_hgcalinput_ref_h

#include "DataFormats/L1TParticleFlow/interface/layer1_emulator.h"

namespace edm {
  class ParameterSet;
}

namespace l1ct {
  class HgcalClusterDecoderEmulator {
    bool slim_;

  public:
    HgcalClusterDecoderEmulator(bool slim = false) : slim_{slim} {};
    HgcalClusterDecoderEmulator(const edm::ParameterSet &pset);

    ~HgcalClusterDecoderEmulator();
    l1ct::HadCaloObjEmu decode(const l1ct::PFRegionEmu &sector, const ap_uint<256> &in) const;
  };
}  // namespace l1ct

#endif
