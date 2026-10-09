#ifndef L1Trigger_Phase2L1ParticleFlow_newfirmware_gctcommoninput_ref_h
#define L1Trigger_Phase2L1ParticleFlow_newfirmware_gctcommoninput_ref_h

#include <vector>

#include "DataFormats/L1TParticleFlow/interface/layer1_emulator.h"

#include "L1Trigger/Phase2L1ParticleFlow/interface/l1-converters/gcteminput_ref.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/l1-converters/gcthadinput_ref.h"

namespace l1ct {
  // Unpacks the raw GCT cluster words of the TMUX18 links into common calo objects (EM and hadronic clusters in the
  // same format), expressed in the frame of the link sector, as they are seen by the regionizer.
  //
  // There are 3 links (one per TMUX18 sector) and 162 clocks per link. In each half of 81 clocks: clock 0 is a
  // technical word, then 2x16 EM words and 2x24 hadronic words, each half-block coming from a different GCT SLR
  // (a 60 degrees sector on one eta side). The SLR of a word is fixed by its link and position in the frame.
  class GctCommonCaloDecoderEmulator {
  public:
    static constexpr unsigned int NLINKS = 3;
    static constexpr unsigned int NCLOCKS = 162;

    // The regions are those set up by the caller (in the producer: the decoded calo sectors and the raw GCT sectors):
    // the 12 decoded sectors (60 degrees x 2 eta sides) that the raw words refer to, and the 3 link sectors.
    // The decoders are not copied, they must outlive this object.
    GctCommonCaloDecoderEmulator(const GctEmClusterDecoderEmulator &emDecoder,
                                 const GctHadClusterDecoderEmulator &hadDecoder,
                                 const std::vector<l1ct::PFRegionEmu> &decodedSectors,
                                 const std::vector<l1ct::PFRegionEmu> &linkSectors);

    // one word of a link; technical words and words with pt == 0 give an empty object
    l1ct::CommonCaloObjEmu decode(unsigned int link, unsigned int iclock, const ap_uint<64> &word) const;

    // all the words: one sector per link, with the link region, and one object per clock
    void decode(const std::vector<l1ct::DetectorSector<ap_uint<64>>> &raw,
                std::vector<l1ct::DetectorSector<l1ct::CommonCaloObjEmu>> &gctcommon) const;

  private:
    const GctEmClusterDecoderEmulator &emDecoder_;
    const GctHadClusterDecoderEmulator &hadDecoder_;
    std::vector<l1ct::PFRegionEmu> decodedSectors_, linkSectors_;
  };
}  // namespace l1ct

#endif
