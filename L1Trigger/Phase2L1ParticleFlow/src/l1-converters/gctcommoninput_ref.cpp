#include "L1Trigger/Phase2L1ParticleFlow/interface/l1-converters/gctcommoninput_ref.h"

#include <cassert>

#include "DataFormats/L1TCalorimeterPhase2/interface/GCTEmDigiCluster.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/GCTHadDigiCluster.h"

namespace {
  // the object is expressed in the frame of the "from" sector (as the decoders return it): move it to the "to" sector
  template <typename T>
  void rebase(T &cl, const l1ct::PFRegionEmu &from, const l1ct::PFRegionEmu &to) {
    l1ct::glbeta_t glbEta = from.hwGlbEtaOf(cl);
    l1ct::glbphi_t glbPhi = from.hwGlbPhiOf(cl);
    assert(to.containsHw(glbEta, glbPhi) && "GCT cluster out of TMUX18 sector bounds!");
    cl.hwEta = glbEta - to.hwEtaCenter;
    ap_int<l1ct::glbphi_t::width + 1> locPhi = glbPhi - to.hwPhiCenter;
    if (locPhi > l1ct::Scales::INTPHI_PI)
      locPhi -= l1ct::Scales::INTPHI_TWOPI;
    else if (locPhi <= -l1ct::Scales::INTPHI_PI)
      locPhi += l1ct::Scales::INTPHI_TWOPI;
    cl.hwPhi = locPhi;
  }
}  // namespace

l1ct::GctCommonCaloDecoderEmulator::GctCommonCaloDecoderEmulator(const GctEmClusterDecoderEmulator &emDecoder,
                                                                 const GctHadClusterDecoderEmulator &hadDecoder,
                                                                 const std::vector<l1ct::PFRegionEmu> &decodedSectors,
                                                                 const std::vector<l1ct::PFRegionEmu> &linkSectors)
    : emDecoder_(emDecoder), hadDecoder_(hadDecoder), decodedSectors_(decodedSectors), linkSectors_(linkSectors) {
  assert(decodedSectors_.size() == 4 * NLINKS && linkSectors_.size() == NLINKS);
}

l1ct::CommonCaloObjEmu l1ct::GctCommonCaloDecoderEmulator::decode(unsigned int link,
                                                                  unsigned int iclock,
                                                                  const ap_uint<64> &word) const {
  constexpr unsigned int NEM_WORDS = 16;
  constexpr unsigned int NHAD_WORDS = 24;
  // the decoded sector (SLR) of each block of a link is slr_order_per_link[block] + 2 * link; all the SLRs of a link
  // lie in the link sector
  static constexpr unsigned int slr_order_per_link[4] = {7, 1, 6, 0};

  assert(link < NLINKS && iclock < NCLOCKS);
  l1ct::CommonCaloObjEmu ret;
  ret.clear();
  unsigned int rel_pos = iclock % 81;
  if (rel_pos == 0)
    return ret;  // technical words -> ignored
  unsigned int e = rel_pos - 1;
  const bool isEm = e < 2 * NEM_WORDS;
  if (!isEm && e >= 2 * NEM_WORDS + 2 * NHAD_WORDS)
    return ret;

  unsigned int block = isEm ? (e < NEM_WORDS ? 0 : 1) : ((e - 2 * NEM_WORDS) < NHAD_WORDS ? 0 : 1);
  if (iclock > 81)
    block += 2;
  const unsigned int isec = slr_order_per_link[block] + 2 * link;
  const l1ct::PFRegionEmu &decoded = decodedSectors_[isec];
  const l1ct::PFRegionEmu &linkSector = linkSectors_[link];

  if (isEm) {
    if (l1tp2::GCTEmDigiCluster(word).pt() == 0)
      return ret;  // empty slot
    l1ct::EmCaloObjEmu cl = emDecoder_.decode(decoded, word);
    rebase(cl, decoded, linkSector);
    ret.convertFrom(cl);
    ret.src = cl.src;
  } else {
    if (l1tp2::GCTHadDigiCluster(word).pt() == 0)
      return ret;  // empty slot
    l1ct::HadCaloObjEmu cl = hadDecoder_.decode(decoded, word);
    rebase(cl, decoded, linkSector);
    ret.convertFrom(cl);
    ret.src = cl.src;
  }
  return ret;
}

void l1ct::GctCommonCaloDecoderEmulator::decode(
    const std::vector<l1ct::DetectorSector<ap_uint<64>>> &raw,
    std::vector<l1ct::DetectorSector<l1ct::CommonCaloObjEmu>> &gctcommon) const {
  assert(raw.size() == NLINKS);
  gctcommon.resize(raw.size());
  for (unsigned int link = 0; link < NLINKS; ++link) {
    gctcommon[link].region = linkSectors_[link];
    gctcommon[link].clear();
    for (unsigned int iclock = 0, n = raw[link].size(); iclock < n; ++iclock) {
      gctcommon[link].obj.push_back(decode(link, iclock, raw[link][iclock]));
    }
  }
}
