#include "L1Trigger/Phase2L1ParticleFlow/interface/l1-converters/hgcalinput_ref.h"

#ifdef CMSSW_GIT_HASH
#include "FWCore/ParameterSet/interface/ParameterSet.h"
l1ct::HgcalClusterDecoderEmulator::HgcalClusterDecoderEmulator(const edm::ParameterSet &pset)
    : slim_(pset.getParameter<bool>("slim")) {}

#endif

l1ct::HgcalClusterDecoderEmulator::~HgcalClusterDecoderEmulator() {}

l1ct::HadCaloObjEmu l1ct::HgcalClusterDecoderEmulator::decode(const l1ct::PFRegionEmu &sector,
                                                              const ap_uint<256> &in) const {
  ap_uint<14> w_pt = in(13, 0);
  ap_uint<14> w_empt = in(27, 14);
  ap_uint<10> w_abseta = in(73, 64);
  ap_int<9> w_eta = l1ct::glbeta_t(w_abseta.to_int() * (sector.floatEtaCenter() > 0 ? +1 : -1)) - sector.hwEtaCenter;
  ap_int<9> w_phi = in(82, 74);
  ap_uint<4> w_gctqual = in(31, 28);  // GCT quality
  ap_uint<7> w_srrtot = in(191, 185);
  ap_uint<12> w_meanz = in(94, 83);
  ap_uint<8> w_emf = in(39, 32);
  // FIXME: we use a spare space in the word for hoe which is not in the current interface
  //ap_uint<12> w_hoe = in(127, 116);
  // Compute an H/E value: 1/emf - 1
  ap_ufixed<10, 5, AP_RND_CONV, AP_SAT> w_hoe = 256.0 / (w_emf.to_int() + 0.5) - 1;

  l1ct::HadCaloObjEmu out;
  out.clear();
  if (w_pt == 0)
    return out;
  out.hwPt = w_pt * l1ct::pt_t(l1ct::Scales::INTPT_LSB);
  out.hwEta = w_eta;
  out.hwPhi = w_phi;  // relative to the region center, at calo
  out.hwEmPt = w_empt * l1ct::pt_t(l1ct::Scales::INTPT_LSB);

  out.hwEmID = 0;
  out.hwEmID[0] = w_gctqual[0] && (out.floatEmPt() > 1.0) && w_emf > int(0.5 * 256);  // PF ID
  out.hwEmID[2] = w_gctqual[0] && (out.floatEmPt() > 3.0) && w_emf > int(0.7 * 256);  // LOOSE ID
  out.hwEmID[1] = w_gctqual[1] && (out.floatEmPt() > 5.0) && w_emf > int(0.8 * 256);  // TIGHT ID
  if (!slim_) {
    constexpr float HGCAL_SRRTOT_LSB = 0.024584 / (1 << 7);
    // -- integer version ---
    //float SRRTOT_FACTOR_FLT = HGCAL_SRRTOT_LSB*l1ct::Scales::SRRTOT_SCALE/l1ct::Scales::SRRTOT_LSB;
    //ap_uint<10> our_srrtot_raw = (w_srrtot.to_int() * int(round(SRRTOT_FACTOR_FLT*(1<<4))))>>4;;
    //out.hwSrrTot(9,0) = our_srrtot_raw(9,0);
    // --- ap_fixed version (note that our LSB doesn't appear) ---
    const ap_ufixed<8, -6> SRRTOT_FACTOR = HGCAL_SRRTOT_LSB * l1ct::Scales::SRRTOT_SCALE;
    out.hwSrrTot = w_srrtot * SRRTOT_FACTOR;
    out.hwMeanZ = l1ct::meanz_t(std::min(w_meanz.to_int() + 4, (1 << 12) - 1) >> 3);
    out.hwHoe = w_hoe;
  }
  return out;
}
