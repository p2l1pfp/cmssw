#include "L1Trigger/Phase2L1ParticleFlow/interface/egamma/pftkegalgo_ref.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/dbgPrintf.h"

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <memory>
#include <iostream>
#include <bitset>

#include "L1Trigger/Phase2L1ParticleFlow/src/dbgPrintf.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
using namespace l1ct;

#ifdef CMSSW_GIT_HASH
#include "FWCore/ParameterSet/interface/ParameterSet.h"

l1ct::PFTkEGAlgoEmuConfig::PFTkEGAlgoEmuConfig(const edm::ParameterSet &pset)
    : nTRACK(pset.getParameter<uint32_t>("nTRACK")),
      nTRACK_EGIN(pset.getParameter<uint32_t>("nTRACK_EGIN")),
      nEMCALO_EGIN(pset.getParameter<uint32_t>("nEMCALO_EGIN")),
      nEM_EGOUT(pset.getParameter<uint32_t>("nEM_EGOUT")),
      filterHwQuality(pset.getParameter<bool>("filterHwQuality")),
      doBremRecovery(pset.getParameter<bool>("doBremRecovery")),
      writeBeforeBremRecovery(pset.getParameter<bool>("writeBeforeBremRecovery")),
      caloHwQual(pset.getParameter<int>("caloHwQual")),
      doEndcapHwQual(pset.getParameter<bool>("doEndcapHwQual")),
      emClusterPtMin(pset.getParameter<double>("caloEtMin")),
      dEtaMaxBrem(pset.getParameter<double>("dEtaMaxBrem")),
      dPhiMaxBrem(pset.getParameter<double>("dPhiMaxBrem")),
      absEtaBoundaries(pset.getParameter<std::vector<double>>("absEtaBoundaries")),
      dEtaValues(pset.getParameter<std::vector<double>>("dEtaValues")),
      dPhiValues(pset.getParameter<std::vector<double>>("dPhiValues")),
      trkQualityPtMin(pset.getParameter<double>("trkQualityPtMin")),
      writeEgSta(pset.getParameter<bool>("writeEGSta")),
      tkIsoParams_tkEle(pset.getParameter<edm::ParameterSet>("tkIsoParametersTkEle")),
      tkIsoParams_tkEm(pset.getParameter<edm::ParameterSet>("tkIsoParametersTkEm")),
      pfIsoParams_tkEle(pset.getParameter<edm::ParameterSet>("pfIsoParametersTkEle")),
      pfIsoParams_tkEm(pset.getParameter<edm::ParameterSet>("pfIsoParametersTkEm")),
      doTkIso(pset.getParameter<bool>("doTkIso")),
      doPfIso(pset.getParameter<bool>("doPfIso")),
      hwIsoTypeTkEle(static_cast<EGIsoEleObjEmu::IsoType>(pset.getParameter<uint32_t>("hwIsoTypeTkEle"))),
      hwIsoTypeTkEm(static_cast<EGIsoObjEmu::IsoType>(pset.getParameter<uint32_t>("hwIsoTypeTkEm"))),
      doCompositeTkEle(pset.getParameter<bool>("doCompositeTkEle")),
      myCompIDparams(pset.getParameter<edm::ParameterSet>("compositeParametersTkEle")),
      debug(pset.getUntrackedParameter<uint32_t>("debug", 0)) {}

l1ct::PFTkEGAlgoEmuConfig::IsoParameters::IsoParameters(const edm::ParameterSet &pset)
    : IsoParameters(pset.getParameter<double>("tkQualityPtMin"),
                    pset.getParameter<double>("dZ"),
                    pset.getParameter<double>("dRMin"),
                    pset.getParameter<double>("dRMax")) {}

l1ct::PFTkEGAlgoEmuConfig::CompIDParameters::CompIDParameters(const edm::ParameterSet &pset)
    : CompIDParameters(pset.getParameter<double>("hoeMin"),
                       pset.getParameter<double>("hoeMax"),
                       pset.getParameter<double>("tkptMin"),
                       pset.getParameter<double>("tkptMax"),
                       pset.getParameter<double>("srrtotMin"),
                       pset.getParameter<double>("srrtotMax"),
                       pset.getParameter<double>("detaMin"),
                       pset.getParameter<double>("detaMax"),
                       pset.getParameter<double>("dptMin"),
                       pset.getParameter<double>("dptMax"),
                       pset.getParameter<double>("meanzMin"),
                       pset.getParameter<double>("meanzMax"),
                       pset.getParameter<double>("dphiMin"),
                       pset.getParameter<double>("dphiMax"),
                       pset.getParameter<double>("tkchi2Min"),
                       pset.getParameter<double>("tkchi2Max"),
                       pset.getParameter<double>("tkz0Min"),
                       pset.getParameter<double>("tkz0Max"),
                       pset.getParameter<double>("tknstubsMin"),
                       pset.getParameter<double>("tknstubsMax")) {}

#endif

PFTkEGAlgoEmulator::PFTkEGAlgoEmulator(const PFTkEGAlgoEmuConfig &config) : cfg(config), 
composite_bdt_(nullptr), 
debug_(cfg.debug) {
  if(cfg.doCompositeTkEle) {
    //FIXME: make the name of the file configurable
    auto resolvedFileName = edm::FileInPath("L1Trigger/Phase2L1ParticleFlow/src/newfirmware/egamma/compositeID.json").fullPath();
    composite_bdt_ = std::make_unique<conifer::BDT<ap_fixed<22,3,AP_RND_CONV,AP_SAT>,ap_fixed<22,3,AP_RND_CONV,AP_SAT>,0>> (resolvedFileName);
  }
}

void PFTkEGAlgoEmulator::toFirmware(const PFInputRegion &in,
                                    PFRegion &region,
                                    EmCaloObj emcalo[/*nCALO*/],
                                    TkObj track[/*nTRACK*/]) const {
  region = in.region;
  l1ct::toFirmware(in.track, cfg.nTRACK_EGIN, track);
  l1ct::toFirmware(in.emcalo, cfg.nEMCALO_EGIN, emcalo);
  if (debug_ > 0)
    dbgCout() << "# of inpput tracks: " << in.track.size() << " (max: " << cfg.nTRACK_EGIN << ")"
              << " emcalo: " << in.emcalo.size() << "(" << cfg.nEMCALO_EGIN << ")" << std::endl;
}

void PFTkEGAlgoEmulator::toFirmware(const OutputRegion &out, EGIsoObj out_egphs[], EGIsoEleObj out_egeles[]) const {
  l1ct::toFirmware(out.egphoton, cfg.nEM_EGOUT, out_egphs);
  l1ct::toFirmware(out.egelectron, cfg.nEM_EGOUT, out_egeles);
  if (debug_ > 0)
    dbgCout() << "# output photons: " << out.egphoton.size() << " electrons: " << out.egelectron.size() << std::endl;
}

void PFTkEGAlgoEmulator::toFirmware(
    const PFInputRegion &in, const l1ct::PVObjEmu &pvin, PFRegion &region, TkObj track[/*nTRACK*/], PVObj &pv) const {
  region = in.region;
  l1ct::toFirmware(in.track, cfg.nTRACK, track);
  pv = pvin;
  if (debug_ > 0)
    dbgCout() << "# of inpput tracks: " << in.track.size() << " (max: " << cfg.nTRACK << ")" << std::endl;
}

float PFTkEGAlgoEmulator::deltaPhi(float phi1, float phi2) const {
  // reduce to [-pi,pi]
  float x = phi1 - phi2;
  float o2pi = 1. / (2. * M_PI);
  if (std::abs(x) <= float(M_PI))
    return x;
  float n = std::round(x * o2pi);
  return x - n * float(2. * M_PI);
}

void PFTkEGAlgoEmulator::link_emCalo2emCalo(const std::vector<EmCaloObjEmu> &emcalo,
                                            std::vector<int> &emCalo2emCalo) const {
  // NOTE: we assume the input to be sorted!!!
  for (int ic = 0, nc = emcalo.size(); ic < nc; ++ic) {
    auto &calo = emcalo[ic];
    if (emCalo2emCalo[ic] != -1)
      continue;

    for (int jc = ic + 1; jc < nc; ++jc) {
      if (emCalo2emCalo[jc] != -1)
        continue;

      auto &otherCalo = emcalo[jc];

      if (fabs(otherCalo.floatEta() - calo.floatEta()) < cfg.dEtaMaxBrem &&
          fabs(deltaPhi(otherCalo.floatPhi(), calo.floatPhi())) < cfg.dPhiMaxBrem) {
        emCalo2emCalo[jc] = ic;
      }
    }
  }
}


void PFTkEGAlgoEmulator::link_emCalo2tk_elliptic(const PFRegionEmu &r,
                                                 const std::vector<EmCaloObjEmu> &emcalo,
                                                 const std::vector<TkObjEmu> &track,
                                                 std::vector<int> &emCalo2tk) const {
  unsigned int nTrackMax = std::min<unsigned>(track.size(), cfg.nTRACK_EGIN);
  for (int ic = 0, nc = emcalo.size(); ic < nc; ++ic) {
    auto &calo = emcalo[ic];

    float dPtMin = 999;
    for (unsigned int itk = 0; itk < nTrackMax; ++itk) {
      const auto &tk = track[itk];
      if (tk.floatPt() < cfg.trkQualityPtMin)
        continue;

      float d_phi = deltaPhi(tk.floatPhi(), calo.floatPhi());
      float d_eta = tk.floatEta() - calo.floatEta();  // We only use it squared

      auto eta_index =
          std::distance(cfg.absEtaBoundaries.begin(),
                        std::lower_bound(
                            cfg.absEtaBoundaries.begin(), cfg.absEtaBoundaries.end(), abs(r.floatGlbEta(calo.hwEta)))) -
          1;

      float dEtaMax = cfg.dEtaValues[eta_index];
      float dPhiMax = cfg.dPhiValues[eta_index];

      if ((((d_phi / dPhiMax) * (d_phi / dPhiMax)) + ((d_eta / dEtaMax) * (d_eta / dEtaMax))) < 1.) {
        // NOTE: for now we implement only best pt match. This is NOT what is done in the L1TkElectronTrackProducer
        if (fabs(tk.floatPt() - calo.floatPt()) < dPtMin) {
          emCalo2tk[ic] = itk;
          dPtMin = fabs(tk.floatPt() - calo.floatPt());
        }
      }
    }
  }
}


void PFTkEGAlgoEmulator::link_emCalo2tk_composite(const PFRegionEmu &r,
                                        const std::vector<EmCaloObjEmu> &emcalo,
                                        const std::vector<TkObjEmu> &track,
                                        std::vector<int> &emCalo2tk,
                                        const PFTkEGAlgoEmuConfig::CompIDParameters &params) const {
  //FIXME: should be configurable
  const int nCAND_PER_CLUSTER = 4;
  unsigned int nTrackMax = std::min<unsigned>(track.size(), cfg.nTRACK_EGIN);
  for (int ic = 0, nc = emcalo.size(); ic < nc; ++ic) {
    auto &calo = emcalo[ic];

    std::vector<CompositeCandidate> candidates;

    for (unsigned int itk = 0; itk < nTrackMax; ++itk) {
          const auto &tk = track[itk];
          if (tk.floatPt() < cfg.trkQualityPtMin)
            continue;

      // float d_phi = deltaPhi(tk.floatPhi(), calo.floatPhi());
      // float d_eta = tk.floatEta() - calo.floatEta();  // We only use it squared
      float clu_eta=calo.floatEta();
      float clu_phi=calo.floatPhi();
      float trk_eta=tk.floatEta();
      float trk_phi=tk.floatPhi();
      float dR = deltaR(clu_eta,clu_phi,trk_eta,trk_phi);

      // if (((d_phi * d_phi ) + (d_eta * d_eta )) < 0.2 * 0.2) {
      if (dR<0.2){
          CompositeCandidate cand;
          cand.cluster_idx = ic;
          cand.track_idx = itk;
          // cand.dR = sqrt((d_phi * d_phi ) + (d_eta * d_eta ));
          cand.dR = dR;
          cand.dpt_double = tk.floatPt()/calo.floatPt();
          // Normalize feature values
          cand.hoe = (calo.src->hOverE()-params.hoeMin)/(params.hoeMax-params.hoeMin);
          cand.tkpt = (tk.floatPt()-params.tkptMin)/(params.tkptMax-params.tkptMin);
          cand.srrtot = (dynamic_cast<const l1t::HGCalMulticluster *>(calo.src->constituentsAndFractions()[0].first.get())->sigmaRRTot()-params.srrtotMin)/(params.srrtotMax-params.srrtotMin); //FIXME: Get the srrtot from the cluster
          std::cout << "SRRTOT: " << dynamic_cast<const l1t::HGCalMulticluster *>(calo.src->constituentsAndFractions()[0].first.get())->sigmaRRTot() << std::endl;
          cand.deta = (tk.src->caloEta()-calo.floatEta()-params.detaMin)/(params.detaMax-params.detaMin);
          cand.dpt = ((tk.floatPt()/calo.floatPt())-params.dptMin)/(params.dptMax-params.dptMin);
          cand.meanz = (dynamic_cast<const l1t::HGCalMulticluster *>(calo.src->constituentsAndFractions()[0].first.get())->zBarycenter()-params.meanzMin)/(params.meanzMax-params.meanzMin); //FIXME: Get the meanz from the cluster
          cand.dphi = (tk.src->caloPhi()- calo.floatPhi() -params.dphiMin)/(params.dphiMax-params.dphiMin);
          cand.chi2 = (tk.src->chi2()-params.tkchi2Min)/(params.tkchi2Max-params.tkchi2Min);
          cand.tkz0 = (tk.floatZ0()-params.tkz0Min)/(params.tkz0Max-params.tkz0Min);
          cand.nstubs = (tk.src->nStubs()-params.tknstubsMin)/(params.tknstubsMax-params.tknstubsMin);
          candidates.push_back(cand);
      }
    }
    // FIXME: find best sort criteria, for now we use dpt
    std::sort(candidates.begin(), candidates.end(), 
              [](const CompositeCandidate & a, const CompositeCandidate & b) -> bool
                { return abs(a.dpt_double-1) < abs(b.dpt_double-1); });
    unsigned int nCandPerCluster = std::min<unsigned int>(candidates.size(), nCAND_PER_CLUSTER);
    std::cout << "# composit candidates: " << nCandPerCluster << std::endl;
    if(nCandPerCluster == 0) continue;

    float bdtWP = 10;
    float minScore = 999;
    int ibest = -1;
    for(unsigned int icand = 0; icand < nCandPerCluster; icand++) {
      auto &cand = candidates[icand];
      float score = compute_composite_score(cand);
      if((score < bdtWP) && (score < minScore)) {
        minScore = score;
        ibest = icand;
      }
    }    
    if(ibest != -1) emCalo2tk[ic] = candidates[ibest].track_idx;
  }
}


float PFTkEGAlgoEmulator::compute_composite_score(const CompositeCandidate& cand) const {
  vector<ap_fixed<22,3,AP_RND_CONV,AP_SAT>> inputs = { cand.hoe, cand.tkpt, cand.srrtot, cand.deta, cand.dpt, cand.meanz, cand.dphi, cand.chi2, cand.tkz0, cand.nstubs } ;
  auto bdt_score = composite_bdt_->decision_function(inputs);
  std::cout<<"BDT score or composite candidate = "<<bdt_score[0]<<std::endl;
  return bdt_score[0];
}


void PFTkEGAlgoEmulator::sel_emCalo(unsigned int nmax_sel,
                                    const std::vector<EmCaloObjEmu> &emcalo,
                                    std::vector<EmCaloObjEmu> &emcalo_sel) const {
  for (int ic = 0, nc = emcalo.size(); ic < nc; ++ic) {
    const auto &calo = emcalo[ic];
    if ((calo.hwPt == 0) || (cfg.filterHwQuality && calo.hwEmID != cfg.caloHwQual) ||
        (calo.floatPt() < cfg.emClusterPtMin))
      continue;
    emcalo_sel.push_back(calo);
    if (emcalo_sel.size() >= nmax_sel)
      break;
  }
}

void PFTkEGAlgoEmulator::run(const PFInputRegion &in, OutputRegion &out) const {
  if (debug_ > 1) {
    for (int ic = 0, nc = in.emcalo.size(); ic < nc; ++ic) {
      const auto &calo = in.emcalo[ic];
      if (calo.hwPt > 0)
        dbgCout() << "[REF] IN calo[" << ic << "] pt: " << calo.hwPt << " eta: " << calo.hwEta
                  << " (glb eta: " << in.region.floatGlbEta(calo.hwEta) << ") phi: " << calo.hwPhi
                  << "(glb phi: " << in.region.floatGlbPhi(calo.hwPhi) << ") qual: " << std::bitset<4>(calo.hwEmID)
                  << std::endl;
    }
  }

  // FIXME: can be removed in the endcap since now running with the "interceptor".
  // Might still be needed in barrel
  // filter and select first N elements of input clusters
  std::vector<EmCaloObjEmu> emcalo_sel;
  sel_emCalo(cfg.nEMCALO_EGIN, in.emcalo, emcalo_sel);

  std::vector<int> emCalo2emCalo(emcalo_sel.size(), -1);
  if (cfg.doBremRecovery)
    link_emCalo2emCalo(emcalo_sel, emCalo2emCalo);

  std::vector<int> emCalo2tk(emcalo_sel.size(), -1);
  if(cfg.doCompositeTkEle) {
    link_emCalo2tk_composite(in.region, emcalo_sel, in.track, emCalo2tk, cfg.myCompIDparams);
  } else {
    link_emCalo2tk_elliptic(in.region, emcalo_sel, in.track, emCalo2tk);
  }
  
  out.egsta.clear();
  std::vector<EGIsoObjEmu> egobjs;
  std::vector<EGIsoEleObjEmu> egeleobjs;
  eg_algo(in.region, emcalo_sel, in.track, emCalo2emCalo, emCalo2tk, out.egsta, egobjs, egeleobjs);

  unsigned int nEGOut = std::min<unsigned>(cfg.nEM_EGOUT, egobjs.size());
  unsigned int nEGEleOut = std::min<unsigned>(cfg.nEM_EGOUT, egeleobjs.size());

  // init output containers
  out.egphoton.clear();
  out.egelectron.clear();
  ptsort_ref(egobjs.size(), nEGOut, egobjs, out.egphoton);
  ptsort_ref(egeleobjs.size(), nEGEleOut, egeleobjs, out.egelectron);
}

void PFTkEGAlgoEmulator::eg_algo(const PFRegionEmu &region,
                                 const std::vector<EmCaloObjEmu> &emcalo,
                                 const std::vector<TkObjEmu> &track,
                                 const std::vector<int> &emCalo2emCalo,
                                 const std::vector<int> &emCalo2tk,
                                 std::vector<EGObjEmu> &egstas,
                                 std::vector<EGIsoObjEmu> &egobjs,
                                 std::vector<EGIsoEleObjEmu> &egeleobjs) const {
  for (int ic = 0, nc = emcalo.size(); ic < nc; ++ic) {
    auto &calo = emcalo[ic];

    // discard immediately EG objects that would not fall in the fiducial eta-phi region
    if (!region.isFiducial(calo))
      continue;

    if (debug_ > 3)
      dbgCout() << "[REF] SEL emcalo with pt: " << calo.hwPt << " qual: " << calo.hwEmID << " eta: " << calo.hwEta
                << " phi " << calo.hwPhi << std::endl;

    int itk = emCalo2tk[ic];

    // check if brem recovery is on
    if (!cfg.doBremRecovery || cfg.writeBeforeBremRecovery) {
      // 1. create EG objects before brem recovery
      unsigned int egQual = calo.hwEmID;
      // If we write both objects with and without brem-recovery
      // bit 3 is used for the brem-recovery bit: if set = no recovery
      // (for consistency with the barrel hwQual where by default the brem recovery is done upstream)
      if (cfg.writeBeforeBremRecovery && cfg.doBremRecovery) {
        egQual = calo.hwEmID | 0x8;
      }

      addEgObjsToPF(egstas, egobjs, egeleobjs, emcalo, track, ic, egQual, calo.hwPt, itk);
    }

    if (!cfg.doBremRecovery)
      continue;

    // check if the cluster has already been used in a brem reclustering
    if (emCalo2emCalo[ic] != -1)
      continue;

    pt_t ptBremReco = calo.hwPt;
    std::vector<unsigned int> components;

    for (int jc = ic; jc < nc; ++jc) {
      if (emCalo2emCalo[jc] == ic) {
        auto &otherCalo = emcalo[jc];
        ptBremReco += otherCalo.hwPt;
        components.push_back(jc);
      }
    }

    // 2. create EG objects with brem recovery
    // NOTE: duplicating the object is suboptimal but this is done for keeping things as in TDR code...
    addEgObjsToPF(egstas, egobjs, egeleobjs, emcalo, track, ic, calo.hwEmID, ptBremReco, itk, components);
  }
}

EGObjEmu &PFTkEGAlgoEmulator::addEGStaToPF(std::vector<EGObjEmu> &egobjs,
                                           const EmCaloObjEmu &calo,
                                           const unsigned int hwQual,
                                           const pt_t ptCorr,
                                           const std::vector<unsigned int> &components) const {
  EGObjEmu egsta;
  egsta.clear();
  egsta.hwPt = ptCorr;
  egsta.hwEta = calo.hwEta;
  egsta.hwPhi = calo.hwPhi;
  egsta.hwQual = hwQual;
  egobjs.push_back(egsta);

  if (debug_ > 2)
    dbgCout() << "[REF] EGSta pt: " << egsta.hwPt << " eta: " << egsta.hwEta << " phi: " << egsta.hwPhi
              << " qual: " << std::bitset<4>(egsta.hwQual) << " packed: " << egsta.pack().to_string(16) << std::endl;

  return egobjs.back();
}

EGIsoObjEmu &PFTkEGAlgoEmulator::addEGIsoToPF(std::vector<EGIsoObjEmu> &egobjs,
                                              const EmCaloObjEmu &calo,
                                              const unsigned int hwQual,
                                              const pt_t ptCorr) const {
  EGIsoObjEmu egiso;
  egiso.clear();
  egiso.hwPt = ptCorr;
  egiso.hwEta = calo.hwEta;
  egiso.hwPhi = calo.hwPhi;
  unsigned int egHwQual = hwQual;
  if (cfg.doEndcapHwQual) {
    // 1. zero-suppress the loose EG-ID (bit 1)
    // 2. for now use the standalone tight definition (bit 0) to set the tight point for photons (bit 2)
    egHwQual = (hwQual & 0x9) | ((hwQual & 0x1) << 2);
  }
  egiso.hwQual = egHwQual;
  egiso.srcCluster = calo.src;
  egobjs.push_back(egiso);

  if (debug_ > 2)
    dbgCout() << "[REF] EGIsoObjEmu pt: " << egiso.hwPt << " eta: " << egiso.hwEta << " phi: " << egiso.hwPhi
              << " qual: " << std::bitset<4>(egiso.hwQual) << " packed: " << egiso.pack().to_string(16) << std::endl;

  return egobjs.back();
}

EGIsoEleObjEmu &PFTkEGAlgoEmulator::addEGIsoEleToPF(std::vector<EGIsoEleObjEmu> &egobjs,
                                                    const EmCaloObjEmu &calo,
                                                    const TkObjEmu &track,
                                                    const unsigned int hwQual,
                                                    const pt_t ptCorr) const {
  EGIsoEleObjEmu egiso;
  egiso.clear();
  egiso.hwPt = ptCorr;
  egiso.hwEta = calo.hwEta;
  egiso.hwPhi = calo.hwPhi;
  unsigned int egHwQual = hwQual;
  if (cfg.doEndcapHwQual) {
    // 1. zero-suppress the loose EG-ID (bit 1)
    // 2. for now use the standalone tight definition (bit 0) to set the tight point for eles (bit 1)
    egHwQual = (hwQual & 0x9) | ((hwQual & 0x1) << 1);
  }
  egiso.hwQual = egHwQual;
  egiso.hwDEta = track.hwVtxEta() - egiso.hwEta;
  egiso.hwDPhi = abs(track.hwVtxPhi() - egiso.hwPhi);
  egiso.hwZ0 = track.hwZ0;
  egiso.hwCharge = track.hwCharge;
  egiso.srcCluster = calo.src;
  egiso.srcTrack = track.src;
  egobjs.push_back(egiso);

  if (debug_ > 2)
    dbgCout() << "[REF] EGIsoEleObjEmu pt: " << egiso.hwPt << " eta: " << egiso.hwEta << " phi: " << egiso.hwPhi
              << " qual: " << std::bitset<4>(egiso.hwQual) << " packed: " << egiso.pack().to_string(16) << std::endl;

  return egobjs.back();
}

void PFTkEGAlgoEmulator::addEgObjsToPF(std::vector<EGObjEmu> &egstas,
                                       std::vector<EGIsoObjEmu> &egobjs,
                                       std::vector<EGIsoEleObjEmu> &egeleobjs,
                                       const std::vector<EmCaloObjEmu> &emcalo,
                                       const std::vector<TkObjEmu> &track,
                                       const int calo_idx,
                                       const unsigned int hwQual,
                                       const pt_t ptCorr,
                                       const int tk_idx,
                                       const std::vector<unsigned int> &components) const {
  int sta_idx = -1;
  if (writeEgSta()) {
    addEGStaToPF(egstas, emcalo[calo_idx], hwQual, ptCorr, components);
    sta_idx = egstas.size() - 1;
  }
  EGIsoObjEmu &egobj = addEGIsoToPF(egobjs, emcalo[calo_idx], hwQual, ptCorr);
  egobj.sta_idx = sta_idx;
  if (tk_idx != -1) {
    EGIsoEleObjEmu &eleobj = addEGIsoEleToPF(egeleobjs, emcalo[calo_idx], track[tk_idx], hwQual, ptCorr);
    eleobj.sta_idx = sta_idx;
  }
}

void PFTkEGAlgoEmulator::runIso(const PFInputRegion &in,
                                const std::vector<l1ct::PVObjEmu> &pvs,
                                OutputRegion &out) const {
  if (cfg.doTkIso) {
    compute_isolation(out.egelectron, in.track, cfg.tkIsoParams_tkEle, pvs[0].hwZ0);
    compute_isolation(out.egphoton, in.track, cfg.tkIsoParams_tkEm, pvs[0].hwZ0);
  }
  if (cfg.doPfIso) {
    compute_isolation(out.egelectron, out.pfcharged, out.pfneutral, cfg.pfIsoParams_tkEle, pvs[0].hwZ0);
    compute_isolation(out.egphoton, out.pfcharged, out.pfneutral, cfg.pfIsoParams_tkEm, pvs[0].hwZ0);
  }

  std::for_each(out.egelectron.begin(), out.egelectron.end(), [&](EGIsoEleObjEmu &obj) {
    obj.hwIso = obj.hwIsoVar(cfg.hwIsoTypeTkEle);
  });
  std::for_each(
      out.egphoton.begin(), out.egphoton.end(), [&](EGIsoObjEmu &obj) { obj.hwIso = obj.hwIsoVar(cfg.hwIsoTypeTkEm); });
}

void PFTkEGAlgoEmulator::compute_isolation(std::vector<EGIsoObjEmu> &egobjs,
                                           const std::vector<TkObjEmu> &objects,
                                           const PFTkEGAlgoEmuConfig::IsoParameters &params,
                                           z0_t z0) const {
  for (int ic = 0, nc = egobjs.size(); ic < nc; ++ic) {
    auto &egphoton = egobjs[ic];
    iso_t sumPt = 0.;
    iso_t sumPtPV = 0.;
    compute_sumPt(sumPt, sumPtPV, objects, cfg.nTRACK, egphoton, params, z0);
    egphoton.setHwIso(EGIsoObjEmu::IsoType::TkIso, sumPt);
    egphoton.setHwIso(EGIsoObjEmu::IsoType::TkIsoPV, sumPtPV);
  }
}

void PFTkEGAlgoEmulator::compute_isolation(std::vector<EGIsoEleObjEmu> &egobjs,
                                           const std::vector<TkObjEmu> &objects,
                                           const PFTkEGAlgoEmuConfig::IsoParameters &params,
                                           z0_t z0) const {
  for (int ic = 0, nc = egobjs.size(); ic < nc; ++ic) {
    auto &egele = egobjs[ic];
    iso_t sumPt = 0.;
    iso_t sumPtPV = 0.;
    compute_sumPt(sumPt, sumPtPV, objects, cfg.nTRACK, egele, params, z0);
    egele.setHwIso(EGIsoEleObjEmu::IsoType::TkIso, sumPtPV);
  }
}

void PFTkEGAlgoEmulator::compute_isolation(std::vector<EGIsoObjEmu> &egobjs,
                                           const std::vector<PFChargedObjEmu> &charged,
                                           const std::vector<PFNeutralObjEmu> &neutrals,
                                           const PFTkEGAlgoEmuConfig::IsoParameters &params,
                                           z0_t z0) const {
  for (int ic = 0, nc = egobjs.size(); ic < nc; ++ic) {
    auto &egphoton = egobjs[ic];
    iso_t sumPt = 0.;
    iso_t sumPtPV = 0.;
    // FIXME: set max # of PF objects for iso
    compute_sumPt(sumPt, sumPtPV, charged, charged.size(), egphoton, params, z0);
    compute_sumPt(sumPt, sumPtPV, neutrals, neutrals.size(), egphoton, params, z0);
    egphoton.setHwIso(EGIsoObjEmu::IsoType::PfIso, sumPt);
    egphoton.setHwIso(EGIsoObjEmu::IsoType::PfIsoPV, sumPtPV);
  }
}

void PFTkEGAlgoEmulator::compute_isolation(std::vector<EGIsoEleObjEmu> &egobjs,
                                           const std::vector<PFChargedObjEmu> &charged,
                                           const std::vector<PFNeutralObjEmu> &neutrals,
                                           const PFTkEGAlgoEmuConfig::IsoParameters &params,
                                           z0_t z0) const {
  for (int ic = 0, nc = egobjs.size(); ic < nc; ++ic) {
    auto &egele = egobjs[ic];
    iso_t sumPt = 0.;
    iso_t sumPtPV = 0.;
    compute_sumPt(sumPt, sumPtPV, charged, charged.size(), egele, params, z0);
    compute_sumPt(sumPt, sumPtPV, neutrals, neutrals.size(), egele, params, z0);
    egele.setHwIso(EGIsoEleObjEmu::IsoType::PfIso, sumPtPV);
  }
}
