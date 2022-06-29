#include "pftkegalgo_ref.h"

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <memory>
#include <iostream>
#include <bitset>

#include "L1Trigger/Phase2L1ParticleFlow/src/dbgPrintf.h"

using namespace l1ct;

#ifdef CMSSW_GIT_HASH
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

l1ct::CompositeIDAlgoEmuConfig::CompositeIDAlgoEmuConfig(const edm::ParameterSet &pset)
    : nTRACK(pset.getParameter<uint32_t>("nTRACK")),
      nTRACK_EGIN(pset.getParameter<uint32_t>("nTRACK_EGIN")),
      nEMCALO_EGIN(pset.getParameter<uint32_t>("nEMCALO_EGIN")),
      nEM_EGOUT(pset.getParameter<uint32_t>("nEM_EGOUT")),
      filterHwQuality(pset.getParameter<bool>("filterHwQuality")),//yes 
      doBremRecovery(pset.getParameter<bool>("doBremRecovery")),//yes
      writeBeforeBremRecovery(pset.getParameter<bool>("writeBeforeBremRecovery")),
      caloHwQual(pset.getParameter<int>("caloHwQual")),//yes
      doEndcapHwQual(pset.getParameter<bool>("doEndcapHwQual")),
      emClusterPtMin(pset.getParameter<double>("caloEtMin")),//yes
      dEtaMaxBrem(pset.getParameter<double>("dEtaMaxBrem")),//yes
      dPhiMaxBrem(pset.getParameter<double>("dPhiMaxBrem")),//yes
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
      debug(pset.getUntrackedParameter<uint32_t>("debug", 0)) {}

l1ct::CompositeIDAlgoEmuConfig::IsoParameters::IsoParameters(const edm::ParameterSet &pset)
    : IsoParameters(pset.getParameter<double>("tkQualityPtMin"),
                    pset.getParameter<double>("dZ"),
                    pset.getParameter<double>("dRMin"),
                    pset.getParameter<double>("dRMax")) {}

#endif


float CompositeIDAlgoEmulator::deltaPhi(float phi1, float phi2) const {
  // reduce to [-pi,pi]
  float x = phi1 - phi2;
  float o2pi = 1. / (2. * M_PI);
  if (std::abs(x) <= float(M_PI))
    return x;
  float n = std::round(x * o2pi);
  return x - n * float(2. * M_PI);
}

void CompositeIDAlgoEmulator::link_emCalo2emCalo(const std::vector<EmCaloObjEmu> &emcalo,
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

void CompositeIDAlgoEmulator::link_emCalo2tk(const PFRegionEmu &r,
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

void CompositeIDAlgoEmulator::sel_emCalo(unsigned int nmax_sel,
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

void CompositeIDAlgoEmulator::run(const PFInputRegion &in, OutputRegion &out) const {
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
  link_emCalo2tk(in.region, emcalo_sel, in.track, emCalo2tk);

  out.egsta.clear();
  std::vector<EGIsoObjEmu> egobjs;
  std::vector<EGIsoEleObjEmu> egeleobjs;
  comp_algo(in.region, emcalo_sel, in.track, emCalo2emCalo, emCalo2tk, out.egsta, egobjs, egeleobjs);

  // unsigned int nEGOut = std::min<unsigned>(cfg.nEM_EGOUT, egobjs.size());
  // unsigned int nEGEleOut = std::min<unsigned>(cfg.nEM_EGOUT, egeleobjs.size());

  // // init output containers
  // out.egphoton.clear();
  // out.egelectron.clear();
  // ptsort_ref(egobjs.size(), nEGOut, egobjs, out.egphoton);
  // ptsort_ref(egeleobjs.size(), nEGEleOut, egeleobjs, out.egelectron);
}

void CompositeIDAlgoEmulator::comp_algo(const PFRegionEmu &region,
                                 const std::vector<EmCaloObjEmu> &emcalo,
                                 const std::vector<TkObjEmu> &track,
                                 const std::vector<int> &emCalo2emCalo,
                                 const std::vector<int> &emCalo2tk,
                                 std::vector<EGObjEmu> &egstas,
                                 std::vector<EGIsoObjEmu> &egobjs,
                                 std::vector<EGIsoEleObjEmu> &egeleobjs) const {

  // LOAD THE COMPOSITE ID BDT
  auto resolvedFileName = edm::FileInPath("compositeID.json").fullPath();
  conifer::BDT<double,double,0> bdt(resolvedFileName);

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

      for (int it = 0, ntk = track.size(); it < ntk; ++it) {
        auto &trk = track[it];

      // FORM THE COMPOSITE CANDIDATE
        float clu_eta=calo.floatEta();
        float clu_phi=calo.floatPhi();
        // float trk_eta=trk.floatPhi();  #Check difference with caloPhi/Eta
        // float trk_phi=trk.floatPhi();
        float trk_eta=trk.src->caloEta();
        float trk_phi=trk.src->caloPhi();

        float dR = deltaR(clu_eta,clu_phi,trk_eta,trk_phi);
        if (dR<0.2){
          // FEATURE SCALING (TO 0-1) BASED ON PARAMETERS FROM BDT TRAINING
          // features=['hoe','tkpt','srrtot','deta','dpt','meanz','dphi','tkchi2','tkz0','tknstubs']
          float hoe = (calo.src->hOverE()+1.0)/(1566.547607421875+1);
          float tkpt = (trk.floatPt()-1.9501149654388428)/(11102.0048828125-1.9501149654388428);
          float srrtot = (0.-0.0)/(0.01274710614234209-0.0);
          float deta = (trk_eta-clu_eta+0.24224889278411865)/(0.23079538345336914+0.24224889278411865);
          float dpt = ((trk.floatPt()/calo.floatPt())-0.010325592942535877)/(184.92538452148438-0.010325592942535877);
          float meanz = (0.-325.0653991699219)/(499.6089782714844-325.0653991699219);
          float dphi = (trk_phi-clu_phi+6.281332015991211)/(6.280326843261719+6.281332015991211);
          float tkchi2 =(trk.src->chi2()-0.024048099294304848)/(1258.37158203125-0.024048099294304848);
          float tkz0 = (trk.floatZ0()+14.94140625)/(14.94140625+14.94140625);
          float tknstubs = (trk.src->nStubs()-4.0)/(6.0-4.0);

          std::cout<<hoe<<"\t"<<tkpt<<"\t"<<srrtot<<"\t"<<deta<<"\t"<<dpt<<"\t"<<meanz<<"\t"<<dphi<<"\t"<<tkchi2<<"\t"<<tkz0<<"\t"<<tknstubs<<std::endl;

          // BDT INFERENCE
          vector<double> inputs = { hoe, tkpt, srrtot, deta, dpt, meanz, dphi, tkchi2, tkz0, tknstubs } ;
          auto outputs = bdt.decision_function(inputs);
          for (int iout = 0, nout = outputs.size(); iout < nout; ++iout) { 
            std::cout<<outputs[iout]<<std::endl;
          }

        }

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

EGObjEmu &CompositeIDAlgoEmulator::addEGStaToPF(std::vector<EGObjEmu> &egobjs,
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

void CompositeIDAlgoEmulator::addEgObjsToPF(std::vector<EGObjEmu> &egstas,
                                       std::vector<EGIsoObjEmu> &egobjs,
                                       std::vector<EGIsoEleObjEmu> &egeleobjs,
                                       const std::vector<EmCaloObjEmu> &emcalo,
                                       const std::vector<TkObjEmu> &track,
                                       const int calo_idx,
                                       const unsigned int hwQual,
                                       const pt_t ptCorr,
                                       const int tk_idx,
                                       const std::vector<unsigned int> &components) const {
  if (writeEgSta()) {
    addEGStaToPF(egstas, emcalo[calo_idx], hwQual, ptCorr, components);
  }
}
