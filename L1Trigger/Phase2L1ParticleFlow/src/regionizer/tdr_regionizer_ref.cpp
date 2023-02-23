#include "L1Trigger/Phase2L1ParticleFlow/interface/regionizer/tdr_regionizer_ref.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/regionizer/tdr_regionizer_elements_ref.icc"

#include <iostream>

#ifdef CMSSW_GIT_HASH
#include "FWCore/ParameterSet/interface/ParameterSet.h"

l1ct::TDRRegionizerEmulator::TDRRegionizerEmulator(const edm::ParameterSet& iConfig)
    : TDRRegionizerEmulator(
          iConfig.getParameter<uint32_t>("nTrack"),
          iConfig.getParameter<uint32_t>("nCalo"),
          iConfig.getParameter<uint32_t>("nEmCalo"),
          iConfig.getParameter<uint32_t>("nMu"),
          iConfig.getParameter<int32_t>("nClocks"),
          iConfig.getParameter<std::vector<int32_t> >("bigRegionEdges"),
          iConfig.getParameter<bool>("doSort")) {
  debug_ = iConfig.getUntrackedParameter<bool>("debug", false);
}
#endif

l1ct::TDRRegionizerEmulator::TDRRegionizerEmulator(uint32_t ntk,
                                                   uint32_t ncalo,
                                                   uint32_t nem,
                                                   uint32_t nmu,
                                                   int32_t nclocks,
                                                   std::vector<int32_t> bigRegionEdges,
                                                   bool dosort)
    : RegionizerEmulator(),
      ntk_(ntk),
      ncalo_(ncalo),
      nem_(nem),
      nmu_(nmu),
      nclocks_(nclocks),
      bigRegionEdges_(bigRegionEdges),
      dosort_(dosort),
      netaInBR_(6),
      nphiInBR_(3),
      init_(false)
{
  nBigRegions_ = bigRegionEdges_.size() - 1;
  MAX_TK_OBJ_ = nclocks_ * 2 / 3;  // 2 objects per 3 clocks
  MAX_EMCALO_OBJ_ = nclocks_;  // 1 object per clock
  MAX_HADCALO_OBJ_ = nclocks_;  // 1 object per clock
  MAX_MU_OBJ_ = nclocks_;  // 1 object per clock  
}

l1ct::TDRRegionizerEmulator::~TDRRegionizerEmulator() {}

void l1ct::TDRRegionizerEmulator::initSectorsAndRegions(const RegionizerDecodedInputs& in,
                                                        const std::vector<PFInputRegion>& out) {
  if (debug_) {
    dbgCout() << "doing init, out_size = " << out.size() << std::endl;
  }
  assert(!init_);
  init_ = true;
  nregions_ = out.size();


  for (unsigned int i = 0; i < nBigRegions_; i++) {
    tkRegionizers_.emplace_back(
        netaInBR_, nphiInBR_, nregions_, ntk_, bigRegionEdges_[i], bigRegionEdges_[i+1], nclocks_);
    hadCaloRegionizers_.emplace_back(
        netaInBR_, nphiInBR_, nregions_, ncalo_, bigRegionEdges_[i], bigRegionEdges_[i+1], nclocks_);
    emCaloRegionizers_.emplace_back(
        netaInBR_, nphiInBR_, nregions_, nem_, bigRegionEdges_[i], bigRegionEdges_[i+1], nclocks_);
    muRegionizers_.emplace_back(
        netaInBR_, nphiInBR_, nregions_, nmu_, bigRegionEdges_[i], bigRegionEdges_[i+1], nclocks_);
  }

  dbgCout() << "in.track.size() = " << in.track.size() << std::endl;
  dbgCout() << "in.hadcalo.size() = " << in.hadcalo.size() << std::endl;
  dbgCout() << "in.emcalo.size() = " << in.emcalo.size() << std::endl;

  if (ntk_) {
    for (unsigned int i = 0; i < nBigRegions_; i++) {
      tkRegionizers_[i].initSectors(in.track);
      tkRegionizers_[i].initRegions(out);
    }
  }
  if (ncalo_) {
    for (unsigned int i = 0; i < nBigRegions_; i++) {
      hadCaloRegionizers_[i].initSectors(in.hadcalo);
      hadCaloRegionizers_[i].initRegions(out);
    }
  }
  if (nem_) {
    for (unsigned int i = 0; i < nBigRegions_; i++) {
      emCaloRegionizers_[i].initSectors(in.emcalo);
      emCaloRegionizers_[i].initRegions(out);
    }
  }
  if (nmu_) {
    for (unsigned int i = 0; i < nBigRegions_; i++) {
      muRegionizers_[i].initSectors(in.muon);
      muRegionizers_[i].initRegions(out);
    }
  }
}

void l1ct::TDRRegionizerEmulator::fillLinks(const l1ct::RegionizerDecodedInputs& in,
                                            std::vector<std::vector<l1ct::TkObjEmu>>& links) {
  if (ntk_ == 0)
    return;

  links.clear();
  links.resize(in.track.size());

  //one link per sector
  for (unsigned int il = 0; il < in.track.size(); il++) {
    const l1ct::DetectorSector<l1ct::TkObjEmu>& sec = in.track[il];
    for (unsigned int io = 0; io < sec.size(); io++) {
      links[il].push_back(sec[io]);
      if (links[il].size() == MAX_TK_OBJ_) {
        break;
      }
    }
  }
}

void l1ct::TDRRegionizerEmulator::fillLinks(const l1ct::RegionizerDecodedInputs& in,
                                            std::vector<std::vector<l1ct::HadCaloObjEmu>>& links) {
  if (ncalo_ == 0)
    return;

  links.clear();
  links.resize(in.hadcalo.size());

  //one link per sector
  for (unsigned int il = 0; il < in.hadcalo.size(); il++) {
    const l1ct::DetectorSector<l1ct::HadCaloObjEmu>& sec = in.hadcalo[il];
    for (unsigned int io = 0; io < sec.size(); io++) {
      links[il].push_back(sec[io]);
      if (links[il].size() == MAX_HADCALO_OBJ_) {
        break;
      }
    }
  }
}

void l1ct::TDRRegionizerEmulator::fillLinks(const l1ct::RegionizerDecodedInputs& in,
                                            std::vector<std::vector<l1ct::EmCaloObjEmu>>& links) {
  if (nem_ == 0)
    return;

  links.clear();
  links.resize(in.emcalo.size());

  //one link per sector
  for (unsigned int il = 0; il < in.emcalo.size(); il++) {
    const l1ct::DetectorSector<l1ct::EmCaloObjEmu>& sec = in.emcalo[il];
    for (unsigned int io = 0; io < sec.size(); io++) {
      links[il].push_back(sec[io]);
      if (links[il].size() == MAX_EMCALO_OBJ_) {
        break;
      }
    }
  }
}

void l1ct::TDRRegionizerEmulator::fillLinks(const l1ct::RegionizerDecodedInputs& in,
                                            std::vector<std::vector<l1ct::MuObjEmu>>& links) {
  if (nmu_ == 0)
    return;

  links.clear();
  links.resize(1);  //muons are global

  const l1ct::DetectorSector<l1ct::MuObjEmu>& sec = in.muon;
  for (unsigned int io = 0; io < sec.size(); io++) {
    links[0].push_back(sec[io]);
    if (links[0].size() == MAX_MU_OBJ_) {
      break;
    }
  }
}

void l1ct::TDRRegionizerEmulator::toFirmware(const std::vector<l1ct::TkObjEmu>& emu, TkObj fw[NTK_SECTORS][NTK_LINKS]) {
  if (ntk_ == 0)
    return;
  assert(emu.size() == NTK_SECTORS * NTK_LINKS * nBigRegions_);
  for (unsigned int is = 0, idx = 0; is < NTK_SECTORS * nBigRegions_; ++is) {  // tf sectors
    for (unsigned int il = 0; il < NTK_LINKS; ++il, ++idx) {
      fw[is][il] = emu[idx];
    }
  }
}
void l1ct::TDRRegionizerEmulator::toFirmware(const std::vector<l1ct::HadCaloObjEmu>& emu,
                                             HadCaloObj fw[NCALO_SECTORS][NCALO_LINKS]) {
  if (ncalo_ == 0)
    return;
  assert(emu.size() == NCALO_SECTORS * NCALO_LINKS * nBigRegions_);
  for (unsigned int is = 0, idx = 0; is < NCALO_SECTORS * nBigRegions_; ++is) {  // tf sectors
    for (unsigned int il = 0; il < NCALO_LINKS; ++il, ++idx) {
      fw[is][il] = emu[idx];
    }
  }
}

void l1ct::TDRRegionizerEmulator::toFirmware(const std::vector<l1ct::EmCaloObjEmu>& emu,
                                             EmCaloObj fw[NCALO_SECTORS][NCALO_LINKS]) {
  if (nem_ == 0)
    return;
  assert(emu.size() == NCALO_SECTORS * NCALO_LINKS * nBigRegions_);
  for (unsigned int is = 0, idx = 0; is < NCALO_SECTORS * nBigRegions_; ++is) {  // tf sectors
    for (unsigned int il = 0; il < NCALO_LINKS; ++il, ++idx) {
      fw[is][il] = emu[idx];
    }
  }
}

void l1ct::TDRRegionizerEmulator::toFirmware(const std::vector<l1ct::MuObjEmu>& emu, MuObj fw[NMU_LINKS]) {
  if (nmu_ == 0)
    return;
  assert(emu.size() == NMU_LINKS);
  for (unsigned int il = 0, idx = 0; il < NMU_LINKS; ++il, ++idx) {
    fw[il] = emu[idx];
  }
}

void l1ct::TDRRegionizerEmulator::run(const RegionizerDecodedInputs& in, std::vector<PFInputRegion>& out) {
  if (debug_) {
    dbgCout() << "TDRRegionizerEmulator::run called, out.size =  " << out.size() << std::endl;
  }

  if (!init_) {
    initSectorsAndRegions(in, out);
  }

  std::vector<std::vector<l1ct::TkObjEmu>> tk_links_in;
  std::vector<std::vector<l1ct::EmCaloObjEmu>> em_links_in;
  std::vector<std::vector<l1ct::HadCaloObjEmu>> calo_links_in;
  std::vector<std::vector<l1ct::MuObjEmu>> mu_links_in;

  // read the inputs
  fillLinks(in, tk_links_in);
  fillLinks(in, em_links_in);
  fillLinks(in, calo_links_in);
  fillLinks(in, mu_links_in);
  //this is overkill and could be improved, for now its ok (the sectors outside each board just wont do anything)

  for (unsigned int ie = 0; ie < nBigRegions_; ie++) {
    //add objects from link
    tkRegionizers_[ie].reset();
    tkRegionizers_[ie].setPipes(tk_links_in);
    tkRegionizers_[ie].initTimes();
    if (debug_) {
      dbgCout() << ie << "SECTORS/LINKS " << ie << std::endl;
      for (unsigned int i = 0; i < tk_links_in.size(); i++) {
        for (unsigned int j = 0; j < tk_links_in[i].size(); j++) {
          dbgCout() << "\t" << i << " " << j << "\t" << tk_links_in[i][j].hwPt.to_int() << "\t"
                    << tk_links_in[i][j].hwEta.to_int() << "\t" << tk_links_in[i][j].hwPhi.to_int() << std::endl;
        }
        dbgCout() << "-------------------------------" << std::endl;
      }
    }
    tkRegionizers_[ie].run(debug_);

    emCaloRegionizers_[ie].reset();
    emCaloRegionizers_[ie].setPipes(em_links_in);
    emCaloRegionizers_[ie].initTimes();
    emCaloRegionizers_[ie].run();

    hadCaloRegionizers_[ie].reset();
    hadCaloRegionizers_[ie].setPipes(calo_links_in);
    hadCaloRegionizers_[ie].initTimes();
    hadCaloRegionizers_[ie].run();

    muRegionizers_[ie].reset();
    muRegionizers_[ie].setPipes(mu_links_in);
    muRegionizers_[ie].initTimes();
    muRegionizers_[ie].run();
  }

  for (unsigned int ie = 0; ie < nBigRegions_; ie++) {
    for (unsigned int ireg = 0; ireg < nregions_; ireg++) {
      std::vector<l1ct::TkObjEmu> out_tks = tkRegionizers_[ie].getSmallRegion(ireg);
      if (!out_tks.empty()) {
        if (dosort_) {
          std::sort(
              out_tks.begin(), out_tks.end(), [](const l1ct::TkObjEmu a, const l1ct::TkObjEmu b) { return a > b; });
        }
        out[ireg].track = out_tks;
      }
      std::vector<l1ct::EmCaloObjEmu> out_emcalos = emCaloRegionizers_[ie].getSmallRegion(ireg);
      if (!out_emcalos.empty()) {
        if (dosort_) {
          std::sort(out_emcalos.begin(), out_emcalos.end(), [](const l1ct::EmCaloObjEmu a, const l1ct::EmCaloObjEmu b) {
            return a > b;
          });
        }
        out[ireg].emcalo = out_emcalos;
      }
      std::vector<l1ct::HadCaloObjEmu> out_hadcalos = hadCaloRegionizers_[ie].getSmallRegion(ireg);
      if (!out_hadcalos.empty()) {
        if (dosort_) {
          std::sort(out_hadcalos.begin(),
                    out_hadcalos.end(),
                    [](const l1ct::HadCaloObjEmu a, const l1ct::HadCaloObjEmu b) { return a > b; });
        }
        out[ireg].hadcalo = out_hadcalos;
      }
      std::vector<l1ct::MuObjEmu> out_mus = muRegionizers_[ie].getSmallRegion(ireg);
      if (!out_mus.empty()) {
        if (dosort_) {
          std::sort(
              out_mus.begin(), out_mus.end(), [](const l1ct::MuObjEmu a, const l1ct::MuObjEmu b) { return a > b; });
        }
        out[ireg].muon = out_mus;
      }
    }
  }
}
