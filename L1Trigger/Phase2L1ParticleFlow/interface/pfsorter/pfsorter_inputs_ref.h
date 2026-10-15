#ifndef pfsorter_inputs_ref_h
#define pfsorter_inputs_ref_h

#include "DataFormats/L1TParticleFlow/interface/layer1_emulator.h"

#include <vector>

// Shape of the input links of the sorter. The links are split in two groups, the first
// NLINKCHARGED carrying the charged PF candidates and the following NLINKNEUTRAL the
// neutral ones (photons and neutral hadrons); each region occupies a fixed slot of
// NOBJPERLINK clock cycles on every link, padded with empty objects.
//
// The links are the ones of the tmux18 demonstrator firmware: the PF algorithm outputs
// NSLOTCHARGED charged and NSLOTNEUTRAL neutral slots per region (NTKSORTED and
// NPFNEUTRAL there, empty slots included), which are serialised at 240 MHz over
// NWRITE240 cycles per region and then re-packed at 360 MHz over NOBJPERLINK cycles; see
// PFSorterEmulator::firmwareSlotLayout for the exact mapping.
#ifndef NLINKCHARGED
#define NLINKCHARGED 4
#endif
#ifndef NLINKNEUTRAL
#define NLINKNEUTRAL 3
#endif
#ifndef NOBJPERLINK
#define NOBJPERLINK 9
#endif
#ifndef NSLOTCHARGED
#define NSLOTCHARGED 30
#endif
#ifndef NSLOTNEUTRAL
#define NSLOTNEUTRAL 20
#endif
#ifndef NWRITE240
#define NWRITE240 6
#endif

namespace l1ct {
  namespace pfsorter {

    constexpr unsigned int nLinksCharged = NLINKCHARGED;
    constexpr unsigned int nLinksNeutral = NLINKNEUTRAL;
    constexpr unsigned int nLinksTotal = NLINKCHARGED + NLINKNEUTRAL;
    constexpr unsigned int nObjPerLink = NOBJPERLINK;
    constexpr unsigned int nSlotsCharged = NSLOTCHARGED;
    constexpr unsigned int nSlotsNeutral = NSLOTNEUTRAL;
    constexpr unsigned int nWrite240 = NWRITE240;

    // One particle coming out of the ParticleFlow algorithm, as seen on one input link of
    // the sorter. PF particles come in two flavours, charged (l1ct::PFChargedObjEmu, also
    // used for muons) and neutral (l1ct::PFNeutralObjEmu, also used for photons), and their
    // coordinates are local to the region they were reconstructed in; the sorter tree
    // instead stores l1ct::PuppiObjEmu, in global coordinates, so every particle has to be
    // converted with toPuppi() and the region it comes from before entering the tree.
    struct PFParticleEmu {
      enum Kind { None = 0, Charged = 1, Neutral = 2 };

      Kind kind;
      l1ct::PFChargedObjEmu charged;
      l1ct::PFNeutralObjEmu neutral;
      l1ct::puppiWgt_t puppiWgt;  // weight of a neutral particle (1 = no puppi correction)
      l1ct::pt_t puppiPt;         // pt of a neutral particle: the same as the PF one, the
                                  // weight is only recorded, not applied
      l1ct::eta_t hwEta;
      l1ct::phi_t hwPhi;

      PFParticleEmu() { clear(); }
      explicit PFParticleEmu(const l1ct::PFChargedObjEmu &src) { set(src); }
      explicit PFParticleEmu(const l1ct::PFNeutralObjEmu &src, l1ct::puppiWgt_t wgt = 1.0) { set(src, wgt); }

      void clear() {
        kind = None;
        charged.clear();
        neutral.clear();
        puppiWgt = 0;
        puppiPt = 0;
        hwEta = 0;
        hwPhi = 0;
      }

      void set(const l1ct::PFChargedObjEmu &src) {
        clear();
        if (src.hwPt == 0)
          return;
        kind = Charged;
        charged = src;
        hwEta = charged.hwEta;
        hwPhi = charged.hwPhi;
      }

      void set(const l1ct::PFNeutralObjEmu &src, l1ct::puppiWgt_t wgt = 1.0) {
        clear();
        if (src.hwPt == 0)
          return;
        kind = Neutral;
        neutral = src;
        puppiWgt = wgt;
        // the pt is taken over unchanged from the PF candidate: eta and phi are the only
        // values the conversion modifies, and the weight is only recorded in the payload
        puppiPt = src.hwPt;
        hwEta = neutral.hwEta;
        hwPhi = neutral.hwPhi;
      }

      // pt the particle will have once converted (0 for an empty link)
      l1ct::pt_t hwPt() const {
        if (kind == Charged)
          return charged.hwPt;
        else if (kind == Neutral)
          return puppiPt;
        return l1ct::pt_t(0);
      }
      bool valid() const { return kind != None && hwPt() != 0; }
      int intPt() const { return l1ct::Scales::intPt(hwPt()); }
    };

    // convert one PF particle into a PuppiObj with global coordinates. An invalid particle
    // (an empty link slot, or one with no pt) gives an empty object, and so does one that
    // is not fiducial in its region: regions are padded by etaExtra/phiExtra, so a particle
    // sitting in the pad is reconstructed a second time by the neighbouring region that
    // holds it as fiducial, and keeping both copies would duplicate it in the sorter
    // output. This is the same cut that fetchPF() and linpuppi_ref() apply to these very
    // same PF collections. The cut lives here, and not where the links are built, because
    // it is part of what the conversion does in firmware: it only needs the local
    // coordinates of the particle and the half-widths of the region, both of which
    // pf_to_puppi already has in hand.
    inline void toPuppi(const l1ct::PFRegionEmu &region, const PFParticleEmu &in, l1ct::PuppiObjEmu &out) {
      if (!in.valid() || !region.isFiducial(in)) {
        out.clear();
      } else if (in.kind == PFParticleEmu::Charged) {
        out.fill(region, in.charged);
      } else {
        out.fill(region, in.neutral, in.puppiPt, in.puppiWgt);
      }
    }
    inline l1ct::PuppiObjEmu toPuppi(const l1ct::PFRegionEmu &region, const PFParticleEmu &in) {
      l1ct::PuppiObjEmu ret;
      toPuppi(region, in, ret);
      return ret;
    }

    // convert a whole PF region: charged particles, photons and neutral hadrons, in this
    // order, skipping the empty slots. Returns the number of particles appended.
    //
    // Note that pf.pfmuon is deliberately NOT read: both PFAlgo2HGC and PFAlgo3 write a
    // muon-matched track into pfcharged (with hwId.isMuon() set) as well as into pfmuon, so
    // reading both would emit every muon twice. This matches what the rest of the chain
    // does -- fetchPF() and linpuppi_ref() also take muons from pfcharged only.
    inline unsigned int toPFParticlesCharged(const l1ct::OutputRegion &pf, std::vector<PFParticleEmu> &out) {
      unsigned int n = 0;
      for (const auto &c : pf.pfcharged) {
        if (c.hwPt != 0) {
          out.emplace_back(c);
          n++;
        }
      }
      return n;
    }

    // the neutral particles of a region: photons first, then neutral hadrons
    inline unsigned int toPFParticlesNeutral(const l1ct::OutputRegion &pf,
                                             std::vector<PFParticleEmu> &out,
                                             l1ct::puppiWgt_t wgt = 1.0) {
      unsigned int n = 0;
      for (const auto &p : pf.pfphoton) {
        if (p.hwPt != 0) {
          out.emplace_back(p, wgt);
          n++;
        }
      }
      for (const auto &h : pf.pfneutral) {
        if (h.hwPt != 0) {
          out.emplace_back(h, wgt);
          n++;
        }
      }
      return n;
    }

    inline unsigned int toPFParticles(const l1ct::OutputRegion &pf,
                                      std::vector<PFParticleEmu> &out,
                                      l1ct::puppiWgt_t wgt = 1.0) {
      unsigned int n = toPFParticlesCharged(pf, out);
      return n + toPFParticlesNeutral(pf, out, wgt);
    }

    // same, converting the particles into PuppiObj at the same time
    inline unsigned int toPuppi(const l1ct::PFRegionEmu &region,
                                const l1ct::OutputRegion &pf,
                                std::vector<l1ct::PuppiObjEmu> &out,
                                l1ct::puppiWgt_t wgt = 1.0) {
      std::vector<PFParticleEmu> particles;
      unsigned int n = toPFParticles(pf, particles, wgt);
      for (const auto &p : particles)
        out.push_back(toPuppi(region, p));
      return n;
    }

  }  // namespace pfsorter
}  // namespace l1ct

#endif
