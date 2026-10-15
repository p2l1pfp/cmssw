#ifndef pfsorter_ref_h
#define pfsorter_ref_h

#include "DataFormats/L1TParticleFlow/interface/layer1_emulator.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/pfsorter/pfsorter_tree_elements_ref.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/pfsorter/pfsorter_inputs_ref.h"

#include <algorithm>
#include <vector>

#ifdef CMSSW_GIT_HASH
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#endif

namespace l1ct {

  // Clock-cycle emulation of the pfsorter merge tree, built on the same fifo / staging
  // area / queue structure as the multififo regionizer (see pfsorter_tree_elements_ref.h).
  //
  //   link 0 --> fifo 0 --.
  //                        >-- o --.
  //   link 1 --> fifo 1 --'         \
  //   link 2 --> fifo 2 --.          >-- o --> one object per clock cycle
  //                        >-- o --'
  //   link 3 --> fifo 3 --'
  //
  // One buffer per link (1:1:1 link / buffer / fifo of the first layer), any number of
  // links, at most 3 output nodes, one object per node.
  //
  // The link fifos are double buffered: there are two sets of fifos, with a read pointer
  // on the set the links are written into and a write pointer on the set the merge tree
  // is drained from. On every new event the two pointers swap sets, so an event is read
  // in while the previous one is still being written out. In the event-by-event mode
  // (continuous = false, the CMSSW default) both pointers are kept on the same set, so
  // that each event is self contained.
  //
  // Usage:
  //    l1ct::PFSorterEmulator sorter(/*nlinks=*/12, /*noutputs=*/3);
  //    l1ct::PuppiObjEmu out;
  //    for (int iclock = 0; iclock < nclocks; ++iclock) {
  //      bool newEvent = (iclock % nclocksPerEvent == 0);
  //      sorter.step(newEvent, links[iclock], out);  // out is empty if nothing came out
  //    }
  template <typename T>
  class PFSorterEmulatorT {
  public:
    PFSorterEmulatorT() : nevt_(0) {}
    PFSorterEmulatorT(unsigned int nlinks, unsigned int noutputs = 3, unsigned int nclocks = 162) : nevt_(0), nclocks_(nclocks) {
      buffer_.initFifos(nlinks, noutputs);
    }

    void initFifos(unsigned int nlinks, unsigned int noutputs = 3) {
      buffer_.initFifos(nlinks, noutputs);
    }

    void reset() {
      buffer_.reset();
      nevt_ = 0;
    }

    // single clock emulation; on a new event the tree rolls over as the firmware does.
    // 'inputs' must have one (possibly empty, i.e. hwPt == 0) object per link, and they go
    // into the fifo with the same index. 'out' is the object leaving the tree in this clock
    // cycle; the return value says whether it is a valid one.
    bool step(bool newEvent, const std::vector<T>& inputs, T& out) {
      if (newEvent) {
        buffer_.roll();
        nevt_++;
      }
      assert(inputs.size() == buffer_.nfifos());
      for (unsigned int i = 0, n = inputs.size(); i < n; ++i) {
        buffer_.maybe_push(i, inputs[i]);
      }
      out = buffer_.pop();
      return out.hwPt != 0;
    }

    // keep clocking with no new inputs until the tree is empty, appending to 'out'
    // everything that comes out. Returns the number of clock cycles used.
    unsigned int drain(std::vector<T>& out, unsigned int maxclocks = 10000) {
      unsigned int nclocks = 0;
      while (buffer_.inFlight() != 0 && nclocks < maxclocks) {
        T obj = buffer_.pop();
        nclocks++;
        if (obj.hwPt != 0)
          out.push_back(obj);
      }
      return nclocks;
    }

    const l1ct::pfsorter::RegionBuffer<T>& buffer() const { return buffer_; }
    l1ct::pfsorter::RegionBuffer<T>& buffer() { return buffer_; }
    unsigned int nEvents() const { return nevt_; }
    unsigned int nClocks() const { return nclocks_; }
    // continuous = true runs events back to back as the firmware does, the read and write
    // pointers swapping to opposite fifo sets on every new event; false is the
    // event-by-event CMSSW mode, where both pointers stay on the same set
    bool continuous() const { return continuous_; }
    void setContinuous(bool continuous) { continuous_ = continuous; }

  private:
    l1ct::pfsorter::RegionBuffer<T> buffer_;
    unsigned int nevt_;
    unsigned int nclocks_ = 162;  // also covers the default constructor
    bool continuous_ = false;
  };

  // The sorter as it is used in practice: it stores PuppiObj, and it can be fed either
  // directly with PuppiObj or with the particles coming out of the ParticleFlow algorithm,
  // which are converted into PuppiObj on the way in (see pfsorter_inputs_ref.h).
  class PFSorterEmulator : public PFSorterEmulatorT<l1ct::PuppiObjEmu> {
  public:
    typedef l1ct::pfsorter::PFParticleEmu PFParticle;
    using PFSorterEmulatorT<l1ct::PuppiObjEmu>::PFSorterEmulatorT;
    using PFSorterEmulatorT<l1ct::PuppiObjEmu>::step;

    // The links are split in two groups: the first 'nlinksCharged' carry the charged PF
    // candidates and the following 'nlinksNeutral' the neutral ones. Each region occupies
    // a fixed slot of 'nobjperlink' clock cycles on every link, into which the
    // 'nslotsCharged' charged and 'nslotsNeutral' neutral slots the PF algorithm outputs
    // for the region are laid out as the firmware does (see makePFLinks).
    PFSorterEmulator(unsigned int nlinksCharged,
                     unsigned int nlinksNeutral,
                     unsigned int noutputs,
                     unsigned int nclocks,
                     unsigned int nobjperlink,
                     unsigned int nslotsCharged = l1ct::pfsorter::nSlotsCharged,
                     unsigned int nslotsNeutral = l1ct::pfsorter::nSlotsNeutral)
        : PFSorterEmulatorT<l1ct::PuppiObjEmu>(nlinksCharged + nlinksNeutral, noutputs, nclocks),
          nLinksCharged_(nlinksCharged),
          nLinksNeutral_(nlinksNeutral),
          nObjPerLink_(nobjperlink),
          nSlotsCharged_(nslotsCharged),
          nSlotsNeutral_(nslotsNeutral) {}

#ifdef CMSSW_GIT_HASH
    PFSorterEmulator(const edm::ParameterSet& iConfig)
        : PFSorterEmulator(iConfig.getParameter<uint32_t>("nLinksCharged"),
                           iConfig.getParameter<uint32_t>("nLinksNeutral"),
                           iConfig.getParameter<uint32_t>("nOutputs"),
                           iConfig.getParameter<uint32_t>("nClocks"),
                           iConfig.getParameter<uint32_t>("nObjPerLink"),
                           iConfig.getParameter<uint32_t>("nSlotsCharged"),
                           iConfig.getParameter<uint32_t>("nSlotsNeutral")) {}

    static edm::ParameterSetDescription getParameterSetDescription() {
      edm::ParameterSetDescription description;
      description.add<uint32_t>("nLinksCharged", l1ct::pfsorter::nLinksCharged);
      description.add<uint32_t>("nLinksNeutral", l1ct::pfsorter::nLinksNeutral);
      description.add<uint32_t>("nObjPerLink", l1ct::pfsorter::nObjPerLink);
      description.add<uint32_t>("nSlotsCharged", l1ct::pfsorter::nSlotsCharged);
      description.add<uint32_t>("nSlotsNeutral", l1ct::pfsorter::nSlotsNeutral);
      description.add<uint32_t>("nOutputs", 3);
      description.add<uint32_t>("nClocks", 0);  // 0 = nObjPerLink * (number of regions)
      return description;
    }
#endif

    unsigned int nLinksCharged() const { return nLinksCharged_; }
    unsigned int nLinksNeutral() const { return nLinksNeutral_; }
    unsigned int nLinks() const { return nLinksCharged_ + nLinksNeutral_; }
    unsigned int nObjPerLink() const { return nObjPerLink_; }
    unsigned int nSlotsCharged() const { return nSlotsCharged_; }
    unsigned int nSlotsNeutral() const { return nSlotsNeutral_; }

    // Where the firmware puts the PF slots of one region on the links of one group
    // (charged or neutral): layout[link][cycle] is the index of the PF slot that link
    // carries in that cycle of the region, or -1 for an empty cycle. It models the two
    // stages of the tmux18 demonstrator between the PF algorithm and the sorter:
    //
    //  1. parallel2serial at 240 MHz, NWRITE = nwrite240 cycles per region: stream s
    //     carries slots s*nwrite240 ... s*nwrite240+nwrite240-1, in this order;
    //  2. stream_compress_9to6 at 360 MHz, after cdc_and_deserializer has replayed those
    //     nwrite240 (6) cycles at the start of each 9-cycle region: every group of three
    //     streams 3g, 3g+1, 3g+2 becomes two links 2g, 2g+1, the first two streams going
    //     straight through in cycles 0-5 and the third one being queued and sent two per
    //     cycle in cycles 6-8 (slots 0, 2, 4 on link 2g+1 and 1, 3, 5 on link 2g).
    //
    // With 30 charged slots that is link 0 = slots 0-5,13,15,17; link 1 = 6-11,12,14,16;
    // link 2 = 18-23; link 3 = 24-29. The number of links it returns, ceil(nslots/9), must
    // be the number of links of the group.
    static std::vector<std::vector<int>> firmwareSlotLayout(unsigned int nslots,
                                                            unsigned int nobjperlink = l1ct::pfsorter::nObjPerLink,
                                                            unsigned int nwrite240 = l1ct::pfsorter::nWrite240) {
      // stream_compress_9to6 is hard wired for 6 cycles in and 9 out
      assert(nwrite240 == 6 && nobjperlink == 9);
      const unsigned int nstream240 = (nslots + nwrite240 - 1) / nwrite240;
      auto in = [&](unsigned int stream, unsigned int cycle) -> int {
        unsigned int slot = stream * nwrite240 + cycle;
        return (stream < nstream240 && slot < nslots) ? int(slot) : -1;
      };
      const unsigned int ngroups = (nslots + 17) / 18;
      const unsigned int nlinks = (nslots + nobjperlink - 1) / nobjperlink;
      std::vector<std::vector<int>> out(2 * ngroups, std::vector<int>(nobjperlink, -1));
      for (unsigned int g = 0; g < ngroups; ++g) {
        for (unsigned int c = 0; c < nwrite240; ++c) {
          out[2 * g][c] = in(3 * g, c);
          out[2 * g + 1][c] = in(3 * g + 1, c);
        }
        for (unsigned int k = 0; k < (nobjperlink - nwrite240); ++k) {
          out[2 * g + 1][nwrite240 + k] = in(3 * g + 2, 2 * k);
          out[2 * g][nwrite240 + k] = in(3 * g + 2, 2 * k + 1);
        }
      }
      out.resize(nlinks);  // the compressor only drives ceil(nslots/9) links
      return out;
    }

    // flush the tree and the conversion register
    void reset() {
      PFSorterEmulatorT<l1ct::PuppiObjEmu>::reset();
      clearConverted_();
    }

    // single clock cycle taking one PF particle per link, each with the region it was
    // reconstructed in (PF coordinates are local to the region): the particles are
    // converted into PuppiObj and then pushed into the fifos.
    //
    // The conversion costs one clock cycle, as it does in firmware (pf_to_puppi registers
    // its output), so what is pushed into the fifos in this clock cycle is what was
    // converted in the previous one, and newEvent is delayed by the same cycle.
    bool step(bool newEvent,
              const std::vector<l1ct::PFRegionEmu>& regions,
              const std::vector<PFParticle>& inputs,
              l1ct::PuppiObjEmu& out) {
      assert(regions.size() == inputs.size());
      std::vector<l1ct::PuppiObjEmu> converted(inputs.size());
      for (unsigned int i = 0, n = inputs.size(); i < n; ++i) {
        l1ct::pfsorter::toPuppi(regions[i], inputs[i], converted[i]);
      }
      if (converted_.size() != inputs.size()) {
        converted_.resize(inputs.size());
        clearConverted_();
      }
      bool ret = step(convertedNewEvent_, converted_, out);
      converted_.swap(converted);
      convertedNewEvent_ = newEvent;
      return ret;
    }

    // same, when all the links carry particles of the same region
    bool step(bool newEvent,
              const l1ct::PFRegionEmu& region,
              const std::vector<PFParticle>& inputs,
              l1ct::PuppiObjEmu& out) {
      return step(newEvent, std::vector<l1ct::PFRegionEmu>(inputs.size(), region), inputs, out);
    }

    // keep clocking with no new inputs until the tree is empty. The first of those clock
    // cycles is the one that lets the last converted particles into the fifos.
    unsigned int drain(std::vector<l1ct::PuppiObjEmu>& out, unsigned int maxclocks = 10000) {
      unsigned int nclocks = 0;
      if (!converted_.empty()) {
        l1ct::PuppiObjEmu obj;
        if (PFSorterEmulatorT<l1ct::PuppiObjEmu>::step(convertedNewEvent_, converted_, obj))
          out.push_back(obj);
        clearConverted_();
        nclocks++;
      }
      return nclocks + PFSorterEmulatorT<l1ct::PuppiObjEmu>::drain(out, maxclocks);
    }

    struct PFLinkItem {
      l1ct::pfsorter::PFParticleEmu particle;
      l1ct::PFRegionEmu region;
    };

    // Collect the puppi candidates of an event and distribute them over 'nlinks' input
    // links, one link per region, round robin if there are more regions than links.
    // Returns the number of objects. If 'regionIndex' is given, it is filled in parallel
    // to 'links' with the index of the region each object came from.
    static unsigned int makeLinks(const l1ct::Event &event,
                                  unsigned int nlinks,
                                  std::vector<std::vector<l1ct::PuppiObjEmu>> &links,
                                  std::vector<std::vector<unsigned int>> *regionIndex = nullptr) {
      links.clear();
      links.resize(nlinks);
      if (regionIndex) {
        regionIndex->clear();
        regionIndex->resize(nlinks);
      }
      unsigned int nobj = 0;
      for (unsigned int ireg = 0, nreg = event.out.size(); ireg < nreg; ++ireg) {
        for (const auto &p : event.out[ireg].puppi) {
          if (p.hwPt == 0)
            continue;
          links[ireg % nlinks].push_back(p);
          if (regionIndex)
            (*regionIndex)[ireg % nlinks].push_back(ireg);
          nobj++;
        }
      }
      return nobj;
    }

    // Distribute the particles produced by the PF algorithm over the input links, laid out
    // exactly as the tmux18 demonstrator firmware lays them out. They are handed to the
    // sorter as PF particles and converted into PuppiObj on the way in.
    //
    // The links are split in two groups, 'nlinksCharged' for the charged candidates
    // followed by 'nlinksNeutral' for the neutral ones, and the event is laid out region
    // by region: each region gets a slot of 'nobjperlink' consecutive clock cycles on every
    // link. Within it, the PF slots of the region keep their position: charged slot i is
    // pfcharged[i] and neutral slot i is pfneutral[i], empty slots included, and they go
    // where firmwareSlotLayout() puts them. An empty slot or cycle is an empty object.
    //
    // A region with more than 'nslotsCharged' / 'nslotsNeutral' slots cannot be laid out:
    // the extra slots are dropped and, if 'ndropped' is given, the non-empty ones among
    // them are counted there (the PF algorithm is configured with exactly as many slots as
    // the firmware carries, so this is expected to stay 0).
    //
    // Returns the number of non-empty particles written into the links.
    static unsigned int makePFLinks(const std::vector<l1ct::PFInputRegion> &pfinputs,
                                    const std::vector<l1ct::OutputRegion> &pfouts,
                                    unsigned int nlinksCharged,
                                    unsigned int nlinksNeutral,
                                    unsigned int nobjperlink,
                                    unsigned int nslotsCharged,
                                    unsigned int nslotsNeutral,
                                    std::vector<std::vector<PFLinkItem>> &links,
                                    l1ct::puppiWgt_t wgt = 1.0,
                                    std::vector<std::vector<unsigned int>> *regionIndex = nullptr,
                                    unsigned int *ndropped = nullptr) {
      const std::vector<std::vector<int>> layoutCharged = firmwareSlotLayout(nslotsCharged, nobjperlink);
      const std::vector<std::vector<int>> layoutNeutral = firmwareSlotLayout(nslotsNeutral, nobjperlink);
      // the link counts must be the ones the firmware drives
      assert(layoutCharged.size() == nlinksCharged);
      assert(layoutNeutral.size() == nlinksNeutral);
      const unsigned int nlinks = nlinksCharged + nlinksNeutral;
      const unsigned int nreg = pfouts.size();
      links.clear();
      links.resize(nlinks);
      if (regionIndex) {
        regionIndex->clear();
        regionIndex->resize(nlinks);
      }
      // every output region must come with its input region, or the conversion has no
      // centre to move the particles to global coordinates
      assert(pfouts.size() == pfinputs.size());
      unsigned int nobj = 0, ndrop = 0;
      for (unsigned int ireg = 0; ireg < nreg; ++ireg) {
        const l1ct::PFRegionEmu &region = pfinputs[ireg].region;
        const l1ct::OutputRegion &pf = pfouts[ireg];
        const unsigned int base = ireg * nobjperlink;
        for (unsigned int l = 0; l < nlinks; ++l) {
          links[l].resize(base + nobjperlink);  // pad with empty objects
          if (regionIndex)
            (*regionIndex)[l].resize(base + nobjperlink, ireg);
        }
        for (unsigned int l = 0; l < nlinksCharged; ++l) {
          for (unsigned int c = 0; c < nobjperlink; ++c) {
            int slot = layoutCharged[l][c];
            if (slot < 0 || unsigned(slot) >= pf.pfcharged.size())
              continue;
            PFLinkItem item{l1ct::pfsorter::PFParticleEmu(pf.pfcharged[slot]), region};
            if (item.particle.valid())
              nobj++;
            links[l][base + c] = item;
          }
        }
        for (unsigned int l = 0; l < nlinksNeutral; ++l) {
          for (unsigned int c = 0; c < nobjperlink; ++c) {
            int slot = layoutNeutral[l][c];
            if (slot < 0 || unsigned(slot) >= pf.pfneutral.size())
              continue;
            PFLinkItem item{l1ct::pfsorter::PFParticleEmu(pf.pfneutral[slot], wgt), region};
            if (item.particle.valid())
              nobj++;
            links[nlinksCharged + l][base + c] = item;
          }
        }
        // slots beyond what the links carry
        for (unsigned int i = nslotsCharged; i < pf.pfcharged.size(); ++i)
          ndrop += (pf.pfcharged[i].hwPt != 0);
        for (unsigned int i = nslotsNeutral; i < pf.pfneutral.size(); ++i)
          ndrop += (pf.pfneutral[i].hwPt != 0);
      }
      if (ndropped)
        *ndropped = ndrop;
      return nobj;
    }

    // same, from a full event (the two vectors above are event.pfinputs and event.out)
    static unsigned int makePFLinks(const l1ct::Event &event,
                                    unsigned int nlinksCharged,
                                    unsigned int nlinksNeutral,
                                    unsigned int nobjperlink,
                                    unsigned int nslotsCharged,
                                    unsigned int nslotsNeutral,
                                    std::vector<std::vector<PFLinkItem>> &links,
                                    l1ct::puppiWgt_t wgt = 1.0,
                                    std::vector<std::vector<unsigned int>> *regionIndex = nullptr,
                                    unsigned int *ndropped = nullptr) {
      return makePFLinks(event.pfinputs,
                         event.out,
                         nlinksCharged,
                         nlinksNeutral,
                         nobjperlink,
                         nslotsCharged,
                         nslotsNeutral,
                         links,
                         wgt,
                         regionIndex,
                         ndropped);
    }

    // same, using the link shape this emulator was configured with
    unsigned int makePFLinks(const std::vector<l1ct::PFInputRegion> &pfinputs,
                             const std::vector<l1ct::OutputRegion> &pfouts,
                             std::vector<std::vector<PFLinkItem>> &links,
                             l1ct::puppiWgt_t wgt = 1.0,
                             std::vector<std::vector<unsigned int>> *regionIndex = nullptr,
                             unsigned int *ndropped = nullptr) const {
      return makePFLinks(pfinputs,
                         pfouts,
                         nLinksCharged_,
                         nLinksNeutral_,
                         nObjPerLink_,
                         nSlotsCharged_,
                         nSlotsNeutral_,
                         links,
                         wgt,
                         regionIndex,
                         ndropped);
    }

    // Run the sorter over one event's worth of PF output, one link per region, and collect
    // the sorted puppi candidates. This is the whole of the algorithm; callers that do not
    // have an l1ct::Event at hand (a standalone testbench working region by region) can
    // pass the two vectors directly.
    //
    // The tree is sized on the number of regions if it has not been sized already, or if it
    // was sized for a different number of links, so that a default constructed emulator can
    // be used as is.
    void run(const std::vector<l1ct::PFInputRegion> &pfinputs,
             const std::vector<l1ct::OutputRegion> &pfouts,
             std::vector<l1ct::PuppiObjEmu> &out) {
      std::vector<l1ct::PuppiObjEmu> stream;
      runStream(pfinputs, pfouts, stream);
      for (const auto &o : stream) {
        if (o.hwPt != 0)
          out.push_back(o);
      }
    }

    // Same, but keeping the output link as it is in firmware: one entry per clock cycle,
    // empty (hwPt == 0) in the cycles where nothing leaves the tree; run() above is this
    // one with the empty cycles dropped.
    //
    // By default each event is self contained: the tree starts empty and is drained at the
    // end, so 'stream' has the event's clock cycles plus however many it took to empty the
    // tree. That is what CMSSW needs, where nothing may outlive the event it came from.
    // The read and write pointers are kept on the same set of fifos, so the objects of
    // the event are read out of the very fifos they were written into.
    //
    // With 'continuous' the sorter runs as the firmware does, events back to back: no
    // reset and no drain, exactly nClocks() entries per event with newEvent on the first
    // one, the read and write pointers swapping to opposite fifo sets on every new event,
    // and whatever is still in the tree at the end of an event handled by the next one.
    // Use it to write pattern files that line up frame by frame with the firmware.
    void runStream(const std::vector<l1ct::PFInputRegion> &pfinputs,
                   const std::vector<l1ct::OutputRegion> &pfouts,
                   std::vector<l1ct::PuppiObjEmu> &stream,
                   bool continuous = false) {
      const unsigned int nlinks = nLinks();
      if (buffer().nfifos() != nlinks)
        initFifos(nlinks, buffer().noutputs() ? buffer().noutputs() : 3);
      this->setContinuous(continuous);
      if (continuous) {
        runContinuous_(pfinputs, pfouts, stream);
        return;
      }
      // Start from an empty tree AND an empty conversion register. The newEvent flag alone
      // is not enough: it only reaches the tree one clock cycle later, because it travels
      // with the particles through the conversion, so without this reset the objects left
      // in the conversion register by the previous event would be pushed into the fifos in
      // the first clock cycle of this one and could be read out before the flush.
      reset();
      // event-by-event mode: the read and write pointers point to the same set of fifos,
      // so that the event is read out of the same fifos it is being written into
      buffer().setWriteToRead();
      std::vector<std::vector<PFLinkItem>> links;
      makePFLinks(pfinputs, pfouts, links);
      // the links are nObjPerLink clock cycles per region long; nClocks() is only a floor,
      // so that a pattern file can be padded to a fixed frame count
      const unsigned int nclk = std::max<unsigned int>(nClocks(), links.empty() ? 0 : links[0].size());
      stream.clear();
      stream.reserve(nclk);
      l1ct::PuppiObjEmu outObj;
      for (unsigned int iclock = 0; iclock < nclk; ++iclock) {
        bool newevt = (iclock == 0);
        std::vector<l1ct::pfsorter::PFParticleEmu> pfin(nlinks);
        std::vector<l1ct::PFRegionEmu> pfregions(nlinks);
        for (unsigned int l = 0; l < nlinks; ++l) {
          if (iclock < links[l].size()) {
            pfin[l] = links[l][iclock].particle;
            pfregions[l] = links[l][iclock].region;
          }
        }
        step(newevt, pfregions, pfin, outObj);
        stream.push_back(outObj);
      }
      // The links are over, but the conversion register and the tree still hold objects:
      // keep clocking with empty links until they are all out, so that nothing of this
      // event is lost. These cycles are appended to the stream like any other, one entry
      // per clock cycle, empty where nothing came out.
      const std::vector<l1ct::pfsorter::PFParticleEmu> nopf(nlinks);
      const std::vector<l1ct::PFRegionEmu> noregions(nlinks);
      // the first of these cycles is the one that lets the last converted particles into
      // the fifos, so it always has to be run
      step(false, noregions, nopf, outObj);
      stream.push_back(outObj);
      while (buffer().inFlight() != 0) {
        step(false, noregions, nopf, outObj);
        stream.push_back(outObj);
      }
    }

#ifdef CMSSW_GIT_HASH
    // writes the sorted candidates into event.sortedpf, so 'event' can't be const here
    void run(l1ct::Event &event) {
      run(event.pfinputs, event.out, event.sortedpf);
    }
#else
    // standalone version: run from raw links (one vector of particles per link)
    void run(const std::vector<l1ct::PFRegionEmu>& regions,
             const std::vector<std::vector<PFParticle>>& links,
             std::vector<l1ct::PuppiObjEmu>& out) {
      l1ct::PuppiObjEmu outObj;
      for (unsigned int iclock = 0; iclock < nClocks(); ++iclock) {
        bool newevt = (iclock == 0);
        std::vector<l1ct::pfsorter::PFParticleEmu> pfin(regions.size());
        std::vector<l1ct::PFRegionEmu> pfregions(regions.size());
        for (unsigned int l = 0; l < regions.size(); ++l) {
          if (iclock < links[l].size()) {
            pfin[l] = links[l][iclock];
            pfregions[l] = regions[l];
          }
        }
        step(newevt, pfregions, pfin, outObj);
        if (outObj.hwPt != 0)
          out.push_back(outObj);
      }
    }
#endif

  private:
    unsigned int nLinksCharged_ = l1ct::pfsorter::nLinksCharged;
    unsigned int nLinksNeutral_ = l1ct::pfsorter::nLinksNeutral;
    unsigned int nObjPerLink_ = l1ct::pfsorter::nObjPerLink;
    unsigned int nSlotsCharged_ = l1ct::pfsorter::nSlotsCharged;
    unsigned int nSlotsNeutral_ = l1ct::pfsorter::nSlotsNeutral;

    // one event of the continuous mode of runStream(): exactly nClocks() clock cycles,
    // newEvent on the first, links shorter than that padded with empty objects
    void runContinuous_(const std::vector<l1ct::PFInputRegion> &pfinputs,
                        const std::vector<l1ct::OutputRegion> &pfouts,
                        std::vector<l1ct::PuppiObjEmu> &stream) {
      const unsigned int nlinks = nLinks();
      std::vector<std::vector<PFLinkItem>> links;
      makePFLinks(pfinputs, pfouts, links);
      const unsigned int nclk = nClocks();
      // an event longer than the period would overlap the next one
      assert(links.empty() || links[0].size() <= nclk);
      stream.clear();
      stream.reserve(nclk);
      l1ct::PuppiObjEmu outObj;
      for (unsigned int iclock = 0; iclock < nclk; ++iclock) {
        std::vector<l1ct::pfsorter::PFParticleEmu> pfin(nlinks);
        std::vector<l1ct::PFRegionEmu> pfregions(nlinks);
        for (unsigned int l = 0; l < nlinks; ++l) {
          if (iclock < links[l].size()) {
            pfin[l] = links[l][iclock].particle;
            pfregions[l] = links[l][iclock].region;
          }
        }
        step(iclock == 0, pfregions, pfin, outObj);
        stream.push_back(outObj);
      }
    }

    // the output register of the conversion: what it holds enters the fifos on the next
    // clock cycle, together with the newEvent that came with it
    std::vector<l1ct::PuppiObjEmu> converted_;
    bool convertedNewEvent_ = false;

    void clearConverted_() {
      for (auto& o : converted_)
        o.clear();
      convertedNewEvent_ = false;
    }
  };
}  // namespace l1ct

#endif
