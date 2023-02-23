#ifndef tdr_regionizer_elements_ref_h
#define tdr_regionizer_elements_ref_h

#include "DataFormats/L1TParticleFlow/interface/layer1_emulator.h"

#include <list>
#include <vector>
#include <cassert>
#include <algorithm>

#ifdef CMSSW_GIT_HASH
#include "L1Trigger/Phase2L1ParticleFlow/interface/dbgPrintf.h"
#else
#include "../../utils/dbgPrintf.h"
#endif

namespace l1ct {
  namespace tdr_regionizer {

    inline int dphi_wrap(int local_phi) {
      if (local_phi > l1ct::Scales::INTPHI_PI)
        local_phi -= l1ct::Scales::INTPHI_TWOPI;
      else if (local_phi <= -l1ct::Scales::INTPHI_PI)
        local_phi += l1ct::Scales::INTPHI_TWOPI;
      return local_phi;
    }


    template <typename T>
    class PipeObject {
    public:
      PipeObject() {}
      PipeObject(const T& obj,
                 unsigned int phiindex,
                 unsigned int etaindex,
                 bool phioverlap,
                 bool etaoverlap,
                 int glbphi,
                 int glbeta,
                 unsigned int clk);

      const unsigned int getClock() { return linkobjclk_; }
      void setClock(unsigned int clock) { linkobjclk_ = clock; }
      const unsigned int getPhi() { return phiindex_; }
      const unsigned int getEta() { return etaindex_; }
      const bool getPhiOverlap() { return phioverlap_; }
      const bool getEtaOverlap() { return etaoverlap_; }
      const unsigned int getCount() { return objcount_; }
      unsigned int getCountAndInc() { return objcount_++; }
      void incCount() { objcount_++; }
      const int getPt() { return obj_.hwPt.to_int(); }
      const int getGlbPhi() { return glbphi_; }
      const int getGlbEta() { return glbeta_; }

      T getObj() { return obj_; }

    private:
      T obj_;
      unsigned int phiindex_, etaindex_;
      bool phioverlap_, etaoverlap_;
      int glbphi_, glbeta_;
      unsigned int linkobjclk_, objcount_;
    };

    template <typename T>
    class Pipe {
    public:
      Pipe(unsigned int nphi = 9) : clkindex_(0), nphi_(nphi) {}

      void addObj(
          T obj, unsigned int phiindex, unsigned int etaindex, bool phioverlap, bool etaoverlap, int glbphi, int glbeta);
      PipeObject<T>& getObj(unsigned int index) { return data_[index]; }
      T getRawObj(unsigned int index) { return data_[index].getObj(); }

      unsigned int getClock(unsigned int index = 0) { return getObj(index).getClock(); }
      void setClock(unsigned int clock, unsigned int index = 0) { return getObj(index).setClock(clock); }
      unsigned int getPhi(unsigned int index = 0) { return getObj(index).getPhi(); }
      unsigned int getEta(unsigned int index = 0) { return getObj(index).getEta(); }
      bool getPhiOverlap(unsigned int index = 0) { return getObj(index).getPhiOverlap(); }
      bool getEtaOverlap(unsigned int index = 0) { return getObj(index).getEtaOverlap(); }
      unsigned int getCount(unsigned int index = 0) { return getObj(index).getCount(); }
      unsigned int getCountAndInc(unsigned int index = 0) { return getObj(index).getCountAndInc(); }
      void incCount(unsigned int index = 0) { getObj(index).incCount(); }
      void erase(unsigned int index = 0) { data_.erase(data_.begin() + index); }
      int getPt(unsigned int index = 0) { return getObj(index).getPt(); }
      int getGlbPhi(unsigned int index = 0) { return getObj(index).getGlbPhi(); }
      int getGlbEta(unsigned int index = 0) { return getObj(index).getGlbEta(); }

      int getClosedIndexForObject(unsigned int index = 0);
      int getPipeIndexForObject(unsigned int index = 0);

      unsigned int getSize() { return data_.size(); }

      void reset() {
        clkindex_ = 0;
        data_.clear();
      }

    private:
      unsigned int clkindex_, nphi_;
      std::vector<PipeObject<T>> data_;
    };

    template <typename T>
    class Regionizer {
    public:
      Regionizer() {}
      Regionizer(
          unsigned int neta, unsigned int nphi,  //the number of eta and phi SRs in a big region (board)
          unsigned int nregions, // The total number of small regions in the full barrel
          unsigned int maxobjects,
          int bigRegionMin, int bigRegionMax,  // the phi range covered by this board
          int nclocks);

      void initSectors(const std::vector<DetectorSector<T>>& sectors);
      void initSectors(const DetectorSector<T>& sector);
      void initRegions(const std::vector<PFInputRegion>& regions);

      // is the given small region in the big region
      bool isInBigRegion(const PFRegionEmu& reg) const;

      unsigned int getSize() { return pipes_.size(); }
      unsigned int getPipeSize(unsigned int index) { return getPipe(index).getSize(); }

      bool setIndicesOverlaps(const T& obj,
                              unsigned int& phiindex,
                              unsigned int& etaindex,
                              bool& phioverlap,
                              bool& etaoverlap,
                              int& glbphi,
                              int& glbeta,
                              unsigned int index);

      void addToPipe(const T& obj, unsigned int index);
      void setPipe(const std::vector<T>& objvec, unsigned int index);
      void setPipes(const std::vector<std::vector<T>>& objvecvec);
      Pipe<T>& getPipe(unsigned int index) { return pipes_[index]; }

      int getPipeTime(int linkIndex, int linkTimeOfObject, int linkAlgoClockRunningTime);
      int popLinkObject(int linkIndex, int currentTimeOfObject);
      int timeNextFromIndex(unsigned int index, int time) { return getPipeTime(index, pipes_[index].getClock(), time); }

      void initTimes();

      int getClosedIndexForObject(unsigned int linknum, unsigned int index = 0) {
        return pipes_[linknum].getClosedIndexForObject(index);
      }
      int getPipeIndexForObject(unsigned int linknum, unsigned int index = 0) {
        return pipes_[linknum].getPipeIndexForObject(index);
      }
      void addToSmallRegion(unsigned int linkNum, unsigned int index = 0);

      void run(bool debug = false);

      void reset();

      std::vector<T> getSmallRegion(unsigned int index);

      void printDebug(int count) {
        dbgCout() << count << "\tindex\tpt\teta\tphi" << std::endl;
        dbgCout() << "PIPES" << std::endl;
        for (unsigned int i = 0; i < getSize(); i++) {
          for (unsigned int j = 0; j < getPipeSize(i); j++) {
            dbgCout() << "\t" << i << " " << j << "\t" << getPipe(i).getPt(j) << "\t" << getPipe(i).getGlbEta(j) << "\t"
                      << getPipe(i).getGlbPhi(j) << std::endl;
          }
          dbgCout() << "-------------------------------" << std::endl;
        }
        dbgCout() << "SMALL REGIONS" << std::endl;
        for (unsigned int i = 0; i < nregions_; i++) {
          for (unsigned int j = 0; j < smallRegionObjects_[i].size(); j++) {
            dbgCout() << "\t" << i << " " << j << "\t" << smallRegionObjects_[i][j].hwPt.to_int() << "\t"
                      // << smallRegionObjects_[i][j].hwEta.to_int() + regionmap_[i].eta << "\t"
                      // << smallRegionObjects_[i][j].hwPhi.to_int() + regionmap_[i].phi 
                       << std::endl;
          }
          dbgCout() << "-------------------------------" << std::endl;
        }
        dbgCout() << "TIMES" << std::endl;
        for (unsigned int i = 0; i < timeOfNextObject_.size(); i++) {
          dbgCout() << "  " << timeOfNextObject_[i];
        }
        dbgCout() << "\n-------------------------------" << std::endl;
      }

    private:

      // this function is for sorting small regions first in phi and then in eta.
      // It takes regions_ indices
      bool sortRegionsRegular(size_t a, size_t b) const;

      /// The numbers of eta and phi in a big region (board)
      unsigned int neta_, nphi_;
      /// The total number of small regions in the barrel (not just in the board)
      unsigned int nregions_;
      /// The maximum number of objects to output per small region
      unsigned int maxobjects_;
      /// The number of input sectors for this type of device
      unsigned int nsectors_;
      /// the minimumum phi of this board
      int bigRegionMin_;
      /// the maximum phi of this board
      int bigRegionMax_;
      /// the number of clocks to receive one event
      int nclocks_;

      /// the region information assopciated with each input sector
      std::vector<l1ct::PFRegionEmu> sectors_;

      /// the region information associated with each SR
      std::vector<l1ct::PFRegionEmu> regions_;

      /// indices of regions that are in the big region (board)
      std::vector<size_t> regionmap_;

      std::vector<Pipe<T>> pipes_;
      std::vector<int> timeOfNextObject_;
      std::vector<std::vector<T>> smallRegionObjects_;  //keep count to see if small region is full
    };

  }  // namespace  tdr_regionizer
}  // namespace l1ct

#endif
