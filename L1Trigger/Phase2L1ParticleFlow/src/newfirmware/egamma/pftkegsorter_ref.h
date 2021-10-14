#ifndef FIRMWARE_PFTKEGSORTER_REF_H
#define FIRMWARE_PFTKEGSORTER_REF_H

#include <cstdio>
#include <vector>

#ifdef CMSSW_GIT_HASH
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "../dataformats/layer1_emulator.h"
#include "L1Trigger/Phase2L1ParticleFlow/src/dbgPrintf.h"
#else
#include "../../../dataformats/layer1_emulator.h"
#include "../../utils/dbgPrintf.h"
#endif

namespace edm {
  class ParameterSet;
}

namespace l1ct {
  class PFTkEGSorterEmulator {
    public:
      PFTkEGSorterEmulator(bool splitEGObjectsToBoards=false,
                           const unsigned int nObjToSort=6,
                           const unsigned int nObjSorted=54,
                           const unsigned int nBoards=5) //FIXME: Numbers to be revised
        : splitEGObjectsToBoards_(splitEGObjectsToBoards),
          nObjToSort_(nObjToSort),
          nObjSorted_(nObjSorted),
          nBoards_(nBoards),
          debug_(false) {}

      ~PFTkEGSorterEmulator() {};
      
      void setDebug(bool debug = true) { debug_ = debug; };
      
#ifdef CMSSW_GIT_HASH
      PFTkEGSorterEmulator(const edm::ParameterSet &iConfig)
        : PFTkEGSorterEmulator(
          iConfig.getUntrackedParameter<bool>("splitEGObjectsToBoards"),
          iConfig.getParameter<uint32_t>("nObjToSort"),
          iConfig.getParameter<uint32_t>("nObjSorted"),
          iConfig.getParameter<uint32_t>("nBoards")) {
        setDebug(iConfig.getUntrackedParameter<bool>("debug", false));
      }
#endif

      template <typename T>
      void run(const std::vector<l1ct::PFInputRegion> & pfregions, 
               const std::vector<l1ct::OutputRegion> & outregions, 
               std::vector<std::vector<T> > & eg_sorted_inBoard){

        std::vector<std::vector<T> >  eg_unsorted_inBoard = eg_sorted_inBoard;
        mergeEGObjFromRegions<T>(pfregions, outregions, eg_unsorted_inBoard);

        if(debug_) {
          dbgCout() << "\nUNSORTED\n";
          for (int i=0, ni=eg_unsorted_inBoard.size(); i<ni; i++) {
            dbgCout() << "\nBoard "<<i<<"\n";
            for (int j=0, nj=eg_unsorted_inBoard[i].size(); j<nj; j++) dbgCout() << "EG["<<j<<"]: pt = " << eg_unsorted_inBoard[i][j].hwPt << ",\t eta = " << eg_unsorted_inBoard[i][j].hwEta << ",\t phi = " << eg_unsorted_inBoard[i][j].hwPhi << "\n";
          }
        }

        if(debug_) dbgCout() << "\nSORTED\n";

        for (int i=0, ni=eg_unsorted_inBoard.size(); i<ni; i++) {
          eg_sorted_inBoard[i] = eg_unsorted_inBoard[i];
          std::reverse(eg_sorted_inBoard[i].begin(), eg_sorted_inBoard[i].end());
          std::stable_sort(eg_sorted_inBoard[i].begin(), eg_sorted_inBoard[i].end(), comparePt<T>);
          eg_sorted_inBoard[i].resize(nObjSorted_);
          while (!eg_sorted_inBoard[i].empty() && eg_sorted_inBoard[i].back().hwPt == 0) eg_sorted_inBoard[i].pop_back(); // Zero suppression

        
          if(debug_) {
            dbgCout() << "\nBoard "<<i<<"\n";
            for (int j=0, nj=eg_sorted_inBoard[i].size(); j<nj; j++) dbgCout() << "EG["<<j<<"]: pt = " << eg_sorted_inBoard[i][j].hwPt << ",\t eta = " << eg_sorted_inBoard[i][j].hwEta << ",\t phi = " << eg_sorted_inBoard[i][j].hwPhi << "\n";
          }
        }
      }


    private:
      bool splitEGObjectsToBoards_;
      unsigned int nObjToSort_, nObjSorted_, nBoards_;
      bool debug_;

      void extractEGObjEmu(const l1ct::OutputRegion outregion, std::vector<l1ct::EGIsoObjEmu> & eg) { eg = outregion.egphoton; }
      void extractEGObjEmu(const l1ct::OutputRegion outregion, std::vector<l1ct::EGIsoEleObjEmu> & eg) { eg = outregion.egelectron; }

      template <typename T>
      static bool comparePt(T obj1, T obj2) { return (obj1.hwPt > obj2.hwPt); }

      template <typename T>
      void mergeEGObjFromRegions(const std::vector<l1ct::PFInputRegion> & pfregions,
                                 const std::vector<l1ct::OutputRegion> & outregions,
                                 std::vector<std::vector<T> > & eg_unsorted_inBoard){
        
        for (int i=0, ni=outregions.size(); i<ni; i++) {
          unsigned int x;
          if(splitEGObjectsToBoards_) {
            l1ct::glbeta_t eta = pfregions[i].region.hwEtaCenter;

            if ( abs(eta) < l1ct::Scales::makeGlbEta(0.5) ) x = 0;
            else if ( abs(eta) < l1ct::Scales::makeGlbEta(1.5) ) x = (eta < 0 ? 1 : 2);
            else if ( abs(eta) < l1ct::Scales::makeGlbEta(2.5) ) x = (eta < 0 ? 3 : 4);
            //else /*if ( fabs(eta) < 3.0 )*/x = 5;
            /*else x = 6;*/ // HF
          }
          else x=0;

          if(debug_) dbgCout() << "\nOutput Region "<<i<<": eta = "<<pfregions[i].region.floatEtaCenter()<<" and phi = "<<pfregions[i].region.floatPhiCenter()<<" \n";
          
          std::vector<T> eg_tmp;
          extractEGObjEmu(outregions[i], eg_tmp);
          for (int j=0, nj=( eg_tmp.size()>nObjToSort_ ? nObjToSort_ : eg_tmp.size() ); j<nj; j++) {
            if(debug_) dbgCout() << "EG["<<j<<"] pt = " << eg_tmp[j].hwPt << ",\t eta = " << eg_tmp[j].hwEta << ",\t phi = " << eg_tmp[j].hwPhi << "\n";
            eg_unsorted_inBoard[x].push_back(eg_tmp[j]);
          }
          if(debug_) dbgCout() << "\n";
        }
      }
  };
}

#endif

