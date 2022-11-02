#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTCDistribution_SA.h"

using namespace std;
using namespace l1thgcfirmware;

HGCalTCDistribution::HGCalTCDistribution(const ClusterAlgoConfig& config) : config_(config) {}

// Rearranges/distributes input streams of TCs so that each output stream contains TCs corresponding to one histogram bin
void HGCalTCDistribution::runTriggerCellDistribution(const HGCalTriggerCellSAPtrCollections& triggerCellsIn,
                                                     HGCalTriggerCellSAPtrCollection& triggerCellsOut) const {
  HGCalTriggerCellSAShrPtrCollection triggerCellsWork;
  triggerCellInput(triggerCellsIn, triggerCellsWork);
  triggerCellDistribution0(triggerCellsWork);
  HGCalTriggerCellSAShrPtrCollections tcDistGrid1;
  triggerCellDistribution1(triggerCellsWork, tcDistGrid1);
  HGCalTriggerCellSAShrPtrCollections tcDistGrid2;
  triggerCellDistribution2(tcDistGrid1, triggerCellsWork, tcDistGrid2);
  HGCalTriggerCellSAShrPtrCollections tcDistGrid3;
  triggerCellDistribution3(tcDistGrid2, triggerCellsWork, tcDistGrid3);
  triggerCellDistribution4(triggerCellsWork);
  triggerCellDistribution5(tcDistGrid3, triggerCellsOut);
}

// Unrolls input TCs (2D vector, frame vs input/link) to 1D vector
void HGCalTCDistribution::triggerCellInput(const HGCalTriggerCellSAPtrCollections& inputTCs,
                                           HGCalTriggerCellSAShrPtrCollection& outputTCs) const {
  outputTCs.clear();
  for (unsigned int iFrame = 0; iFrame < inputTCs.size(); ++iFrame) {
    for (unsigned int iInput = 0; iInput < inputTCs[iFrame].size(); ++iInput) {
      auto& tc = inputTCs[iFrame][iInput];

      tc->setIndex(iInput);
      tc->setClock(iFrame + 1);

      if (tc->dataValid()) {
        // Copying input TCs and converting to shared pointer to avoid further copies in TCDist
        outputTCs.push_back(make_shared<HGCalTriggerCell>(*tc));
      }
    }
  }
}

// Adds empty input to each group of 24, making them groups of 25
void HGCalTCDistribution::triggerCellDistribution0(const HGCalTriggerCellSAShrPtrCollection& triggerCellsIn) const {
  for (auto& tc : triggerCellsIn) {
    unsigned int newIndex = tc->index() + int(tc->index() / (config_.cInputs() / config_.n60Sectors()));
    tc->setIndex(newIndex);
  }
}

// In each group of 25 streams, assign coarse phi region index to sort key of each TC (0-5 within each large phi region)
void HGCalTCDistribution::triggerCellDistribution1(
    const HGCalTriggerCellSAShrPtrCollection& triggerCellsIn,
    HGCalTriggerCellSAShrPtrCollections& outTriggerCellDistributionGrid) const {
  initializeTriggerCellDistGrid(outTriggerCellDistributionGrid, config_.cClocks(), config_.cInputs2());

  const unsigned int stepLatency = config_.getStepLatency(Dist1);
  for (auto& tc : triggerCellsIn) {
    tc->addLatency(stepLatency);
    unsigned int sector = int(tc->index() / (config_.cInputs2() / config_.n60Sectors()));
    outTriggerCellDistributionGrid[tc->clock() - 2][tc->index()] = tc;
    for (int iSortKey = config_.nCoarsePhiDist1() - 1; iSortKey >= 0; --iSortKey) {
      if (int(tc->phi() % config_.phiNValues()) >
          int(config_.cColumns() * iSortKey + config_.phiNValues() / config_.n60Sectors() * sector)) {
        tc->setSortKey(iSortKey);
        break;
      }
    }
  }
}

// First distribution layer
// Use 15 5-to-6 distribution servers to merge multiple streams into 15x6=90 streams
// Within each 3x25 groups of input streams, 5 distribution servers organise the TCs into 5 groups of 6 streams.  The 6 streams within one group correspond to TCs with the same coarse phi index
// Not enough documentation on distribution server vhdl to give more meaningful comments
void HGCalTCDistribution::triggerCellDistribution2(
    const HGCalTriggerCellSAShrPtrCollections& inTriggerCellDistributionGrid,
    HGCalTriggerCellSAShrPtrCollection& triggerCellsOut,
    HGCalTriggerCellSAShrPtrCollections& outTriggerCellDistributionGrid) const {
  const unsigned int latency = config_.getLatencyUpToAndIncluding(Dist2);

  triggerCellsOut.clear();
  initializeTriggerCellDistGrid(outTriggerCellDistributionGrid, config_.cClocks(), config_.cInt());

  runDistServers(inTriggerCellDistributionGrid,
                 outTriggerCellDistributionGrid,
                 triggerCellsOut,
                 latency,
                 config_.nDistServers1(),
                 config_.distServer1_nIn(),
                 config_.distServer1_nOut(),
                 config_.distServer1_nInterleave(),
                 true);
}

// Arrange streams with the same coarse phi index within each group of 30 streams to be next to each other
// Streams then correspond to 18 groups of 5, which each group corresponding to TCs with similar phi (same coarse phi index)
void HGCalTCDistribution::triggerCellDistribution3(
    const HGCalTriggerCellSAShrPtrCollections& inTriggerCellDistributionGrid,
    HGCalTriggerCellSAShrPtrCollection& triggerCellsOut,
    HGCalTriggerCellSAShrPtrCollections& outTriggerCellDistributionGrid) const {
  triggerCellsOut.clear();
  initializeTriggerCellDistGrid(outTriggerCellDistributionGrid, config_.cClocks(), config_.cInt());
  for (unsigned int iClock = 0; iClock < config_.cClocks(); ++iClock) {
    for (unsigned int i = 0; i < config_.n60Sectors(); ++i) {
      for (unsigned int k = 0; k < config_.distServer1_nOut(); ++k) {
        for (unsigned int j = 0; j < config_.distServer1_nIn(); ++j) {
          auto& tc = inTriggerCellDistributionGrid[iClock][config_.distServer1_nOut() * config_.distServer1_nIn() * i +
                                                           config_.distServer1_nOut() * j + k];
          if (tc->dataValid()) {
            tc->setIndex((config_.distServer1_nOut() * config_.distServer1_nIn() * i) +
                         (config_.distServer1_nIn() * k) + j);
            triggerCellsOut.push_back(tc);
            outTriggerCellDistributionGrid[iClock - 2][tc->index()] = tc;
          }
        }
      }
    }
  }
}

// In each group of 5 streams, assign fine phi index to each TC (0-5)
void HGCalTCDistribution::triggerCellDistribution4(const HGCalTriggerCellSAShrPtrCollection& triggerCellsIn) const {
  const unsigned int stepLatency = config_.getStepLatency(Dist4);
  for (auto& lCell : triggerCellsIn) {
    lCell->addLatency(stepLatency);
    unsigned int sector = int(
        lCell->index() /
        (config_.nDistServers1() /
         config_
             .n60Sectors()));  // Is config_.nDistServers1() / config_.n60Sectors() correct?  Equal to 5, which is what's needed (don't think it's a coincidence...)
    for (int iSortKey = config_.nCoarsePhiDist2() - 1; iSortKey >= 0; --iSortKey) {
      if (int(lCell->phi() % config_.phiNValues()) >
          int((config_.n60Sectors() * config_.distServer1_nOut() * iSortKey) + (config_.cColumns() * sector))) {
        lCell->setSortKey(iSortKey);
        break;
      }
    }
  }
}

// Second distribution layer
// Use 18 5-to-6 distribution servers to merge multiple streams into 18x6=108 streams (corresponding to the number of histogram bins in phi)
// Within each 18x5 groups of streams, 1 distribution server organises the TCs into 6 streams.  The 6 streams correspond to the TCs with the same fine phi index i.e the same histogram bin in phi
// Not enough documentation on distribution server vhdl to give more meaningful comments
void HGCalTCDistribution::triggerCellDistribution5(
    const HGCalTriggerCellSAShrPtrCollections& inTriggerCellDistributionGrid,
    HGCalTriggerCellSAPtrCollection& triggerCellsOut) const {
  const unsigned int latency = config_.getLatencyUpToAndIncluding(Dist5);

  // Dummy distribution grid?  After writing the runDistServers function to avoid duplicating identical code from triggerCellDistribution2, I realised the second possible use case of runDistServers isn't completely identical i.e. the triggerCellDistributionGrid isn't used/set
  HGCalTriggerCellSAShrPtrCollections outTriggerCellDistributionGrid;
  HGCalTriggerCellSAShrPtrCollection triggerCellsOutTemp;
  runDistServers(inTriggerCellDistributionGrid,
                 outTriggerCellDistributionGrid,
                 triggerCellsOutTemp,
                 latency,
                 config_.nDistServers2(),
                 config_.distServer2_nIn(),
                 config_.distServer2_nOut(),
                 config_.distServer2_nInterleave(),
                 false);

  // Convert final output to unique_ptr
  triggerCellsOut.clear();
  for (const auto& tc : triggerCellsOutTemp) {
    triggerCellsOut.push_back(make_unique<HGCalTriggerCell>(*tc));
  }
}

// Initialize distribution grids with invalid TCs
void HGCalTCDistribution::initializeTriggerCellDistGrid(HGCalTriggerCellSAShrPtrCollections& grid,
                                                        unsigned int nX,
                                                        unsigned int nY) const {
  for (unsigned int iX = 0; iX < nX; ++iX) {
    HGCalTriggerCellSAShrPtrCollection temp;
    for (unsigned int iY = 0; iY < nY; ++iY) {
      temp.emplace_back(make_shared<HGCalTriggerCell>());
    }
    grid.push_back(temp);
  }
}

// Not enough documentation on distribution server vhdl to give meaningful comments
void HGCalTCDistribution::runDistServers(const HGCalTriggerCellSAShrPtrCollections& gridIn,
                                         HGCalTriggerCellSAShrPtrCollections& gridOut,
                                         HGCalTriggerCellSAShrPtrCollection& tcsOut,
                                         unsigned int latency,
                                         unsigned int nDistServers,
                                         unsigned int nInputs,
                                         unsigned int nOutputs,
                                         unsigned int nInterleave,
                                         bool setOutputGrid) const {
  vector<DistServer> distServers(nDistServers, DistServer(nInputs, nOutputs, nInterleave));

  for (unsigned int iClock = 0; iClock < config_.cClocks(); ++iClock) {
    for (unsigned int iDistServer = 0; iDistServer < nDistServers; ++iDistServer) {
      auto first = gridIn[iClock].cbegin() + nInputs * iDistServer;
      auto last = gridIn[iClock].cbegin() + nInputs * (iDistServer + 1);
      HGCalTriggerCellSAShrPtrCollection inCells(first, last);
      HGCalTriggerCellSAShrPtrCollection lCells = distServers[iDistServer].clock(inCells);

      for (unsigned int iOutput = 0; iOutput < lCells.size(); ++iOutput) {
        auto& tc = lCells[iOutput];
        if (tc->dataValid()) {
          tc->setIndex(nOutputs * iDistServer + iOutput);
          tc->setClock(iClock + latency);

          tcsOut.push_back(tc);
          if (setOutputGrid) {
            gridOut[iClock][tc->index()] = tc;
          }
        }
      }
    }
  }
}
