#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClustering_SA.h"

#include <cmath>
#include <algorithm>

using namespace std;
using namespace l1thgcfirmware;

HGCalHistoClustering::HGCalHistoClustering(const ClusterAlgoConfig& config) : config_(config) {}

void HGCalHistoClustering::runClustering(const HGCalTriggerCellSAPtrCollection& triggerCellsIn,
                                         const HGCalHistogramCellSAPtrCollection& histogramIn,
                                         HGCalTriggerCellSAShrPtrCollection& clusteredTriggerCellsOut,
                                         CentroidHelperPtrCollection& readoutFlagsOut,
                                         HGCalClusterSAPtrCollection& protoClusters) const {
  HGCalTriggerCellSAShrPtrCollection unclusteredTriggerCells;
  clusterizer(triggerCellsIn, histogramIn, clusteredTriggerCellsOut, unclusteredTriggerCells, readoutFlagsOut);
  triggerCellToCluster(clusteredTriggerCellsOut, protoClusters);
}

// Main implementation of clustering
// Takes histogram containing seeds and streams of TCs (each stream corresponds to one column of the histogram)
// Outputs clustered and unclustered TCs, and readoutFlags which contain info on which clustered TCs belong to each cluster
// Require more comments on firmware to provide more meanginful comments here
void HGCalHistoClustering::clusterizer(const HGCalTriggerCellSAPtrCollection& triggerCellsIn,
                                       const HGCalHistogramCellSAPtrCollection& histogram,
                                       HGCalTriggerCellSAShrPtrCollection& clusteredTriggerCellsOut,
                                       HGCalTriggerCellSAShrPtrCollection& unclusteredTriggerCellsOut,
                                       CentroidHelperPtrCollection& readoutFlagsOut) const {
  unsigned int seedCounter = 0;
  CentroidHelperPtrCollections fifos(config_.nFifos());
  vector<unsigned int> clock(config_.cColumns(), config_.clusterizerMagicTime());
  const unsigned dummy_entries_latched = 1;
  CentroidHelperShrPtrCollection latched(
      config_.nFifos() + dummy_entries_latched,
      make_shared<
          CentroidHelper>());  // 1 extra (dummy) entry compared to fifos, to match firmware behaviour (avoids issues with index wrap-around)

  HGCalTriggerCellSAShrPtrCollections clusteredTriggerCells(config_.cColumns(), HGCalTriggerCellSAShrPtrCollection());
  HGCalTriggerCellSAShrPtrCollections unclusteredTriggerCells(config_.cColumns(), HGCalTriggerCellSAShrPtrCollection());
  CentroidHelperPtrCollections readoutFlags(config_.cColumns());

  HGCalTriggerCellSAShrPtrCollectionss triggerCellBuffers(
      config_.cColumns(), HGCalTriggerCellSAShrPtrCollections(config_.cRows(), HGCalTriggerCellSAShrPtrCollection()));
  for (const auto& tc : triggerCellsIn) {
    // Temp copy of tc whilst moving from shared to unique ptr
    triggerCellBuffers.at(tc->index()).at(tc->sortKey()).push_back(make_shared<HGCalTriggerCell>(*tc));
  }

  for (unsigned int iRow = 0; iRow < config_.cRows(); ++iRow) {
    for (unsigned int j = 0; j < config_.nColumnsPerFifo(); ++j) {
      for (unsigned int k = 0; k < config_.nFifos(); ++k) {
        unsigned int col = config_.firstSeedBin() + (config_.nColumnsPerFifo() * k) + j;
        const auto& cell = histogram.at(config_.cColumns() * iRow + col);
        if (cell->S() > 0) {
          auto ch = make_unique<CentroidHelper>(cell->clock() + 1 + j,
                                                config_.nColumnsPerFifo() * k + j,
                                                cell->index(),
                                                cell->sortKey(),
                                                cell->S(),
                                                cell->X(),
                                                cell->Y(),
                                                true);
          fifos[k].push_back(move(ch));
          ++seedCounter;
        }
      }
    }
  }
  while (seedCounter > 0) {
    for (unsigned int i = 0; i < config_.nFifos(); ++i) {
      if (!latched[i]->dataValid()) {
        if (!fifos[i].empty()) {
          latched[i] = move(fifos[i][0]);
          fifos[i].erase(fifos.at(i).begin());
        }
      }
    }

    const unsigned dummy_entries_accepted = 2;
    CentroidHelperShrPtrCollection accepted(
        config_.nFifos() + dummy_entries_accepted,
        make_shared<
            CentroidHelper>());  // 2 extra (dummy) entry compared to fifos (1 compared to latched), to match firmware behaviour (avoids issues with index wrap-around)
    CentroidHelperShrPtrCollection lastLatched(latched);

    for (unsigned int i = 0; i < config_.nFifos(); ++i) {
      // Different implementation to python emulator
      // For i=0, i-1=-1, which would give the last element of lastLatched in python, but is out of bounds in C++
      // Similar for i=17 (==config_.nFifos()-1)
      // Need to find out intended behaviour
      bool deltaMinus =
          (i > 0) ? (lastLatched[i]->column() - lastLatched[i - 1]->column()) > config_.nColumnFifoVeto() : true;
      bool deltaPlus = (i < config_.nFifos() - 1)
                           ? (lastLatched[i + 1]->column() - lastLatched[i]->column()) > config_.nColumnFifoVeto()
                           : true;

      bool compareEMinus = (i > 0) ? (lastLatched[i]->energy() > lastLatched[i - 1]->energy()) : true;
      bool compareEPlus =
          (i < config_.nFifos() - 1) ? (lastLatched[i]->energy() >= lastLatched[i + 1]->energy()) : true;

      if (lastLatched[i]->dataValid()) {
        // Similar out of bounds issue here
        // if ( ( !lastLatched[i+1]->dataValid() || compareEPlus || deltaPlus ) && ( !lastLatched[i-1]->dataValid() || compareEMinus || deltaMinus ) ) {

        bool accept = true;
        if (lastLatched.size() > i + 1) {
          if (!lastLatched[i + 1]->dataValid() || compareEPlus || deltaPlus) {
            accept = true;
          } else {
            accept = false;
          }
        }

        if (i > 0) {
          if (!lastLatched[i - 1]->dataValid() || compareEMinus || deltaMinus) {
            accept = true;
          } else {
            accept = false;
          }
        }

        if (accept) {
          accepted[i] = latched[i];
          latched[i] = make_shared<CentroidHelper>();
          --seedCounter;
        }
      }
    }

    for (const auto& a : accepted) {
      if (a->dataValid()) {
        for (unsigned int iCol = a->column() - config_.nColumnsForClustering();
             iCol < a->column() + config_.nColumnsForClustering() + 1;
             ++iCol) {
          clock[iCol] = clock[a->column()];
        }
      }
    }

    for (const auto& a : accepted) {
      if (a->dataValid()) {
        // Firmware does support layer-dependent delta R threshold
        // But all thresholds currently set to same value
        // And python emulator just has one threshold, which is why it's not configurable here yet
        // i.e. I directly translated what is implemented in the python emulator, not the actual firmware
        unsigned int dR2Cut = config_.deltaR2Cut();
        unsigned int T = 0;

        for (unsigned int iCol = a->column() - config_.nColumnsForClustering();
             iCol < a->column() + config_.nColumnsForClustering() + 1;
             ++iCol) {
          const unsigned stepLatency = 8;
          clock[iCol] += stepLatency;
          for (int k = -1 * config_.nRowsForClustering(); k < int(config_.nRowsForClustering()) + 1; ++k) {
            int row = a->row() + k;
            if (row < 0)
              continue;
            if (row >= int(config_.cRows()))
              continue;  // Not in python emulator, but required to avoid out of bounds access
            if (triggerCellBuffers[iCol][row].empty()) {
              clock[iCol] += 1;
              continue;
            }
            for (auto& tc : triggerCellBuffers[iCol][row]) {
              clock[iCol] += 1;

              unsigned int r1 = tc->rOverZ();
              unsigned int r2 = a->Y();
              int dR = r1 - r2;
              unsigned int absDPhi = abs(int(tc->phi()) - int(a->X()));
              unsigned int dR2 = dR * dR;
              unsigned int cosTerm = (absDPhi > config_.nBinsCosLUT()) ? 2047 : config_.cosLUT(absDPhi);

              const unsigned a = 128;   // 2^7
              const unsigned b = 1024;  // 2^10
              dR2 += int(r1 * r2 / a) * cosTerm / b;
              tc->setClock(clock[iCol] + 1);
              if (clock[iCol] > T)
                T = clock[iCol];

              if (dR2 < dR2Cut) {
                clusteredTriggerCells[iCol].push_back(tc);
              } else {
                unclusteredTriggerCells[iCol].push_back(tc);
              }
            }
          }

          for (const auto& tc : clusteredTriggerCells[iCol]) {
            auto tcMatch = std::find_if(
                triggerCellBuffers[iCol][tc->sortKey()].begin(),
                triggerCellBuffers[iCol][tc->sortKey()].end(),
                [&](const HGCalTriggerCellSAShrPtr tcToMatch) {
                  bool isMatch = tc->index() == tcToMatch->index() && tc->rOverZ() == tcToMatch->rOverZ() &&
                                 tc->layer() == tcToMatch->layer() && tc->energy() == tcToMatch->energy() &&
                                 tc->phi() == tcToMatch->phi() && tc->sortKey() == tcToMatch->sortKey() &&
                                 tc->deltaR2() == tcToMatch->deltaR2() && tc->dX() == tcToMatch->dX() &&
                                 tc->Y() == tcToMatch->Y() && tc->frameValid() == tcToMatch->frameValid() &&
                                 tc->dataValid() == tcToMatch->dataValid() && tc->clock() == tcToMatch->clock();
                  return isMatch;
                });

            if (tcMatch != triggerCellBuffers[iCol][tc->sortKey()].end()) {
              triggerCellBuffers[iCol][tc->sortKey()].erase(tcMatch);
            }
          }
        }

        for (unsigned int iCol = a->column() - config_.nColumnsForClustering();
             iCol < a->column() + config_.nColumnsForClustering() + 1;
             ++iCol) {
          clock[iCol] = T + 1;

          CentroidHelperPtr readoutFlag = make_unique<CentroidHelper>(T - 2, iCol, true);
          const unsigned stepLatency = 14;
          if (readoutFlag->clock() == config_.clusterizerMagicTime() + stepLatency) {
            readoutFlag->setClock(readoutFlag->clock() + 1);
          }

          readoutFlags[iCol].push_back(move(readoutFlag));
        }
      }
    }
  }

  const unsigned largeReadoutTime = 1000;
  for (unsigned int i = 0; i < largeReadoutTime;
       ++i) {  // Magic numbers - a large number to ensure we read out all clustered trigger cells etc.?
    for (unsigned int iCol = 0; iCol < config_.cColumns(); ++iCol) {
      for (const auto& clustered : clusteredTriggerCells[iCol]) {
        if (clustered->clock() == config_.clusterizerMagicTime() + i) {
          clusteredTriggerCellsOut.push_back(clustered);
        }
      }

      for (const auto& unclustered : unclusteredTriggerCells[iCol]) {
        if (unclustered->clock() == config_.clusterizerMagicTime() + i) {
          unclusteredTriggerCellsOut.push_back(unclustered);
        }
      }

      for (auto& readoutFlag : readoutFlags[iCol]) {
        if (readoutFlag) {
          if (readoutFlag->clock() == config_.clusterizerMagicTime() + i) {
            // TODO : Check if we can move the readoutFlag and leave a nullptr
            // Or if the readoutFlag could be used again later on
            readoutFlagsOut.push_back(move(readoutFlag));
          }
        }
      }
    }
  }
}

// Converts clustered TCs into cluster object (one for each TC) ready for accumulation
void HGCalHistoClustering::triggerCellToCluster(const HGCalTriggerCellSAShrPtrCollection& clusteredTriggerCells,
                                                HGCalClusterSAPtrCollection& clustersOut) const {
  const unsigned int stepLatency = config_.getStepLatency(TriggerCellToCluster);

  clustersOut.clear();
  for (const auto& tc : clusteredTriggerCells) {
    auto cluster = make_unique<HGCalCluster>(tc->clock() + stepLatency, tc->index(), true, true);

    // Cluster from single TC
    // Does this ever happen?
    // Removed from newer versions of firmware in any case, but leave for now
    const unsigned singleTCDeltaR2Thresh = 25000;
    if (tc->deltaR2() >= singleTCDeltaR2Thresh) {
      clustersOut.push_back(move(cluster));
      continue;
    }

    const unsigned weightFactor = 4;
    unsigned long int s_TC_W = (int(tc->energy() / weightFactor) == 0) ? 1 : tc->energy() / weightFactor;
    unsigned long int s_TC_Z = config_.depth(tc->layer());

    unsigned int triggerLayer = config_.triggerLayer(tc->layer());
    const unsigned nBitsESums = 18;  // Need to double check this is correct description of constant
    unsigned int s_E_EM =
        ((((unsigned long int)tc->energy() * config_.layerWeight_E_EM(triggerLayer)) + config_.correction()) >>
         nBitsESums);
    if (s_E_EM > config_.saturation())
      s_E_EM = config_.saturation();

    unsigned int s_E_EM_core =
        (((unsigned long int)tc->energy() * config_.layerWeight_E_EM_core(triggerLayer) + config_.correction()) >>
         nBitsESums);
    if (s_E_EM_core > config_.saturation())
      s_E_EM_core = config_.saturation();

    // Alternative constructor perhaps?
    cluster->set_n_tc(1);
    cluster->set_n_tc_w(1);

    cluster->set_e((config_.layerWeight_E(triggerLayer) == 1) ? tc->energy() : 0);
    cluster->set_e_h_early((config_.layerWeight_E_H_early(triggerLayer) == 1) ? tc->energy() : 0);

    cluster->set_e_em(s_E_EM);
    cluster->set_e_em_core(s_E_EM_core);

    cluster->set_w(s_TC_W);
    cluster->set_w2(s_TC_W * s_TC_W);

    cluster->set_wz(s_TC_W * s_TC_Z);
    cluster->set_weta(0);
    cluster->set_wphi(s_TC_W * tc->phi());
    cluster->set_wroz(s_TC_W * tc->rOverZ());

    cluster->set_wz2(s_TC_W * s_TC_Z * s_TC_Z);
    cluster->set_weta2(0);
    cluster->set_wphi2(s_TC_W * tc->phi() * tc->phi());
    cluster->set_wroz2(s_TC_W * tc->rOverZ() * tc->rOverZ());

    const unsigned nTriggerLayers = 36;  // Should get from config/elsewhere in CMSSW
    cluster->set_layerbits(cluster->layerbits() | (((unsigned long int)1) << (nTriggerLayers - triggerLayer)));
    cluster->set_sat_tc(cluster->e() == config_.saturation() || cluster->e_em() == config_.saturation());
    cluster->set_shapeq(1);

    // Temp copy of TC whilst reducing use of shared ptr
    cluster->add_constituent(make_shared<HGCalTriggerCell>(*tc));
    clustersOut.push_back(move(cluster));
  }
}