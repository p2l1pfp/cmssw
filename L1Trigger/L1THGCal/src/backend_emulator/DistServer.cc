#include "L1Trigger/L1THGCal/interface/backend_emulator/DistServer.h"

#include <algorithm>
#include <iostream>

using namespace std;
using namespace l1thgcfirmware;

DistServer::DistServer(unsigned int nInputs, unsigned int nOutputs, unsigned int nInterleaving)
    : nOutputs_(nOutputs), nInterleaving_(nInterleaving), inputs_(nInputs) {
  for (unsigned int iInput = 0; iInput < inputs_.size(); ++iInput) {
    addr_.emplace_back(this->nInterleaving(), 0);
    for (unsigned int iInterleave = 0; iInterleave < this->nInterleaving(); ++iInterleave) {
      addr_[iInput][iInterleave] = iInterleave;
    }
  }
}

// Main method of distribution server, executed on each clock
// Not enough documentation on distribution server vhdl to give meaningful comments
// Other than stating it's very complicated
HGCalTriggerCellSAShrPtrCollection DistServer::clock(HGCalTriggerCellSAShrPtrCollection& data) {
  for (unsigned int iInput = 0; iInput < inputs_.size(); ++iInput) {
    if (data[iInput]->dataValid()) {
      inputs_[iInput].push_back(data[iInput]);
    }
  }

  vector<vector<bool> > lMap(inputs_.size(), vector<bool>(nOutputs()));

  HGCalTriggerCellSAShrPtrCollection lInputs(inputs_.size(), std::make_shared<HGCalTriggerCell>());

  for (unsigned int iInput = 0; iInput < inputs_.size(); ++iInput) {
    unsigned int lAddr = addr_[iInput][0];
    if (lAddr < inputs_[iInput].size()) {
      lInputs[iInput] = inputs_[iInput][lAddr];
      lMap[iInput][lInputs[iInput]->sortKey()] = true;
    }
  }

  for (unsigned int iInput = 0; iInput < inputs_.size(); ++iInput) {
    vector<unsigned int>& toRotate = addr_[iInput];
    std::rotate(toRotate.begin(), toRotate.begin() + 1, toRotate.end());
  }

  HGCalTriggerCellSAShrPtrCollection lOutputs(nOutputs(), std::make_shared<HGCalTriggerCell>());

  unsigned int nOutputs = 0;
  for (unsigned int iOutput = 0; iOutput < lOutputs.size(); ++iOutput) {
    for (unsigned int iInput = 0; iInput < inputs_.size(); ++iInput) {
      if (lMap[iInput][iOutput]) {
        lOutputs[iOutput] = lInputs[iInput];
        addr_[iInput].back() += this->nInterleaving();
        nOutputs++;
        break;
      }
    }
  }

  return lOutputs;
}
