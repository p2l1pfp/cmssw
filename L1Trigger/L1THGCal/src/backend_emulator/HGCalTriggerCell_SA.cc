#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTriggerCell_SA.h"

using namespace l1thgcfirmware;

bool HGCalTriggerCell::operator==(const HGCalTriggerCell& rhs) const {
  bool isEqual = (this->index() == rhs.index());
  return isEqual;
}

bool HGCalTriggerCell::operator==(const std::shared_ptr<HGCalTriggerCell>& rhs) const {
  bool isEqual = (this->index() == rhs->index());
  return isEqual;
}
