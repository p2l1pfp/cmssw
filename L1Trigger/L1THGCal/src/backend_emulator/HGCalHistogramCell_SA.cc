#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistogramCell_SA.h"

using namespace l1thgcfirmware;

const HGCalHistogramCell& HGCalHistogramCell::operator+=(const HGCalHistogramCell& hc) {
  S_ += hc.S();
  X_ += hc.X();
  Y_ += hc.Y();
  N_ += hc.N();

  return *this;
}

const HGCalHistogramCell HGCalHistogramCell::operator/(const unsigned int factor) const {
  HGCalHistogramCell hc(*this);

  if (factor != 0) {
    hc.S_ /= factor;
    hc.X_ = 0;
    hc.Y_ = 0;
    hc.N_ = 0;
  }
  return hc;
}

const HGCalHistogramCell HGCalHistogramCell::operator+(const HGCalHistogramCell& hc) const {
  HGCalHistogramCell sum(*this);
  sum += hc;
  return sum;
}
