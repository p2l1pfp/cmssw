#ifndef L1TkTauEtComparator_HH
#define L1TkTauEtComparator_HH
#include "DataFormats/L1TCorrelator/interface/TkTau.h"
#include "DataFormats/L1TCorrelator/interface/TkTauFwd.h"

namespace L1TkTau{
  class EtComparator {
  public:
    bool operator()(const l1t::TkTau& a, const l1t::TkTau& b) const {
      double et_a = a.et();
      double et_b = b.et();
      return et_a > et_b;
    }
  };

  class PtComparator {
  public:
    bool operator()(const l1t::TkTau& a, const l1t::TkTau& b) const {
      double et_a = a.pt();
      double et_b = b.pt();
      return et_a > et_b;
    }
  };

}
#endif

