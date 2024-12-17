#ifndef L1TRIGGER_PHASE2L1PARTICLEFLOWS_NNVtx_H
#define L1TRIGGER_PHASE2L1PARTICLEFLOWS_NNVtx_H

#include <string>
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/puppi/linpuppi_ref.h"

using namespace l1ct;


class NNVtx {
public:
  NNVtx(tensorflow::Session* AssociationSesh,
        const double AssociationThreshold,
        const std::vector<double>& AssociationNetworkZ0binning,
        const std::vector<double>& AssociationNetworkEtaBounds,
        const std::vector<double>& AssociationNetworkZ0ResBins);
  ~NNVtx(){};

  template <typename T>
  bool TTTrackNetworkSelector(T& t, const l1ct::PVObjEmu& v);

private:
  tensorflow::Session* AssociationSesh_;
  double AssociationThreshold_;
  std::vector<double> z0_binning_;
  std::vector<double> eta_bins_;
  std::vector<double> res_bins_; 
 };
 #endif


