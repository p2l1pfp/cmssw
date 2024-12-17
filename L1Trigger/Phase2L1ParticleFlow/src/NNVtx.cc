#include "L1Trigger/Phase2L1ParticleFlow/interface/NNVtx.h"

NNVtx::NNVtx(tensorflow::Session* AssociationSesh,
            const double AssociationThreshold,
            const std::vector<double>& AssociationNetworkZ0binning,
            const std::vector<double>& AssociationNetworkEtaBounds,
            const std::vector<double>& AssociationNetworkZ0ResBins)
        : AssociationSesh_(AssociationSesh),
          AssociationThreshold_(AssociationThreshold),
          z0_binning_(AssociationNetworkZ0binning),
          eta_bins_(AssociationNetworkEtaBounds),
          res_bins_(AssociationNetworkZ0ResBins) {}


template <typename T>
bool NNVtx::TTTrackNetworkSelector(T& t, const l1ct::PVObjEmu& v){
      tensorflow::Tensor inputAssoc(tensorflow::DT_FLOAT, {1, 4});
      std::vector<tensorflow::Tensor> outputAssoc;

      TTTrack_TrackWord::tanl_t etaEmulationBits = t.etaEmulationBits;
      ap_fixed<16, 3> etaEmulation;
      etaEmulation.V = (etaEmulationBits.range());

      auto lower = std::lower_bound(eta_bins_.begin(), eta_bins_.end(), etaEmulation.to_double());

      int resbin = std::distance(eta_bins_.begin(), lower);
      float binWidth = z0_binning_[2];
      // Calculate integer dZ from track z0 and vertex z0 (use floating point version and convert internally allowing use of both emulator and simulator vertex and track)
      float dZ =
          abs(floor(((t.Z0 + z0_binning_[1]) / (binWidth))) - floor(((v.Z0 + z0_binning_[1]) / (binWidth))));

      // The following constants <22, 9> are defined by the quantisation of the Neural Network

      ap_ufixed<22, 9> ptEmulation_rescale = t.hwPt;
      ap_ufixed<22, 9> resBinEmulation_rescale = res_bins_[resbin];
      ap_ufixed<22, 9> MVAEmulation_rescale = t.MVAQualityBits;
      ap_ufixed<22, 9> dZEmulation_rescale = dZ;

      inputAssoc.tensor<float, 2>()(0, 0) = ptEmulation_rescale.to_double();
      inputAssoc.tensor<float, 2>()(0, 1) = MVAEmulation_rescale.to_double();
      inputAssoc.tensor<float, 2>()(0, 2) = resBinEmulation_rescale.to_double() / 16.0;
      inputAssoc.tensor<float, 2>()(0, 3) = dZEmulation_rescale.to_double();

      // Run Association Network:
      tensorflow::run(AssociationSesh_, {{"assoc:0", inputAssoc}}, {"Identity:0"}, &outputAssoc);

      double NNOutput = (double)outputAssoc[0].tensor<float, 2>()(0, 0);
      double NNOutput_exp = 1.0 / (1.0 + exp(-1.0 * (NNOutput)));
      
      return NNOutput_exp >= AssociationThreshold_;
    }

template bool NNVtx::TTTrackNetworkSelector<const l1ct::TkObjEmu>(const l1ct::TkObjEmu&, const l1ct::PVObjEmu&);
template bool NNVtx::TTTrackNetworkSelector<const l1ct::PFChargedObjEmu>(const l1ct::PFChargedObjEmu&, const l1ct::PVObjEmu&);
