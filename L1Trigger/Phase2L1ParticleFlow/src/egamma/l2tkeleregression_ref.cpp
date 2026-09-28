#include "L1Trigger/Phase2L1ParticleFlow/interface/egamma/l2tkeleregression_ref.h"

#ifdef CMSSW_GIT_HASH
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

l1ct::L2TkEleRegressionEmulator::L2TkEleRegressionEmulator(const edm::ParameterSet& iConfig)
    : L2TkEleRegressionEmulator(iConfig.getParameter<std::vector<double>>("eta_bins"),
                                iConfig.getParameter<std::vector<unsigned int>>("model_types"),
                                iConfig.getParameter<std::vector<std::string>>("model_paths"),
                                iConfig.getUntrackedParameter<int>("debug", 0)) {}

edm::ParameterSetDescription l1ct::L2TkEleRegressionEmulator::getParameterSetDescription() {
  edm::ParameterSetDescription desc;
  desc.add<std::vector<double>>("eta_bins", std::vector<double>{1.5});
  desc.add<std::vector<unsigned int>>("model_types", std::vector<unsigned int>{1});
  desc.add<std::vector<std::string>>(
      "model_paths",
      std::vector<std::string>{"L1Trigger/Phase2L1ParticleFlow/data/egamma/eb_v0/eta_bin0_emf_bin0/model.json"});
  desc.addUntracked<int>("debug", 0);
  return desc;
}

#endif

l1ct::L2TkEleRegressionEmulator::L2TkEleRegressionEmulator(const std::vector<double>& eta_bins,
                                                           const std::vector<unsigned int>& model_types,
                                                           const std::vector<std::string>& model_paths,
                                                           int debug) {
  if (model_types.size() != model_paths.size() || model_types.size() != eta_bins.size()) {
    throw std::invalid_argument("Size of model_types and model_paths must be the same");
  }
  // check the eta bins are in increasing order
  for (size_t i = 1; i < eta_bins.size(); ++i) {
    if (eta_bins[i] <= eta_bins[i - 1]) {
      throw std::invalid_argument("Eta bins must be in increasing order");
    }
  }
  for (unsigned int i = 0; i < model_types.size(); ++i) {
    auto type = ModelType(model_types[i]);
    if (type == ModelType::EB_v0) {
      models_.push_back(std::make_unique<Model_EB_v0>(model_paths[models_.size()], debug));
    } else if (type == ModelType::null) {
      models_.push_back(nullptr);
    } else {
      throw std::invalid_argument("Unsupported regression algorithm");
    }
    eta_bins_.push_back(eta_bins[i]);
  }
}

l1ct::L2TkEleRegressionEmulator::Model_EB_v0::Model_EB_v0(const std::string& model_path, int debug) {
#ifdef CMSSW_GIT_HASH
  auto resolvedFileName = edm::FileInPath(model_path).fullPath();
#else
  auto resolvedFileName = model_path;
#endif
  model_ = std::make_unique<conifer::BDT<bdt_feature_t, bdt_out_t, false>>(resolvedFileName);
}

l1ct::pt_t l1ct::L2TkEleRegressionEmulator::Model_EB_v0::compute_ptCorr(const EGIsoEleObjEmu& ele) const {
  bdt_feature_t scaled_ID = bdt_feature_t(ele.hwIDScore.to_float());
  bdt_feature_t scaled_cl_eta = scale(fabs(ele.hwEta.to_float()), 0., 9);
  bdt_feature_t scaled_cltk_absDphi = scale(ele.hwTkCaloDphi.to_float(), 0., 5);
  bdt_feature_t scaled_tk_chi2RPhi = scale(ele.hwTkRedChi2RPhi.to_float(), 0., 3);
  bdt_feature_t scaled_cl_pt = scale(ele.hwPt.to_float(), 0., 5);
  bdt_feature_t scaled_cl_ss = scale(ele.hwCaloShowerShape.to_float(), 0., -1);
  bdt_feature_t scaled_cltk_ptRatio = scale(ele.hwCaloTkPtRatio.to_float(), 0., 0);

#ifdef L2TKELEREG_DEBUG
  if (ele.hwPt > 0) {
    std::cout << "[EMULATOR] BDT inputs: " << std::endl;
    std::cout << "  scaled_ID: " << scaled_ID << std::endl;
    std::cout << "  scaled_cl_eta: " << scaled_cl_eta << std::endl;
    std::cout << "  scaled_cltk_absDphi: " << scaled_cltk_absDphi << std::endl;
    std::cout << "  scaled_tk_chi2RPhi: " << scaled_tk_chi2RPhi << std::endl;
    std::cout << "  scaled_cl_pt: " << scaled_cl_pt << std::endl;
    std::cout << "  scaled_cl_ss: " << scaled_cl_ss << std::endl;
    std::cout << "  scaled_cltk_ptRatio: " << scaled_cltk_ptRatio << std::endl;
  }
#endif

  // Run BDT inference
  std::vector<bdt_feature_t> inputs = {scaled_ID,
                                       scaled_cl_eta,
                                       scaled_cltk_absDphi,
                                       scaled_tk_chi2RPhi,
                                       scaled_cl_pt,
                                       scaled_cl_ss,
                                       scaled_cltk_ptRatio};

  std::vector<bdt_out_t> bdt_output = model_->decision_function(inputs);

  bdt_out_t corr_factor = bdt_out_t(bdt_output[0]);
  float corr_pt = ele.hwPt.to_float() * (1. + corr_factor.to_float());

#ifdef L2TKELEREG_DEBUG
  if (ele.hwPt > 0) {
    std::cout << "[EMULATOR] BDT output: " << 1. + corr_factor.to_float() << std::endl;
  }
#endif

  return pt_t(corr_pt);
}

void l1ct::L2TkEleRegressionEmulator::run(const std::vector<EGIsoEleObjEmu>& in_eles,
                                          std::vector<EGIsoEleObjEmu>& out_eles) const {
  out_eles.clear();
  for (const auto& ele : in_eles) {
    EGIsoEleObjEmu corrected_ele = ele;
    // find the eta bin the electron falls into and apply the corresponding regression model
    for (size_t i = 0; i < eta_bins_.size(); ++i) {
      if (std::abs(ele.floatEta()) < eta_bins_[i]) {
        if (models_[i]) {
          corrected_ele.hwPt = models_[i]->compute_ptCorr(ele);
        }
        break;
      }
    }
    out_eles.push_back(corrected_ele);
  }
}
