#ifndef DataFormats_L1Trigger_HGCalClusterT_h
#define DataFormats_L1Trigger_HGCalClusterT_h

/* CMSSW */
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "DataFormats/L1THGCal/interface/HGCalTriggerCell.h"
#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerDetId.h"
#include "DataFormats/ForwardDetId/interface/HFNoseTriggerDetId.h"

/* ROOT */
#include "Math/Vector3D.h"

#include <unordered_map>

namespace l1t {
  template <class C>
  class HGCalClusterT : public L1Candidate {
  public:
    typedef typename std::unordered_map<uint32_t, edm::Ptr<C>>::const_iterator const_iterator;

  public:
    HGCalClusterT() {}
    HGCalClusterT(const LorentzVector p4, int pt = 0, int eta = 0, int phi = 0)
        : L1Candidate(p4, pt, eta, phi),
          valid_(true),
          detId_(0),
          centre_(0, 0, 0),
          centreProj_(0., 0., 0.),
          mipPt_(0),
          seedMipPt_(0),
          sumPt_(0) {}

    HGCalClusterT(const edm::Ptr<C>& c, float fraction = 1.)
        : valid_(true),
          detId_(c->detId()),
          centre_(0., 0., 0.),
          centreProj_(0., 0., 0.),
          mipPt_(0.),
          seedMipPt_(0.),
          sumPt_(0.) {
      addConstituent(c, true, fraction);
    }

    ~HGCalClusterT() override{};

    const std::unordered_map<uint32_t, edm::Ptr<C>>& constituents() const { return constituents_; }
    const_iterator constituents_begin() const { return constituents_.begin(); }
    const_iterator constituents_end() const { return constituents_.end(); }
    unsigned size() const { return constituents_.size(); }

    void addConstituent(const edm::Ptr<C>& c, bool updateCentre = true, float fraction = 1.) {
      double cMipt = c->mipPt() * fraction;

      if (constituents_.empty()) {
        detId_ = DetId(c->detId());
        seedMipPt_ = cMipt;
        /* if the centre will not be dynamically calculated
             the seed centre is considere as cluster centre */
        if (!updateCentre) {
          centre_ = c->position();
        }
      }
      updateP4AndPosition(c, updateCentre, fraction);

      constituents_.emplace(c->detId(), c);
      constituentsFraction_.emplace(c->detId(), fraction);
    }

    void removeConstituent(const edm::Ptr<C>& c, bool updateCentre = true) {
      /* remove the pointer to c from the std::vector */
      double fraction = 0;
      const auto& constituent_itr = constituents_.find(c->detId());
      const auto& fraction_itr = constituentsFraction_.find(c->detId());
      if (constituent_itr != constituents_.end()) {
        // remove constituent and get its fraction in the cluster
        fraction = fraction_itr->second;
        constituents_.erase(constituent_itr);
        constituentsFraction_.erase(fraction_itr);

        updateP4AndPosition(c, updateCentre, -fraction);
      }
    }

    bool valid() const { return valid_; }
    void setValid(bool valid) { valid_ = valid; }

    double mipPt() const { return mipPt_; }
    double seedMipPt() const { return seedMipPt_; }
    uint32_t detId() const { return detId_.rawId(); }
    void setDetId(uint32_t id) { detId_ = id; }
    void setPt(double pt) { setP4(math::PtEtaPhiMLorentzVector(pt, eta(), phi(), mass())); }
    double sumPt() const { return sumPt_; }
    /* distance in 'cm' */
    double distance(const l1t::HGCalTriggerCell& tc) const { return (tc.position() - centre_).mag(); }

    const GlobalPoint& position() const { return centre_; }
    const GlobalPoint& centre() const { return centre_; }
    const GlobalPoint& centreProj() const { return centreProj_; }

    double hOverE() const {
      double pt_em = 0.;
      double pt_had = 0.;
      double hOe = 0.;

      for (const auto& id_constituent : constituents()) {
        DetId id(id_constituent.first);
        auto id_fraction = constituentsFraction_.find(id_constituent.first);
        double fraction = (id_fraction != constituentsFraction_.end() ? id_fraction->second : 1.);
        if ((id.det() == DetId::HGCalEE) ||
            (id.det() == DetId::HGCalTrigger &&
             HGCalTriggerDetId(id).subdet() == HGCalTriggerSubdetector::HGCalEETrigger) ||
            (id.det() == DetId::Forward && id.subdetId() == ForwardSubdetector::HFNose && HFNoseDetId(id).isEE()) ||
            (id.det() == DetId::HGCalTrigger &&
             HGCalTriggerDetId(id).subdet() == HGCalTriggerSubdetector::HFNoseTrigger &&
             HFNoseTriggerDetId(id).isEE())) {
          pt_em += id_constituent.second->pt() * fraction;
        } else if ((id.det() == DetId::HGCalHSi) || (id.det() == DetId::HGCalHSc) ||
                   (id.det() == DetId::HGCalTrigger &&
                    HGCalTriggerDetId(id).subdet() == HGCalTriggerSubdetector::HGCalHSiTrigger) ||
                   (id.det() == DetId::Forward && id.subdetId() == ForwardSubdetector::HFNose &&
                    HFNoseDetId(id).isHE()) ||
                   (id.det() == DetId::HGCalTrigger &&
                    HGCalTriggerDetId(id).subdet() == HGCalTriggerSubdetector::HFNoseTrigger &&
                    HFNoseTriggerDetId(id).isHSilicon())) {
          pt_had += id_constituent.second->pt() * fraction;
        }
      }
      if (pt_em > 0)
        hOe = pt_had / pt_em;
      else
        hOe = -1.;
      return hOe;
    }

    uint32_t subdetId() const { return detId_.subdetId(); }

    //shower shape

    int showerLength() const { return showerLength_; }
    int coreShowerLength() const { return coreShowerLength_; }
    int firstLayer() const { return firstLayer_; }
    int maxLayer() const { return maxLayer_; }
    float eMax() const { return eMax_; }
    float sigmaEtaEtaMax() const { return sigmaEtaEtaMax_; }
    float sigmaPhiPhiMax() const { return sigmaPhiPhiMax_; }
    float sigmaEtaEtaTot() const { return sigmaEtaEtaTot_; }
    float sigmaPhiPhiTot() const { return sigmaPhiPhiTot_; }
    float sigmaZZ() const { return sigmaZZ_; }
    float sigmaRRTot() const { return sigmaRRTot_; }
    float sigmaRRMax() const { return sigmaRRMax_; }
    float sigmaRRMean() const { return sigmaRRMean_; }
    float zBarycenter() const { return zBarycenter_; }
    float layer10percent() const { return layer10percent_; }
    float layer50percent() const { return layer50percent_; }
    float layer90percent() const { return layer90percent_; }
    float triggerCells67percent() const { return triggerCells67percent_; }
    float triggerCells90percent() const { return triggerCells90percent_; }

    void showerLength(int showerLength) { showerLength_ = showerLength; }
    void coreShowerLength(int coreShowerLength) { coreShowerLength_ = coreShowerLength; }
    void firstLayer(int firstLayer) { firstLayer_ = firstLayer; }
    void maxLayer(int maxLayer) { maxLayer_ = maxLayer; }
    void eMax(float eMax) { eMax_ = eMax; }
    void sigmaEtaEtaMax(float sigmaEtaEtaMax) { sigmaEtaEtaMax_ = sigmaEtaEtaMax; }
    void sigmaEtaEtaTot(float sigmaEtaEtaTot) { sigmaEtaEtaTot_ = sigmaEtaEtaTot; }
    void sigmaPhiPhiMax(float sigmaPhiPhiMax) { sigmaPhiPhiMax_ = sigmaPhiPhiMax; }
    void sigmaPhiPhiTot(float sigmaPhiPhiTot) { sigmaPhiPhiTot_ = sigmaPhiPhiTot; }
    void sigmaRRMax(float sigmaRRMax) { sigmaRRMax_ = sigmaRRMax; }
    void sigmaRRTot(float sigmaRRTot) { sigmaRRTot_ = sigmaRRTot; }
    void sigmaRRMean(float sigmaRRMean) { sigmaRRMean_ = sigmaRRMean; }
    void sigmaZZ(float sigmaZZ) { sigmaZZ_ = sigmaZZ; }
    void zBarycenter(float zBarycenter) { zBarycenter_ = zBarycenter; }
    void layer10percent(float layer10percent) { layer10percent_ = layer10percent; }
    void layer50percent(float layer50percent) { layer50percent_ = layer50percent; }
    void layer90percent(float layer90percent) { layer90percent_ = layer90percent; }
    void triggerCells67percent(float triggerCells67percent) { triggerCells67percent_ = triggerCells67percent; }
    void triggerCells90percent(float triggerCells90percent) { triggerCells90percent_ = triggerCells90percent; }

    // Firmware-specific cluster properties
    unsigned long int hw_sigma_e_quotient() const { return hw_sigma_e_quotient_; }
    unsigned long int hw_sigma_e_fraction() const { return hw_sigma_e_fraction_; }
    unsigned long int hw_mean_z_quotient() const { return hw_mean_z_quotient_; }
    unsigned long int hw_mean_z_fraction() const { return hw_mean_z_fraction_; }
    unsigned long int hw_mean_phi_quotient() const { return hw_mean_phi_quotient_; }
    unsigned long int hw_mean_phi_fraction() const { return hw_mean_phi_fraction_; }
    unsigned long int hw_mean_eta_quotient() const { return hw_mean_eta_quotient_; }
    unsigned long int hw_mean_eta_fraction() const { return hw_mean_eta_fraction_; }
    unsigned long int hw_mean_roz_quotient() const { return hw_mean_roz_quotient_; }
    unsigned long int hw_mean_roz_fraction() const { return hw_mean_roz_fraction_; }
    unsigned long int hw_sigma_z_quotient() const { return hw_sigma_z_quotient_; }
    unsigned long int hw_sigma_z_fraction() const { return hw_sigma_z_fraction_; }
    unsigned long int hw_sigma_phi_quotient() const { return hw_sigma_phi_quotient_; }
    unsigned long int hw_sigma_phi_fraction() const { return hw_sigma_phi_fraction_; }
    unsigned long int hw_sigma_eta_quotient() const { return hw_sigma_eta_quotient_; }
    unsigned long int hw_sigma_eta_fraction() const { return hw_sigma_eta_fraction_; }
    unsigned long int hw_sigma_roz_quotient() const { return hw_sigma_roz_quotient_; }
    unsigned long int hw_sigma_roz_fraction() const { return hw_sigma_roz_fraction_; }
    unsigned long int hw_e_em_over_e_quotient() const { return hw_e_em_over_e_quotient_; }
    unsigned long int hw_e_em_over_e_fraction() const { return hw_e_em_over_e_fraction_; }
    unsigned long int hw_e_em_core_over_e_em_quotient() const { return hw_e_em_core_over_e_em_quotient_; }
    unsigned long int hw_e_em_core_over_e_em_fraction() const { return hw_e_em_core_over_e_em_fraction_; }
    unsigned long int hw_e_h_early_over_e_quotient() const { return hw_e_h_early_over_e_quotient_; }
    unsigned long int hw_e_h_early_over_e_fraction() const { return hw_e_h_early_over_e_fraction_; }

    void hw_sigma_e_quotient(unsigned long int sigma_e_quotient) { hw_sigma_e_quotient_ = sigma_e_quotient; }
    void hw_sigma_e_fraction(unsigned long int sigma_e_fraction) { hw_sigma_e_fraction_ = sigma_e_fraction; }
    void hw_mean_z_quotient(unsigned long int mean_z_quotient) { hw_mean_z_quotient_ = mean_z_quotient; }
    void hw_mean_z_fraction(unsigned long int mean_z_fraction) { hw_mean_z_fraction_ = mean_z_fraction; }
    void hw_mean_phi_quotient(unsigned long int mean_phi_quotient) { hw_mean_phi_quotient_ = mean_phi_quotient; }
    void hw_mean_phi_fraction(unsigned long int mean_phi_fraction) { hw_mean_phi_fraction_ = mean_phi_fraction; }
    void hw_mean_eta_quotient(unsigned long int mean_eta_quotient) { hw_mean_eta_quotient_ = mean_eta_quotient; }
    void hw_mean_eta_fraction(unsigned long int mean_eta_fraction) { hw_mean_eta_fraction_ = mean_eta_fraction; }
    void hw_mean_roz_quotient(unsigned long int mean_roz_quotient) { hw_mean_roz_quotient_ = mean_roz_quotient; }
    void hw_mean_roz_fraction(unsigned long int mean_roz_fraction) { hw_mean_roz_fraction_ = mean_roz_fraction; }
    void hw_sigma_z_quotient(unsigned long int sigma_z_quotient) { hw_sigma_z_quotient_ = sigma_z_quotient; }
    void hw_sigma_z_fraction(unsigned long int sigma_z_fraction) { hw_sigma_z_fraction_ = sigma_z_fraction; }
    void hw_sigma_phi_quotient(unsigned long int sigma_phi_quotient) { hw_sigma_phi_quotient_ = sigma_phi_quotient; }
    void hw_sigma_phi_fraction(unsigned long int sigma_phi_fraction) { hw_sigma_phi_fraction_ = sigma_phi_fraction; }
    void hw_sigma_eta_quotient(unsigned long int sigma_eta_quotient) { hw_sigma_eta_quotient_ = sigma_eta_quotient; }
    void hw_sigma_eta_fraction(unsigned long int sigma_eta_fraction) { hw_sigma_eta_fraction_ = sigma_eta_fraction; }
    void hw_sigma_roz_quotient(unsigned long int sigma_roz_quotient) { hw_sigma_roz_quotient_ = sigma_roz_quotient; }
    void hw_sigma_roz_fraction(unsigned long int sigma_roz_fraction) { hw_sigma_roz_fraction_ = sigma_roz_fraction; }
    void hw_e_em_over_e_quotient(unsigned long int e_em_over_e_quotient) {
      hw_e_em_over_e_quotient_ = e_em_over_e_quotient;
    }
    void hw_e_em_over_e_fraction(unsigned long int e_em_over_e_fraction) {
      hw_e_em_over_e_fraction_ = e_em_over_e_fraction;
    }
    void hw_e_em_core_over_e_em_quotient(unsigned long int e_em_core_over_e_em_quotient) {
      hw_e_em_core_over_e_em_quotient_ = e_em_core_over_e_em_quotient;
    }
    void hw_e_em_core_over_e_em_fraction(unsigned long int e_em_core_over_e_em_fraction) {
      hw_e_em_core_over_e_em_fraction_ = e_em_core_over_e_em_fraction;
    }
    void hw_e_h_early_over_e_quotient(unsigned long int e_h_early_over_e_quotient) {
      hw_e_h_early_over_e_quotient_ = e_h_early_over_e_quotient;
    }
    void hw_e_h_early_over_e_fraction(unsigned long int e_h_early_over_e_fraction) {
      hw_e_h_early_over_e_fraction_ = e_h_early_over_e_fraction;
    }

    /* operators */
    bool operator<(const HGCalClusterT<C>& cl) const { return mipPt() < cl.mipPt(); }
    bool operator>(const HGCalClusterT<C>& cl) const { return cl < *this; }
    bool operator<=(const HGCalClusterT<C>& cl) const { return !(cl > *this); }
    bool operator>=(const HGCalClusterT<C>& cl) const { return !(cl < *this); }

  private:
    bool valid_ = false;
    DetId detId_;

    std::unordered_map<uint32_t, edm::Ptr<C>> constituents_;
    std::unordered_map<uint32_t, double> constituentsFraction_;

    GlobalPoint centre_;
    GlobalPoint centreProj_;  // centre projected onto the first HGCal layer

    double mipPt_ = 0.;
    double seedMipPt_ = 0.;
    double sumPt_ = 0.;

    //shower shape

    int showerLength_ = 0;
    int coreShowerLength_ = 0;
    int firstLayer_ = 0;
    int maxLayer_ = 0;
    float eMax_ = 0.;
    float sigmaEtaEtaMax_ = 0.;
    float sigmaPhiPhiMax_ = 0.;
    float sigmaRRMax_ = 0.;
    float sigmaEtaEtaTot_ = 0.;
    float sigmaPhiPhiTot_ = 0.;
    float sigmaRRTot_ = 0.;
    float sigmaRRMean_ = 0.;
    float sigmaZZ_ = 0.;
    float zBarycenter_ = 0.;
    float layer10percent_ = 0.;
    float layer50percent_ = 0.;
    float layer90percent_ = 0.;
    float triggerCells67percent_ = 0.;
    float triggerCells90percent_ = 0.;

    // firmware-specific cluster properties
    unsigned long int hw_sigma_e_quotient_ = 0;
    unsigned long int hw_sigma_e_fraction_ = 0;
    unsigned long int hw_mean_z_quotient_ = 0;
    unsigned long int hw_mean_z_fraction_ = 0;
    unsigned long int hw_mean_phi_quotient_ = 0;
    unsigned long int hw_mean_phi_fraction_ = 0;
    unsigned long int hw_mean_eta_quotient_ = 0;
    unsigned long int hw_mean_eta_fraction_ = 0;
    unsigned long int hw_mean_roz_quotient_ = 0;
    unsigned long int hw_mean_roz_fraction_ = 0;
    unsigned long int hw_sigma_z_quotient_ = 0;
    unsigned long int hw_sigma_z_fraction_ = 0;
    unsigned long int hw_sigma_phi_quotient_ = 0;
    unsigned long int hw_sigma_phi_fraction_ = 0;
    unsigned long int hw_sigma_eta_quotient_ = 0;
    unsigned long int hw_sigma_eta_fraction_ = 0;
    unsigned long int hw_sigma_roz_quotient_ = 0;
    unsigned long int hw_sigma_roz_fraction_ = 0;
    unsigned long int hw_e_em_over_e_quotient_ = 0;
    unsigned long int hw_e_em_over_e_fraction_ = 0;
    unsigned long int hw_e_em_core_over_e_em_quotient_ = 0;
    unsigned long int hw_e_em_core_over_e_em_fraction_ = 0;
    unsigned long int hw_e_h_early_over_e_quotient_ = 0;
    unsigned long int hw_e_h_early_over_e_fraction_ = 0;

    void updateP4AndPosition(const edm::Ptr<C>& c, bool updateCentre = true, float fraction = 1.) {
      double cMipt = c->mipPt() * fraction;
      double cPt = c->pt() * fraction;
      /* update cluster positions (IF requested) */
      if (updateCentre) {
        Basic3DVector<float> constituentCentre(c->position());
        Basic3DVector<float> clusterCentre(centre_);

        clusterCentre = clusterCentre * mipPt_ + constituentCentre * cMipt;
        if ((mipPt_ + cMipt) > 0) {
          clusterCentre /= (mipPt_ + cMipt);
        }
        centre_ = GlobalPoint(clusterCentre);

        if (clusterCentre.z() != 0) {
          centreProj_ = GlobalPoint(clusterCentre / std::abs(clusterCentre.z()));
        }
      }

      /* update cluster energies */
      mipPt_ += cMipt;
      sumPt_ += cPt;
      int updatedPt = hwPt() + (int)(c->hwPt() * fraction);
      setHwPt(updatedPt);

      math::PtEtaPhiMLorentzVector updatedP4(p4());
      updatedP4 += (c->p4() * fraction);
      setP4(updatedP4);
    }
  };

}  // namespace l1t

#endif
