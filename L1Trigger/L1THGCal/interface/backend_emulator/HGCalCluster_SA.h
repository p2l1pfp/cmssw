#ifndef L1Trigger_L1THGCal_HGCalCluster_SA_h
#define L1Trigger_L1THGCal_HGCalCluster_SA_h

#include "HGCalTriggerCell_SA.h"

#include <vector>
#include <memory>
#include <iostream>
namespace l1thgcfirmware {

  class HGCalCluster {
  public:
    HGCalCluster(unsigned int clock, unsigned int index, bool frameValid, bool dataValid)
        : clock_(clock),
          index_(index),
          layerBits_(0),
          frameValid_(frameValid),
          dataValid_(dataValid),
          n_tc_(0),
          e_(0),
          e_em_(0),
          e_em_core_(0),
          e_h_early_(0),
          w_(0),
          n_tc_w_(0),
          w2_(0),
          wz_(0),
          weta_(0),
          wphi_(0),
          wroz_(0),
          wz2_(0),
          weta2_(0),
          wphi2_(0),
          wroz2_(0),
          layerbits_(0),
          sat_tc_(false),
          shapeq_(1),
          sigma_e_quotient_(0),
          sigma_e_fraction_(0),
          mean_z_quotient_(0),
          mean_z_fraction_(0),
          mean_phi_quotient_(0),
          mean_phi_fraction_(0),
          mean_eta_quotient_(0),
          mean_eta_fraction_(0),
          mean_roz_quotient_(0),
          mean_roz_fraction_(0),
          sigma_z_quotient_(0),
          sigma_z_fraction_(0),
          sigma_phi_quotient_(0),
          sigma_phi_fraction_(0),
          sigma_eta_quotient_(0),
          sigma_eta_fraction_(0),
          sigma_roz_quotient_(0),
          sigma_roz_fraction_(0),
          first_layer_(0),
          last_layer_(0),
          shower_len_(0),
          core_shower_len_(0),
          e_em_over_e_quotient_(0),
          e_em_over_e_fraction_(0),
          e_em_core_over_e_em_quotient_(0),
          e_em_core_over_e_em_fraction_(0),
          e_h_early_over_e_quotient_(0),
          e_h_early_over_e_fraction_(0) {}

    ~HGCalCluster() {}

    std::pair<unsigned int, unsigned int> sigma_energy(unsigned int N_TC_W,
                                                       unsigned long int Sum_W2,
                                                       unsigned int Sum_W);
    std::pair<unsigned int, unsigned int> mean_coordinate(unsigned int Sum_Wc, unsigned int Sum_W);
    std::pair<unsigned int, unsigned int> sigma_coordinate(unsigned int Sum_W,
                                                           unsigned long int Sum_Wc2,
                                                           unsigned int Sum_Wc);
    std::pair<unsigned int, unsigned int> energy_ratio(unsigned int E_N, unsigned int E_D);
    std::vector<int> showerLengthProperties(unsigned long int layerBits);

    // Setters
    void setClock(const unsigned int clock) { clock_ = clock; }
    void setIndex(const unsigned int index) { index_ = index; }
    void setDataValid(const bool dataValid) { dataValid_ = dataValid; }

    void set_n_tc(unsigned int n_tc) { n_tc_ = n_tc; }
    void set_e(unsigned int e) { e_ = e; }
    void set_e_em(unsigned int e_em) { e_em_ = e_em; }
    void set_e_em_core(unsigned int e_em_core) { e_em_core_ = e_em_core; }
    void set_e_h_early(unsigned int e_h_early) { e_h_early_ = e_h_early; }
    void set_w(unsigned int w) { w_ = w; }
    void set_n_tc_w(unsigned int n_tc_w) { n_tc_w_ = n_tc_w; }
    void set_w2(unsigned long int w2) { w2_ = w2; }
    void set_wz(unsigned int wz) { wz_ = wz; }
    void set_weta(unsigned int weta) { weta_ = weta; }
    void set_wphi(unsigned int wphi) { wphi_ = wphi; }
    void set_wroz(unsigned int wroz) { wroz_ = wroz; }
    void set_wz2(unsigned long int wz2) { wz2_ = wz2; }
    void set_weta2(unsigned long int weta2) { weta2_ = weta2; }
    void set_wphi2(unsigned long int wphi2) { wphi2_ = wphi2; }
    void set_wroz2(unsigned long int wroz2) { wroz2_ = wroz2; }
    void set_layerbits(unsigned long int layerbits) { layerbits_ = layerbits; }
    void set_sat_tc(bool sat_tc) { sat_tc_ = sat_tc; }
    void set_shapeq(unsigned int shapeq) { shapeq_ = shapeq; }

    void set_sigma_e_quotient(unsigned long int sigma_e_quotient) { sigma_e_quotient_ = sigma_e_quotient; }
    void set_sigma_e_fraction(unsigned long int sigma_e_fraction) { sigma_e_fraction_ = sigma_e_fraction; }
    void set_mean_z_quotient(unsigned long int mean_z_quotient) { mean_z_quotient_ = mean_z_quotient; }
    void set_mean_z_fraction(unsigned long int mean_z_fraction) { mean_z_fraction_ = mean_z_fraction; }
    void set_mean_phi_quotient(unsigned long int mean_phi_quotient) { mean_phi_quotient_ = mean_phi_quotient; }
    void set_mean_phi_fraction(unsigned long int mean_phi_fraction) { mean_phi_fraction_ = mean_phi_fraction; }
    void set_mean_eta_quotient(unsigned long int mean_eta_quotient) { mean_eta_quotient_ = mean_eta_quotient; }
    void set_mean_eta_fraction(unsigned long int mean_eta_fraction) { mean_eta_fraction_ = mean_eta_fraction; }
    void set_mean_roz_quotient(unsigned long int mean_roz_quotient) { mean_roz_quotient_ = mean_roz_quotient; }
    void set_mean_roz_fraction(unsigned long int mean_roz_fraction) { mean_roz_fraction_ = mean_roz_fraction; }
    void set_sigma_z_quotient(unsigned long int sigma_z_quotient) { sigma_z_quotient_ = sigma_z_quotient; }
    void set_sigma_z_fraction(unsigned long int sigma_z_fraction) { sigma_z_fraction_ = sigma_z_fraction; }
    void set_sigma_phi_quotient(unsigned long int sigma_phi_quotient) { sigma_phi_quotient_ = sigma_phi_quotient; }
    void set_sigma_phi_fraction(unsigned long int sigma_phi_fraction) { sigma_phi_fraction_ = sigma_phi_fraction; }
    void set_sigma_eta_quotient(unsigned long int sigma_eta_quotient) { sigma_eta_quotient_ = sigma_eta_quotient; }
    void set_sigma_eta_fraction(unsigned long int sigma_eta_fraction) { sigma_eta_fraction_ = sigma_eta_fraction; }
    void set_sigma_roz_quotient(unsigned long int sigma_roz_quotient) { sigma_roz_quotient_ = sigma_roz_quotient; }
    void set_sigma_roz_fraction(unsigned long int sigma_roz_fraction) { sigma_roz_fraction_ = sigma_roz_fraction; }
    void set_firstLayer(unsigned long int FirstLayer) { first_layer_ = FirstLayer; }
    void set_lastLayer(unsigned long int LastLayer) { last_layer_ = LastLayer; }
    void set_showerLen(unsigned long int ShowerLen) { shower_len_ = ShowerLen; }
    void set_coreShowerLen(unsigned long int CoreShowerLen) { core_shower_len_ = CoreShowerLen; }
    void set_e_em_over_e_quotient(unsigned long int quotient) { e_em_over_e_quotient_ = quotient; }
    void set_e_em_over_e_fraction(unsigned long int fraction) { e_em_over_e_fraction_ = fraction; }
    void set_e_em_core_over_e_em_quotient(unsigned long int quotient) { e_em_core_over_e_em_quotient_ = quotient; }
    void set_e_em_core_over_e_em_fraction(unsigned long int fraction) { e_em_core_over_e_em_fraction_ = fraction; }
    void set_e_h_early_over_e_quotient(unsigned long int quotient) { e_h_early_over_e_quotient_ = quotient; }
    void set_e_h_early_over_e_fraction(unsigned long int fraction) { e_h_early_over_e_fraction_ = fraction; }

    void add_constituent(HGCalTriggerCellSAShrPtr constituent) { constituents_.emplace_back(constituent); }

    // Getters
    unsigned int clock() const { return clock_; }
    unsigned int index() const { return index_; }
    bool frameValid() const { return frameValid_; }
    bool dataValid() const { return dataValid_; }

    unsigned int n_tc() const { return n_tc_; }
    unsigned int e() const { return e_; }
    unsigned int e_em() const { return e_em_; }
    unsigned int e_em_core() const { return e_em_core_; }
    unsigned int e_h_early() const { return e_h_early_; }
    unsigned int w() const { return w_; }
    unsigned int n_tc_w() const { return n_tc_w_; }
    unsigned long int w2() const { return w2_; }
    unsigned int wz() const { return wz_; }
    unsigned int weta() const { return weta_; }
    unsigned int wphi() const { return wphi_; }
    unsigned int wroz() const { return wroz_; }
    unsigned long int wz2() const { return wz2_; }
    unsigned long int weta2() const { return weta2_; }
    unsigned long int wphi2() const { return wphi2_; }
    unsigned long int wroz2() const { return wroz2_; }
    unsigned long int layerbits() const { return layerbits_; }
    bool sat_tc() const { return sat_tc_; }
    unsigned int shapeq() const { return shapeq_; }

    unsigned long int sigma_e_quotient() const { return sigma_e_quotient_; }
    unsigned long int sigma_e_fraction() const { return sigma_e_fraction_; }
    unsigned long int mean_z_quotient() const { return mean_z_quotient_; }
    unsigned long int mean_z_fraction() const { return mean_z_fraction_; }
    unsigned long int mean_phi_quotient() const { return mean_phi_quotient_; }
    unsigned long int mean_phi_fraction() const { return mean_phi_fraction_; }
    unsigned long int mean_eta_quotient() const { return mean_eta_quotient_; }
    unsigned long int mean_eta_fraction() const { return mean_eta_fraction_; }
    unsigned long int mean_roz_quotient() const { return mean_roz_quotient_; }
    unsigned long int mean_roz_fraction() const { return mean_roz_fraction_; }
    unsigned long int sigma_z_quotient() const { return sigma_z_quotient_; }
    unsigned long int sigma_z_fraction() const { return sigma_z_fraction_; }
    unsigned long int sigma_phi_quotient() const { return sigma_phi_quotient_; }
    unsigned long int sigma_phi_fraction() const { return sigma_phi_fraction_; }
    unsigned long int sigma_eta_quotient() const { return sigma_eta_quotient_; }
    unsigned long int sigma_eta_fraction() const { return sigma_eta_fraction_; }
    unsigned long int sigma_roz_quotient() const { return sigma_roz_quotient_; }
    unsigned long int sigma_roz_fraction() const { return sigma_roz_fraction_; }
    unsigned long int firstLayer() const { return first_layer_; }
    unsigned long int lastLayer() const { return last_layer_; }
    unsigned long int showerLen() const { return shower_len_; }
    unsigned long int coreShowerLen() const { return core_shower_len_; }
    unsigned long int e_em_over_e_quotient() const { return e_em_over_e_quotient_; }
    unsigned long int e_em_over_e_fraction() const { return e_em_over_e_fraction_; }
    unsigned long int e_em_core_over_e_em_quotient() const { return e_em_core_over_e_em_quotient_; }
    unsigned long int e_em_core_over_e_em_fraction() const { return e_em_core_over_e_em_fraction_; }
    unsigned long int e_h_early_over_e_quotient() const { return e_h_early_over_e_quotient_; }
    unsigned long int e_h_early_over_e_fraction() const { return e_h_early_over_e_fraction_; }

    HGCalTriggerCellSAShrPtrCollection& constituents() { return constituents_; }

    // Operators
    const HGCalCluster& operator+=(HGCalCluster& hc);

  private:
    unsigned int clock_;
    unsigned int index_;
    unsigned int layerBits_;
    bool frameValid_;
    bool dataValid_;

    unsigned int n_tc_;
    unsigned int e_;
    unsigned int e_em_;
    unsigned int e_em_core_;
    unsigned int e_h_early_;
    unsigned int w_;
    unsigned int n_tc_w_;
    unsigned long int w2_;
    unsigned int wz_;
    unsigned int weta_;
    unsigned int wphi_;
    unsigned int wroz_;
    unsigned long int wz2_;
    unsigned long int weta2_;
    unsigned long int wphi2_;
    unsigned long int wroz2_;
    unsigned long int layerbits_;
    bool sat_tc_;
    unsigned int shapeq_;
    unsigned long int sigma_e_quotient_;
    unsigned long int sigma_e_fraction_;
    unsigned long int mean_z_quotient_;
    unsigned long int mean_z_fraction_;
    unsigned long int mean_phi_quotient_;
    unsigned long int mean_phi_fraction_;
    unsigned long int mean_eta_quotient_;
    unsigned long int mean_eta_fraction_;
    unsigned long int mean_roz_quotient_;
    unsigned long int mean_roz_fraction_;
    unsigned long int sigma_z_quotient_;
    unsigned long int sigma_z_fraction_;
    unsigned long int sigma_phi_quotient_;
    unsigned long int sigma_phi_fraction_;
    unsigned long int sigma_eta_quotient_;
    unsigned long int sigma_eta_fraction_;
    unsigned long int sigma_roz_quotient_;
    unsigned long int sigma_roz_fraction_;
    unsigned long int first_layer_;
    unsigned long int last_layer_;
    unsigned long int shower_len_;
    unsigned long int core_shower_len_;
    unsigned long int e_em_over_e_quotient_;
    unsigned long int e_em_over_e_fraction_;
    unsigned long int e_em_core_over_e_em_quotient_;
    unsigned long int e_em_core_over_e_em_fraction_;
    unsigned long int e_h_early_over_e_quotient_;
    unsigned long int e_h_early_over_e_fraction_;

    // Extra variables, not available in firmware
    // Perhaps move to separate "extra" class?
    HGCalTriggerCellSAShrPtrCollection constituents_;
  };

  typedef std::vector<HGCalCluster> HGCalClusterSACollection;
  typedef std::unique_ptr<HGCalCluster> HGCalClusterSAPtr;
  typedef std::vector<HGCalClusterSAPtr> HGCalClusterSAPtrCollection;
  typedef std::vector<HGCalClusterSAPtrCollection> HGCalClusterSAPtrCollections;

}  // namespace l1thgcfirmware

#endif
