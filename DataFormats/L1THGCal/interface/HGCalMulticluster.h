#ifndef DataFormats_L1Trigger_HGCalMulticluster_h
#define DataFormats_L1Trigger_HGCalMulticluster_h

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1THGCal/interface/HGCalClusterT.h"
#include "DataFormats/L1THGCal/interface/HGCalCluster.h"
#include <boost/iterator/transform_iterator.hpp>
#include <functional>
#include <bitset>

namespace l1t {

  class HGCalMulticluster : public HGCalClusterT<l1t::HGCalCluster> {
  public:
    HGCalMulticluster() : hOverEValid_(false) {}
    HGCalMulticluster(const LorentzVector p4, int pt = 0, int eta = 0, int phi = 0);

    HGCalMulticluster(const edm::Ptr<l1t::HGCalCluster>& tc, float fraction = 1);

    ~HGCalMulticluster() override;

    float hOverE() const {
      // --- this below would be faster when reading old objects, as HoE will only be computed once,
      // --- but it may not be allowed by CMS rules because of the const_cast
      // --- and could potentially cause a data race
      // if (!hOverEValid_) (const_cast<HGCalMulticluster*>(this))->saveHOverE();
      // --- this below is safe in any case
      return hOverEValid_ ? hOverE_ : l1t::HGCalClusterT<l1t::HGCalCluster>::hOverE();
    }

    void saveHOverE() {
      hOverE_ = l1t::HGCalClusterT<l1t::HGCalCluster>::hOverE();
      hOverEValid_ = true;
    }

    enum EnergyInterpretation { EM = 0, EM_CORE = 1, H_EARLY = 2 };

    void saveEnergyInterpretation(const HGCalMulticluster::EnergyInterpretation eInt, double energy);

    double iEnergy(const HGCalMulticluster::EnergyInterpretation eInt) const {
      return energy() * interpretationFraction(eInt);
    }

    double iPt(const HGCalMulticluster::EnergyInterpretation eInt) const { return pt() * interpretationFraction(eInt); }

    math::XYZTLorentzVector iP4(const HGCalMulticluster::EnergyInterpretation eInt) const {
      return p4() * interpretationFraction(eInt);
    }

    math::PtEtaPhiMLorentzVector iPolarP4(const HGCalMulticluster::EnergyInterpretation eInt) const {
      return math::PtEtaPhiMLorentzVector(pt() * interpretationFraction(eInt), eta(), phi(), 0.);
    }

  private:
    template <typename Iter>
    struct KeyGetter : std::unary_function<typename Iter::value_type, typename Iter::value_type::first_type> {
      const typename Iter::value_type::first_type& operator()(const typename Iter::value_type& p) const {
        return p.first;
      }
    };

    template <typename Iter>
    boost::transform_iterator<KeyGetter<Iter>, Iter> key_iterator(Iter itr) const {
      return boost::make_transform_iterator<KeyGetter<Iter>, Iter>(itr, KeyGetter<Iter>());
    }

  public:
    typedef boost::transform_iterator<KeyGetter<std::map<EnergyInterpretation, double>::const_iterator>,
                                      std::map<EnergyInterpretation, double>::const_iterator>
        EnergyInterpretation_const_iterator;

    std::pair<EnergyInterpretation_const_iterator, EnergyInterpretation_const_iterator> energyInterpretations() const {
      return std::make_pair(key_iterator(energyInterpretationFractions_.cbegin()),
                            key_iterator(energyInterpretationFractions_.cend()));
    }

    EnergyInterpretation_const_iterator interpretations_begin() const {
      return key_iterator(energyInterpretationFractions_.cbegin());
    }

    EnergyInterpretation_const_iterator interpretations_end() const {
      return key_iterator(energyInterpretationFractions_.cend());
    }

    size_type interpretations_size() const { return energyInterpretationFractions_.size(); }

  private:
    double interpretationFraction(const HGCalMulticluster::EnergyInterpretation eInt) const;

    float hOverE_;
    bool hOverEValid_;
    std::map<EnergyInterpretation, double> energyInterpretationFractions_;

  public:
    // Types for firmware representation of cluster data sent to L1T
    // These are duplicated with those in interface/backend_emulator/HGCalCluster_SA.h
    // Can we define in one place?
    static constexpr int wordLength = 64;
    static constexpr int nWordsPerCluster = 4;
    typedef uint64_t ClusterWord;
    typedef std::array<ClusterWord, nWordsPerCluster> ClusterWords;
    void setHwData( ClusterWords clusterWords ) { hwData_ = clusterWords; }
    ClusterWords getHwData() const { return hwData_; }
    void setHwSector( unsigned int sector ) { hwSector_ = sector; }
    unsigned int getHwSector() const { return hwSector_; }
    void setHwZSide( int zside ) { hwZSide_ = zside; }
    int getHwZSide() const { return hwZSide_; }

    // Types for firmware representation of cluster sum data, input to cluster properties step
    static constexpr int clusterSumWordLength = 64;
    static constexpr int nWordsPerClusterSum = 8;
    typedef uint64_t ClusterSumWord;
    typedef std::array<ClusterSumWord, nWordsPerClusterSum> ClusterSumWords;
    void setHwClusterSumData( ClusterSumWords clusterSumWords ) { hwClusterSumData_ = clusterSumWords; }
    ClusterSumWords getHwClusterSumData() const { return hwClusterSumData_; }

  private:
    ClusterWords hwData_;
    ClusterSumWords hwClusterSumData_;
    // Store sector number for pattern files
    // TODO check if detID of emulated cluster always gives same sector
    // Emulated cluster takes detID of first TC added to it, which is not necessarily a TC in the seed histo bin
    unsigned int hwSector_;
    // Safe to take zside from the detID of the emualted cluster
    // But keep consistent with sector for now
    int hwZSide_;
  };

  typedef BXVector<HGCalMulticluster> HGCalMulticlusterBxCollection;

}  // namespace l1t

#endif
