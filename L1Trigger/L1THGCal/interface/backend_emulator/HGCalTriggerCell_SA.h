#ifndef L1Trigger_L1THGCal_HGCalTriggerCell_SA_h
#define L1Trigger_L1THGCal_HGCalTriggerCell_SA_h

#include <vector>
#include <memory>
#include <iostream>
namespace l1thgcfirmware {

  class HGCalTriggerCell {
  public:
    HGCalTriggerCell() : HGCalTriggerCell(false, false, 0, 0, 0, 0) {}
    HGCalTriggerCell(
        bool frameValid, bool dataValid, unsigned int rOverZ, unsigned int phi, unsigned int layer, unsigned int energy)
        : clock_(0),
          index_(0),
          rOverZ_(rOverZ),
          layer_(layer),
          energy_(energy),
          phi_(phi),
          sortKey_(0),
          deltaR2_(0),
          dX_(0),
          Y_(0),
          frameValid_(frameValid),
          dataValid_(dataValid),
          index_cmssw_(std::pair{-1, 01}) {}

    ~HGCalTriggerCell() {}

    // Setters
    void setClock(const unsigned int clock) { clock_ = clock; }
    void addLatency(const unsigned int latency) { clock_ += latency; }
    void setIndex(const unsigned int index) { index_ = index; }
    void setSortKey(const unsigned int sortKey) { sortKey_ = sortKey; }
    void setDX(const int dX) { dX_ = dX; }
    void setY(const unsigned int Y) { Y_ = Y; }
    void setDeltaR2(const unsigned int deltaR2) { deltaR2_ = deltaR2; }
    void setCmsswIndex(const std::pair<int, int> index) { index_cmssw_ = index; }

    // Getters
    unsigned int clock() const { return clock_; }
    unsigned int index() const { return index_; }
    unsigned int rOverZ() const { return rOverZ_; }
    unsigned int layer() const { return layer_; }
    unsigned int energy() const { return energy_; }
    unsigned int phi() const { return phi_; }
    unsigned int sortKey() const { return sortKey_; }
    unsigned int deltaR2() const { return deltaR2_; }
    int dX() const { return dX_; }
    unsigned int Y() const { return Y_; }
    bool frameValid() const { return frameValid_; }
    bool dataValid() const { return dataValid_; }
    const std::pair<int, int>& cmsswIndex() const { return index_cmssw_; }

    // Operators
    bool operator==(const HGCalTriggerCell& rhs) const;
    bool operator==(const std::shared_ptr<HGCalTriggerCell>& rhs) const;

  private:
    unsigned int clock_;
    unsigned int index_;
    unsigned int rOverZ_;
    unsigned int layer_;
    unsigned int energy_;
    unsigned int phi_;
    unsigned int sortKey_;
    unsigned int deltaR2_;
    int dX_;
    unsigned int Y_;
    bool frameValid_;
    bool dataValid_;
    std::pair<int, int> index_cmssw_;
  };

  typedef std::vector<HGCalTriggerCell> HGCalTriggerCellSACollection;

  typedef std::unique_ptr<HGCalTriggerCell> HGCalTriggerCellSAPtr;
  typedef std::vector<HGCalTriggerCellSAPtr> HGCalTriggerCellSAPtrCollection;
  typedef std::vector<std::vector<std::unique_ptr<HGCalTriggerCell> > > HGCalTriggerCellSAPtrCollections;
  typedef std::vector<std::vector<std::vector<std::unique_ptr<HGCalTriggerCell> > > > HGCalTriggerCellSAPtrCollectionss;

  typedef std::shared_ptr<HGCalTriggerCell> HGCalTriggerCellSAShrPtr;
  typedef std::vector<HGCalTriggerCellSAShrPtr> HGCalTriggerCellSAShrPtrCollection;
  typedef std::vector<std::vector<std::shared_ptr<HGCalTriggerCell> > > HGCalTriggerCellSAShrPtrCollections;
  typedef std::vector<std::vector<std::vector<std::shared_ptr<HGCalTriggerCell> > > >
      HGCalTriggerCellSAShrPtrCollectionss;
}  // namespace l1thgcfirmware

#endif
