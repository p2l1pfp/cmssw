#ifndef __L1Trigger_L1THGCal_CentroidHelper_h__
#define __L1Trigger_L1THGCal_CentroidHelper_h__

#include <vector>

namespace l1thgcfirmware {

  class CentroidHelper {
  public:
    CentroidHelper() : CentroidHelper(0, 0, 0, 0, 0, 0, 0, false) {}
    CentroidHelper(unsigned int clock, unsigned int index, unsigned int dataValid)
        : CentroidHelper(clock, index, 0, 0, 0, 0, 0, dataValid) {}
    CentroidHelper(unsigned int clock,
                   unsigned int index,
                   unsigned int column,
                   unsigned int row,
                   unsigned int energy,
                   unsigned int X,
                   unsigned int Y,
                   unsigned int dataValid)
        : clock_(clock),
          index_(index),
          column_(column),
          row_(row),
          energy_(energy),
          x_(X),
          y_(Y),
          dataValid_(dataValid) {}
    ~CentroidHelper() {}

    // Getters
    unsigned int clock() const { return clock_; }
    unsigned int index() const { return index_; }
    unsigned int column() const { return column_; }
    unsigned int row() const { return row_; }
    unsigned int energy() const { return energy_; }
    unsigned int X() const { return x_; }
    unsigned int Y() const { return y_; }
    bool dataValid() const { return dataValid_; }

    // Setters
    void setClock(const unsigned int clock) { clock_ = clock; }
    void setIndex(const unsigned int index) { index_ = index; }

  private:
    unsigned int clock_;
    unsigned int index_;
    unsigned int column_;
    unsigned int row_;
    unsigned int energy_;
    unsigned int x_;
    unsigned int y_;
    bool dataValid_;
  };

  typedef std::unique_ptr<CentroidHelper> CentroidHelperPtr;
  typedef std::vector<CentroidHelperPtr> CentroidHelperPtrCollection;
  typedef std::vector<CentroidHelperPtrCollection> CentroidHelperPtrCollections;

  typedef std::shared_ptr<CentroidHelper> CentroidHelperShrPtr;
  typedef std::vector<CentroidHelperShrPtr> CentroidHelperShrPtrCollection;
}  // namespace l1thgcfirmware

#endif
