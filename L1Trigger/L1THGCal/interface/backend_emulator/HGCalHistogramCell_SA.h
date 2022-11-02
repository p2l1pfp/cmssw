#ifndef L1Trigger_L1THGCal_HGCalHistogramCell_SA_h
#define L1Trigger_L1THGCal_HGCalHistogramCell_SA_h

#include <vector>
#include <memory>
#include <iostream>
namespace l1thgcfirmware {

  class HGCalHistogramCell {
  public:
    HGCalHistogramCell(unsigned int clock,
                       unsigned int index,
                       unsigned int S,
                       unsigned int X,
                       unsigned int Y,
                       unsigned int N,
                       unsigned int sortKey)
        : clock_(clock),
          index_(index),
          S_(S),
          X_(X),
          Y_(Y),
          N_(N),
          sortKey_(sortKey),
          frameValid_(true),
          dataValid_(true) {}

    HGCalHistogramCell()
        : clock_(0), index_(0), S_(0), X_(0), Y_(0), N_(0), sortKey_(0), frameValid_(false), dataValid_(false) {}

    HGCalHistogramCell(unsigned int clock, unsigned int index, unsigned int sortKey, unsigned latencyOffset)
        : HGCalHistogramCell(latencyOffset * sortKey + clock, index, 0, 0, 0, 0, sortKey) {}

    ~HGCalHistogramCell() {}

    // Setters
    void setClock(const unsigned int clock) { clock_ = clock; }
    void addLatency(const unsigned int latency) { clock_ += latency; }
    void setIndex(const unsigned int index) { index_ = index; }
    void setSortKey(const unsigned int sortKey) { sortKey_ = sortKey; }
    void setS(const unsigned int S) { S_ = S; }
    void setX(const unsigned int X) { X_ = X; }
    void setY(const unsigned int Y) { Y_ = Y; }
    void setN(const unsigned int N) { N_ = N; }

    // Getters
    unsigned int clock() const { return clock_; }
    unsigned int index() const { return index_; }
    unsigned int S() const { return S_; }
    unsigned int X() const { return X_; }
    unsigned int Y() const { return Y_; }
    unsigned int N() const { return N_; }
    unsigned int sortKey() const { return sortKey_; }
    bool frameValid() const { return frameValid_; }
    bool dataValid() const { return dataValid_; }

    // Operators
    const HGCalHistogramCell& operator+=(const HGCalHistogramCell& hc);
    const HGCalHistogramCell operator/(const unsigned int factor) const;
    const HGCalHistogramCell operator+(const HGCalHistogramCell& hc) const;

  private:
    unsigned int clock_;
    unsigned int index_;
    unsigned int S_;  // Sum of cell contents i.e. total p_T in histogram bin
    unsigned int
        X_;  // Average X coordinate of trigger cells in a cell.  Will contain sum of X coordinates whilst trigger cells are being added to this cell
    unsigned int
        Y_;  // Average X coordinate of trigger cells in a cell.  Will contain sum of Y coordinates whilst trigger cells are being added to this cell
    unsigned int N_;  // Number of entries (trigger cells) in a cell
    unsigned int sortKey_;
    bool frameValid_;
    bool dataValid_;
  };

  typedef std::vector<HGCalHistogramCell> HGCalHistogramCellSACollection;
  typedef std::vector<std::unique_ptr<HGCalHistogramCell>> HGCalHistogramCellSAPtrCollection;
}  // namespace l1thgcfirmware

#endif
