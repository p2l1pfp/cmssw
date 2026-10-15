#ifndef pfsorter_tree_elements_ref_h
#define pfsorter_tree_elements_ref_h

// Merge tree of the pfsorter, adapted from the multififo regionizer
// (L1Trigger/Phase2L1ParticleFlow/interface/regionizer/multififo_regionizer_elements_ref.h)
// keeping the same structure, object and function names.
//
// Differences with respect to the regionizer:
//  * there is a 1:1:1 correspondence between input link, buffer and fifo of the first
//    layer, so an object pushed on link i always goes into fifo i: there is no routing and
//    no eta/phi selection, hence no sector or region geometry;
//  * the number of links is free (the regionizer only supports 1, 2, 3, 4, 6, 8 and 12),
//    so a layer with an odd number of nodes merges some pairs and carries the remaining
//    node through;
//  * the tree stops at the requested number of output nodes, at most 3, as in the
//    regionizer.
// As in the regionizer, every node of the tree (staging areas and queues alike) holds one
// object at most, and pop() returns one object per clock cycle: the first output node that
// holds one, the others waiting for their turn.
//
// There is no migration between nodes of the same layer: the sources of the nodes of a
// layer partition the nodes of the previous one, so the object held by a node has exactly
// one possible destination and can never end up in one of its siblings. The only movement
// inside a node is the one from its staging area to its queue register.

#include <algorithm>
#include <list>
#include <vector>
#include <cassert>
#include <cstdio>

namespace l1ct {
  namespace pfsorter {

    // same helpers as the multififo regionizer: an object is "empty" if hwPt == 0
    template <typename T>
    inline void shift(T& from, T& to) {
      to = from;
      from.clear();
    }
    template <typename TL, typename T>
    inline void pop_back(TL& from, T& to) {
      assert(!from.empty());
      to = from.back();
      from.pop_back();
    }

    // the fifos are filled from the front and read from the back
    template <typename T>
    inline void push_to_fifo(const T& t, std::list<T>& fifo) {
      fifo.push_front(t);
    }

    // no region to select on, so an object is pushed as soon as it is not empty
    template <typename T>
    inline void maybe_push(const T& t, std::list<T>& fifo) {
      if (t.hwPt != 0)
        push_to_fifo(t, fifo);
    }

    // trace() identifies the objects with the id returned by this kind of function; a test
    // typically returns the barcode it stored in one of the 'src' pointers of the object.
    // Without one, the objects are identified by their pt.
    template <typename T>
    using ObjId = unsigned int (*)(const T&);

    // one object in a fixed-width column of a trace line: its id, or '.' if empty
    template <typename T>
    inline void trace_cell(FILE* f, const T& t, ObjId<T> objId) {
      if (t.hwPt == 0)
        fprintf(f, "   . ");
      else
        fprintf(f, "%4u ", objId ? (*objId)(t) : unsigned(t.intPt()));
    }

    template <typename T>
    class RegionBuffer {
    public:
      RegionBuffer() : nfifos_(0), noutputs_(0)  {}
      // one fifo per link in each of the two banks, merged down to 'noutputs' output
      // nodes (at most 3)
      void initFifos(unsigned int nfifos, unsigned int noutputs = 3);
      void flush();
      // roll to a new event, as the firmware does: the link fifos are double buffered, so
      // the read and write pointers swap banks and the contents of the new write bank
      // survive and drain into the new event; the merge nodes discard their staging areas
      // and queues on roll. So only the queues are cleared, not the fifos.
      void roll();
      // flush() + reset the dropped-object and postponed-readout counters
      void reset();
      // point the write pointer at the bank the read pointer is on, so that pushes and
      // pops share the same set of fifos. This is the event-by-event mode CMSSW uses;
      // the firmware instead keeps the two on opposite banks (see roll()).
      void setWriteToRead() { iwrite_ = iread_; }
      void maybe_push(int fifo, const T& t);
      T pop();

      // --- tracing and bookkeeping (not in the regionizer) ---
      unsigned int nfifos() const { return nfifos_; }
      unsigned int noutputs() const { return noutputs_; }
      unsigned int nstages() const { return queues_.size(); }
      // the banks the read and write pointers are on: inputs go into bank iread(), the
      // output is popped out of bank iwrite()
      unsigned int iread() const { return iread_; }
      unsigned int iwrite() const { return iwrite_; }
      unsigned int fifoSize(unsigned int i) const { return fifos_[iread_][i].size(); }
      // the object that will be popped first out of fifo i
      const T* fifoHead(unsigned int i) const {
        return fifos_[iread_][i].empty() ? nullptr : &fifos_[iread_][i].back();
      }
      const std::vector<T>& stage(unsigned int i) const { return queues_[i].first; }
      const std::vector<T>& queue(unsigned int i) const { return queues_[i].second; }
      // objects still inside the fifos and the merge stages
      unsigned int inFlight() const;
      // objects dropped because their fifo was full, and read-outs postponed because
      // another output node was served first
      void trace(FILE* f, ObjId<T> objId = nullptr) const;

    private:
      unsigned int nfifos_, noutputs_;
      // two banks of link fifos, one per event in flight: the read pointer iread_ is on
      // the bank the links are written into, the write pointer iwrite_ on the one the
      // merge tree is drained from, and the two swap on every roll()
      std::vector<std::list<T>> fifos_[2];
      unsigned int iread_ = 0, iwrite_ = 1;
      std::vector<std::pair<std::vector<T>, std::vector<T>>> queues_;

      T pop_next_trivial_();
      void fifos_to_stage_(std::vector<T>& staging_area);
      void queue_to_stage_(std::vector<T>& queue, std::vector<T>& staging_area);
      void stage_to_queue_(std::vector<T>& staging_area, std::vector<T>& queue);
      T pop_queue_(std::vector<T>& queue);

      // first and last of the 'nin' sources merged into node j of a layer of 'nout' nodes;
      // in the regionizer this is always the pair (2j, 2j+1), here a layer with an odd
      // number of nodes has some nodes with a single source
      static unsigned int src_first_(unsigned int j, unsigned int nin, unsigned int nout) { return (j * nin) / nout; }
      static unsigned int src_last_(unsigned int j, unsigned int nin, unsigned int nout) {
        return ((j + 1) * nin) / nout - 1;
      }
    };

  }  // namespace pfsorter
}  // namespace l1ct

// the implementation is a template one, so it is included here instead of being
// explicitly instantiated in a .cpp file
#include "L1Trigger/Phase2L1ParticleFlow/interface/pfsorter/pfsorter_tree_elements_ref.icc"

#endif
