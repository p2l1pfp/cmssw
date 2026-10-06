#ifndef L1Trigger_Phase2L1ParticleFlow_bitonic_vhdl_ref_h
#define L1Trigger_Phase2L1ParticleFlow_bitonic_vhdl_ref_h

#include <algorithm>
#include <cassert>
#include <vector>

// Bit-exact model of the bitonic-sort-VHDL network (bitonic_sort_ii2 and its sub-entities), including the
// placement of objects with equal keys. Array index i corresponds to VHDL index i; after sorting, index 0 holds
// the highest key. The helpers in l1ct::bitonic_vhdl work on an internal buffer; use bitonic_vhdl_sort_and_crop_ref.
// Only hwPt is used for the comparison, as in the firmware (COMPARISON_WIDTH = PT_BIT_WIDTH).
namespace l1ct {
  namespace bitonic_vhdl {

    // bitonic_split: compare in_a(i) with in_b(i); a = upper-index half, b = lower-index half
    template <typename T>
    void split(T a[], T b[], unsigned int n, bool plus) {
      for (unsigned int i = 0; i < n; ++i) {
        if (!(a[i].hwPt <= b[i].hwPt))  // otherwise lo = a, hi = b
          std::swap(a[i], b[i]);
        // now a = lo, b = hi
        if (!plus)
          std::swap(a[i], b[i]);
      }
    }

    // bitonic_merge (and bitonic_merge_ii2, which only time-multiplexes the same comparisons)
    template <typename T>
    void merge(T a[], T b[], unsigned int n, bool plus) {
      split(a, b, n, plus);
      if (n > 1) {
        merge(a + n / 2, a, n / 2, plus);
        merge(b + n / 2, b, n / 2, plus);
      }
    }

    // bitonic_sort: upper half sorted with PLUS = '1', lower half with PLUS = '0', then merged
    template <typename T>
    void sort(T x[], unsigned int n, bool plus) {
      if (n > 1) {
        sort(x + n / 2, n / 2, true);
        sort(x, n / 2, false);
        merge(x + n / 2, x, n / 2, plus);
      }
    }

    // bitonic_sort_ii2 (PLUS = '1'): both halves sorted with the same sorter, the lower one reversed before the merge
    template <typename T>
    void sort_ii2(T x[], unsigned int n) {
      assert(n > 0 && (n & (n - 1)) == 0);
      if (n > 1) {
        sort(x + n / 2, n / 2, true);
        sort(x, n / 2, true);
        std::reverse(x, x + n / 2);
        merge(x + n / 2, x, n / 2, true);
      }
    }

  }  // namespace bitonic_vhdl
}  // namespace l1ct

// The VHDL sorter always takes power-of-2 inputs, so the input is padded with cleared objects
// (nominally it produces the same size output, though things may get optimized away in the implementation)
template <typename T>
void bitonic_vhdl_sort_and_crop_ref(unsigned int nIn, unsigned int nOut, const T in[], T out[]) {
  unsigned int nPadded = 1;
  while (nPadded < nIn)
    nPadded <<= 1;
  assert(nOut <= nPadded);
  std::vector<T> work(nPadded);
  for (unsigned int i = 0; i < nPadded; ++i) {
    if (i < nIn)
      work[i] = in[i];
    else
      work[i].clear();
  }
  l1ct::bitonic_vhdl::sort_ii2(work.data(), nPadded);
  for (unsigned int i = 0; i < nOut; ++i)
    out[i] = work[i];
}

#endif
