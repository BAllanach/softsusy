// fastalloc.cpp  --  optional drop-in fast allocator for SOFTSUSY executables.
//
// SOFTSUSY spends a large fraction of its time allocating and freeing many
// short-lived DoubleVector / DoubleMatrix temporaries (each backed by a
// std::valarray<double>, which always heap-allocates). These appear in the
// RGE beta functions, the loop / Passarino-Veltman functions, the self
// energies and the matrix diagonalisation routines.
//
// This translation unit replaces the global operator new / delete with
// size-segregated free lists ("pools"). Freed small blocks are recycled
// instead of returned to the system allocator, which removes almost all of
// that malloc/free overhead.
//
//   * It changes NO numerical results (verified bit-identical to 17 s.f.).
//   * It is single-threaded, matching SOFTSUSY's execution model. Do NOT use
//     it if you parallelise a single spectrum calculation across threads.
//   * Just compile it and link it into the executable (see softsusy.amk patch).
//
// Correctness notes:
//   - The minimum pooled payload is 16 bytes, which safely holds the free-list
//     "next" pointer even for a zero-byte request (e.g. valarray<double>(0)).
//   - malloc() returns memory aligned for max_align_t; the fixed 16-byte header
//     preserves 16-byte alignment of the returned pointer.
//   - Over-aligned (alignas) allocations use the C++17 aligned new/delete,
//     which are intentionally NOT overridden here and fall through to the
//     default implementation, so alignment guarantees are never weakened.

#include <cstdlib>
#include <new>

namespace {
  const int         NB  = 65;   // free lists for blocks up to 65*16 = 1040 bytes
  const std::size_t HDR = 16;   // per-block header, keeps 16-byte alignment
  void* freelist[NB] = {0};
  // Bucket index for a payload of n bytes. Clamped to >= 1 so that every
  // pooled block is at least 16 bytes and can hold the free-list pointer.
  inline int bucketOf(std::size_t n){ int b = int((n + 15) / 16); return b < 1 ? 1 : b; }
}

void* operator new(std::size_t n){
  int b = bucketOf(n);
  if (b < NB){
    void* p;
    if (freelist[b]){ p = freelist[b]; freelist[b] = *reinterpret_cast<void**>(p); }
    else {
      char* raw = (char*)std::malloc(std::size_t(b)*16 + HDR);
      if(!raw) throw std::bad_alloc();
      *reinterpret_cast<int*>(raw) = b;
      p = raw + HDR;
    }
    return p;
  }
  char* raw = (char*)std::malloc(n + HDR);
  if(!raw) throw std::bad_alloc();
  *reinterpret_cast<int*>(raw) = -1;
  return raw + HDR;
}
void* operator new[](std::size_t n){ return operator new(n); }

void operator delete(void* p) noexcept {
  if(!p) return;
  char* raw = (char*)p - HDR;
  int b = *reinterpret_cast<int*>(raw);
  if (b >= 0 && b < NB){ *reinterpret_cast<void**>(p) = freelist[b]; freelist[b] = p; }
  else std::free(raw);
}
void operator delete[](void* p) noexcept { operator delete(p); }
void operator delete(void* p, std::size_t) noexcept { operator delete(p); }
void operator delete[](void* p, std::size_t) noexcept { operator delete(p); }
