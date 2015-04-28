#ifndef MTL_MACROS_H
//A Macro to get the version of gcc. if we're not compiling on gcc, all those defines are zero. and thus GCC_VERSION is zero. If GCC reaches a version number of 100 this macro gets invalid...
#define GCC_VER(x,y,z)  ((x) * 10000 + (y) * 100 + (z))

#ifndef GCC_VERSION
#define GCC_VERSION GCC_VER(__GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__)
#endif
//Helper Macros
#define Scalar(MatrixType) typename MatrixType::Elementtype

//the last version I heard of that supports this seems to be 3.4.1
#if (GCC_VERSION > GCC_VER(3,4,1))
#define MTL_RESTRICT __restrict__
#endif

#if GCC_VERSION >  GCC_VER(3,0,0)
#define MTL_NEVER_NEEDED_INLINE __attribute__ ((always_inline, flatten, visibility("internal")))
#define MTL_PURE_FUNCTION __attribute__ ((pure))
#define MTL_CONST_FUNCTION __attribute__ ((const))
#endif

//_MSC_VER <= 16 means that we require Visual Studio 2010
#if (_MSC_VER >=1600) || (GCC_VERSION >= GCC_VER(4,3,0) && __GXX_EXPERIMENTAL_CXX0X__) || (__INTEL_COMPILER >= 1110) || defined(__clang__ )
#define HAS_RVALUE_REFERENCES
#endif

#if (_MSC_VER >=1600) || (GCC_VERSION >= GCC_VER(4,4,0) && __GXX_EXPERIMENTAL_CXX0X__) || (__INTEL_COMPILER >= 1100) || defined(__clang__ ) || (__IBMCPP__ >= 1110)
#define HAS_AUTO_KEYWORD
#endif

#if (_MSC_VER >=1600) || (GCC_VERSION >= GCC_VER(4,3,0) && __GXX_EXPERIMENTAL_CXX0X__) || (__INTEL_COMPILER >= 1100) || (__clang_major__ >= 3) || (__IBMCPP__ >= 1110)
#define HAS_DECL_TYPE
#endif

#ifdef __GNUC__
#ifndef DEBUGPRED
#define unlikely(expr) __builtin_expect(!!(expr), 0)
#define likely(expr) __builtin_expect(!!(expr), 1)
#else
asm (".section predict_data, \"aw\"; .previous\n"
     ".section predict_line, \"a\"; .previous\n"
     ".section predict_file, \"a\"; .previous");
#ifdef __x86_64__
#define debugpred__(e, E) \
({ long int _e = !!(e); \
asm volatile (".pushsection predict_data\n" \
              "..predictcnt%=: .quad 0; .quad 0\n" \
              ".section predict_line; .quad %c1\n" \
              ".section predict_file; .quad %c2; .popsection\n" \
              "addq $1,..predictcnt%=(,%0,8)" \
              : : "r" ((unsigned long)(_e == E)), "i" (__LINE__), "i" (__FILE__)); \
__builtin_expect (_e, E); \
})
#else
#error "other architectures not supported"
#endif
#define unlikely(expr) debugpred__((expr), 0)
#define likely(expr) debugpred__((expr), 1)
#endif
#else
#define unlikely(expr) expr
#define likely(expr) expr
#endif

//supported MACROS
#ifndef MTL_NEVER_NEEDED_INLINE
#define MTL_NEVER_NEEDED_INLINE
#endif
#ifndef MTL_PURE_FUNCTION
#define MTL_PURE_FUNCTION
#endif
#ifndef MTL_CONST_FUNCTION
#define MTL_CONST_FUNCTION
#endif
#ifndef MTL_RESTRICT
#define MTL_RESTRICT
#endif

#endif
