#ifndef LR3_TYPES
#define LR3_TYPES


typedef signed short int int16;
typedef signed long int int32;
typedef signed long long int int64;

typedef unsigned short int uint16;
typedef unsigned long int uint32;
typedef unsigned long long int uint64;

#ifndef u32
#	define u32 uint32
#endif

//#ifndef uint
//#	define uint uint32
//#endif

// [McM] ZScript-style pointer:
// [KpH] Sry, but it is being used in Qt :<
//#define self (*this)

// Looking at disasm, compilers have been seen to be stupid about inlining some single-instruction SIMD intrinsics functions, so use this to force.
#ifdef _DEBUG
#define FORCE_INLINE inline
#elif defined(_MSC_VER)
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline __attribute__((always_inline))
#endif

// The CONST_WIN32 is a #define which resolves to 'const' on Windows, and null on other
// platforms. This #define is used on Windows to detect accidental programming errors
// occurring from an expression "const float3 vec; vec[1] = 5;". Trying to return
// const float from operator[] on GCC gives a warning "type qualifiers ignored on function return type",
// so hence this is only enabled on Visual Studio.
#ifdef _MSC_VER
#define CONST_WIN32 const
#else
#define CONST_WIN32
#endif

#endif
