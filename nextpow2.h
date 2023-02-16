#pragma once

/*
 * Macro `nextpow2` for computing the next largest power of 2 using the size_t
 * type defined in <stddef.h>. Example of how to use:
 *
 *     size_t value;
 *
 *     value = 34;
 *     nextpow2(value);
 *     assert(value == 64);
 *
 *     value = 899;
 *     nextpow2(value);
 *     assert(value == 1024);
 *
 *     value = 1;
 *     nextpow2(value);
 *     assert(value == 2);
 *
 * There is undefined behaviour if value <= 0, and if the expression
 * passed to `nextpow2` is not an lvalue. For example:
 *
 *     value = 0;
 *     nextpow2(value); // undefined
 *
 *     value = -2;
 *     nextpow2(value); // undefined
 *
 *     value = 100;
 *     nextpow2(value + 10); // will not compile.
 */

#include <stdint.h>
#include <limits.h>

#if SIZE_MAX == 0xffffffff
#define nextpow2(x) up32(x)
#elif SIZE_MAX == 0xffffffffffffffff
#define nextpow2(x) up64(x)
#else
#error "sizeof(size_t) is abnormal; it is expected to be either 4 or 8 bytes"
#endif

#define up32(x) (((x) = (((x) != 1)? (x)-1 : (x))), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)++)
#define up64(x) (((x) = (((x) != 1)? (x)-1 : (x))), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, (x)++)

