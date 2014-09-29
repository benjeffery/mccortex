/*
Adapted from lookup3.c by Isaac Turner Sept 2014
-------------------------------------------------------------------------------
lookup3.c, by Bob Jenkins, May 2006, Public Domain.

These are functions for producing 32-bit hashes for hash table lookup.
bklk3_hashword(), bklk3_hashlittle(), bklk3_hashlittle2(), hashbig(), bklk3_mix(), and bklk3_final()
are externally useful functions.  Routines to test the hash are included
if SELF_TEST is defined.  You can use this free for any purpose.  It's in
the public domain.  It has no warranty.

You probably want to use bklk3_hashlittle().  bklk3_hashlittle() and hashbig()
hash byte arrays.  bklk3_hashlittle() is is faster than hashbig() on
little-endian machines.  Intel and AMD are little-endian machines.
On second thought, you probably want bklk3_hashlittle2(), which is identical to
bklk3_hashlittle() except it returns two 32-bit hashes for the price of one.
You could implement hashbig2() if you wanted but I haven't bothered here.

If you want to find a hash of, say, exactly 7 integers, do
  a = i1;  b = i2;  c = i3;
  bklk3_mix(a,b,c);
  a += i4; b += i5; c += i6;
  bklk3_mix(a,b,c);
  a += i7;
  bklk3_final(a,b,c);
then use c as the hash value.  If you have a variable length array of
4-byte integers to hash, use bklk3_hashword().  If you have a byte array (like
a character string), use bklk3_hashlittle().  If you have several byte arrays, or
a mix of things, see the comments above bklk3_hashlittle().

Why is this so big?  I read 12 bytes at a time into 3 4-byte integers,
then mix those integers.  This is fast (you can do a lot more thorough
mixing with 12*3 instructions on 3 integers than you can with 3 instructions
on 1 byte), but shoehorning those bytes into integers efficiently is messy.
-------------------------------------------------------------------------------
*/

#ifndef BK_LOOKUP3_H_
#define BK_LOOKUP3_H_

#include <stdint.h>     /* defines uint32_t etc */

#define bklk3_rot(x,k) (((x)<<(k)) | ((x)>>(32-(k))))

/*
-------------------------------------------------------------------------------
mix -- mix 3 32-bit values reversibly.

This is reversible, so any information in (a,b,c) before bklk3_mix() is
still in (a,b,c) after bklk3_mix().

If four pairs of (a,b,c) inputs are run through bklk3_mix(), or through
bklk3_mix() in reverse, there are at least 32 bits of the output that
are sometimes the same for one pair and different for another pair.
This was tested for:
* pairs that differed by one bit, by two bits, in any combination
  of top bits of (a,b,c), or in any combination of bottom bits of
  (a,b,c).
* "differ" is defined as +, -, ^, or ~^.  For + and -, I transformed
  the output delta to a Gray code (a^(a>>1)) so a string of 1's (as
  is commonly produced by subtraction) look like a single 1-bit
  difference.
* the base values were pseudorandom, all zero but one bit set, or
  all zero plus a counter that starts at zero.

Some k values for my "a-=c; a^=bklk3_rot(c,k); c+=b;" arrangement that
satisfy this are
    4  6  8 16 19  4
    9 15  3 18 27 15
   14  9  3  7 17  3
Well, "9 15 3 18 27 15" didn't quite get 32 bits diffing
for "differ" defined as + with a one-bit base and a two-bit delta.  I
used http://burtleburtle.net/bob/hash/avalanche.html to choose
the operations, constants, and arrangements of the variables.

This does not achieve avalanche.  There are input bits of (a,b,c)
that fail to affect some output bits of (a,b,c), especially of a.  The
most thoroughly mixed value is c, but it doesn't really even achieve
avalanche in c.

This allows some parallelism.  Read-after-writes are good at doubling
the number of bits affected, so the goal of mixing pulls in the opposite
direction as the goal of parallelism.  I did what I could.  Rotates
seem to cost as much as shifts on every machine I could lay my hands
on, and rotates are much kinder to the top and bottom bits, so I used
rotates.
-------------------------------------------------------------------------------
*/
#define bklk3_mix(a,b,c) \
{ \
  a -= c;  a ^= bklk3_rot(c, 4);  c += b; \
  b -= a;  b ^= bklk3_rot(a, 6);  a += c; \
  c -= b;  c ^= bklk3_rot(b, 8);  b += a; \
  a -= c;  a ^= bklk3_rot(c,16);  c += b; \
  b -= a;  b ^= bklk3_rot(a,19);  a += c; \
  c -= b;  c ^= bklk3_rot(b, 4);  b += a; \
}

/*
-------------------------------------------------------------------------------
final -- final mixing of 3 32-bit values (a,b,c) into c

Pairs of (a,b,c) values differing in only a few bits will usually
produce values of c that look totally different.  This was tested for
* pairs that differed by one bit, by two bits, in any combination
  of top bits of (a,b,c), or in any combination of bottom bits of
  (a,b,c).
* "differ" is defined as +, -, ^, or ~^.  For + and -, I transformed
  the output delta to a Gray code (a^(a>>1)) so a string of 1's (as
  is commonly produced by subtraction) look like a single 1-bit
  difference.
* the base values were pseudorandom, all zero but one bit set, or
  all zero plus a counter that starts at zero.

These constants passed:
 14 11 25 16 4 14 24
 12 14 25 16 4 14 24
and these came close:
  4  8 15 26 3 22 24
 10  8 15 26 3 22 24
 11  8 15 26 3 22 24
-------------------------------------------------------------------------------
*/
#define bklk3_final(a,b,c) \
{ \
  c ^= b; c -= bklk3_rot(b,14); \
  a ^= c; a -= bklk3_rot(c,11); \
  b ^= a; b -= bklk3_rot(a,25); \
  c ^= b; c -= bklk3_rot(b,16); \
  a ^= c; a -= bklk3_rot(c,4);  \
  b ^= a; b -= bklk3_rot(a,14); \
  c ^= b; c -= bklk3_rot(b,24); \
}

/*
-------------------------------------------------------------------------------
bklk3_hashlittle() -- hash a variable-length key into a 32-bit value
  k       : the key (the unaligned variable-length array of bytes)
  length  : the length of the key, counting by bytes
  initval : can be any 4-byte value
Returns a 32-bit value.  Every bit of the key affects every bit of
the return value.  Two keys differing by one or two bits will have
totally different hash values.

The best hash table sizes are powers of 2.  There is no need to do
mod a prime (mod is sooo slow!).  If you need less than 32 bits,
use a bitmask.  For example, if you need only 10 bits, do
  h = (h & bklk3_hashmask(10));
In which case, the hash table should have bklk3_hashsize(10) elements.

If you are hashing n strings (uint8_t **)k, do it like this:
  for (i=0, h=0; i<n; ++i) h = bklk3_hashlittle( k[i], len[i], h);

By Bob Jenkins, 2006.  bob_jenkins@burtleburtle.net.  You may use this
code any way you wish, private, educational, or commercial.  It's free.

Use for hash table lookup, or anything where one collision in 2^^32 is
acceptable.  Do NOT use for cryptographic purposes.
-------------------------------------------------------------------------------
*/

static inline uint32_t bklk3_hashlittle(const BinaryKmer bkmer, uint32_t initval)
{
  uint32_t a, b, c; /* internal state */

  /* Set up the internal state */
  a = b = c = 0xdeadbeef + ((uint32_t)BKMER_BYTES) + initval;

  const uint32_t *k = (const uint32_t *)bkmer.b;    /* read 32-bit chunks */

  #if BKMER_BYTES > 12
    /*------ all but last block: aligned reads and affect 32 bits of (a,b,c) */
    size_t i;
    for(i = 0; i+12 < BKMER_BYTES; i += 12)
    {
      a += k[0];
      b += k[1];
      c += k[2];
      bklk3_mix(a, b, c);
      k += 3;
    }
  #endif

  #if (BKMER_BYTES % 12 == 0) && (BKMER_BYTES > 0)
    a += k[0];
    b += k[1];
    c += k[2];
  #elif BKMER_BYTES % 12 == 8
    a += k[0];
    b += k[1];
  #elif BKMER_BYTES % 12 == 4
    a += k[0];
  #else
    #error BKMER_BYTES is not a multiple of 4 and greater than zero
  #endif

  bklk3_final(a, b, c);
  return c;
}

#endif /* BK_LOOKUP3_H_ */
