#ifndef _SX_ENDIAN_H_
#define _SX_ENDIAN_H_

#include <SxError.h>
#include <SxTypeDefs.h>

/** Endianness.

   - little-endian: x86
   - middle-endian: not supported
   - big-endian: PowerPC

   - network-endian: big-endian
   - host-endian: little/big-endian (automatic detection)

   Standard htons(), htonl(), ntohs(), ntohl() are probably faster
   since they can use special instructions for conversions.

   \par Features:
   - Support for little endian.
   - (64-bit integers, not yet).
   - Memory block swapping. A convenience for some image file formats.

\code
    // --- Big-Endian (MSB first) : define HAS_BIG_ENDIAN

    uint32_t crc = 0x504c5445;
    unsigned char *p = (unsigned char*)&crc;

    SX_CHECK( p[0] == 0x50 ); // Most Significant Byte
    SX_CHECK( p[1] == 0x4c );
    SX_CHECK( p[2] == 0x54 );
    SX_CHECK( p[3] == 0x45 );
\endcode

    \ingroup group_crypto
    \author Vaclav Bubnik, bubnik@mpie.de
*/
class SxEndian
{
   public:

   static bool isBig ();
   static bool isLittle ();

   static uint32_t bitReverse (uint32_t value, int nBits);

   static uint16_t swap16 (uint16_t value);
   static uint32_t swap32 (uint32_t value);

   // --- unaligned swap
   static void swap16 (const uint8_t *src, uint8_t *dst, ssize_t nbytes);
   static void swap32 (const uint8_t *src, uint8_t *dst, ssize_t nbytes);
   static void swap64 (const uint8_t *src, uint8_t *dst, ssize_t nbytes);

   // --- big-endian
   static uint16_t hostToBig16 (uint16_t value);
   static uint32_t hostToBig32 (uint32_t value);
   static void hostToBig16 (uint16_t src, uint8_t *dst);
   static void hostToBig32 (uint32_t src, uint8_t *dst);
   static void hostToBig64 (uint64_t src, uint8_t *dst);
   static void hostToBig16 (const uint8_t *src, uint8_t *dst, ssize_t nbytes);
   static void hostToBig32 (const uint8_t *src, uint8_t *dst, ssize_t nbytes);
   static void hostToBig64 (const uint8_t *src, uint8_t *dst, ssize_t nbytes);

   static uint16_t bigToHost16 (uint16_t value);
   static uint32_t bigToHost32 (uint32_t value);
   static uint16_t bigToHost16 (const uint8_t *src);
   static uint32_t bigToHost32 (const uint8_t *src);
   static uint64_t bigToHost64 (const uint8_t *src);

   // --- little-endian
   static uint16_t hostToLittle16 (uint16_t value);
   static uint32_t hostToLittle32 (uint32_t value);
   static void hostToLittle16 (uint16_t src, uint8_t *dst);
   static void hostToLittle32 (uint32_t src, uint8_t *dst);
   static void hostToLittle16 (const uint8_t *src, uint8_t *dst,ssize_t nbytes);
   static void hostToLittle32 (const uint8_t *src, uint8_t *dst,ssize_t nbytes);

   static uint16_t littleToHost16 (uint16_t value);
   static uint32_t littleToHost32 (uint32_t value);
   static uint16_t littleToHost16 (const uint8_t *src);
   static uint32_t littleToHost32 (const uint8_t *src);
};

#include <SxEndian.hpp>

#endif /* _SX_ENDIAN_H_ */
