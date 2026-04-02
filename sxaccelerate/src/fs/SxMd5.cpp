// ---------------------------------------------------------------------------
//
//           The general purpose cross platform C/C++ framework
//
//                       S x A c c e l e r a t e
//
//           Home:       https://www.sxlib.de
//           License:    Apache 2
//           Authors:    see src/AUTHORS
//
// ---------------------------------------------------------------------------


#include <SxMd5.h>
#include <SxFileIO.h>
#include <SxEndian.h>
#include <SxException.h>

static void transform (uint32_t *digest_, uint32_t *chunk_)
{
   SX_TRACE ();
   SX_CHECK (digest_);
   SX_CHECK (chunk_);

#  define SX_MD5_F(x, y, z) (((x) & (y)) | ((~x) & (z)))
#  define SX_MD5_G(x, y, z) (((x) & (z)) | ((y) & (~z)))
#  define SX_MD5_H(x, y, z) ((x) ^ (y) ^ (z))
#  define SX_MD5_I(x, y, z) ((y) ^ ((x) | (~z)))

#  define SX_MD5_ROTATE_LEFT(x, n) (((x) << (n)) | ((x) >> (32-(n))))

#  define SX_MD5_FF(a, b, c, d, x, s, ac)  {                         \
           (a) += SX_MD5_F ((b), (c), (d)) + (x) + (uint32_t)(ac); \
           (a) = SX_MD5_ROTATE_LEFT ((a), (s));                    \
           (a) += (b);                                             \
        }
#  define SX_MD5_GG(a, b, c, d, x, s, ac)  {                         \
           (a) += SX_MD5_G ((b), (c), (d)) + (x) + (uint32_t)(ac); \
           (a) = SX_MD5_ROTATE_LEFT ((a), (s));                    \
           (a) += (b);                                             \
        }
#  define SX_MD5_HH(a, b, c, d, x, s, ac)  {                         \
           (a) += SX_MD5_H ((b), (c), (d)) + (x) + (uint32_t)(ac); \
           (a) = SX_MD5_ROTATE_LEFT ((a), (s));                    \
           (a) += (b);                                             \
        }
#  define SX_MD5_II(a, b, c, d, x, s, ac)  {                         \
           (a) += SX_MD5_I ((b), (c), (d)) + (x) + (uint32_t)(ac); \
           (a) = SX_MD5_ROTATE_LEFT ((a), (s));                    \
           (a) += (b);                                             \
        }

   uint32_t a = digest_[0];
   uint32_t b = digest_[1];
   uint32_t c = digest_[2];
   uint32_t d = digest_[3];

   // --- Round 1
   SX_MD5_FF (a, b, c, d, chunk_[0], 7, 0xd76aa478); // 1
   SX_MD5_FF (d, a, b, c, chunk_[1], 12, 0xe8c7b756); // 2
   SX_MD5_FF (c, d, a, b, chunk_[2], 17, 0x242070db); // 3
   SX_MD5_FF (b, c, d, a, chunk_[3], 22, 0xc1bdceee); // 4
   SX_MD5_FF (a, b, c, d, chunk_[4], 7, 0xf57c0faf); // 5
   SX_MD5_FF (d, a, b, c, chunk_[5], 12, 0x4787c62a); // 6
   SX_MD5_FF (c, d, a, b, chunk_[6], 17, 0xa8304613); // 7
   SX_MD5_FF (b, c, d, a, chunk_[7], 22, 0xfd469501); // 8
   SX_MD5_FF (a, b, c, d, chunk_[8], 7, 0x698098d8); // 9
   SX_MD5_FF (d, a, b, c, chunk_[9], 12, 0x8b44f7af); // 10
   SX_MD5_FF (c, d, a, b, chunk_[10], 17, 0xffff5bb1); // 11
   SX_MD5_FF (b, c, d, a, chunk_[11], 22, 0x895cd7be); // 12
   SX_MD5_FF (a, b, c, d, chunk_[12], 7, 0x6b901122); // 13
   SX_MD5_FF (d, a, b, c, chunk_[13], 12, 0xfd987193); // 14
   SX_MD5_FF (c, d, a, b, chunk_[14], 17, 0xa679438e); // 15
   SX_MD5_FF (b, c, d, a, chunk_[15], 22, 0x49b40821); // 16

   // --- Round 2
   SX_MD5_GG (a, b, c, d, chunk_[1], 5, 0xf61e2562); // 17
   SX_MD5_GG (d, a, b, c, chunk_[6], 9, 0xc040b340); // 18
   SX_MD5_GG (c, d, a, b, chunk_[11], 14, 0x265e5a51); // 19
   SX_MD5_GG (b, c, d, a, chunk_[0], 20, 0xe9b6c7aa); // 20
   SX_MD5_GG (a, b, c, d, chunk_[5], 5, 0xd62f105d); // 21
   SX_MD5_GG (d, a, b, c, chunk_[10], 9, 0x2441453); // 22
   SX_MD5_GG (c, d, a, b, chunk_[15], 14, 0xd8a1e681); // 23
   SX_MD5_GG (b, c, d, a, chunk_[4], 20, 0xe7d3fbc8); // 24
   SX_MD5_GG (a, b, c, d, chunk_[9], 5, 0x21e1cde6); // 25
   SX_MD5_GG (d, a, b, c, chunk_[14], 9, 0xc33707d6); // 26
   SX_MD5_GG (c, d, a, b, chunk_[3], 14, 0xf4d50d87); // 27
   SX_MD5_GG (b, c, d, a, chunk_[8], 20, 0x455a14ed); // 28
   SX_MD5_GG (a, b, c, d, chunk_[13], 5, 0xa9e3e905); // 29
   SX_MD5_GG (d, a, b, c, chunk_[2], 9, 0xfcefa3f8); // 30
   SX_MD5_GG (c, d, a, b, chunk_[7], 14, 0x676f02d9); // 31
   SX_MD5_GG (b, c, d, a, chunk_[12], 20, 0x8d2a4c8a); // 32

   // --- Round 3
   SX_MD5_HH (a, b, c, d, chunk_[5], 4, 0xfffa3942); // 33
   SX_MD5_HH (d, a, b, c, chunk_[8], 11, 0x8771f681); // 34
   SX_MD5_HH (c, d, a, b, chunk_[11], 16, 0x6d9d6122); // 35
   SX_MD5_HH (b, c, d, a, chunk_[14], 23, 0xfde5380c); // 36
   SX_MD5_HH (a, b, c, d, chunk_[1], 4, 0xa4beea44); // 37
   SX_MD5_HH (d, a, b, c, chunk_[4], 11, 0x4bdecfa9); // 38
   SX_MD5_HH (c, d, a, b, chunk_[7], 16, 0xf6bb4b60); // 39
   SX_MD5_HH (b, c, d, a, chunk_[10], 23, 0xbebfbc70); // 40
   SX_MD5_HH (a, b, c, d, chunk_[13], 4, 0x289b7ec6); // 41
   SX_MD5_HH (d, a, b, c, chunk_[0], 11, 0xeaa127fa); // 42
   SX_MD5_HH (c, d, a, b, chunk_[3], 16, 0xd4ef3085); // 43
   SX_MD5_HH (b, c, d, a, chunk_[6], 23, 0x4881d05); // 44
   SX_MD5_HH (a, b, c, d, chunk_[9], 4, 0xd9d4d039); // 45
   SX_MD5_HH (d, a, b, c, chunk_[12], 11, 0xe6db99e5); // 46
   SX_MD5_HH (c, d, a, b, chunk_[15], 16, 0x1fa27cf8); // 47
   SX_MD5_HH (b, c, d, a, chunk_[2], 23, 0xc4ac5665); // 48

   // --- Round 4
   SX_MD5_II (a, b, c, d, chunk_[0], 6, 0xf4292244); // 49
   SX_MD5_II (d, a, b, c, chunk_[7], 10, 0x432aff97); // 50
   SX_MD5_II (c, d, a, b, chunk_[14], 15, 0xab9423a7); // 51
   SX_MD5_II (b, c, d, a, chunk_[5], 21, 0xfc93a039); // 52
   SX_MD5_II (a, b, c, d, chunk_[12], 6, 0x655b59c3); // 53
   SX_MD5_II (d, a, b, c, chunk_[3], 10, 0x8f0ccc92); // 54
   SX_MD5_II (c, d, a, b, chunk_[10], 15, 0xffeff47d); // 55
   SX_MD5_II (b, c, d, a, chunk_[1], 21, 0x85845dd1); // 56
   SX_MD5_II (a, b, c, d, chunk_[8], 6, 0x6fa87e4f); // 57
   SX_MD5_II (d, a, b, c, chunk_[15], 10, 0xfe2ce6e0); // 58
   SX_MD5_II (c, d, a, b, chunk_[6], 15, 0xa3014314); // 59
   SX_MD5_II (b, c, d, a, chunk_[13], 21, 0x4e0811a1); // 60
   SX_MD5_II (a, b, c, d, chunk_[4], 6, 0xf7537e82); // 61
   SX_MD5_II (d, a, b, c, chunk_[11], 10, 0xbd3af235); // 62
   SX_MD5_II (c, d, a, b, chunk_[2], 15, 0x2ad7d2bb); // 63
   SX_MD5_II (b, c, d, a, chunk_[9], 21, 0xeb86d391); // 64

   digest_[0] += a;
   digest_[1] += b;
   digest_[2] += c;
   digest_[3] += d;
}

SxMd5::SxMd5 ()
{
   SX_TRACE ();
   memset (digest, 0, 16);
   initialize ();
}

SxMd5::SxMd5 (const SxArray<uint32_t> &digest_)
{
   SX_TRACE ();
   SX_CHECK (digest_.getSize () == 4);
   memset (digest, 0, 16);
   initialize ();

   digest[0] = digest_(0);
   digest[1] = digest_(1);
   digest[2] = digest_(2);
   digest[3] = digest_(3);
}

SxMd5::SxMd5 (const SxString &checksum_)
{
   SX_TRACE ();
   SX_CHECK (checksum_.getSize () == 32);

   memset (digest, 0, 16);
   initialize ();

   bool error = false;
   bool anyError = false;

   digest [0] = SxEndian::swap32 (checksum_.subString (0, 7)
                                  .toNumber<uint32_t> (&error, 16));
   if (error) anyError = true;
   digest [1] = SxEndian::swap32 (checksum_.subString (8, 15)
                                  .toNumber<uint32_t> (&error, 16));
   if (error) anyError = true;
   digest [2] = SxEndian::swap32 (checksum_.subString (16, 23)
                                  .toNumber<uint32_t> (&error, 16));
   if (error) anyError = true;
   digest [3] = SxEndian::swap32 (checksum_.subString (24)
                                  .toNumber<uint32_t> (&error, 16));
   if (error) anyError = true;
   if (anyError) SX_THROW ("Cannot convert " + checksum_ + " to an Md5 digest.");
}


SxMd5::~SxMd5 ()
{
   SX_TRACE ();
}

void SxMd5::initialize ()
{
   SX_TRACE ();

   digest[0] = 0x67452301L;
   digest[1] = 0xefcdab89L;
   digest[2] = 0x98badcfeL;
   digest[3] = 0x10325476L;
}

void SxMd5::update (const uint8_t *buffer, ssize_t len, int64_t totalSize)
{
   SX_TRACE ();

   SX_CHECK (buffer);
   SX_CHECK (len >= 0, len);

   uint32_t nFullChunks = static_cast<uint32_t>(len >> 6);
   uint32_t chunk[16];
   uint8_t  *chunk8 = (uint8_t*)&chunk[0];

   for (size_t iChunk = 0; iChunk < nFullChunks; ++iChunk) {
      SxEndian::hostToLittle32 (buffer + ((ssize_t)iChunk << 6), chunk8, 64);
      transform (digest, chunk);
   }

   if (totalSize < 0LL) {
      SX_CHECK (len % 64 == 0 && "no bytes left", len);
   } else {
      // --- terminate last block of input data
      size_t nBytesLeft = (size_t)(len - ((size_t)nFullChunks << 6));
      const uint8_t *pSrcLeft = buffer + ((size_t)nFullChunks << 6);

      SX_CHECK (nBytesLeft < 64, nBytesLeft);
      memcpy (chunk8, pSrcLeft, nBytesLeft);
      memset (chunk8 + nBytesLeft, 0, 64 - nBytesLeft);

      // --- insert stop bit
      chunk8[nBytesLeft] = 0x80;

      size_t n = ((size_t)nBytesLeft >> 2) + (size_t)1;
      SX_CHECK (n <= 16, n);
      for (size_t i = 0; i < n; ++i) {
         chunk[i] = SxEndian::hostToLittle32 (chunk[i]);
      }

      // --- insert size of the message in bits
      if (nBytesLeft <= 55) {
         uint64_t nBits = (uint64_t)totalSize << 3;
         chunk[14] = (uint32_t)(nBits & 0xffffffff);
         chunk[15] = (uint32_t)((nBits >> 32) & 0xffffffff);
      }
      transform (digest, chunk);

      if (nBytesLeft > 55) {
         // --- extra block, align size of the message to 512 bits
         memset (chunk8, 0, 64);
         uint64_t nBits = (uint64_t)totalSize << 3;
         chunk[14] = (uint32_t)(nBits & 0xffffffff);
         chunk[15] = (uint32_t)((nBits >> 32) & 0xffffffff);
         transform (digest, chunk);
      }
   }

   // --- reset memory, can be a copy of password from input buffer
   memset (chunk, 0, 64);
}

SxString SxMd5::toString () const
{
   SX_TRACE ();

   uint32_t buffer[4];

   buffer[0] = SxEndian::hostToLittle32 (digest[0]);
   buffer[1] = SxEndian::hostToLittle32 (digest[1]);
   buffer[2] = SxEndian::hostToLittle32 (digest[2]);
   buffer[3] = SxEndian::hostToLittle32 (digest[3]);

   char resultASCII[33];
   uint8_t *p = (uint8_t *)(buffer);
   for (ssize_t i = 0; i < 16; i++) {
      sprintf (resultASCII + 2 * i, "%02x", p[i]);
   }
   resultASCII[32] = '\0';

   return SxString (resultASCII);
}

SxArray<uint32_t> SxMd5::toDigest () const
{
   SxArray<uint32_t> res(4);
   for (int i = 0; i < 4; ++i)  res(i) = digest[i];
   return res;
}

SxString SxMd5::get (const uint8_t *buffer_, ssize_t len_)
{
   SX_TRACE ();

   SX_CHECK (len_ >= 0, len_);

   SxMd5 md5;
   md5.update (buffer_, len_, len_);

   return md5.toString ();
}

SxString SxMd5::get (const SxString &string_)
{
   SX_TRACE ();

   return SxMd5::get ((const uint8_t*)string_.getElems (),
                      string_.getNBytes ());
}

SxString SxMd5::file (const SxString &filename_, uint64_t *size_)
{
   SX_TRACE ();
   SX_CHECK (filename_ != "");

   SxMd5 md5;

   // 4 Megabytes of buffer to increase read performance
   // A multiple of DX_BLOCK_SIZE
   const ssize_t len = 4 * 1024 * 1024;
   SxArray<uint8_t> buffer (len);
   int64_t nReadTl = 0;

   if (size_) {
      *size_ = 0;
   }

   SxFileIO io;

   try {
      io.open (filename_, "rb");

      while (nReadTl < (int64_t)io.getSize ()) {
         uint64_t nRead = io.read (&buffer, len);
         nReadTl += nRead;
         if (nReadTl != (int64_t)io.getSize ()) {
            md5.update (buffer.elements, (ssize_t)nRead, -1LL);
         } else {
            md5.update (buffer.elements, (ssize_t)nRead, nReadTl);
         }
      }

      io.close ();
   } catch (SxException e) {
      SX_THROW (e, "Cannot hash file " + filename_);
   }
   if (size_) {
      *size_ = nReadTl;
   }
   return md5.toString ();
}

SxArray<uint32_t> SxMd5::stringToDigest (const SxString &md5sum)
{
   SX_TRACE ();
   SX_CHECK (md5sum.getSize () == 32);

   SxArray<uint32_t> digest (4);
   bool error, anyError;
   error = false;
   anyError = false;
   digest (0) = SxEndian::swap32 (md5sum.subString (0, 7)
                                  .toNumber<uint32_t> (&error, 16));
   if (error) anyError = true;
   digest (1) = SxEndian::swap32 (md5sum.subString (8, 15)
                                  .toNumber<uint32_t> (&error, 16));
   if (error) anyError = true;
   digest (2) = SxEndian::swap32 (md5sum.subString (16, 23)
                                  .toNumber<uint32_t> (&error, 16));
   if (error) anyError = true;
   digest (3) = SxEndian::swap32 (md5sum.subString (24)
                                  .toNumber<uint32_t> (&error, 16));
   if (error) anyError = true;
   if (anyError) SX_THROW ("Cannot convert " + md5sum + " to an Md5 digest.");
   return digest;
}

SxString SxMd5::digestToString (const uint32_t *digest_)
{
   uint32_t digest[4];

   digest[0] = SxEndian::hostToLittle32 (digest_[0]);
   digest[1] = SxEndian::hostToLittle32 (digest_[1]);
   digest[2] = SxEndian::hostToLittle32 (digest_[2]);
   digest[3] = SxEndian::hostToLittle32 (digest_[3]);

   char resultASCII[33];
   uint8_t* p = (uint8_t*)&digest[0];
   for (ssize_t i = 0; i < 16; i++) {
      sprintf (resultASCII + 2 * i, "%02x", p[i]);
   }
   resultASCII[32] = '\0';

   return SxString (resultASCII);
}
