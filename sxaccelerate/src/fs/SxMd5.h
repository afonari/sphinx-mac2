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

#ifndef _SX_MD5_H_
#define _SX_MD5_H_

#include <SxFS.h>
#include <SxString.h>
#include <SxHashFunction.h>

/** MD5 checksum.

    The MD5 Message-Digest Algorithm (RFC1321).

    \par Features
    - MD5 is considered to be weak for cryptography.
    - 80 MB/s on 1.33 GHz PowerPC G4, 512 MB RAM
    - 550 MB/s on 2.93 GHz Core i7

\par Simple interface
\code
   cout << SxMd5::get ("A simple string.") << endl;
\endcode


\par Advanced interface
\code
   // --- Accumulate md5 from n fragments
   //     Non-terminating blocks of data must be multiple of 64 bytes.
   uint8_t  buffer1[128]; // 128 bytes used
   uint8_t  buffer2[128]; // 10 bytes used

   SxMd5 md5;
   md5.update (digest, buffer1  128,  -1);
   md5.update (digest, buffer2,  10, 138);

   cout << md5.toString () << endl;
\endcode


    \par Testing
    - md5.cpp Tested on 32-bit big endian hardware (5.12.2009).
    - md5.cpp Tested on 64-bit little endian hardware (2.12.2011).


    \ingroup group_crypto
    \ingroup group_checksum

    \author Vaclav Bubnik, bubnik@mpie.de
    \author Bjoern Lange, bjoern.lange@dacs-labs.com */
class SX_EXPORT_FS SxMd5
{
   public:

      SxMd5 ();
      SxMd5 (const SxArray<uint32_t> &digest_);
      SxMd5 (const SxString &checksum_);
     ~SxMd5 ();

      void initialize ();

      /** Accumulate checksum from fragmented input data.
          \param buffer Raw data in memory.
          \param len The number of bytes in buffer. Must be multiple of 64 bytes
                    with exception for the last update which can use any number.
          \param totalSize Set to -1 if there are more data to update.
                           The number of bytes in complete message. */
      void update (const uint8_t *buffer, ssize_t len, int64_t totalSize);

      /** Return md5 digest as a string.
          example: "d41d8cd98f00b204e9800998ecf8427e" */
      SxString toString () const;

      /** Return md5 digest as an array of 4 32-bit integers. */
      SxArray<uint32_t> toDigest () const;

      bool operator== (const SxMd5 &in) const
      {
         return    (digest[0] == in.digest[0]) && (digest[1] == in.digest[1])
                && (digest[2] == in.digest[2]) && (digest[3] == in.digest[3]);
      }

      bool operator!= (const SxMd5 &in) const
      {
         return    (digest[0] != in.digest[0]) || (digest[1] != in.digest[1])
                || (digest[2] != in.digest[2]) || (digest[3] != in.digest[3]);
      }

      static size_t hash (const SxMd5 &in)
      {
         return SxHashFunction::murmur (in.digest, 16, 12345);
      }

      // Convenient functions
      static SxString get (const uint8_t *buffer, ssize_t len);
      static SxString get (const SxString &string);
      static SxString file (const SxString &filename, uint64_t *size = NULL);

      // Converter functions
      // Recreate digest array from md5 string
      static SxArray<uint32_t> stringToDigest (const SxString &md5sum);
      // Recreate md5 string from digest array
      static SxString digestToString (const uint32_t *digest);

   private:

      // Md5 digest storage (128 bit / 16 byte)
      uint32_t digest[4];
};

#endif /* _SX_MD5_H_ */
