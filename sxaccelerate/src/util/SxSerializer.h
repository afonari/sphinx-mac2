#ifndef _SX_SERIALIZER_H_
#define _SX_SERIALIZER_H_

#include <SxTypeDefs.h>
#include <SxString.h>
#include <SxEndian.h>

namespace SxSerializer
{
   class Stream
   {
      public:
         Stream ();
         Stream (const SxArray<char> &data);
         ~Stream ();

         const char *get () const;
         ssize_t getSize () const;

         bool seek (ssize_t offset);

      protected:
         SxArray<char> data;
         ssize_t size;
         ssize_t pos;
   };

   SxString binToHex (const SxString &);
   SxString hexToBin (const SxString &);

   template<class T>
   SxString pack (const T &)
   {
      SX_EXIT; // not implemented
      return SxString ();
   }

   SxArray<char> pack (int64_t value);
   bool unpack (Stream *stream, ssize_t len, int64_t *value);

   SxArray<char> pack (double value);
   bool unpack (Stream *stream, ssize_t len, double *value);

   SxArray<char> pack (int64_t major, int64_t minor);
   bool unpack (Stream *stream, int64_t *major, int64_t *minor);

   SxArray<char> pack (const SxString &str);
   ssize_t unpack (const char *data, ssize_t len, SxString *str);

   uint8_t getSize (int64_t value);
   void pack (int64_t value, uint8_t *buffer);
   ssize_t unpack (const char *data, ssize_t len, int64_t *n);

   uint8_t getSize (double value);
   void pack (double value, uint8_t *buffer);
   ssize_t unpack (const char *data, ssize_t len, double *n);

}

#include <SxSerializer.hpp>

#endif /* _SX_SERIALIZER_H_ */
