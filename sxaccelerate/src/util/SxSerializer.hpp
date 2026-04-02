#include <SxException.h>

namespace SxSerializer
{

   inline Stream::Stream ()
      : size(0),
        pos(0)
   {
      // empty
   }

   inline Stream::Stream (const SxArray<char> &data_)
      : data(data_),
        size(data_.getSize ()),
        pos(0)
   {
      // empty
   }

   inline Stream::~Stream ()
   {
      // empty
   }

   inline const char *Stream::get () const
   {
      return data.elements + pos;
   }

   inline ssize_t Stream::getSize () const
   {
      return size;
   }

   inline bool Stream::seek (ssize_t offset)
   {
      if (offset < 0 || offset > size)  {
         return false;
      } else {
         pos += offset;
         size -= offset;
         return true;
      }
   }

   inline SxString binToHex (const SxString &str_)
   {
      SxString res;
      char token[3];
      const unsigned char *buf = (const unsigned char *)str_.elements;
      for (ssize_t i=0; i < str_.SxArray<char>::getSize (); ++i)  {
         sprintf (token, "%02x", buf[i]);
         res += token;
      }
      return res;
   }

   inline SxString hexToBin (const SxString &str_)
   {
      SxString res;
      res.resize (str_.SxArray<char>::getSize()/2);
      for (ssize_t i=0; i < str_.SxArray<char>::getSize()/2; ++i)  {
         SxString token = str_.subString (2*i, 2*i+1);
         unsigned int val = 0;
         int nElem = sscanf (token.ascii(), "%x", &val);
         if (nElem != 1)  {
            SX_THROW ("Invalid hex on ascii input");
         }
         res(i) = (char)val;
      }
      return res;
   }

   inline SxArray<char> pack (int64_t value)
   {
      SX_TRACE ();
      SxArray<char> res;

      uint8_t valueNBytes = getSize (value);
      res.resize (valueNBytes + (ssize_t)(2));
      uint8_t *buffer = (uint8_t*)res.elements;

      buffer[0] = valueNBytes;
      pack (value, buffer + 1);

      return res;
   }

   inline bool unpack (Stream *stream, ssize_t len, int64_t *value)
   {
      SX_TRACE ();
      SX_CHECK (stream);

      if (len < 2)
         return false;

      ssize_t nBytes = 0;
      int64_t data = 0;

      nBytes = unpack (stream->get (), len, &data);
      if (!stream->seek (nBytes))
         return false;

      if (value)
         *value = data;

      return true;
   }

   inline SxArray<char> pack (double value)
   {
      SX_TRACE ();
      SxArray<char> res;

      uint8_t valueNBytes = getSize (value);
      res.resize (valueNBytes + (ssize_t)(2));
      uint8_t *buffer = (uint8_t*)res.elements;

      buffer[0] = valueNBytes;
      pack (value, buffer + 1);

      return res;
   }

   inline bool unpack (Stream *stream, ssize_t len, double *value)
   {
      SX_TRACE ();
      SX_CHECK (stream);

      if (len < 2)
         return false;

      ssize_t nBytes = 0;
      double data = 0.f;

      nBytes = unpack (stream->get (), len, &data);
      if (!stream->seek (nBytes))
         return false;

      if (value)
         *value = data;

      return true;
   }

   inline SxArray<char> pack (int64_t major, int64_t minor)
   {
      SxArray<char> res;

      uint8_t majorLen = getSize (major);
      uint8_t minorLen = getSize (minor);
      res.resize ((ssize_t)(majorLen) + (ssize_t)(minorLen) + (ssize_t)(2));

      ssize_t pos = 0;
      uint8_t *buffer = (uint8_t*)res.elements;

      // --- major
      buffer[pos] = majorLen;
      pos += 1;
      pack (major, buffer + pos);
      pos += majorLen;

      // --- minor
      buffer[pos] = minorLen;
      pos += 1;
      pack (minor, buffer + pos);
      pos += minorLen;
      //buffer[pos] = static_cast<uint8_t>(0);
      return res;
   }

   inline bool unpack (Stream *stream, int64_t *major, int64_t *minor)
   {
      SX_CHECK (stream);
      if (stream->getSize () < 2)
         return false;

      ssize_t nBytes = 0;
      int64_t data = 0;

      // --- major
      nBytes = unpack (stream->get (), stream->getSize (), &data);
      if (!stream->seek (nBytes))
         return false;

      if (major)
         *major = data;//static_cast<int>(data);

      // --- minor
      nBytes = unpack (stream->get (), stream->getSize (), &data);
      if (!stream->seek (nBytes))
         return false;

      if (minor)
         *minor = data;//static_cast<int>(data);

      return true;
   }

   inline SxArray<char> pack (const SxString &str)
   {
      SxArray<char> res;
      ssize_t pos = 0;
      int64_t n = str.getSize ();
      uint8_t nLen = getSize(n);
      res.resize (n + nLen + 1);
      uint8_t *buffer = (uint8_t*)res.elements;
      // --- size
      //GxEndian::hostToBig64 (static_cast<uint64_t>(n), buffer);
      buffer[0] = nLen;
      pos += 1;
      pack (n, buffer + pos);
      pos += nLen;
      // --- data
      memcpy (buffer + pos, str.elements, (size_t)n);
      return res;
   }

   inline ssize_t unpack (const char *data, ssize_t len, SxString *str)
   {
      SX_CHECK (str);
      ssize_t pos = 0;
      if (len < 2)  {
         return -1;
      }
      SX_CHECK (data);
      //const uint8_t *buffer = (const uint8_t*)(data);
      // --- size
      int64_t n = 0;
      pos = unpack (data, len, &n);
      if (pos < 0)  {
         return -1;
      }
      // --- data
      if (len < pos + n)  {
         return -1;
      }
      *str = SxString (data + pos, n);
      pos += n;
      return pos;
   }


   inline uint8_t getSize (int64_t value)
   {
      if (value <= UINT8_MAX)        return 1;
      else if (value <= UINT16_MAX)  return 2;
      else if (value <= UINT32_MAX)  return 4;
      else                           return 8;
   }

   inline void pack (int64_t value, uint8_t *buffer)
   {
      SX_CHECK (buffer);
      if (value <= UINT8_MAX)  {
         buffer[0] = static_cast<uint8_t>(value);
      } else if (value <= UINT16_MAX)  {
         SxEndian::hostToBig16 (static_cast<uint16_t>(value), buffer);
      } else if (value <= UINT32_MAX)  {
         SxEndian::hostToBig32 (static_cast<uint32_t>(value), buffer);
      } else  {
         SxEndian::hostToBig64 (static_cast<uint64_t>(value), buffer);
      }
   }

   inline ssize_t unpack (const char *data, ssize_t len, int64_t *n)
   {
      SX_CHECK (n);
      ssize_t pos = 0;
      if (len < 2)  {
         return -1;
      }
      SX_CHECK (data);
      const uint8_t *buffer = (const uint8_t*)data;
      ssize_t nLen = buffer[0];
      pos += 1;
      if (len < pos + nLen)  {
         return -1;
      }
      uint64_t size = 0;
      if (nLen == 1)  {
         size = buffer[pos];
      } else if (nLen == 2)  {
         size = SxEndian::bigToHost16 (buffer + pos);
      } else if (nLen == 4)  {
         size = SxEndian::bigToHost32 (buffer + pos);
      } else if (nLen == 8)  {
         size = SxEndian::bigToHost64 (buffer + pos);
      }
      *n = static_cast<int64_t>(size);
      pos += nLen;
      return pos;
   }

   inline uint8_t getSize (double value)
   {
      SX_UNUSED (value);
      return 8;
   }

   inline void pack (double value, uint8_t *buffer)
   {
      SX_CHECK (buffer);
      SxEndian::hostToBig64 (((const uint8_t *)(&value)), buffer, 8);
   }

   inline ssize_t unpack (const char *data, ssize_t len, double *n)
   {
      SX_CHECK (n);
      ssize_t pos = 0;
      if (len < 2)  {
         return -1;
      }
      SX_CHECK (data);
      const uint8_t *buffer = (const uint8_t*)data;
      ssize_t nLen = buffer[0];
      pos += 1;
      if (len < pos + nLen)  {
         return -1;
      }
      uint64_t size = 0;
      if (nLen == 1)  {
         size = buffer[pos];
      } else if (nLen == 2)  {
         size = SxEndian::bigToHost16 (buffer + pos);
      } else if (nLen == 4)  {
         size = SxEndian::bigToHost32 (buffer + pos);
      } else if (nLen == 8)  {
         size = SxEndian::bigToHost64 (buffer + pos);
      }
      void *tmp = (void *)&size;
      *n = *((double *)(tmp));
      pos += nLen;
      return pos;
   }

}
