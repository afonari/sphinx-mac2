#ifndef _SX_SERIALIZER_UTIL_H_
#define _SX_SERIALIZER_UTIL_H_

#include <SxSerializer.h>
#include <SxVariant.h>

namespace SxSerializer
{
   inline SxArray<char> pack (const SxVariant &var)
   {
      return var.pack ();
   }

   inline ssize_t unpack (const char *data, ssize_t len, SxVariant *var)
   {
      SX_CHECK (var);
      return var->unpack (data, len);
   }

   template<class T>
   SxString pack (const SxList<T> &in);
   template<class T>
   ssize_t unpack (const char *data, ssize_t len, SxList<T> *in);

   template<class K, class V>
   SxString pack (const SxMap<K,V> &in);
   template<class K, class V>
   ssize_t unpack (const char *data, ssize_t len, SxMap<K,V> *in);
}

#include <SxSerializerUtil.hpp>

#endif /* _SX_SERIALIZER_UTIL_H_ */
