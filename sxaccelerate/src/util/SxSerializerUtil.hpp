
namespace SxSerializer
{

   template<class T>
   inline SxString pack (const SxList<T> &in)
   {
      ssize_t pos = 0;
      SxList<SxString> list;
      typename SxList<T>::ConstIterator it;
      for (it = in.begin(); it != in.end(); ++it)  {
         list << pack (*it);
         pos += list.last().getSize ();
      }
      int64_t n = list.getSize ();
      uint8_t nLen = getSize (n);

      SxString res;
      res.SxArray<char>::resize (pos + nLen + 2);

      pos = 0;
      uint8_t *buffer = (uint8_t*)res.elements;
      // --- size
      buffer[0] = static_cast<uint8_t>(nLen);
      pos += 1;
      pack (n, buffer + pos);
      pos += nLen;

      // --- data
      typename SxList<SxString>::ConstIterator it2;
      for (it2 = list.begin(); it2 != list.end(); ++it2)  {
         memcpy (buffer + pos, it2->elements, (size_t)it2->getSize ());
         pos += it2->getSize ();
      }

      return res;
   }

   template<class T>
   inline ssize_t unpack (const char *data, ssize_t len, SxList<T> *list)
   {
      SX_CHECK (list);
      ssize_t pos = 0;
      if (len < 2)  {
         return -1;
      }
      SX_CHECK (data);

      // --- size
      int64_t n = 0;
      pos = unpack (data, len, &n);
      if (pos < 0)  {
         return -1;
      }
      // --- data
      for (ssize_t i=0; i < n; ++i)  {
         T item;
         ssize_t nBytes = unpack (data + pos, len - pos, &item);
         if (nBytes < 0)  {
            return -1;
         }
         list->append (item);
         pos += nBytes;
      }

      return pos;
   }

   template<class K, class V>
   inline SxString pack (const SxMap<K,V> &in)
   {
      ssize_t pos = 0;
      SxList<SxArray<char>> list;
      typename SxMap<K,V >::ConstIterator m;
      for (m = in.begin(); m != in.end(); ++m)  {
         list << pack (m.getKey ());
         pos += list.last().getSize ();
         list << pack (m.getValue ());
         pos += list.last().getSize ();
      }
      int64_t n = in.getSize ();
      uint8_t nLen = getSize (n);
      SxString res;
      res.SxArray<char>::resize (pos + nLen + 2);

      pos = 0;
      uint8_t *buffer = (uint8_t*)res.elements;
      // --- size
      buffer[0] = static_cast<uint8_t>(nLen);
      pos += 1;
      pack (n, buffer + pos);
      pos += nLen;

      // --- data
      typename SxList<SxArray<char>>::ConstIterator it2;
      for (it2 = list.begin(); it2 != list.end(); ++it2)  {
         memcpy (buffer + pos, it2->elements, (size_t)it2->getSize ());
         pos += it2->getSize ();
      }
      //buffer[pos] = static_cast<uint8_t>(0);
      return res;
   }

   template<class K, class V>
   inline ssize_t unpack (const char *data, ssize_t len, SxMap<K,V> *in)
   {
      SX_CHECK (in);
      ssize_t pos = 0;
      if (len < 2)  {
         return -1;
      }
      SX_CHECK (data);
      // --- size
      int64_t n = 0;
      pos = unpack (data, len, &n);
      if (pos < 0)  {
         return -1;
      }
      // --- data
      for (ssize_t i=0; i < n; ++i)  {
         // --- key
         K key;
         ssize_t nBytes = unpack (data + pos, len - pos, &key);
         if (nBytes < 0)  {
           return -1;
         }
         pos += nBytes;
         // --- value
         V val;
         nBytes = unpack (data + pos, len - pos, &val);
         if (nBytes < 0)  {
            return -1;
         }
         pos += nBytes;
         (*in)(key) = val;
      }
      return pos;
   }

}
