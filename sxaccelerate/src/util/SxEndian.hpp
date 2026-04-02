
/** Detect Big endian.
    \return true if the architecture has big endian. */
inline bool SxEndian::isBig ()
{
   // --- it is possible that endian is not big neither little but rare
   uint32_t test = 1;
   char *ptr=(char*)&test;
   return ptr[sizeof(uint32_t)-1] == 1; // LSB last
}

/** Detect Little endian.
    \return true if the architecture has little endian. */
inline bool SxEndian::isLittle ()
{
   uint32_t test = 1;
   char *ptr=(char*)&test;
   return ptr[0] == 1; // LSB first
}

/** Reverse last nBits_ bits in the input value.
    \return reversed bits. */
inline uint32_t SxEndian::bitReverse (const uint32_t value_, const int nBits_)
{
   int i;
   uint32_t result=0;

   for (i=0; i < nBits_; i++)  {
      result = (result << 1) | ((value_ >> i) & 0x01);
   }

   return result;
}

inline uint16_t SxEndian::swap16 (const uint16_t value_)
{
   return (uint16_t)(((value_ >> 8) & 0xff)
                   + ((value_ & 0xff) << 8));
}

inline uint32_t SxEndian::swap32 (const uint32_t value_)
{
   return (uint32_t)(((value_ >> 24) &   0xff)
                   + ((value_ >>  8) & 0xff00)
                   + ((value_ & 0xff00) <<  8)
                   + ((value_ &   0xff) << 24));
}

/** Swap LSB and MSB by 16-bits.
    \param dst_ Destination.
    \param src_ Source.
    \param nbytes_ The number of bytes. */
inline void SxEndian::swap16 (const uint8_t *src_,
                              uint8_t       *dst_,
                              const ssize_t nbytes_)
{
   ssize_t i;
   ssize_t n = nbytes_ >> 1;
   for (i=0; i < n; i++)  {
      dst_[(i << 1)    ] = src_[(i << 1) + 1];
      dst_[(i << 1) + 1] = src_[(i << 1)];
   }
}

/** Swap LSB and MSB by 32-bits.
    \param dst_ Destination.
    \param src_ Source.
    \param nbytes_ The number of bytes. */
inline void SxEndian::swap32 (const uint8_t *src_,
                              uint8_t       *dst_,
                              const ssize_t nbytes_)
{
   ssize_t i;
   ssize_t n = nbytes_ >> 2;
   for (i=0; i < n; i++)  {
      dst_[(i << 2)    ] = src_[(i << 2) + 3];
      dst_[(i << 2) + 1] = src_[(i << 2) + 2];
      dst_[(i << 2) + 2] = src_[(i << 2) + 1];
      dst_[(i << 2) + 3] = src_[(i << 2)    ];
   }
}

inline void SxEndian::swap64 (const uint8_t *src_,
                              uint8_t       *dst_,
                              const ssize_t nbytes_)
{
   ssize_t i;
   ssize_t n = nbytes_ >> 3;
   for (i=0; i < n; i++)  {
      dst_[(i << 3)    ] = src_[(i << 3) + 7];
      dst_[(i << 3) + 1] = src_[(i << 3) + 6];
      dst_[(i << 3) + 2] = src_[(i << 3) + 5];
      dst_[(i << 3) + 3] = src_[(i << 3) + 4];
      dst_[(i << 3) + 4] = src_[(i << 3) + 3];
      dst_[(i << 3) + 5] = src_[(i << 3) + 2];
      dst_[(i << 3) + 6] = src_[(i << 3) + 1];
      dst_[(i << 3) + 7] = src_[(i << 3)    ];
   }
}

// ----------------------------------------------------------------------------
// hardware dependent code, big versus little endian

/** Convert host endian to big endian.
    \param value_ Specifies the input value in host endian.
    \return swapped bytes in big endian. */
inline uint16_t SxEndian::hostToBig16 (const uint16_t value_)
{
#ifdef HAS_BIG_ENDIAN
   return value_;
#else
   return swap16 (value_);
#endif
}

/** Convert host endian to big endian.
    \param value_ Specifies the input value in host endian.
    \return swapped bytes in big endian. */
inline uint32_t SxEndian::hostToBig32 (const uint32_t value_)
{
#ifdef HAS_BIG_ENDIAN
   return value_;
#else
   return swap32 (value_);
#endif
}

/** Copy and convert host endian to big endian.
    16-bit integers do not need to be memory aligned.
    \param dst_ Destination.
    \param src_ Source.
    \param nbytes_ The number of bytes. */
inline void SxEndian::hostToBig16 (uint16_t src_, uint8_t *dst_)
{
   dst_[0] = static_cast<uint8_t>(src_ >> 8);
   dst_[1] = static_cast<uint8_t>(src_);
}

inline void SxEndian::hostToBig32 (uint32_t src_, uint8_t *dst_)
{
   dst_[0] = static_cast<uint8_t>(src_ >> 24);
   dst_[1] = static_cast<uint8_t>(src_ >> 16);
   dst_[2] = static_cast<uint8_t>(src_ >> 8);
   dst_[3] = static_cast<uint8_t>(src_);
}

inline void SxEndian::hostToBig64 (uint64_t src_, uint8_t *dst_)
{
   dst_[0] = static_cast<uint8_t>(src_ >> 56);
   dst_[1] = static_cast<uint8_t>(src_ >> 48);
   dst_[2] = static_cast<uint8_t>(src_ >> 40);
   dst_[3] = static_cast<uint8_t>(src_ >> 32);
   dst_[4] = static_cast<uint8_t>(src_ >> 24);
   dst_[5] = static_cast<uint8_t>(src_ >> 16);
   dst_[6] = static_cast<uint8_t>(src_ >> 8);
   dst_[7] = static_cast<uint8_t>(src_);
}

/** Copy and convert host endian to big endian.
    16-bit integers do not need to be memory aligned.
    \param dst_ Destination.
    \param src_ Source.
    \param nbytes_ The number of bytes. */
inline void SxEndian::hostToBig16 (const uint8_t *src_,
                                   uint8_t       *dst_,
                                   const ssize_t nbytes_)
{
#ifdef HAS_BIG_ENDIAN
   memcpy (dst_, src_, nbytes_);
#else
   swap16 (src_, dst_, nbytes_);
#endif
}

inline void SxEndian::hostToBig32 (const uint8_t *src_,
                                   uint8_t       *dst_,
                                   const ssize_t nbytes_)
{
#ifdef HAS_BIG_ENDIAN
   memcpy (dst_, src_, nbytes_);
#else
   swap32 (src_, dst_, nbytes_);
#endif
}

inline void SxEndian::hostToBig64 (const uint8_t *src_,
                                   uint8_t       *dst_,
                                   const ssize_t nbytes_)
{
#ifdef HAS_BIG_ENDIAN
   memcpy (dst_, src_, nbytes_);
#else
   swap64 (src_, dst_, nbytes_);
#endif
}

/** Convert big endian to host endian.
    \param value_ Specifies the input value in big endian.
    \return swapped bytes in host endian. */
inline uint16_t SxEndian::bigToHost16 (const uint16_t value_)
{
#ifdef HAS_BIG_ENDIAN
   return value_;
#else
   return swap16 (value_);
#endif
}

/** Convert big endian to host endian.
    \param value_ Specifies the input value in big endian.
    \return swapped bytes in host endian. */
inline uint32_t SxEndian::bigToHost32 (const uint32_t value_)
{
#ifdef HAS_BIG_ENDIAN
   return value_;
#else
   return swap32 (value_);
#endif
}


inline uint16_t SxEndian::bigToHost16 (const uint8_t *src_)
{
   // ((uint16_t)src_[0] << 8) | src_[1]; result is int?
   return (uint16_t)(((uint16_t)src_[0] << 8) | src_[1]);
}


inline uint32_t SxEndian::bigToHost32 (const uint8_t *src_)
{
   return ((uint32_t)src_[0] << 24)
        | ((uint32_t)src_[1] << 16)
        | ((uint32_t)src_[2] << 8)
        | src_[3];
}

inline uint64_t SxEndian::bigToHost64 (const uint8_t *src_)
{
   return ((uint64_t)src_[0] << 56)
        | ((uint64_t)src_[1] << 48)
        | ((uint64_t)src_[2] << 40)
        | ((uint64_t)src_[3] << 32)
        | ((uint64_t)src_[4] << 24)
        | ((uint64_t)src_[5] << 16)
        | ((uint64_t)src_[6] << 8)
        | src_[7];
}

/** Convert host endian to little endian.
    \param value_ Specifies the input value in host endian.
    \return swapped bytes in little endian. */
inline uint16_t SxEndian::hostToLittle16 (const uint16_t value_)
{
#ifdef HAS_BIG_ENDIAN
   return swap16 (value_);
#else
   return value_;
#endif
}

/** Convert host endian to little endian.
    \param value_ Specifies the input value in host endian.
    \return swapped bytes in little endian. */
inline uint32_t SxEndian::hostToLittle32 (const uint32_t value_)
{
#ifdef HAS_BIG_ENDIAN
   return swap32 (value_);
#else
   return value_;
#endif
}


inline void SxEndian::hostToLittle16 (uint16_t src_, uint8_t *dst_)
{
   dst_[0] = static_cast<uint8_t>(src_);
   dst_[1] = static_cast<uint8_t>(src_ >> 8);
}


inline void SxEndian::hostToLittle32 (uint32_t src_, uint8_t *dst_)
{
   dst_[0] = static_cast<uint8_t>(src_);
   dst_[1] = static_cast<uint8_t>(src_ >> 8);
   dst_[2] = static_cast<uint8_t>(src_ >> 16);
   dst_[3] = static_cast<uint8_t>(src_ >> 24);
}


inline void SxEndian::hostToLittle16 (const uint8_t *src_,
                                      uint8_t       *dst_,
                                      const ssize_t nbytes_)
{
#ifdef HAS_BIG_ENDIAN
   swap16 (src_, dst_, nbytes_);
#else
   memcpy (dst_, src_, (size_t)nbytes_);
#endif
}

/** Copy and convert host endian to little endian.
    32-bit integers do not need to be memory aligned.
    \param dst_ Destination.
    \param src_ Source.
    \param nbytes_ The number of bytes. */
inline void SxEndian::hostToLittle32 (const uint8_t *src_,
                                      uint8_t       *dst_,
                                      const ssize_t nbytes_)
{
#ifdef HAS_BIG_ENDIAN
   swap32 (src_, dst_, nbytes_);
#else
   memcpy (dst_, src_, (size_t)nbytes_);
#endif
}

/** Convert little endian to host endian.
    \param value_ Specifies the input value in little endian.
    \return swapped bytes in host endian. */
inline uint16_t SxEndian::littleToHost16 (const uint16_t value_)
{
#ifdef HAS_BIG_ENDIAN
   return swap16 (value_);
#else
   return value_;
#endif
}

/** Convert little endian to host endian.
    \param value_ Specifies the input value in little endian.
    \return swapped bytes in host endian. */
inline uint32_t SxEndian::littleToHost32 (const uint32_t value_)
{
#ifdef HAS_BIG_ENDIAN
   return swap32 (value_);
#else
   return value_;
#endif
}

inline uint16_t SxEndian::littleToHost16 (const uint8_t *src_)
{
   // ((uint16_t)src_[1] << 8) | src_[0]; result is int?
   return (uint16_t)(((uint16_t)src_[1] << 8) | src_[0]);
}

inline uint32_t SxEndian::littleToHost32 (const uint8_t *src_)
{
   return ((uint32_t)src_[3] << 24)
        | ((uint32_t)src_[2] << 16)
        | ((uint32_t)src_[1] << 8)
        | src_[0];
}
