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

#ifndef _SX_EXCEPTION_BASE_H_
#define _SX_EXCEPTION_BASE_H_

#include <SxUtil.h>
#include <SxMacroLib.h>
#include <SxString.h>
#include <SxVariant.h>
#include <SxMap.h>
#include <atomic>

extern atomic<bool> dumpCoreFile;

/** \brief SxAccelerate Exception class

    \b SxException = SxAccelerate Exception Handling Class

    This is a standard exception class for all SxAccelerate exceptions.
    Exception handling in C++ can cause spurious bus errors when they are
    not caught properly in a try-catch block. If SxAccelerate aborts with a bus
    error which could be related to such a missing try-catch block this class
    allows to handle the exception with an SxAccelerate SX_EXIT statement
    instead. In the DEBUG mode the resulting core file can be analyzed
    (\ref dbg_analysis).

    \author Sixten Boeck
  */
class SX_EXPORT_UTIL SxException
{
   public:

      // print mode for toString()
      enum PrintMode {
         Simple=0x00,      ///< top-most exception w/o sender info
         Stack=0x01,       ///< exception stack w/o sender info
         DebugSimple=0x02, ///< top-most exception w/ sender info
         DebugStack=0x03,  ///< exception stack w/ sender info
      };

      /// \brief Standard constructor
      SxException ();
      /// \brief Copy constructor
      SxException (const SxException &);
      SxException (const SxString &message,
                   const SxString &senderFile="unknown",
                   int senderLine=0);
      SxException (const SxString &tag,
                   const SxString &message,
                   const SxString &file="unknown", int senderLine=0);
      SxException (const SxException &chain, const SxString &message,
                   const SxString &senderFile="unknown", int senderLine=0);
      SxException (const SxException &chain,
                   const SxString &tag, const SxString &message,
                   const SxString &senderFile="unknown", int senderLine=0);

      // constructors for throwing with specific 'type'
      SxException (const SxString &tag,    uint32_t tagHash,
                   const SxString &subTag, uint32_t subtagHash,
                   const SxList<SxVariant> &args,
                   const SxString &senderFile="unknown", int senderLine=0);
      SxException (const SxException &chain,
                   const SxString &tag,    uint32_t tagHash,
                   const SxString &subTag, uint32_t subtagHash,
                   const SxList<SxVariant> &args,
                   const SxString &senderFile="unknown", int senderLine=0);
      /// \brief Destructor
      ~SxException ();

      /// \brief assignment operator
      SxException &operator= (const SxException &);

      /** \brief Print the message body

          if \em printSender is true additionally the file and the line
          from where this exception was trown will be printed. */
      void print (bool printSender=false) const;

      SxString toString (int mode = (int)Simple,
                         const SxString &delimiter="\n") const;

      void printStack (const SxString &delimiter="\n") const;

      bool hasTag (const SxString &tag_) const;   // current instance only
      bool findTag (const SxString &tag_) const;  // search in chain recusively


      /// \brief returns the message body.
      const SxString &getMessage () const  { return message; }

      /// \brief returns the sender file.
      const SxString &getFileName () const  { return senderFile; }

      const SxString &getTag () const { return tag; }
      const SxString &getSubTag () const { return subTag; }

      template<uint32_t tagHash, uint32_t subtagHash>
      bool is () const;
      template<uint32_t tagHash, uint32_t subtagHash>
      bool has () const;
      template<uint32_t tagHash> bool isCategory () const;
      template<uint32_t tagHash> bool hasCategory () const;
      template<uint32_t subtagHash> bool isTag () const;
      template<uint32_t subtagHash> bool hasTag () const;

      /// \brief returns the sender line.
      int getLineNumber () const { return senderLine; }

      /// \brief returns next exception in the chain.
      const SxException *getNext () const { return next; }


      /** \brief Treat all exception as SX_EXIT statements.

          Usually an uncaught exception causes an "Abort" signal. So no
          core file will be prompted and backtracing in the debug mode
          becomes impossible.
          By calling this static function every exception violates the
          memory explicily in order to dump a qualified core file.
          For debugging only
          \sa SxException::throwExceptions
          \sa SxException::dumpCoreFile */
      static void causeSegFault   ();
      /** \brief Treat all exceptions normally (no SX_EXIT emulation)

          For debugging only. Switch off dumping of core files rather
          than thwing exceptions.
          \sa SxException::causeSegFault
          \sa SxException::dumpCodeFile */
      static void throwExceptions ();

      /** Program global variable.

          If SxException::dumpCoreFile is set to true all exception are
          treated with a SFHIngX SX_EXIT statement instead. This allows the
          developer to trace missing try-catch blocks.

          This variable shouldn't be set directly. Use better
          SxException::causeSegFaults or SxException::throwExceptions.

          \code
             SxException::dumpSegFaults();
          \endcode

          \sa \ref page_debug  */
      static atomic<bool> dumpCoreFile;


      static SxMap<SxString, SxMap<SxString,SxList<SxString> > > & getExceptionTable ();

      inline const SxVariant& getArg (ssize_t i) const { return arguments(i); }

   protected:

      /// \brief The exception message body
      SxString message;

      /// \brief Optional, exception key, such as Permission, BrokenPipe
      SxString tag;
      SxString subTag;

      /// \brief Hash of the tag. Hashed with SxHashFunction::jenkins
      uint32_t tagHash;
      uint32_t subtagHash;

      SxList<SxVariant> arguments;

      /// \brief Contains the file name from where the exception was trown.
      SxString senderFile;

      /// \brief Contains the line from where the exception was trown.
      int senderLine;

      /// \brief Next exception in the chain
      SxException *next;

      /// \brief previous exception in the chain
      SxException *prev;

      /// \brief last exception in the chain
      SxException *tail;

      void init (const SxException &);
      void init (const SxString &message,
                 const SxString &tag,
                 const SxString &senderFile);

      // init with 'type'
      void init (const SxString &tag,
                 const SxString &subtag,
                 const SxList<SxVariant> &args,
                 const SxString &senderFile);

      void destroy ();
      void copy (const SxException &);
      bool foundTagHash (uint32_t hashTag) const;
      bool foundSubtagHash (uint32_t subtagHash) const;
};

struct SX_EXPORT_UTIL SxExceptionEntry {
   explicit SxExceptionEntry (const char *,
                              const char *,
                              const SxList<SxString>&);
};

struct SxTrueType  { static constexpr bool value = true;  };
struct SxFalseType { static constexpr bool value = false; };

template<uint32_t tagHash>
struct SX_EXPORT_UTIL SxExceptionTag : SxFalseType { };


template<uint32_t tagHash, uint32_t subtagHash>
struct SX_EXPORT_UTIL SxExceptionType : SxFalseType { };

template<uint32_t tagHash_, uint32_t subtagHash_>
bool SxException::is () const
{
   static_assert (SxExceptionType<tagHash_,subtagHash_>::value,
                 "Undefined exception type");
   return (tagHash_ == tagHash && subtagHash_ == subtagHash);
}

template<uint32_t tagHash_, uint32_t subtagHash_>
bool SxException::has () const
{
   static_assert (SxExceptionType<tagHash_,subtagHash_>::value,
                 "Undefined exception type");
   return (foundTagHash(tagHash_) && foundSubtagHash(subtagHash_));
}

template<uint32_t tagHash_>
bool SxException::isCategory () const
{
   static_assert (SxExceptionTag<tagHash_>::value,
                  "Undefined Exception Catergory");
   return (tagHash_ == tagHash);
}

template<uint32_t tagHash_>
bool SxException::hasCategory () const
{
   static_assert (SxExceptionTag<tagHash_>::value,
                  "Undefined Exception Catergory");
   return foundTagHash (tagHash_);
}

template<uint32_t subtagHash_>
bool SxException::isTag () const
{
   static_assert (SxExceptionTag<subtagHash_>::value, "Undefined Exception Tag");
   return (subtagHash_ == subtagHash);
}

template<uint32_t subtagHash_>
bool SxException::hasTag () const
{
   static_assert (SxExceptionTag<subtagHash_>::value, "Undefined Exception Tag");
   return foundSubtagHash (subtagHash_);
}

// SX_THROW

#define _SxThrow1(msg)                                                         \
   throw SxException(msg, __FILE__, __LINE__)

#define _SxThrow2(var,msg)                                                     \
   throw SxException(var, msg, __FILE__, __LINE__)

#define _SxThrow3(tag,subtag,arg0)                                             \
   static_assert (SxExceptionTag<tag ""_SX>::value, "Unknown Tag: " tag);      \
   static_assert (SxExceptionType<tag ""_SX, subtag ""_SX>::value,             \
                  "Unknown tag and subtag combination: [" tag "," subtag "]"); \
   throw SxException(tag, tag ""_SX, subtag, subtag ""_SX,                     \
                     SxExceptionType<tag ""_SX, subtag ""_SX>()(arg0),         \
                     __FILE__, __LINE__)

#define _SxThrow4(tag,subtag,arg0,arg1)                                        \
   static_assert (SxExceptionTag<tag ""_SX>::value, "Unknown Tag: " tag);      \
   static_assert (SxExceptionType<tag ""_SX, subtag ""_SX>::value,             \
                  "Unknown tag and subtag combination: [" tag "," subtag "]"); \
   throw SxException(tag, tag ""_SX, subtag, subtag ""_SX,                     \
                     SxExceptionType<tag ""_SX, subtag ""_SX>()(arg0,arg1),    \
                     __FILE__, __LINE__)

#define _SxThrow5(tag,subtag,arg0,arg1,arg2)                                   \
   static_assert (SxExceptionTag<tag ""_SX>::value, "Unknown Tag: " tag);      \
   static_assert (SxExceptionType<tag ""_SX, subtag ""_SX>::value,             \
                  "Unknown tag and subtag combination: [" tag "," subtag "]"); \
   throw SxException(tag, tag ""_SX, subtag, subtag ""_SX,                     \
                     SxExceptionType<tag ""_SX, subtag ""_SX>()(arg0,arg1,arg2)\
                     ,__FILE__, __LINE__)

#ifdef MSVC
// disable integral constant overflow warning
#  define SX_THROW(...)                                                       \
      __pragma(warning (suppress:4307))                                       \
      SX_VMACRO(_SxThrow,__VA_ARGS__)
#else
#  define SX_THROW(...) SX_VMACRO(_SxThrow,__VA_ARGS__)
#endif

// SX_RETHROW

#define _SxReThrow2(e,msg)                          \
   throw SxException(e,msg, __FILE__,__LINE__)

#define _SxReThrow3(e,tag,msg)                      \
   throw SxException(e,tag,msg, __FILE__,__LINE__)

#define _SxReThrow4(e,tag,subtag,arg0)                                         \
   static_assert (SxExceptionTag<tag ""_SX>::value, "Unknown Tag: " tag);      \
   static_assert (SxExceptionType<tag ""_SX, subtag ""_SX>::value,             \
                  "Unknown tag and subtag combination: [" tag "," subtag "]"); \
   throw SxException(e, tag, tag ""_SX, subtag, subtag ""_SX,                  \
                     SxExceptionType<tag ""_SX, subtag ""_SX>()(arg0),         \
                     __FILE__, __LINE__)

#define _SxReThrow5(e,tag,subtag,arg0,arg1)                                    \
   static_assert (SxExceptionTag<tag ""_SX>::value, "Unknown Tag: " tag);      \
   static_assert (SxExceptionType<tag ""_SX, subtag ""_SX>::value,             \
                  "Unknown tag and subtag combination: [" tag "," subtag "]"); \
   throw SxException(e, tag, tag ""_SX, subtag, subtag ""_SX,                  \
                     SxExceptionType<tag ""_SX, subtag ""_SX>()(arg0,arg1),    \
                     __FILE__, __LINE__)

#define _SxReThrow6(e,tag,subtag,arg0,arg1,arg2)                               \
   static_assert (SxExceptionTag<tag ""_SX>::value, "Unknown Tag: " tag);      \
   static_assert (SxExceptionType<tag ""_SX, subtag ""_SX>::value,             \
                  "Unknown tag and subtag combination: [" tag "," subtag "]"); \
   throw SxException(e, tag, tag ""_SX, subtag, subtag ""_SX,                  \
                     SxExceptionType<tag ""_SX, subtag ""_SX>()(arg0,arg1,arg2)\
                     ,__FILE__, __LINE__)
#ifdef MSVC
// disable integral constant overflow warning
#  define SX_RETHROW(...)                                                     \
      __pragma(warning (suppress:4307))                                       \
      SX_VMACRO(_SxReThrow,__VA_ARGS__)
#else
#  define SX_RETHROW(...) SX_VMACRO(_SxReThrow,__VA_ARGS__)
#endif

#ifdef MSVC
// disable integral constant overflow warning
#  define SX_EXCEPTION_TAG(tag)                                               \
      static_assert (sizeof (tag) != 0, "Incorrect Tag, size == 0");          \
      __pragma(warning (suppress:4307))                                       \
      template<> struct SxExceptionTag<tag ""_SX> : SxTrueType {}
#else
#  define SX_EXCEPTION_TAG(tag)                                               \
      static_assert (sizeof (tag) != 0, "Incorrect Tag, size == 0");          \
      template<> struct SxExceptionTag<tag ""_SX> : SxTrueType {}
#endif

#ifdef MSVC
// disable integral constant overflow warning
#  define SX_EXCEPTION_CAT(tag)                                               \
      static_assert (sizeof (tag) != 0, "Incorrect Category, size == 0");     \
      __pragma(warning (suppress:4307))                                       \
      template<> struct SxExceptionTag<tag ""_SX> : SxTrueType {}
#else
#  define SX_EXCEPTION_CAT(tag)                                                \
      static_assert (sizeof (tag) != 0, "Incorrect Category, size == 0");      \
      template<> struct SxExceptionTag<tag ""_SX> : SxTrueType {}
#endif

#define _SxThrowType4(tag,subtag,argName0,ArgType0)                            \
   static_assert (SxExceptionTag<tag ""_SX>::value, "Unknown Category: " tag); \
   static_assert (SxExceptionTag<subtag ""_SX>::value, "Unknown Tag: " tag);   \
   static SxExceptionEntry SX_UNIQUE_ID(entry) (                               \
      tag,subtag,{SxString(argName0) + " : " #ArgType0});                      \
   template<> struct SxExceptionType<tag ""_SX, subtag ""_SX> : SxTrueType     \
   {                                                                           \
      SxList<SxVariant> operator() (const ArgType0 &t0) {                      \
         return SxList<SxVariant>{ t0 };                                       \
      }                                                                        \
   }

#define _SxThrowType6(tag,subtag,argName0,ArgType0,argName1,ArgType1)          \
   static_assert (SxExceptionTag<tag ""_SX>::value, "Unknown Category: " tag); \
   static_assert (SxExceptionTag<subtag ""_SX>::value, "Unknown Tag: " tag);   \
   static SxExceptionEntry SX_UNIQUE_ID(entry) (                               \
      tag,subtag,{SxString(argName0) + " : " #ArgType0,                        \
                  SxString(argName1) + " : " #ArgType1});                      \
   template<> struct SxExceptionType<tag ""_SX, subtag ""_SX> : SxTrueType     \
   {                                                                           \
      SxList<SxVariant> operator() (const ArgType0 &t0,                        \
                                    const ArgType1 &t1)                        \
      {                                                                        \
         return SxList<SxVariant>{ t0, t1 };                                   \
      }                                                                        \
   }

#define _SxThrowType8(tag,subtag,                                              \
                      argName0,ArgType0,                                       \
                      argName1,ArgType1,                                       \
                      argName2,ArgType2)                                       \
                                                                               \
   static_assert (SxExceptionTag<tag ""_SX>::value, "Unknown Category: " tag); \
   static_assert (SxExceptionTag<subtag ""_SX>::value, "Unknown Tag: " tag);   \
   static SxExceptionEntry SX_UNIQUE_ID(entry) (                               \
      tag,subtag,{SxString(argName0) + " : " #ArgType0,                        \
                  SxString(argName1) + " : " #ArgType1,                        \
                  SxString(argName2) + " : " #ArgType2});                      \
   template<> struct SxExceptionType<tag ""_SX, subtag ""_SX> : SxTrueType     \
   {                                                                           \
      SxList<SxVariant> operator() (const ArgType0 &t0,                        \
                                    const ArgType1 &t1,                        \
                                    const ArgType2 &t2)                        \
      {                                                                        \
         return SxList<SxVariant>{ t0, t1, t2 };                               \
      }                                                                        \
   }

#ifdef MSVC
// disable integral constant overflow warning
#  define SX_EXCEPTION_TYPE(...)                                              \
      __pragma(warning (suppress:4307))                                       \
      SX_VMACRO_DECL(_SxThrowType,__VA_ARGS__)
#else
#  define SX_EXCEPTION_TYPE(...) SX_VMACRO_DECL(_SxThrowType,__VA_ARGS__)
#endif

#endif /* _SX_EXCEPTION_BASE_H_ */
