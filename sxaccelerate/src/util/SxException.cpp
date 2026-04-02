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


#include <SxException.h>
#include <SxList.h>
#include <SxError.h>
#include <SxConfig.h>
#include <SxHashFunction.h>

atomic<bool> SxException::dumpCoreFile(false);

SxException::SxException ()
   : message(),
     tag(),
     subTag(),
     tagHash(0),
     subtagHash(0),
     senderFile (),
     senderLine (-1),
     next(NULL),
     prev(NULL),
     tail(NULL)
{
   tail = this;
}

SxException::SxException (const SxException &in)
   : message(),
     tag(),
     subTag(),
     tagHash(0),
     subtagHash(0),
     senderFile (),
     senderLine (in.senderLine),
     next(NULL),
     prev(NULL),
     tail(NULL)
{
   try {
      set_new_handler(NULL);
      tail = this;
      copy (in);
   } catch (std::bad_alloc &excep) {
      SX_UNUSED(excep);
      in.print (true);
      sxOutOfMemoryHandler ();
   }
   set_new_handler(sxOutOfMemoryHandler);

   if (dumpCoreFile)  {
      print ();
#     ifndef SX_ANDROID
         fflush (stdout);
         fflush (stderr);
#     endif
      SX_BREAK;
   }
}


SxException::SxException (const SxString &message_,
                          const SxString &senderFile_, int senderLine_)
   : message(),
     tag(),
     subTag(),
     tagHash(0),
     subtagHash(0),
     senderFile (),
     senderLine (senderLine_),
     next(NULL),
     prev(NULL),
     tail(NULL)
{
   SX_DBG_MSG ("throwing (" << message_ << ")" << " from "
                            << senderFile_ << ":" << senderLine_);
   try {
      set_new_handler(NULL);
      tail = this;
      init (message_, "", senderFile_);
   } catch (std::bad_alloc &excep) {
      SX_UNUSED(excep);
      std::cerr << std::endl << message_ << std::endl;
      std::cerr << "This Exception was emitted from "
                << senderFile_
                << " (line " << senderLine_ << ")\n";
      sxOutOfMemoryHandler ();
   }
   set_new_handler(sxOutOfMemoryHandler);

   if (dumpCoreFile)  {
      print ();
#     ifndef SX_ANDROID
         fflush (stdout);
         fflush (stderr);
#     endif
      SX_BREAK;
   }
}


SxException::SxException (const SxString &tag_, const SxString &message_,
                          const SxString &senderFile_, int senderLine_)
   : message(),
     tag(),
     subTag(),
     tagHash(0),
     subtagHash(0),
     senderFile (),
     senderLine (senderLine_),
     next(NULL),
     prev(NULL),
     tail(NULL)
{
   SX_DBG_MSG ("throwing (" << message_ << ")" << " from "
                            << senderFile_ << ":" << senderLine_);
   try {
      set_new_handler(NULL);
      tail = this;
      init (message_, tag_, senderFile_);
   } catch (std::bad_alloc &excep) {
      SX_UNUSED(excep);
      std::cerr << std::endl << message_ << std::endl;
      std::cerr << "This Exception was emitted from "
                << senderFile_
                << " (line " << senderLine_ << ")\n";
      sxOutOfMemoryHandler ();
   }
   set_new_handler(sxOutOfMemoryHandler);

   if (dumpCoreFile)  {
      print ();
#     ifndef SX_ANDROID
         fflush (stdout);
         fflush (stderr);
#     endif
      SX_BREAK;
   }
}



SxException::SxException (const SxException &chain_,
                          const SxString &message_,
                          const SxString &senderFile_, int senderLine_)
   : message(),
     tag(),
     subTag(),
     tagHash(0),
     subtagHash(0),
     senderFile (),
     senderLine (senderLine_),
     next(NULL),
     prev(NULL),
     tail(NULL)
{
   SX_DBG_MSG ("Msg: " << message_ << " from "
                       << senderFile_ << ":" << senderLine_);

   try {
      set_new_handler(NULL);
      next = new SxException();
      next->copy (chain_);
      init (message_, "", senderFile_);
      this->tail = next->tail;
   } catch (std::bad_alloc &excep) {
      SX_UNUSED(excep);
      cerr << "Exception Stack: \n";
      const SxException *ptr = chain_.tail;
      int counter = 1;
      int i = 0;
      int strSize = 0;
      while (ptr) {
         cerr << counter++ << ": ";
         if (ptr->tag.getSize () > 0)
            cerr << "[" << ptr->tag << ", " << ptr->subTag << "] ";
         if (ptr->message.getSize () > 0) cerr << ptr->message << " ";
         strSize = (int)ptr->senderFile.getSize ();
         i = strSize - 1;
         while (i >= 0 && ptr->senderFile(i) != '\\' &&
                ptr->senderFile(i) != '/') i--;
         cerr << "(";
         while (((ssize_t)i+1) < strSize)  { cerr << ptr->senderFile((ssize_t)i+1); i++; }
         cerr << ":" << ptr->senderLine << ")";
         ssize_t argSize = ptr->arguments.getSize ();
         if (argSize > 0) {
            cerr << " ";
            for (const SxVariant &v : ptr->arguments)
               cerr << v << " ";
         }
         cerr << endl;
         ptr = ptr->prev;
      }
      cerr << counter++ << ": ";
      if (message_.getSize () > 0) cerr << message_ << " ";
      strSize = (int)senderFile_.getSize ();
      i = strSize - 1;
      while (i >= 0 && senderFile_(i) != '\\' && senderFile_(i) != '/') i--;
      cerr << "(";
      while (((ssize_t)i+1) < strSize)  { cerr << senderFile_((ssize_t)i+1); i++; }
      cerr << ":" << senderLine_ << ")" << endl;
      sxOutOfMemoryHandler ();
   }
   set_new_handler(sxOutOfMemoryHandler);

   if (dumpCoreFile)  {
      print ();
#     ifndef SX_ANDROID
         fflush (stdout);
         fflush (stderr);
#     endif
      SX_BREAK;
   }
}

SxException::SxException (const SxException &chain_,
                          const SxString &tag_, const SxString &message_,
                          const SxString &senderFile_, int senderLine_)
   : message(),
     tag(),
     subTag(),
     tagHash(0),
     subtagHash(0),
     senderFile (),
     senderLine (senderLine_),
     next(NULL),
     prev(NULL),
     tail(NULL)
{
   SX_DBG_MSG ("Msg: " << message_ << " from "
                       << senderFile_ << ":" << senderLine_);
   try {
      set_new_handler(NULL);
      next = new SxException();
      next->copy (chain_);
      init (message_, tag_, senderFile_);
      this->tail = next->tail;
   } catch (std::bad_alloc &excep) {
      SX_UNUSED(excep);
      cerr << "Exception Stack: \n";
      const SxException *ptr = chain_.tail;
      int counter = 1;
      int i = 0;
      int strSize = 0;
      while (ptr) {
         cerr << counter++ << ": ";
         if (ptr->tag.getSize () > 0)
            cerr << "[" << ptr->tag << ", " << ptr->subTag << "] ";
         if (ptr->message.getSize () > 0) cerr << ptr->message << " ";
         strSize = (int)ptr->senderFile.getSize ();
         i = strSize - 1;
         while (i >= 0 && ptr->senderFile(i) != '\\' &&
                ptr->senderFile(i) != '/') i--;
         cerr << "(";
         while (((ssize_t)i+1) < strSize)  { cerr << ptr->senderFile((ssize_t)i+1); i++; }
         cerr << ":" << ptr->senderLine << ")";
         ssize_t argSize = ptr->arguments.getSize ();
         if (argSize > 0) {
            cerr << " ";
            for (const SxVariant &v : ptr->arguments)
               cerr << v << " ";
         }
         cerr << endl;
         ptr = ptr->prev;
      }
      cerr << counter++ << ": ";
      if (tag_.getSize () > 0)  cerr << "[" << tag_ << "] ";
      if (message_.getSize () > 0) cerr << message_ << " ";
      strSize = (int)senderFile_.getSize ();
      i = strSize - 1;
      while (i >= 0 && senderFile_(i) != '\\' && senderFile_(i) != '/') i--;
      cerr << "(";
      while (((ssize_t)i+1) < strSize)  { cerr << senderFile_((ssize_t)i+1); i++; }
      cerr << ":" << senderLine_ << ")" << endl;

      sxOutOfMemoryHandler ();
   }
   set_new_handler(sxOutOfMemoryHandler);

   if (dumpCoreFile)  {
      print ();
#     ifndef SX_ANDROID
         fflush (stdout);
         fflush (stderr);
#     endif
      SX_BREAK;
   }
}

SxException::SxException (const SxString &tag_,
                          uint32_t tagHash_,
                          const SxString &subTag_,
                          uint32_t subtagHash_,
                          const SxList<SxVariant> &args,
                          const SxString &senderFile_,
                          int senderLine_)
   : message(),
     tag(),
     subTag(),
     tagHash(tagHash_),
     subtagHash(subtagHash_),
     senderFile (),
     senderLine (senderLine_),
     next(NULL),
     prev(NULL),
     tail(NULL)
{
   SX_DBG_MSG ("throwing (" << tag_ << "," << subTag_<< ")" << " from "
                            << senderFile_ << ":" << senderLine_);
   try {
      set_new_handler(NULL);
      tail = this;
      init (tag_, subTag_, args, senderFile_);
   } catch (std::bad_alloc &excep) {
      SX_UNUSED(excep);
      std::cerr << std::endl << tag_ << "," << subTag_ << std::endl;
      std::cerr << args << std::endl;
      std::cerr << "This Exception was emitted from "
                << senderFile_
                << " (line " << senderLine_ << ")\n";
      sxOutOfMemoryHandler ();
   }
   set_new_handler(sxOutOfMemoryHandler);

   if (dumpCoreFile)  {
      print ();
#     ifndef SX_ANDROID
         fflush (stdout);
         fflush (stderr);
#     endif
      SX_BREAK;
   }

}

SxException::SxException (const SxException &chain_,
                          const SxString &tag_,
                          uint32_t tagHash_,
                          const SxString &subTag_,
                          uint32_t subtagHash_,
                          const SxList<SxVariant> &args,
                          const SxString &senderFile_,
                          int senderLine_)
   : message(),
     tag(),
     subTag(),
     tagHash(tagHash_),
     subtagHash(subtagHash_),
     senderFile (),
     senderLine (senderLine_),
     next(NULL),
     prev(NULL),
     tail(NULL)
{
   SX_DBG_MSG ("throwing (" << tag_ << "," << subTag_<< ")" << " from "
                            << senderFile_ << ":" << senderLine_);
   try {
      set_new_handler(NULL);
      next = new SxException();
      next->copy (chain_);
      init (tag_, subTag_, args, senderFile_);
      this->tail = next->tail;
   } catch (std::bad_alloc &excep) {
      SX_UNUSED(excep);
      cerr << "Exception Stack: \n";
      const SxException *ptr = chain_.tail;
      int counter = 1;
      int i = 0;
      int strSize = 0;
      while (ptr) {
         cerr << counter++ << ": ";
         if (ptr->tag.getSize () > 0)
            cerr << "[" << ptr->tag << ", " << ptr->subTag << "] ";
         if (ptr->message.getSize () > 0) cerr << ptr->message << " ";
         strSize = (int)ptr->senderFile.getSize ();
         i = strSize - 1;
         while (i >= 0 && ptr->senderFile(i) != '\\' &&
                ptr->senderFile(i) != '/') i--;
         cerr << "(";
         while (((ssize_t)i+1) < strSize)  { cerr << ptr->senderFile((ssize_t)i+1); i++; }
         cerr << ":" << ptr->senderLine << ")";
         ssize_t argSize = ptr->arguments.getSize ();
         if (argSize > 0) {
            cerr << " ";
            for (const SxVariant &v : ptr->arguments)
               cerr << v << " ";
         }
         cerr << endl;
         ptr = ptr->prev;
      }
      cerr << counter++ << ": ";
      cerr << "[" << tag_, subTag_ << "] ";
      strSize = (int)senderFile_.getSize ();
      i = strSize - 1;
      while (i >= 0 && senderFile_(i) != '\\' && senderFile_(i) != '/') i--;
      cerr << "(";
      while (((ssize_t)i+1) < strSize)  { cerr << senderFile_((ssize_t)i+1); i++; }
      cerr << ":" << senderLine_ << ")";
      cerr << " ";
      for (const SxVariant &v : args) {
         cerr << v << " ";
      }
      cerr << endl;
      sxOutOfMemoryHandler ();
   }
   set_new_handler(sxOutOfMemoryHandler);

   if (dumpCoreFile)  {
      print ();
#     ifndef SX_ANDROID
         fflush (stdout);
         fflush (stderr);
#     endif
      SX_BREAK;
   }
}

SxException::~SxException ()
{
   destroy ();
}

SxException &SxException::operator= (const SxException &in)
{
   if (&in != this)  {
      copy (in);
   }
   return *this;
}

void SxException::init (const SxException &in)
{
   init (in.message, in.tag, in.senderFile);
   subTag     = in.subTag;
   tagHash    = in.tagHash;
   subtagHash = in.subtagHash;
   arguments  = in.arguments;
   senderLine = in.senderLine;
}

void SxException::init (const SxString &message_,
                        const SxString &tag_,
                        const SxString &senderFile_)
{
   message    = message_;
   tag        = tag_;
   senderFile = senderFile_;
}

void SxException::init (const SxString &tag_,
                        const SxString &subTag_,
                        const SxList<SxVariant> &args,
                        const SxString &senderFile_)
{
   tag = tag_;
   subTag = subTag_;
   arguments = args;
   senderFile = senderFile_;
}

void SxException::destroy ()
{
   if (next)  {
      delete next;
      next = NULL;
   }

   if (prev) {
      prev = NULL;
   }
}

void SxException::copy (const SxException &in)
{
   destroy ();
   init (in);
   SxException *dst = this;
   SxException *src = in.next;
   while (src)  {
      dst->next = new SxException();
      dst->next->init (*src);
      dst->next->prev = dst;
      dst = dst->next;
      src = src->next;
   }
   this->tail = dst;
}

void SxException::causeSegFault ()
{
//#  ifdef NDEBUG
//      std::cout ("Exceptions can be traced only in the SFHIngX DEBUG mode.\n");
//      dumpCoreFile = false;
//#  else
      std::cout <<  "Exceptions will cause segmentation faults and "
              "dump a core file.\n";
      dumpCoreFile = true;
//#  endif /* NDEBUG */
}


void SxException::throwExceptions ()
{
   dumpCoreFile = false;
}

void SxException::print (bool printSender) const
{
   const SxException *ptr = this;
   while (ptr)  {
      std::cerr << std::endl;
      if (tag.getSize () > 0) {
         std::cerr << "[" << ptr->tag << ", " << ptr->subTag << "]";
      }
      ssize_t argSize = ptr->arguments.getSize ();
      if (argSize > 0) {
         std::cerr << " ";
         for (const SxVariant &v : ptr->arguments)
             std::cerr << v.toString () << " ";
         std::cerr << std::endl;
      }
      if (ptr->message != "")
         std::cerr << ptr->message << std::endl;

      if (printSender)
         std::cerr << "This Exception was emitted from "
                   << ptr->senderFile
                   << " (line " << ptr->senderLine << ")\n";

      ptr = ptr->next;
   }
}

SxString SxException::toString (int mode,
                                const SxString &delimiter_) const
{
   SX_CHECK (mode >= 0x00 && mode <= 0x03, mode);

   SxList<SxString> list;

   const SxException *ptr = this;
   while (ptr)  {
      SxString str, type;
      if (mode & DebugSimple)  {
         str = ptr->senderFile;
         ssize_t i = str.getSize () - 1;
         while (i >= 0 && str(i) != '\\' && str(i) != '/') i--;
         if (i >= 0) str = str.subString (i + 1);
         str = " (" + str + ":" + SxString(ptr->senderLine) + ")";
      }
      if (mode & Stack && ptr->tag.getSize () > 0)  {
         type = "[" + ptr->tag + ", " + ptr->subTag + "] ";
      }
      str = type + ptr->message + str;
      ssize_t argSize = ptr->arguments.getSize ();
      if (argSize > 0) {
         str += " ";
         for (const SxVariant &v : ptr->arguments)
             str += v.toString () + " ";
      }
      list.prepend (str);
      ptr = ptr->next;
   }

   if (mode & Stack)  {
     int i=1;
     list.foreach ([&i](auto it) { *it = SxString(i++) + ": " + *it; });
     return SxString::join (list, delimiter_);
   } else {
     return list.last ();
   }
}

void SxException::printStack (const SxString &delimiter_) const
{
   const SxException *ptr = this->tail;
   int counter = 1;
   int i = 0;
   int strSize = 0;
   while (ptr) {
      cerr << counter++ << ": ";
      if (ptr->tag.getSize () > 0)
         cerr << "[" << ptr->tag << ", " << ptr->subTag << "] ";
      if (ptr->message.getSize () > 0) cerr << ptr->message << " ";
      strSize = (int)ptr->senderFile.getSize ();
      i = strSize - 1;
      while (i >= 0 && ptr->senderFile(i) != '\\' && ptr->senderFile(i) != '/')
         i--;
      cerr << "(";
      while (((ssize_t)i+1) < strSize)  { cerr << ptr->senderFile((ssize_t)i+1); i++; }
      cerr << ":" << ptr->senderLine << ")";
      ssize_t argSize = ptr->arguments.getSize ();
      if (argSize > 0) {
         cerr << " ";
         for (const SxVariant &v : ptr->arguments)
            cerr << v << " ";
      }
      cerr << delimiter_;
      ptr = ptr->prev;
   }
}

bool SxException::hasTag (const SxString &tag_) const
{
   SX_CHECK (tag_.getSize () > 0);
   return (tag_ == tag);
}

bool SxException::findTag (const SxString &tag_) const
{
   SX_CHECK (tag_ != NULL);
   SxString needle = tag_;
   const SxException *ptr = this;
   while (ptr)  {
      if (ptr->getTag() == needle)  return true;
      ptr = ptr->next;
   }
   return false;
}

bool SxException::foundTagHash (uint32_t tagHash_) const
{
   SxException *it = this->tail;
   SX_CHECK (it);

   do {
      if (it->tagHash == tagHash_) return true;
   } while ((it = it->next));
   return false;
}

bool SxException::foundSubtagHash (uint32_t subtagHash_) const
{
   SxException *it = this->tail;
   SX_CHECK (it);
   do {
      if (it->subtagHash == subtagHash_) return true;
   } while ((it = it->next));
   return false;
}

/* This exception table is populated by static initialization; before
 * entering main (). Currently, accesses to this map are NOT thread safe.
 * Potential races may occur by dynamically loading libraries from threads
 * different from main, if those libraries also use SX_EXCEPTION_TYPE macros.
*/

SxMap<SxString,SxMap<SxString,SxList<SxString> > > &SxException_getExceptionTable ()
{
   static SxMap<SxString, SxMap<SxString, SxList<SxString> > > exceptionTable;
   return exceptionTable;
}

SxMap<SxString,SxMap<SxString,SxList<SxString> > > &SxException::getExceptionTable ()
{
   return SxException_getExceptionTable ();
}


static void SxException_addTableEntry (
   SxString tag,
   SxString subtag,
   SxList<SxString> argNames)
{
   SxMap<SxString, SxMap<SxString,SxList<SxString> > > &exceptionTable
      = SxException_getExceptionTable ();
   (exceptionTable(tag))(subtag) = argNames;
}

SxExceptionEntry::SxExceptionEntry (
   const char *tag_,
   const char *subtag_,
   const SxList<SxString> &argNames)
{
   SxString tag (tag_);
   SxString subtag (subtag_);
   SxException_addTableEntry (tag, subtag, argNames);
}
