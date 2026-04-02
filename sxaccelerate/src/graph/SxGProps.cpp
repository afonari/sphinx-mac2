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

#include <SxGProps.h>

size_t sxHash (const SxGProps &in)
{
   return SxHashFunction::hash (in.getId ());
}


SxGProps::SxGProps () : id(-1)
{
   // empty
}

SxGProps::SxGProps (ssize_t id_) : id(id_)
{
   // empty
}

SxGProps::SxGProps (const SxGProps &in)
{
   SX_TRACE ();
   id    = in.id;
   props = in.props;
}

SxGProps &SxGProps::operator= (const SxGProps &in)
{
   SX_TRACE ();
   if (this == &in)  return *this;
   id    = in.id;
   props = in.props;
   return *this;
}

ssize_t SxGProps::getId () const
{
   SX_TRACE ();
   return id;
}

SxVariant &SxGProps::getProperty (const SxString &key)
{
   SX_TRACE ();
   if(props.containsKey (key))  return props(key);
   else                         SX_CHECK (false, key);
}

const SxVariant &SxGProps::getProperty (const SxString &key) const
{
   if(props.containsKey (key))  return props(key);
   else                         SX_CHECK (false, key);
}

const SxMap<SxString,SxVariant> &
SxGProps::getProperties () const
{
   SX_TRACE ();
   return props;
}

bool SxGProps::operator== (const SxGProps &n1) const
{
   return id == n1.id;
}

bool SxGProps::hasProperty (const SxString &key) const
{
   SX_TRACE ();
   return this->props.containsKey (key);
}

bool SxGProps::removeProperty (const SxString &key)
{
   SX_TRACE ();
   return this->props.removeKey (key);
}

SxGProps::operator ssize_t() const
{
   SX_TRACE ();
   return id;
} 
