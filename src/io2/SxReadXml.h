// ---------------------------------------------------------------------------
//
//      The ab-initio based multiscale library
//
//                  S / P H I / n X
//
//      Copyright:  Max-Planck-Institut fuer Eisenforschung GmbH
//                  40237 Duesseldorf, Germany
//
//      Contact:    https://sxrepo.mpie.de
//      Authors:    see sphinx/AUTHORS
//      License:    see sphinx/LICENSE
//
// ---------------------------------------------------------------------------

#ifndef _SX_READ_XML_H_
#define _SX_READ_XML_H_

#include <SxIO2.h>
#include <SxString.h>
#include <SxVector.h>
#include <rapidxml-1.13/rapidxml.hpp>

// --- Auxiliary functions for xml parsing via rapidxml

/// Output xml node name
std::ostream& SX_EXPORT_IO2
operator<< (std::ostream &out, const rapidxml::xml_node<char> &node);

/// Get a sub-node, or fail with an error
rapidxml::xml_node<char>* SX_EXPORT_IO2
get (const rapidxml::xml_node<char> *parent,
     const char* name,
     const SxString &filename,
     bool verbose = false);

/// Get an attribute, or fail with an error
rapidxml::xml_attribute<char>* SX_EXPORT_IO2
getAttr (const rapidxml::xml_node<char> *node, const char* name);

/// Count number of sub-nodes with a given name
int count (const rapidxml::xml_node<char> *node, const char *name);

/// Convert value into SxString
inline SxString SX_EXPORT_IO2 toString (const rapidxml::xml_base<char> *item)
{
   SX_CHECK (item);
   return SxString (item->value (), item->value_size ());
}

/// Convert value into double
double toDouble SX_EXPORT_IO2 (const rapidxml::xml_base<char> *item);

/// Convert value into int
int toInt SX_EXPORT_IO2 (const rapidxml::xml_base<char> *item);

/// Convert value into vector of doubles
SxVector<double> SX_EXPORT_IO2 toVector (const rapidxml::xml_node<char> *item);

#endif /* _SX_READ_XML_H_ */
