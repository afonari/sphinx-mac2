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

#include <SxReadXml.h>
#include <SxException.h>

using namespace rapidxml;

std::ostream& operator<< (std::ostream &out,
                          const xml_node<char> &node)
{
   out << '<';
   const char *text = node.name ();
   ssize_t n = node.name_size ();
   for (ssize_t i = 0; i < n; ++i) out << text[i];
   out << '>';
   return out;
}

xml_node<char>* get (const xml_node<char> *parent,
                     const char* name,
                     const SxString &filename,
                     bool verbose)
{
   xml_node<char> *node = parent->first_node (name);
   if (!node) {
      cout << "Failed to find <" << name << ">";
      if (parent->parent ())
         cout << " in " << *parent->parent ();
      cout << " in xml file '" << filename << "'." << endl;
      SX_QUIT;
   }
   if (verbose) cout << "Reading <" << name << ">" << endl;
   return node;
}

xml_attribute<char>* getAttr (const xml_node<char> *node, const char* name)
{
   xml_attribute<char> *attr= node->first_attribute (name);
   if (!attr) {
      cout << "Failed to find attribute '" << name << "' in "
           << *node << "." << endl;
      SX_QUIT;
   }
   return attr;
}

double toDouble (const xml_base<char> *item)
{
   SX_CHECK (item);
   double res = -1.;
   SxString val = toString (item);
   try {
      res = val.toDouble ();
   } catch (SxException &e) {
      if (xml_node<char> *parent = item->parent ())
         cout << "in " << *parent << ": " << endl;
      cout << "Error while translating '"
           << SxString (item->name (), item->name_size ())
           << "' = '" << val << "' to double." << endl;
      e.print ();
      SX_QUIT;
   }
   return res;
}

int toInt (const xml_base<char> *item)
{
   SX_CHECK (item);
   int res = -1;
   SxString val = toString (item);
   try {
      res = val.toInt ();
   } catch (SxException &e) {
      if (xml_node<char> *parent = item->parent ())
         cout << "in " << *parent << ": " << endl;
      cout << "Error while translating '"
           << SxString (item->name (), item->name_size ())
           << "' = '" << val << "' to int." << endl;
      e.print ();
      SX_QUIT;
   }
   return res;
}

SxVector<double> toVector (const xml_node<char> *item)
{
   const char *text = item->value ();
   ssize_t nChar = item->value_size ();
   SxStack<double> vals;
   ssize_t i = 0;
   for (; i < nChar; )  {
      int step = 0;
      double val = 0.;
      int nConv = sscanf (text + i, "%lf%n", &val, &step);
      i += step;
      if (nConv != 1 || step == 0) break;
      // --- capture ill-formed numbers of type 1.23456789-100
      //     (after normal 9.87654321E-99 )
      //     that I observed for some JTH potentials and probably arise from
      //     Fortran workaround if number of digits in exponent are exhausted
      if (text[i] == '-')  {
         cout << "Parsing error in <" << SxString (item->name (), item->name_size ());
         for (xml_attribute<char> *attr = item->first_attribute ();
              attr ; attr = attr->next_attribute ())  {
            cout << ' ' << SxString (attr->name (), attr->name_size ())
                 << '=' << toString (attr);
         }
         cout << ">:" << endl;
         cout << "Ill-formed number after reading " << vals.getSize () << " numbers: ";
         while (i > 0 && !::isspace(text[i-1])) --i;
         while (i < nChar && !::isspace(text[i])) cout << text[i++];
         cout << endl;
         SX_QUIT;
      }
      vals << val;
   }
   SX_CHECK (i < nChar, i, nChar);
   // check that rest of text is only spaces
   for ( ; i < nChar; ++i)  {
      if (!::isspace (text[i]))  {
         cout << "In " << *item << ": Cannot transform to vector: '";
         for ( ; i < nChar; ++i)
            cout << text[i];
         cout << "' after reading " << vals.getSize () << " items." << endl;
      }
   }
   return vals;
}

int count (const xml_node<char> *node, const char *name)
{
   int N = 0;
   for (xml_node<char> *it = node->first_node (name); it;
        it = it->next_sibling (name))
   {
      N++;
   }
   return N;
}

