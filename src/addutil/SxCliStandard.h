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

#ifndef _SX_CLI_STANDARD_H_
#define _SX_CLI_STANDARD_H_

#include <SxAddUtil.h>
#include <SxCLI.h>
#include <SxAtomicStructure.h>
#include <SxAOBasis.h>

/// Add "-b|--sxbfile" option, to switch between sx and sxb input
inline bool sxbOption (SxCLI &cli)
{
   return cli.option ("-b|--sxb", "input file is binary waves file "
                      "rather than S/PHI/nX input file").toBool ();
}

/// Add "-i|--input" option (default: input.sx or waves.sxb)
inline SxString inputOption (SxCLI &cli, bool sxbFile)
{
   cli.option ("-i|--input", "input file", "take original input file");
   cli.last ().defaultValue = "default: input.sx, or waves.sxb for --sxb flag";
   return cli.last ().toString (sxbFile ? "waves.sxb" : "input.sx");
}

/** \brief Read structure from input file
    @param inFile input file name
    @param sxbFile if true, input file is treated as binary file
*/
SxAtomicStructure SX_EXPORT_ADD_UTIL
structureFromInput (const SxString &inFile, bool sxbFile = false,
                    bool withMovable = false);

/** \brief Write a structure
    @param structure the structure to write
    @param outFile   filename. If zero, write to stdout (not sxb, of course)
    @param outSxb    if true, write to sxb file
  */
void SX_EXPORT_ADD_UTIL writeStructure (const SxAtomicStructure &structure,
                                        const SxString &outFile,
                                        bool outSxb = false);

class SX_EXPORT_ADD_UTIL SxStructOut {
   public:
      SxString outFile;
      bool outSxb;

      /** \brief Declare -o and --outsxb options
          @param cli the CLI object
          @param sxbFile - if true, enforce sxb file output
       */
      SxStructOut (SxCLI &cli, bool sxbFile = false);

      /// Write structure
      inline void write (const SxAtomicStructure &str)  {
         writeStructure (str, outFile, outSxb);
      }
};

/** \brief Create SxAOBasis from input
    @param table     parsed input file
    @param basisFile filename of basis descriptor file. If given, takes precedence
    @param gkBasis   G+k basis
*/

SxPtr<SxAOBasis> SX_EXPORT_ADD_UTIL
aoBasisFromInput (const SxParser::Table &table,
                  const SxString &basisFile,
                  const SxGkBasis &gkBasis);

#endif /* _SX_CLI_STANDARD_H_ */
