#include <SxCLI.h>
#include <SxJSONParser.h>
#include <SxFileIO.h>
#include <SxSchema.h>
#include <SxGQuery.h>
#include <SxSymbol.h>
#include <SxTimer.h>

typedef typename SxGQExprBase::SelSet SelSet;
using namespace sx;


void readAtoms ()
{

   SxJSONParser demo;
   SxJSONParser schemaParser;
   try  {
      demo.setSearchPath (".");
      schemaParser.setSearchPath (".");
      SxTimer t(1);
      t.start ();

      demo.readFile ("atoms.dat");

      t.stop ();

      cout << "readfile time:" << t.getTime () << endl;

      t.reset ();
      t.start ();
      schemaParser.readFile ("atomSchema.dat");

      SxSchema schema (schemaParser.getAst ());

      t.stop ();
      std::cout << "schema conversion time: " << t.getTime () << std::endl;

      t.reset ();
      t.start ();

      std::cout << "schema res: "
                << schema.validate (demo.getAst ())
                << endl;
      t.stop ();
      std::cout << "validation time: " << t.getTime () << std::endl;

      SxPtr<SxGraph<SxGProps> > ast = demo.getAst ();
      auto it = ast->begin (0);
      ++it; // skip dummy root

      
      // ~1.5 secs, 100k
      // faster because we first search all the
      // species in targeted area and collect the atoms only
      // from those individual species. hence minimizing the 
      // number of nodes to be checked
      SxGQuery q = (N(".key") == "species");

      t.reset ();
      t.start ();
      SelSet sels = q.matchAll (ast, it.neighbors ());
      std::cout << "number of species: " << sels->getSize () << std::endl;
      t.stop ();
      std::cout << "query time: " << t.getTime () << std::endl;

      size_t c = 0;
      t.reset ();
      t.start ();
      SxList<double> lst;
      SxGQuery q2 = (N(".key") == "atom");
      for (auto selIt = sels->begin (); selIt != sels->end (); ++selIt)
      {
         SelSet aSels = q2.matchAll (ast, ast->begin ((**selIt)(0)).neighbors ());

         for (auto aselIt = aSels->begin (); aselIt != aSels->end (); ++aselIt)
         {
            auto coordIt = ast->begin ((**aselIt)(0)).out (0);
            for (ssize_t i = 0; i < coordIt.getSizeOut (); ++i)  {
               coordIt.out (i)->getProperty (".val").toDouble ();
            }
            c++;
         }

      }
      t.stop ();
      std::cout << "number of atoms: " << c << std::endl;
      std::cout << "atom coords: " << lst << std::endl;
      std::cout << "atoms query time: " << t.getTime () << std::endl;

   }  catch (SxException e)  {
      e.print ();
   }

}

void readSymbols ()
{
   using namespace SxParserKit;

   SxJSONParser demo;
   SxJSONParser schemaParser;
   try  {
      demo.setSearchPath (".");
      schemaParser.setSearchPath (".");
      SxTimer t(1);
      t.start ();

      demo.readFile ("atoms.dat");

      t.stop ();

      cout << "readfile time:" << t.getTime () << endl;

      t.reset ();
      t.start ();
      schemaParser.readFile ("atomSchema.dat");

      SxSchema schema (schemaParser.getAst ());

      t.stop ();
      std::cout << "schema conversion time: " << t.getTime () << std::endl;

      t.reset ();
      t.start ();

      std::cout << "schema res: "
                << schema.validate (demo.getAst ())
                << endl;
      t.stop ();
      std::cout << "validation time: " << t.getTime () << std::endl;

      t.reset ();
      t.start ();


      //auto it = demo.begin ();
      SxSymbol it = demo.getRootSymbol ();
      it = it.getElem ("species");
      cout << "name: " << it.getKey () << endl;
      int c = 0;
      SxList<SxSymbol> lst = it.toList ();

      for (auto itt = lst.begin (); itt.isValid (); ++itt)  {
         cout << itt->getElem ("coords").toDoubleList () << endl;
         c++;
      }

      t.stop ();
      std::cout << "access time: " << t.getTime () << std::endl;
      cout << "count: " << c << endl;

   }  catch (SxException e)  {
      e.print ();
   }
}

void prettyPrint (const SxString &inputFile, const SxString &schemaFile,
                  const SxString &outFile)
{

   using namespace SxParserKit;

   SxJSONParser dataParser;
   SxJSONParser schemaParser;

   try  {
      dataParser.setSearchPath (".");
      schemaParser.setSearchPath (".");

      // -- parse input json file
      dataParser.readFile (inputFile);

      // -- parse schema json file
      schemaParser.readFile (schemaFile);
      SxSchema schema (schemaParser.getAst ());

      schema.validate (dataParser.getAst ());

      SxSymbol tr = dataParser.getRootSymbol ();

      SxFileIO f;
      f.open (outFile, "w");

      tr.print (f);
   }  catch (SxException e)  {
      e.printStack ();
   }
}

void readImgSx ()  {

   using namespace SxParserKit;

   SxJSONParser dataParser;
   SxJSONParser schemaParser;

   // --- type of elements
   typedef typename SxParserAst::ElemType Type;

   try  {
      dataParser.setSearchPath (".");
      schemaParser.setSearchPath (".");

      // -- parse input json file
      dataParser.readFile ("img.json");

      // -- parse schema json file
      schemaParser.readFile ("dxf.jtd");
      SxSchema schema (schemaParser.getAst ());

      schema.validate (dataParser.getAst ());

      SxSymbol table = dataParser.getRootSymbol ();

      // --- read elem of type group from current group/level
      SxSymbol unitGrp = table.getElem ("unit");
      auto version = unitGrp.getElem ("version").toString ();
      auto appId = unitGrp.getElem ("appId").toString ();
      auto patchId = unitGrp.getElem ("patchId").toString ();
      bool hasRegistry = false;
      SX_UNUSED (version, appId, patchId, hasRegistry);
      if (unitGrp.hasElem ("hasRegistry"))  hasRegistry = true;

      SxString contact;
      SX_UNUSED (contact);
      if (unitGrp.hasElem ("contact"))  {
         contact = unitGrp.getElem ("contact").toString ();
      }  else  {
         contact = "feedback@rockitlaunched.com";
      }

      SxString guid;
      SX_UNUSED (guid);
      if (unitGrp.hasElem ("guid"))  {
         guid = unitGrp.getElem ("guid").toString ();
      }  else  {
         guid = "";
      }


      SxSymbol imageGrp = unitGrp.getElem ("image");
      if (!imageGrp.isValid ())  {
         SX_THROW ("could not find 'image'");
      }

      auto nBlocksRAW = imageGrp.getElem ("size").toString ().toInt64 ();
      SX_UNUSED (nBlocksRAW);

      int64_t emulatedSectorSize;
      SX_UNUSED (emulatedSectorSize);

      if (imageGrp.hasElem ("emulatedSectorSize"))  {
         emulatedSectorSize = imageGrp.getElem ("emulatedSectorSize").toString ().toInt64 ();
      }  else  {
         emulatedSectorSize = 4096;
      }

      auto fileSize = imageGrp.getElem ("fileSize").toString ().toInt64 ();
      SX_UNUSED (fileSize);
      int64_t fileMaxSize;
      SX_UNUSED (fileMaxSize);
      if (imageGrp.hasElem ("fileMaxSize"))  {
         fileMaxSize = imageGrp.getElem ("fileMaxSize").toString ().toInt64 ();
      }
      int64_t cacheSizeMB;
      SX_UNUSED (cacheSizeMB);
      if (imageGrp.hasElem ("cacheSizeMB"))  {
         cacheSizeMB = imageGrp.getElem ("cacheSizeMB").toString ().toInt64 ();
      }
      int64_t nBlocksMnt = imageGrp.getElem ("mount").toString ().toInt64 ();
      SX_UNUSED (nBlocksMnt);

      int64_t nBlocksExtra = 0;
      SX_UNUSED (nBlocksExtra);
      if (imageGrp.hasElem ("extra"))
         nBlocksExtra = imageGrp.getElem ("extra").toString ().toInt64 ();

      if (imageGrp.hasElem ("mountMode"))  {
         SxString mntMode = imageGrp.getElem ("mountMode").toString ();
         if (mntMode == "ro")  {
            SX_EXIT;
            // mountMode = ...
         }
      }

      // --- unit.image.underflow
      SxSymbol underflowGrp = imageGrp.getElem ("underflow");
      SxString actionStr = underflowGrp.getElem ("action").toString ();

      int uAction = 0; // ...
      double ratio = 0.;
      SX_UNUSED (uAction, ratio);
      if (uAction)  {
         SxSymbol resumeGrp = underflowGrp.getElem ("resume");
         if (resumeGrp.hasElem ("ratio"))  {
            ratio = resumeGrp.getElem ("ratio").toDouble () / 100.;
         }  else if (resumeGrp.hasElem ("minMB"))  {
            ratio = (resumeGrp.getElem ("minMB").toDouble () * 1e6) / (double)fileSize;
         }

         if (resumeGrp.hasElem ("minBandwidth")) {
            int64_t minBandwidth = resumeGrp.getElem ("minBandwidth").toInt ();
            SX_UNUSED (minBandwidth);
         }
         uint64_t nBlocksInit; // ...
         SX_UNUSED (nBlocksInit);
         // ...
         if (resumeGrp.hasElem ("trace"))  {
            bool tracing = resumeGrp.getElem ("trace").toBool ();
            SX_UNUSED (tracing);
         }
         // --- underflow handling
         SxString appData = "";
         SxString configFile = "";
         SxString windowName = "";
         SX_UNUSED (appData, configFile, windowName);
         if (underflowGrp.hasElem ("windowed"))  {
            // uMode = Windowed;
            SxSymbol windowGrp = underflowGrp.getElem ("windowed");
            SX_UNUSED (windowGrp);
            SX_CHECK (windowGrp.getType () == Type::Group);

            if (windowGrp.hasElem ("appData"))  {
               appData = windowGrp.getElem ("appData").toString ();
            }
            if (windowGrp.hasElem ("configFile"))  {
               configFile = windowGrp.getElem ("configFile").toString ();
            }
            if (windowGrp.hasElem ("delayTime")) {
               int64_t delayTime = windowGrp.getElem ("delayTime").toInt ();
               SX_UNUSED (delayTime);
            }
            if (windowGrp.hasElem ("launchDelay")) {
               int64_t launchDelay = windowGrp.getElem ("launchDelay").toInt ();
               SX_UNUSED (launchDelay);
            }
            if (windowGrp.hasElem ("windowName")) {
               windowName = windowGrp.getElem ("windowName").toString ();
            }
            if (windowGrp.hasElem ("noSuspend")) {
               bool noSuspend = windowGrp.getElem ("noSuspend").toBool ();
               SX_UNUSED (noSuspend);
            }

            if (windowGrp.hasElem ("minimizeMethod"))  {
               SxString minModeStr = windowGrp.getElem ("minimizeMethod").toString ();
               SX_UNUSED (minModeStr);
               // ...
            }
         }  else if (underflowGrp.hasElem ("onDemand"))  {
            // uMode = onDemand;
            SxSymbol onDemandGrp = underflowGrp.getElem ("onDemand");
            if (onDemandGrp.hasElem ("minBlocksDemand"))  {
               int64_t minBlocksDemand = onDemandGrp.getElem ("minBlocksDemand").toInt ();
               SX_UNUSED (minBlocksDemand);
            }

            if (onDemandGrp.hasElem ("maxBlockDemand"))  {
               int64_t maxBlocksDemand = onDemandGrp.getElem ("maxBlockDemand").toInt ();
               SX_UNUSED (maxBlocksDemand);
            }

            if (onDemandGrp.hasElem ("nUnderflowMinimize"))  {
               int64_t nUnderflowMinimize = onDemandGrp.getElem ("nUnderflowMinimize").toInt ();
               SX_UNUSED (nUnderflowMinimize);
            }

            if (onDemandGrp.hasElem ("detectInterval"))  {
               uint64_t detectIntervalMS = onDemandGrp.getElem ("detectInterval").toInt ()
                                           * 1000ULL;
               SX_UNUSED (detectIntervalMS);
            }
         }  else  {
            SX_EXIT;
            // uMode = Shutdown;
         }
      }  else  {
         bool zeroOnUnderflow = true;
         int64_t dT = 0;
         SX_UNUSED (zeroOnUnderflow, dT);
      }

      if (unitGrp.hasElem ("loadprofile"))  {
         // --- to be implemented
         SX_EXIT;
         // --- might need to change constructor here
         // loadProfiler = DxLoadProfiler (unitGrp.getElem ("loadprofile"));
      }

      SxSymbol execGrp = unitGrp.getElem ("execute");
      SxString cmd = execGrp.getElem ("cmd").toString ();
      SX_UNUSED (execGrp, cmd);
      if (execGrp.hasElem ("procName")) {
         SxString procName = execGrp.getElem ("procName").toString ();
         SX_UNUSED (procName);
      }

      if (execGrp.hasElem ("whitelist"))  {
         SX_EXIT;
         // SxString whitelist = execGrp ...
      }

      if (unitGrp.hasElem ("fileinfo"))  {
         SxSymbol fileInfoGrp = unitGrp.getElem ("fileinfo");
         SX_UNUSED (fileInfoGrp);

         if (fileInfoGrp.hasElem ("imgdxf"))  {
            SxSymbol imgdxfGrp = fileInfoGrp.getElem ("imgdxf");
            SxString md5DXF = imgdxfGrp.getElem ("md5sum").toString ();
            SX_UNUSED (imgdxfGrp, md5DXF);
         }

         if (fileInfoGrp.hasElem ("imgidx"))  {
            SxSymbol imgIdxGrp = fileInfoGrp.getElem ("imgidx");
            SxString md5Idx = imgIdxGrp.getElem ("md5sum").toString ();
            int64_t sizeIdx = imgIdxGrp.getElem ("size").toString ().toInt64 ();
            SX_UNUSED (imgIdxGrp, md5Idx, sizeIdx);
         }

         if (fileInfoGrp.hasElem ("imgraw"))  {
            SxSymbol imgRAWGrp = fileInfoGrp.getElem ("imgraw");
            SxString md5RAW = imgRAWGrp.getElem ("md5sum").toString ();
            SX_UNUSED (imgRAWGrp, md5RAW);
         }
      }

      // --- rest of the code
      // .....


   }  catch (SxException e)  {
      if (e.isCategory<"ParserKit"_SX>() || e.hasCategory<"ParserKit"_SX> ())
      {
         cout << "Syntax error: " << e.toString () << endl;
      }
   }

}

void fromUtf16 ()
{

   using namespace SxParserKit;

   SxJSONParser dataParser;
   SxJSONParser schemaParser;

   try  {
      dataParser.setSearchPath (".");
      schemaParser.setSearchPath (".");

      SxString data = SxFileIO::readUtf16 ("abc.dat");

      // --- parse input json file
      dataParser.readString (data, "abc.dat");

      // --- parse schema json file
      schemaParser.readFile ("def.dat");
      SxSchema schema (schemaParser.getAst ());
      schema.validate (dataParser.getAst ());

      SxSymbol tr = dataParser.getRootSymbol ();

      cout << tr.getElem ("k", true).toString () << endl;
      cout << tr.getElem ("num").toDouble () << endl;
      cout << tr.getElem ("num2").toInt () << endl;
      cout << tr.getElem ("num3").toInt () << endl;
      cout << tr.getElem ("num4").toInt () << endl;
      cout << tr.getElem ("bol").toBool () << endl;

      SxList<SxSymbol> lists = tr.getArrayList ("arr");

      for (auto lstIt = lists.begin (); lstIt.isValid (); ++lstIt)  {
         SxList<SxSymbol> elems = lstIt->toList ();
         for (auto it = elems.begin (); it.isValid (); ++it)  {
            cout << it->toDoubleList () << endl;
         }
      }

   }  catch (SxException e)  {
      if (e.isCategory<"ParserKit"_SX> ())  {
         cout << "Syntax error:" << e.toString () << endl;
      }  else if (e.isCategory<"SchemaValidation"_SX> ())  {
         cout << "Validation error:" << e.toString () << endl;
      }  else if (e.isCategory<"SxSymbol"_SX> ())  {
         cout << "SymbolIO error: " << e.toString () << endl;
      }  else  {
         e.printStack ();
      }
   }
}

void writerTest ()
{
   using namespace SxParserKit;
   typedef typename SxParserAst::ElemType Type;
   try {
      SxSymbol s;
      s.setType (Type::Group);
      //s.append ("k\\\"hel\noo\"1", "he\\nlπαράδl\too  \"bar\"");
      s.append ("k1", "helloo  \"bar\"");
      s.append ("k2", "ตัวอย่าง");
      s.append ("k3", "ตัวอย่าง");
      s.append ("παράδειγμα", 44.5);
      s.append ("k4", "v4");
      s.append ("k5", "v5");
      s.append ("k6", 11);


      auto s2 = s.append ("myList", Type::List);
      //s.setType (Type::List);
      s2.append ("v1");
      s2.append ("v2");
      s2.append ("ตัวอย่าง");
      s2.append ("v4");
      s2.append ("v5");
      s2.append ("v6");

      auto s3 = s2.append (Type::Group);
      s3.append ("k1", "v1");
      s3.append ("k2", "v2");
      s3.append ("k3", 44.5);
      s3.append ("k4", "v4");
      s3.append ("k5", "v5");
      s3.append ("k6", 11);

      auto lst = s.findAll ("k3");
      for (auto it = lst.begin (); it.isValid (); ++it) {
         cout << "found: " << it->getKey () << ":: val: " << it->getValue () << endl;
      }

      s2.remove (3);
      s3.remove ("k3");
      auto tmpSym = s2.find ("k3");

      if (tmpSym.isValid ())
         cout << "found: " << tmpSym.getKey ()
              << ":: val: " << tmpSym.getValue () << endl;

      s.write ("tmp.json");

   } catch (SxException e) {
      e.printStack ();
   }
}

void updateTest ()
{
   using namespace SxParserKit;
   try {
      SxJSONParser demo;
      demo.readFile ("tmp.json");
      SxSymbol s = demo.getRootSymbol ();

      SxSymbol lst = s.find ("myList");
      if (lst.isValid ()) {
         lst.append ("newElem1");
         lst.append ("newElem2");
         lst.append ("newElem3");
      }

      s.write ("tmp.json");

   } catch (SxException e) {
      e.printStack ();
   }
}

void writeAppJson ()
{
   using namespace SxParserKit;
   typedef typename SxParserAst::ElemType Type;
   try {
      SxSymbol root;
      root.setType (Type::Group);
      root.append ("version", 2);
      root.append ("dxfVersion", 1);
      root.append ("guid", "44924089-1f35-4145-8e07-52692037395d");

      SxSymbol app = root.append ("app", Type::Group);
      app.append ("cmd","{{DXDRIVE}}\\ROCKIT-protected\\rayman origins\\");
      app.append ("delay", 5);

      SxSymbol disk = root.append ("disk", Type::Group);
      disk.append ("sizeMnt", 9035776);
      disk.append ("cacheSize", 0);
      disk.append ("mountMode", "ro");

      SxSymbol underflow = root.append ("underflow", Type::Group);
      underflow.append ("delay", 5);
      underflow.append ("suspendDrive", false);
      underflow.append ("suspendProc", false);
      underflow.append ("minOnDemand", 819200);
      underflow.append ("maxOnDemand", 819200);

      SxSymbol checkSums = root.append ("checksums", Type::Group);
      checkSums.append ("mnt", "587ce23188733ec4821226558d9da7a3");

      root.write ("app.json");
   } catch (SxException e) {
      e.printStack ();
   }
}

void writeRuntimeJson ()
{
   using namespace SxParserKit;
   typedef typename SxParserAst::ElemType Type;
   try {
      SxSymbol root;
      root.setType (Type::Group);
      root.append ("version", 2);
      root.append ("guid", "");
      root.append ("time0", 0);
      root.append ("nStarts", 1);
      root.append ("runTime", 3);
      root.append ("preloadSize", 37202424);
      root.append ("timeToStart", 2);
      root.append ("downloadFromZero", true);
      root.append ("downloadComplete", false);
      root.append ("avgSpeed", 99.768);
      root.append ("sigmaSpeed", 0.0979529);
      root.append ("nUnderflows", 0);
      root.append ("avgUnderflowTime", 0);
      root.append ("avgUnderflowSize", 0);
      root.append ("maxUnderflowTime", 0);
      root.append ("maxUnderflowSize", 0);
      root.append ("nOnDemandFetches", 4);
      root.append ("avgOnDemandSize", 28295168);
      root.append ("avgOnDemandTime", 329);

      root.write ("runtime.json");
   } catch (SxException e) {
      e.printStack ();
   }
}

int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   cli.finalize ();


   using namespace SxParserKit;


   SxJSONParser demo;
   SxJSONParser demo2;
   try  {
      demo.setSearchPath (".");
      demo2.setSearchPath (".");

      // --- parse invalid comment
      //demo.readString ("*/abc");

      // --- parse from file
      //SxString s = SxString::fromUtf8 (SxString::read ("abc.dat").getElems());

      // --- parse UTF-8 string and test file include
      //SxString s = SxString::fromUtf8 (u8"\n#include<./abc.dat> \"παράδειγμα\":\"aad\\n\\nds\xADfe\"");

      // --- parse multilevel comments
      //SxString s = SxString::fromUtf8 (u8"/*co/*mme*/nt*/\"παράδειγμα\":\"aad\\n\\nds\xADfe\"");

      // --- parse UTF-8 string with comments and new lines

      //SxString s = SxString::unicodeFromUtf8 (u8"/*comment*/\"παράδειγμα\":\"ตัว\\n\\nอย่าง\"");
      //SxString s = SxString::unicodeFromUtf8 (u8"/*comment*/\"παράδειγμα\":\"ตัวอย่าง\"");
      //SxString s2 = SxString::unicodeFromUtf8 (u8"ตัวอย่าง");
      //demo.readString (s);

      // --- example json file
      demo.readFile ("abc.txt");

      // --- example json schema file
      demo2.readFile ("def.txt");
      SxSchema schema (demo2.getAst ());

      std::cout << "Schema res: " << schema.validate (demo.getAst ()) << std::endl;

      SxSymbol tr = demo.getRootSymbol ();

      tr.print ();
   }  catch (SxException e)  {
      e.print ();
   }

   return 0;
}
