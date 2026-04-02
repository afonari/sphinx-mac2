#include <SxException.h>
#include <SxCLI.h>
#include <stdio.h>


SX_EXCEPTION_TAG ("NoPermissions");
SX_EXCEPTION_TAG ("FileNotFound");

SX_EXCEPTION_CAT ("FileIO");
SX_EXCEPTION_TYPE ("FileIO", "FileNotFound",  "filepath", SxString);
SX_EXCEPTION_TYPE ("FileIO", "NoPermissions", "filename", SxString,
                                              "mode",     SxString);

SX_EXCEPTION_CAT ("DxComm");
SX_EXCEPTION_TAG ("404");
SX_EXCEPTION_TYPE ("DxComm", "FileNotFound", "filename", SxString);
SX_EXCEPTION_TYPE ("DxComm", "404", "url", SxString, "other", int64_t);

struct File
{
   File (char const *path, char const *mode);
   File (File const &)            = delete;
   File &operator= (File const &) = delete;
  ~File ();

   FILE *fp;
};

File::File (char const *path, char const *mode)
   : fp (NULL)
{
   fp = fopen (path, mode);
   if (!fp) {
      try {
         SX_THROW ("FileIO", "FileNotFound", path);
      } catch (SxException e) {
         SX_RETHROW (e, "FileIO", "NoPermissions", path, mode);
      }
   }
}

File::~File ()
{
   if (fp) fclose (fp);
}


int main (int argc, char **argv) {

   SxCLI cli(argc, argv);
   cli.finalize ();

   // --- printing exception table
   {
      cout << SX_SEPARATOR;
      auto excpTable = SxException::getExceptionTable ();
      for (auto it = excpTable.begin ();
           it != excpTable.end ();
           ++it)
      {
         sxprintf ("%s {\n", it.getKey ().getElems ());
         SxMap<SxString, SxList<SxString> > &entries = *it; //.getValue ();
         for (auto entryIt = entries.begin ();
              entryIt != entries.end ();
               ++entryIt)
         {
            sxprintf ("   %s: ", entryIt.getKey ().getElems ());
            SxList<SxString> &argNames = *entryIt;
            SxString print = "{ " + SxString::join (argNames, ", ") + " }\n";
            //for (SxString &arg : argNames)
            //   sxprintf (" %s", arg.getElems ());
            //sxprintf (" }\n");
            sxprintf ("%s", print.getElems ());
         }
         sxprintf ("}\n");
      }
      cout << SX_SEPARATOR;
   }

   try {

      File f ("not_there.file", "r");

   } catch (SxException e) {

      if (e.isCategory<"FileIO"_SX> ()) {

         e.printStack ();
         cout << SX_SEPARATOR;
         e.print (true);
         SxString s = e.toString (1 | 2, "\n");
         cout << SX_SEPARATOR;
         cout << s;
         cout << endl;

      //} else if (e.hasTag<"SomeUsefulTag"_SX> ()) {
      //} else if (e.hasTag<"SomethingElse"_SX> ()) {

      } else {
         printf ("An unknown error occured\n");
         e.printStack ();
         e.print ();
         SxString s = e.toString ();
         cout << SX_SEPARATOR;
         cout << s;
         cout << endl;
      }
   }

#if 0
   try {
      SxServer s ("0.0.0.0", "3306");
   } catch (SxException e) {

      if (e.is <"Network"_SX, "PortInUse"_SX> ()) {
         //...
      } else {

      }
   }
#endif

   return 0;
}
