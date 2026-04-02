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

#include <SxParallelHierarchy.h>
#include <SxParser.h>
#include <SxLoopMPI.h>

void writeParallelHierarchy(const SxString &filename)
{
   SxParallelHierarchy &ph = getMPIHierarchy ();
   // find the root node
   SxTaskGroup* root = ph.taskGroups.begin ().getValue ();
   SX_CHECK (root);
   while (root->getParentId () != "ROOT") root = root->getParent ();
   if (!root->master ()) return; // only the master node must write this

   // --- printout
   ofstream out(filename.ascii ());
   out << "format parallelHierarchy;" << endl;
   root->write (out);
   out.close ();
}

static
void build_hierarchy_recursively(SxSymbolTable * node, SxSymbolTable * parent,
                                 SxParallelHierarchy &ph)
{
   if (node == NULL)
      return;

   SxString name = "";
   int siblings = 0;
   int members = 0;

   if (node->contains("name", true)) {
      name = node->get("name")->toString();
      siblings = node->get("siblings")->toInt();
      members = node->get("members")->toInt();

      // catch default values, see /src/share/sphinx/parallelHierarchy.sx
      if (siblings < 1)
         siblings = SxLoopMPI::nr ();
      if (members < 1)
         members = SxLoopMPI::nr ();

      SxPtr<SxMpiTaskGroup> mpiTg = SxPtr<SxMpiTaskGroup>::create ();

      SxMpiComm mpiComm;
      SxMpiComm parentIntraComm;
      SxString parentId;

      mpiTg->setName(name);
      mpiTg->setNsiblings(siblings);
      mpiTg->setNmembers(members);

      if (parent->contains("name", true)) {
         SxMpiTaskGroup * parentTg;
         parentId = parent->get("name")->toString();
         parentTg = ph.mpiTaskGroups(parentId).getPtr();
         parentIntraComm = parentTg->getIntraCommunicator();
         mpiComm = SxMpi::divideMpiComm(parentIntraComm, siblings);
      }
      else {
         // --- better set parent pointer of "top" to itself?
         parentId = "ROOT";
         //mpiTg->setParent(NULL);
         // ---
         mpiTg->setParent( mpiTg.getPtr() );
         mpiComm = parentIntraComm = SxMpi::mpiCommWorld();
      }
      mpiTg->setIntraCommunicator(mpiComm);
      mpiTg->setParentId(parentId);

      mpiComm = SxMpi::joinMpiComm(mpiTg->getMasterCommunicator(), parentIntraComm);
      mpiTg->setInterCommunicator(mpiComm);

      mpiTg->setSiblingRank( SxMpi::rank(parentIntraComm) / (SxMpi::size(parentIntraComm)/siblings) );
      mpiTg->setMemberRank( SxMpi::rank(mpiTg->getIntraCommunicator()) );

      // work with strings for now and build the parent-child pointer maps later
      // as the memory addresses inside the map are not yet fixed
      // TODO eliminate childIds
      for (SxList<SxSymbolTable*>::ConstIterator child = node->children.begin();
            child != node->children.end(); child++) {
         if ((*child)->name == "level")
            mpiTg->addChildId( (*child)->get("name")->toString() );
      }

      ph.mpiTaskGroups(name) = mpiTg;
   }

   for (SxList<SxSymbolTable*>::ConstIterator child = node->children.begin();
         child != node->children.end(); child++) {
      if ((*child)->name == "level")
         build_hierarchy_recursively(*child, node, ph);
   }

   return;
}

void buildHierarchy(const SxString &hierarchyFile, SxParallelHierarchy &ph)
{
   SxParser parser;
   SxParser::Table table = parser.read(hierarchyFile, "std/parallelHierarchy.std");
   SxSymbolTable * top = const_cast<SxSymbolTable*> (table.getPtr());

   // (0) brief sloppy consistency check of the symbol table
   try {
      if (top->containsGroup("level")) {
         SxSymbolTable * child = top->getGroup("level");
         if ((child->get("name")->toString() == SxString("top-level")) &&
               (child->get("members")->toInt() > 0) &&
               (child->get("members")->toInt() != SxLoopMPI::nr ())) {
            cout << SX_SEPARATOR;
            cout << "| Error: Parallel hierarchy is incompatible with actual"
                    " MPI resources" << endl;
            cout << "| Number of top-level members in hierarchy: "
                 << child->get("members")->toInt() << endl;
            cout << "| Number of MPI processes:                  "
                 << SxLoopMPI::nr () << endl;
            cout << "| Check hierarchy file '" << hierarchyFile
                 << "' or mpirun parameters." << endl;
            cout << SX_SEPARATOR;
            SX_QUIT;
         }
      }
   } catch (SxException e) {
      e.print ();
      SX_EXIT;
   }

   // (1) recursively parse the symbol table and build the parallel hierarchy accordingly
   build_hierarchy_recursively (top, NULL, ph);

}

