#include <SxMpiComm.h>

#include <SxTaskGroup.h>
#include <SxMpiTaskGroup.h>
#include <SxParallelHierarchy.h>

#include <SxMap.h>
#include <SxString.h>
#include <SxList.h>
#include <SxError.h>
#include <SxLoopMPI.h>


bool (*SxParallelHierarchy::setupRoutine)(SxParallelHierarchy &ph) = NULL;

SxParallelHierarchy::SxParallelHierarchy()
{
   if (setupRoutine && setupRoutine (*this))  {
      // set-up pointers between the task groups
      build_hierarchy_pointers ();

      info ();
   } else {
      // --- build default top-level group
      SxPtr<SxMpiTaskGroup> mpiTg = SxPtr<SxMpiTaskGroup>::create ();

      SxMpiComm mpiComm;
      SxMpiComm parentIntraComm;
      SxString parentId;

      int mpiSize = SxLoopMPI::nr ();
      SxString name = "top";
      int siblings = 1;

      mpiTg->setName(name);
      mpiTg->setNsiblings(siblings);
      mpiTg->setNmembers(mpiSize);

      parentId = "ROOT";
      // --- better set parent pointer of "top" to itself?
//      mpiTg->setParent(NULL);
      mpiTg->setParent( mpiTg.getPtr() );
      // ---
      if (mpiSize > 1) {
         mpiComm = parentIntraComm = SxMpi::mpiCommWorld();
         mpiTg->setIntraCommunicator(mpiComm);
      }
      mpiTg->setParentId(parentId);

      if (mpiSize > 1) {
         mpiComm = SxMpi::joinMpiComm(mpiTg->getMasterCommunicator(), parentIntraComm);
         mpiTg->setInterCommunicator(mpiComm);

         mpiTg->setSiblingRank( SxMpi::rank(parentIntraComm) / (SxMpi::size(parentIntraComm)/siblings) );
         mpiTg->setMemberRank( SxMpi::rank(mpiTg->getIntraCommunicator()) );
      } else {
         mpiTg->setSiblingRank (0);
         mpiTg->setMemberRank (0);
      }

      mpiTaskGroups(name) = mpiTg;

      // complete the setup
      build_hierarchy_pointers ();
   }
}

void SxParallelHierarchy::destroy()
{
   taskGroups.removeAll();
   mpiTaskGroups.removeAll();
   cout << SX_SEPARATOR;
   cout << "| parallel hierarchy cleared by ::destroy()" << endl;
   cout << SX_SEPARATOR;
}


void SxParallelHierarchy::build_hierarchy_pointers ()
{
   for (SxMap< SxString, SxPtr<SxMpiTaskGroup> >::Iterator it = mpiTaskGroups.begin(); it != mpiTaskGroups.end(); ++it)
   {
      taskGroups(it.getKey()) = it.getValue().getPtr();
   }

   // set up SxTaskGroup* connections between parents and children
   for (SxMap<SxString, SxTaskGroup*>::Iterator it = taskGroups.begin(); it != taskGroups.end(); ++it)
   {
      SxTaskGroup * parentTg;
      parentTg = it.getValue();
      SX_CHECK(parentTg);
      // autoLvl routines make the links themselves
      if (parentTg->hasAutoLvlName())
         continue;
      if (parentTg->childId.getSize () > 0)
      {
         for (SxList<SxString>::Iterator it2 = parentTg->childId.begin();
               it2 != parentTg->childId.end(); ++it2)
         {
            SxTaskGroup * childTg = taskGroups(*it2);
            SX_CHECK(childTg);
            parentTg->addChild(childTg);
            childTg->setParent(parentTg);
         }
         for (SxMap<SxString, SxTaskGroup*>::Iterator it2 = parentTg->children.begin(); it2 != parentTg->children.end(); ++it2)
         {
            SX_CHECK( it2.getValue() );
         }
      }
   }
}


void SxParallelHierarchy::info()
{
   cout << SX_SEPARATOR;
   cout << "| SxParallelHierarchy::info() ..." << endl;
   if (taskGroups.getSize() > 0)
   {
      for (SxMap<SxString, SxTaskGroup*>::Iterator it = taskGroups.begin(); it != taskGroups.end(); ++it)
      {
         SxString id = it.getKey();
         SxTaskGroup *tg  = it.getValue();
         SX_CHECK( (id == tg->getName()) || (id == tg->getAutoLvlName()) );
         sxprintf("| %s : members=%d, siblings=%d", tg->getName().ascii(), tg->getNmembers(), tg->getNsiblings());
         if (tg->hasChildren())
         {
            sxprintf("; children :");
            for (SxMap<SxString,SxTaskGroup*>::Iterator it2 = tg->children.begin(); it2 != tg->children.end(); ++it2)
            {
               SxString id2 = it2.getKey();
#ifndef NDEBUG
               SxTaskGroup *tg2 = it2.getValue();
               SX_CHECK( id2 == tg2->getName() );
#endif
               sxprintf(" %s", id2.ascii());
            }
         }
         else
         {
            sxprintf("; no children");
         }
         sxprintf("\n");
      }
   }
   else
   {
      sxprintf("| no entries\n");
   }
   cout << SX_SEPARATOR;
}


SxTaskGroup * SxParallelHierarchy::getTaskGroup(const SxString &name)
{
   if (taskGroups.containsKey(name))
      return taskGroups(name);
   else
   {
      sxprintf("\n%s", SX_SEPARATOR.ascii());
      sxprintf("| In %s :\n", SX_FUNC);
      sxprintf("|    SxTaskGroup \"%s\" is not defined.\n", name.ascii());
      sxprintf("%s", SX_SEPARATOR.ascii());
      SX_EXIT;
      return NULL;
   }
}



