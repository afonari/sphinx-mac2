#ifndef _SX_PARALLELHIERARCHY_H_
#define _SX_PARALLELHIERARCHY_H_

#include <SxTaskGroup.h>
#include <SxMpiTaskGroup.h>
#include <SxMap.h>
#include <SxString.h>
#include <SxMpiComm.h>

class SX_EXPORT_MPI SxParallelHierarchy
{
   protected:
      friend SxParallelHierarchy& getMPIHierarchy ();

      /// Constructor
      SxParallelHierarchy();
   public:

      /// tear down the parallel hierarchy manually
      void destroy();

      /// access a task group specified by the string label
      SxTaskGroup * getTaskGroup(const SxString &label);

      /// print basic information about the hierarchy
      void info();

      /// store the instances of SxMpiTaskGroup via SxPtr to allow to modify the instances later
      SxMap< SxString, SxPtr<SxMpiTaskGroup> > mpiTaskGroups;
      /// store plain C pointers to the SxTaskGroup objects in the map mpiTaskGroups
      SxMap< SxString, SxTaskGroup* > taskGroups;

      /** \brief Setup callback

        If it returns false, ph is not set up and a default initialization is run.
        */
      static bool (*setupRoutine)(SxParallelHierarchy &ph);

   private:
      /// set up the taskGroups pointer map and create child-parent pointer links
      void build_hierarchy_pointers ();

   public:
      /// access an MPI task group specified by the string label
      SxMpiTaskGroup * getMpiTaskGroup(const SxString &label)
      {
         SX_CHECK( mpiTaskGroups.containsKey(label) );
         return mpiTaskGroups(label).getPtr ();
      }
};

inline SxParallelHierarchy& SX_EXPORT_MPI getMPIHierarchy ()
{
   static SxParallelHierarchy ph;
   return ph;
}

#endif
