/****************************************************************************
 *
 *      Module:       DisableUnneededParams.c
 *
 *      Description:  This function is called by task 0 during
 *                    initialization to explicitly mark certain
 *                    control parameters in order to prevent them
 *                    from being written out when creating restart
 *                    files.  This is done because not all the known
 *                    control parameters are applicable to every
 *                    simulation.  For instance, the set of mobility
 *                    parameters used during a given simulation
 *                    is dependent on the mobility law in use, and
 *                    will only be a subset of all the available
 *                    mobility parameters.
 *
 ****************************************************************************/
#include "Home.h"
#include "Parse.h"

/*
 *      Need only be called by task zero
 */
void DisableUnneededParams(Home_t *home) {
  Param_t *param;

  param = home->param;

  /*
   *      Not all mobility parameters are needed by all mobility functions.
   *      The easiest thing to do is first disable all of them and then only
   *      re-enable the ones appropriate to the selected mobility.
   */
  MarkParamDisabled(home->ctrlParamList, "MobScrew");
  MarkParamDisabled(home->ctrlParamList, "MobEdge");
  MarkParamDisabled(home->ctrlParamList, "MobClimb");
  MarkParamDisabled(home->ctrlParamList, "TempK");

  /*
   *      If osmotic forces are not enabled, disable/enable any params
   *      as needed.
   *
   *      Note: TempK was disabled along with mobility params above but
   *            it is needed to do osmotic forces
   */
  if (param->vacancyConcEquilibrium <= 0) {
    MarkParamDisabled(home->ctrlParamList, "vacancyConc");
    MarkParamDisabled(home->ctrlParamList, "vacancyConcEquilibrium");
  } else {
    MarkParamEnabled(home->ctrlParamList, "TempK");
  }

  /*
   *      The selected <loadType> affects which parameters are used.
   *      do that setup now.  As with the mobility parameters above,
   *      for some of the parameters, it's easier to disable a group
   *      of them and only re-enable them as needed.
   */
  MarkParamDisabled(home->ctrlParamList, "cTimeOld");
  MarkParamDisabled(home->ctrlParamList, "dCyclicStrain");
  MarkParamDisabled(home->ctrlParamList, "netCyclicStrain");
  MarkParamDisabled(home->ctrlParamList, "numLoadCycle");
  MarkParamDisabled(home->ctrlParamList, "eAmp");

  switch (param->loadType) {
  case 0:
    MarkParamDisabled(home->ctrlParamList, "indxErate");
    break;
  case 1:
    break;
  case 2:
    break;
  case 3:
    break;
  /*
   *          For loadType's 4 and 5, we use re-enabled the same set
   *          of parameters
   */
  case 4:
  case 5:
    MarkParamEnabled(home->ctrlParamList, "cTimeOld");
    MarkParamEnabled(home->ctrlParamList, "dCyclicStrain");
    MarkParamEnabled(home->ctrlParamList, "netCyclicStrain");
    MarkParamEnabled(home->ctrlParamList, "numLoadCycle");
    MarkParamEnabled(home->ctrlParamList, "eAmp");
    break;
  }

  /*
   *      If no sessile burgers vectors have been specified, don't write
   *      the related arrays out.
   */
  if (param->sessileburgspec[0] <= 0) {
    MarkParamDisabled(home->ctrlParamList, "sessileburgspec");
    MarkParamDisabled(home->ctrlParamList, "sessilelinespec");
  }

  /*
   *      If the geometries are not specified in a user-defined laboratory
   *      frame but use the standard crystalographic frame, don't write
   *      out any axes for a laboratory frame.
   */
  if (param->useLabFrame == 0) {
    MarkParamDisabled(home->ctrlParamList, "labFrameXDir");
    MarkParamDisabled(home->ctrlParamList, "labFrameYDir");
    MarkParamDisabled(home->ctrlParamList, "labFrameZDir");
  }

  /*
   *      Disable writing of the output-related parameters that do not
   *      apply to the types of output actually selected.
   */

  if (param->savepropdt > 0) {
    MarkParamDisabled(home->ctrlParamList, "savepropfreq");
  } else {
    MarkParamDisabled(home->ctrlParamList, "savepropdt");
    MarkParamDisabled(home->ctrlParamList, "saveproptime");
  }

  if (param->tecplot == 0) {
    MarkParamDisabled(home->ctrlParamList, "tecplotdt");
    MarkParamDisabled(home->ctrlParamList, "tecplottime");
    MarkParamDisabled(home->ctrlParamList, "tecplotcounter");
  }

  if (param->tecplotdt <= 0) {
    MarkParamDisabled(home->ctrlParamList, "tecplotdt");
    MarkParamDisabled(home->ctrlParamList, "tecplottime");
  }

  if ((param->savedensityspec[0] == 0) || (param->savedensityspec[1] == 0) ||
      (param->savedensityspec[2] == 0)) {
    MarkParamDisabled(home->ctrlParamList, "savedensityspec");
  }
  if (param->HandleCrossSlipPerSlipSystem == 0) {
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem1BulkCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem2BulkCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem3BulkCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem4BulkCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem5BulkCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem6BulkCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem7BulkCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem8BulkCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem9BulkCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem10BulkCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem11BulkCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem12BulkCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "OtherSlipSystemBulkCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem1SurfaceCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem2SurfaceCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem3SurfaceCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem4SurfaceCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem5SurfaceCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem6SurfaceCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem7SurfaceCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem8SurfaceCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem9SurfaceCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem10SurfaceCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem11SurfaceCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem12SurfaceCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "OtherSlipSystemSurfaceCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem1RepulsiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem2RepulsiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem3RepulsiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem4RepulsiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem5RepulsiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem6RepulsiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem7RepulsiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem8RepulsiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem9RepulsiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem10RepulsiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem11RepulsiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem12RepulsiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "OtherSlipSystemRepulsiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem1AttractiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem2AttractiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem3AttractiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem4AttractiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem5AttractiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem6AttractiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem7AttractiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem8AttractiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem9AttractiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem10AttractiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem11AttractiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "SlipSystem12AttractiveCrossSlipEventsCount");
    MarkParamDisabled(home->ctrlParamList,
                      "OtherSlipSystemAttractiveCrossSlipEventsCount");
  }
}
