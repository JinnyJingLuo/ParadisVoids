/***************************************************************************
 *
 *  Function    : InputSanity
 *  Description : Do some checks on the consistency of the input
 *
 *      Last Modified: 04/08/2008 gh - Added check that FMM is used if
 *                                     PBC is disabled.
 *
 **************************************************************************/

#include "Home.h"
#include "Param.h"
#include "Util.h"
#include "math.h"

void InputSanity(Home_t *home) {
  int i, cellCount, domainCount;
  Param_t *param;

  param = home->param;
  /*
   *	Verify that the specified geometry matches the domain count
   */
  domainCount = param->nXdoms * param->nYdoms * param->nZdoms;
  if (domainCount != home->numDomains) {
    Fatal("%s; geometry (%dX%dX%d) mismatch with domain count %d\n",
          "InputSanity", param->nXdoms, param->nYdoms, param->nZdoms,
          home->numDomains);
  }

  /*
   *	turn off dynamic load balance if uniprocessor run
   */
  if (param->nXdoms * param->nYdoms * param->nZdoms == 1)
    param->DLBfreq = 0;

  /*
   *      If the <loadType> is zero, explicitly set <eRate> to 1 and
   *      if <edotdir> is zeroed, set it to the default so there are
   *      no problems in plastic strain calculations.
   */
  if (param->loadType == 0) {
    if (param->eRate != 1.0) {
      param->eRate = 1.0;
    }
    if ((param->edotdir[0] == 0.0) && (param->edotdir[1] == 0.0) &&
        (param->edotdir[2] == 0.0)) {
      param->edotdir[0] = 1.0;
      param->edotdir[1] = 0.0;
      param->edotdir[2] = 0.0;
    }
  }

  /*
   *      Make sure the frequency at which to do multi-node splits is > 0.
   */
  param->splitMultiNodeFreq = MAX(1, param->splitMultiNodeFreq);
}
