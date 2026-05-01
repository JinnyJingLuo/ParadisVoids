#ifndef _TwinPlaneCrossSlip_h
#define _TwinPlaneCrossSlip_h
#include "ParadisCrossSlipServer.h"
void TwinPlaneCrossSlip(Home_t *home);
void MakeFourArmsNode(Home_t *home);
void DissociateFourArmsNode(Home_t *home);
void ExemptCollisionAfterDissociation(Home_t *home);
bool HingeJointCanHappen(Home_t *home, Node_t *poNode1, Node_t *poNode2,
                         Node_t *poNode3);
bool SegSegCanHappen(Home_t *home, Node_t *poNode1, Node_t *poNode2,
                     Node_t *poNode3, Node_t *poNode4);
void SurfaceNodeDissociation(Home_t *home, Node_t *node);
void RemeshCSplane(Home_t *home);
#endif