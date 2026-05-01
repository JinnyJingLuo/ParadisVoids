// Ahmed M. Hussein

#include "ParadisCrossSlipServer.h"
#include "Comm.h"
#include "Matrix.h"
#include "map"
#include "DSMPI.h"
#include "Point.h"

using namespace std;
using namespace DislocationSystem;
using namespace EZ;

ParadisCrossSlipServer *ParadisCrossSlipServer::m_poServer = NULL; // ask
ParadisCrossSlipServer *ParadisCrossSlipServer::GetInstance() {
  if (m_poServer == NULL) {
    m_poServer = new ParadisCrossSlipServer; // ask
  }
  return m_poServer;
}
void ParadisCrossSlipServer::Free() {
  if (m_poServer != NULL) {
    delete m_poServer; // asl
    m_poServer = NULL;
  }
}
ParadisCrossSlipServer::ParadisCrossSlipServer() { Initialize(); }
void ParadisCrossSlipServer::Initialize() {}
ParadisCrossSlipServer::~ParadisCrossSlipServer() { Reset(); }
void ParadisCrossSlipServer::Reset() { m_oNetwork.Collapse(); }
void ParadisCrossSlipServer::HandleCrossSlip(Home_t *poHome) // ask
{
  if (poHome->param->HandleCrossSlip == 0) {
    return;
  }
  bool bIntersectionOnlyCheck = false;
  if ((poHome->cycle % poHome->param->CrossSlipHandlingFrequency) != 0) {
    // bIntersectionOnlyCheck = true;
    return;
  }
  // before anything, back up all the cross slip statistics
  double dCrossSlipOldStatistics[60];
  dCrossSlipOldStatistics[0] = (double)poHome->param->BulkCrossSlipEventsCount;
  dCrossSlipOldStatistics[1] =
      (double)poHome->param->SurfaceCrossSlipEventsCount;
  dCrossSlipOldStatistics[2] =
      (double)poHome->param->AttractiveCrossSlipEventsCount;
  dCrossSlipOldStatistics[3] =
      (double)poHome->param->RepulsiveCrossSlipEventsCount;
  dCrossSlipOldStatistics[4] = poHome->param->TotalBulkCrossSlippedChainsLength;
  dCrossSlipOldStatistics[5] =
      poHome->param->TotalSurfaceCrossSlippedChainsLength;
  dCrossSlipOldStatistics[6] =
      poHome->param->TotalAttractiveCrossSlippedChainsLength;
  dCrossSlipOldStatistics[7] =
      poHome->param->TotalRepulsiveCrossSlippedChainsLength;
  if (poHome->param->HandleCrossSlipPerSlipSystem == 1) {
    dCrossSlipOldStatistics[8] =
        (double)poHome->param->SlipSystem1BulkCrossSlipEventsCount;
    dCrossSlipOldStatistics[9] =
        (double)poHome->param->SlipSystem2BulkCrossSlipEventsCount;
    dCrossSlipOldStatistics[10] =
        (double)poHome->param->SlipSystem3BulkCrossSlipEventsCount;
    dCrossSlipOldStatistics[11] =
        (double)poHome->param->SlipSystem4BulkCrossSlipEventsCount;
    dCrossSlipOldStatistics[12] =
        (double)poHome->param->SlipSystem5BulkCrossSlipEventsCount;
    dCrossSlipOldStatistics[13] =
        (double)poHome->param->SlipSystem6BulkCrossSlipEventsCount;
    dCrossSlipOldStatistics[14] =
        (double)poHome->param->SlipSystem7BulkCrossSlipEventsCount;
    dCrossSlipOldStatistics[15] =
        (double)poHome->param->SlipSystem8BulkCrossSlipEventsCount;
    dCrossSlipOldStatistics[16] =
        (double)poHome->param->SlipSystem9BulkCrossSlipEventsCount;
    dCrossSlipOldStatistics[17] =
        (double)poHome->param->SlipSystem10BulkCrossSlipEventsCount;
    dCrossSlipOldStatistics[18] =
        (double)poHome->param->SlipSystem11BulkCrossSlipEventsCount;
    dCrossSlipOldStatistics[19] =
        (double)poHome->param->SlipSystem12BulkCrossSlipEventsCount;
    dCrossSlipOldStatistics[20] =
        (double)poHome->param->OtherSlipSystemBulkCrossSlipEventsCount;
    dCrossSlipOldStatistics[21] =
        (double)poHome->param->SlipSystem1SurfaceCrossSlipEventsCount;
    dCrossSlipOldStatistics[22] =
        (double)poHome->param->SlipSystem2SurfaceCrossSlipEventsCount;
    dCrossSlipOldStatistics[23] =
        (double)poHome->param->SlipSystem3SurfaceCrossSlipEventsCount;
    dCrossSlipOldStatistics[24] =
        (double)poHome->param->SlipSystem4SurfaceCrossSlipEventsCount;
    dCrossSlipOldStatistics[25] =
        (double)poHome->param->SlipSystem5SurfaceCrossSlipEventsCount;
    dCrossSlipOldStatistics[26] =
        (double)poHome->param->SlipSystem6SurfaceCrossSlipEventsCount;
    dCrossSlipOldStatistics[27] =
        (double)poHome->param->SlipSystem7SurfaceCrossSlipEventsCount;
    dCrossSlipOldStatistics[28] =
        (double)poHome->param->SlipSystem8SurfaceCrossSlipEventsCount;
    dCrossSlipOldStatistics[29] =
        (double)poHome->param->SlipSystem9SurfaceCrossSlipEventsCount;
    dCrossSlipOldStatistics[30] =
        (double)poHome->param->SlipSystem10SurfaceCrossSlipEventsCount;
    dCrossSlipOldStatistics[31] =
        (double)poHome->param->SlipSystem11SurfaceCrossSlipEventsCount;
    dCrossSlipOldStatistics[32] =
        (double)poHome->param->SlipSystem12SurfaceCrossSlipEventsCount;
    dCrossSlipOldStatistics[33] =
        (double)poHome->param->OtherSlipSystemSurfaceCrossSlipEventsCount;
    dCrossSlipOldStatistics[34] =
        (double)poHome->param->SlipSystem1RepulsiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[35] =
        (double)poHome->param->SlipSystem2RepulsiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[36] =
        (double)poHome->param->SlipSystem3RepulsiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[37] =
        (double)poHome->param->SlipSystem4RepulsiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[38] =
        (double)poHome->param->SlipSystem5RepulsiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[39] =
        (double)poHome->param->SlipSystem6RepulsiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[40] =
        (double)poHome->param->SlipSystem7RepulsiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[41] =
        (double)poHome->param->SlipSystem8RepulsiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[42] =
        (double)poHome->param->SlipSystem9RepulsiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[43] =
        (double)poHome->param->SlipSystem10RepulsiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[44] =
        (double)poHome->param->SlipSystem11RepulsiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[45] =
        (double)poHome->param->SlipSystem12RepulsiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[46] =
        (double)poHome->param->OtherSlipSystemRepulsiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[47] =
        (double)poHome->param->SlipSystem1AttractiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[48] =
        (double)poHome->param->SlipSystem2AttractiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[49] =
        (double)poHome->param->SlipSystem3AttractiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[50] =
        (double)poHome->param->SlipSystem4AttractiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[51] =
        (double)poHome->param->SlipSystem5AttractiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[52] =
        (double)poHome->param->SlipSystem6AttractiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[53] =
        (double)poHome->param->SlipSystem7AttractiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[54] =
        (double)poHome->param->SlipSystem8AttractiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[55] =
        (double)poHome->param->SlipSystem9AttractiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[56] =
        (double)poHome->param->SlipSystem10AttractiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[57] =
        (double)poHome->param->SlipSystem11AttractiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[58] =
        (double)poHome->param->SlipSystem12AttractiveCrossSlipEventsCount;
    dCrossSlipOldStatistics[59] =
        (double)poHome->param->OtherSlipSystemAttractiveCrossSlipEventsCount;
  }

  DSMPI::Barrier();
  ClearOpList(poHome);
  Set(poHome);
  list<DislocationChain *> loChains = GenerateDislocationChains(poHome);
  double dAngle = PI / 180.0 * 15.0;
  SeparateScrewChains(poHome->param, &loChains, dAngle, bIntersectionOnlyCheck);
  RefineScrewChains(&loChains, bIntersectionOnlyCheck);
  HandleChainCrossSlip(poHome, &loChains, bIntersectionOnlyCheck);
  DSMPI::Barrier(); // ask
  CommSendRemesh(poHome);
  DSMPI::Barrier();
  FixRemesh(poHome); // ask
  DisposeChains(&loChains);
  Reset();
  DSMPI::Barrier();
  // the final step is to sync the cross slip stats across all processors
  // first, get the new cross slip events statistics
  double dCrossSlipNewStatistics[60];
  dCrossSlipNewStatistics[0] = (double)poHome->param->BulkCrossSlipEventsCount -
                               dCrossSlipOldStatistics[0];
  dCrossSlipNewStatistics[1] =
      (double)poHome->param->SurfaceCrossSlipEventsCount -
      dCrossSlipOldStatistics[1];
  dCrossSlipNewStatistics[2] =
      (double)poHome->param->AttractiveCrossSlipEventsCount -
      dCrossSlipOldStatistics[2];
  dCrossSlipNewStatistics[3] =
      (double)poHome->param->RepulsiveCrossSlipEventsCount -
      dCrossSlipOldStatistics[3];
  dCrossSlipNewStatistics[4] =
      poHome->param->TotalBulkCrossSlippedChainsLength -
      dCrossSlipOldStatistics[4];
  dCrossSlipNewStatistics[5] =
      poHome->param->TotalSurfaceCrossSlippedChainsLength -
      dCrossSlipOldStatistics[5];
  dCrossSlipNewStatistics[6] =
      poHome->param->TotalAttractiveCrossSlippedChainsLength -
      dCrossSlipOldStatistics[6];
  dCrossSlipNewStatistics[7] =
      poHome->param->TotalRepulsiveCrossSlippedChainsLength -
      dCrossSlipOldStatistics[7];
  if (poHome->param->HandleCrossSlipPerSlipSystem == 1) {
    dCrossSlipNewStatistics[8] =
        (double)poHome->param->SlipSystem1BulkCrossSlipEventsCount -
        dCrossSlipOldStatistics[8];
    dCrossSlipNewStatistics[9] =
        (double)poHome->param->SlipSystem2BulkCrossSlipEventsCount -
        dCrossSlipOldStatistics[9];
    dCrossSlipNewStatistics[10] =
        (double)poHome->param->SlipSystem3BulkCrossSlipEventsCount -
        dCrossSlipOldStatistics[10];
    dCrossSlipNewStatistics[11] =
        (double)poHome->param->SlipSystem4BulkCrossSlipEventsCount -
        dCrossSlipOldStatistics[11];
    dCrossSlipNewStatistics[12] =
        (double)poHome->param->SlipSystem5BulkCrossSlipEventsCount -
        dCrossSlipOldStatistics[12];
    dCrossSlipNewStatistics[13] =
        (double)poHome->param->SlipSystem6BulkCrossSlipEventsCount -
        dCrossSlipOldStatistics[13];
    dCrossSlipNewStatistics[14] =
        (double)poHome->param->SlipSystem7BulkCrossSlipEventsCount -
        dCrossSlipOldStatistics[14];
    dCrossSlipNewStatistics[15] =
        (double)poHome->param->SlipSystem8BulkCrossSlipEventsCount -
        dCrossSlipOldStatistics[15];
    dCrossSlipNewStatistics[16] =
        (double)poHome->param->SlipSystem9BulkCrossSlipEventsCount -
        dCrossSlipOldStatistics[16];
    dCrossSlipNewStatistics[17] =
        (double)poHome->param->SlipSystem10BulkCrossSlipEventsCount -
        dCrossSlipOldStatistics[17];
    dCrossSlipNewStatistics[18] =
        (double)poHome->param->SlipSystem11BulkCrossSlipEventsCount -
        dCrossSlipOldStatistics[18];
    dCrossSlipNewStatistics[19] =
        (double)poHome->param->SlipSystem12BulkCrossSlipEventsCount -
        dCrossSlipOldStatistics[19];
    dCrossSlipNewStatistics[20] =
        (double)poHome->param->OtherSlipSystemBulkCrossSlipEventsCount -
        dCrossSlipOldStatistics[20];
    dCrossSlipNewStatistics[21] =
        (double)poHome->param->SlipSystem1SurfaceCrossSlipEventsCount -
        dCrossSlipOldStatistics[21];
    dCrossSlipNewStatistics[22] =
        (double)poHome->param->SlipSystem2SurfaceCrossSlipEventsCount -
        dCrossSlipOldStatistics[22];
    dCrossSlipNewStatistics[23] =
        (double)poHome->param->SlipSystem3SurfaceCrossSlipEventsCount -
        dCrossSlipOldStatistics[23];
    dCrossSlipNewStatistics[24] =
        (double)poHome->param->SlipSystem4SurfaceCrossSlipEventsCount -
        dCrossSlipOldStatistics[24];
    dCrossSlipNewStatistics[25] =
        (double)poHome->param->SlipSystem5SurfaceCrossSlipEventsCount -
        dCrossSlipOldStatistics[25];
    dCrossSlipNewStatistics[26] =
        (double)poHome->param->SlipSystem6SurfaceCrossSlipEventsCount -
        dCrossSlipOldStatistics[26];
    dCrossSlipNewStatistics[27] =
        (double)poHome->param->SlipSystem7SurfaceCrossSlipEventsCount -
        dCrossSlipOldStatistics[27];
    dCrossSlipNewStatistics[28] =
        (double)poHome->param->SlipSystem8SurfaceCrossSlipEventsCount -
        dCrossSlipOldStatistics[28];
    dCrossSlipNewStatistics[29] =
        (double)poHome->param->SlipSystem9SurfaceCrossSlipEventsCount -
        dCrossSlipOldStatistics[29];
    dCrossSlipNewStatistics[30] =
        (double)poHome->param->SlipSystem10SurfaceCrossSlipEventsCount -
        dCrossSlipOldStatistics[30];
    dCrossSlipNewStatistics[31] =
        (double)poHome->param->SlipSystem11SurfaceCrossSlipEventsCount -
        dCrossSlipOldStatistics[31];
    dCrossSlipNewStatistics[32] =
        (double)poHome->param->SlipSystem12SurfaceCrossSlipEventsCount -
        dCrossSlipOldStatistics[32];
    dCrossSlipNewStatistics[33] =
        (double)poHome->param->OtherSlipSystemSurfaceCrossSlipEventsCount -
        dCrossSlipOldStatistics[33];
    dCrossSlipNewStatistics[34] =
        (double)poHome->param->SlipSystem1RepulsiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[34];
    dCrossSlipNewStatistics[35] =
        (double)poHome->param->SlipSystem2RepulsiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[35];
    dCrossSlipNewStatistics[36] =
        (double)poHome->param->SlipSystem3RepulsiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[36];
    dCrossSlipNewStatistics[37] =
        (double)poHome->param->SlipSystem4RepulsiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[37];
    dCrossSlipNewStatistics[38] =
        (double)poHome->param->SlipSystem5RepulsiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[38];
    dCrossSlipNewStatistics[39] =
        (double)poHome->param->SlipSystem6RepulsiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[39];
    dCrossSlipNewStatistics[40] =
        (double)poHome->param->SlipSystem7RepulsiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[40];
    dCrossSlipNewStatistics[41] =
        (double)poHome->param->SlipSystem8RepulsiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[41];
    dCrossSlipNewStatistics[42] =
        (double)poHome->param->SlipSystem9RepulsiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[42];
    dCrossSlipNewStatistics[43] =
        (double)poHome->param->SlipSystem10RepulsiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[43];
    dCrossSlipNewStatistics[44] =
        (double)poHome->param->SlipSystem11RepulsiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[44];
    dCrossSlipNewStatistics[45] =
        (double)poHome->param->SlipSystem12RepulsiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[45];
    dCrossSlipNewStatistics[46] =
        (double)poHome->param->OtherSlipSystemRepulsiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[46];
    dCrossSlipNewStatistics[47] =
        (double)poHome->param->SlipSystem1AttractiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[47];
    dCrossSlipNewStatistics[48] =
        (double)poHome->param->SlipSystem2AttractiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[48];
    dCrossSlipNewStatistics[49] =
        (double)poHome->param->SlipSystem3AttractiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[49];
    dCrossSlipNewStatistics[50] =
        (double)poHome->param->SlipSystem4AttractiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[50];
    dCrossSlipNewStatistics[51] =
        (double)poHome->param->SlipSystem5AttractiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[51];
    dCrossSlipNewStatistics[52] =
        (double)poHome->param->SlipSystem6AttractiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[52];
    dCrossSlipNewStatistics[53] =
        (double)poHome->param->SlipSystem7AttractiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[53];
    dCrossSlipNewStatistics[54] =
        (double)poHome->param->SlipSystem8AttractiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[54];
    dCrossSlipNewStatistics[55] =
        (double)poHome->param->SlipSystem9AttractiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[55];
    dCrossSlipNewStatistics[56] =
        (double)poHome->param->SlipSystem10AttractiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[56];
    dCrossSlipNewStatistics[57] =
        (double)poHome->param->SlipSystem11AttractiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[57];
    dCrossSlipNewStatistics[58] =
        (double)poHome->param->SlipSystem12AttractiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[58];
    dCrossSlipNewStatistics[59] =
        (double)poHome->param->OtherSlipSystemAttractiveCrossSlipEventsCount -
        dCrossSlipOldStatistics[59];
  }
  // then reduce this data across all processors
  double dTotalCrossSlipNewStatistics[60];
  if (poHome->param->HandleCrossSlipPerSlipSystem == 1) {
    DSMPI::AllReduce(dCrossSlipNewStatistics, dTotalCrossSlipNewStatistics, 60,
                     MPI_DOUBLE, MPI_SUM); // the data count should be updated
    unsigned int i = 0;
    for (i = 0; i < 60; i++) {
      dTotalCrossSlipNewStatistics[i] =
          dTotalCrossSlipNewStatistics[i] + dCrossSlipOldStatistics[i];
    }
  } else if (poHome->param->HandleCrossSlipPerSlipSystem == 0) {
    DSMPI::AllReduce(dCrossSlipNewStatistics, dTotalCrossSlipNewStatistics, 8,
                     MPI_DOUBLE, MPI_SUM); // the data count should be updated
    unsigned int i = 0;
    for (i = 0; i < 8; i++) {
      dTotalCrossSlipNewStatistics[i] =
          dTotalCrossSlipNewStatistics[i] + dCrossSlipOldStatistics[i];
    }
  }
  poHome->param->BulkCrossSlipEventsCount =
      (unsigned int)floor(dTotalCrossSlipNewStatistics[0] + 0.5);
  poHome->param->SurfaceCrossSlipEventsCount =
      (unsigned int)floor(dTotalCrossSlipNewStatistics[1] + 0.5);
  poHome->param->AttractiveCrossSlipEventsCount =
      (unsigned int)floor(dTotalCrossSlipNewStatistics[2] + 0.5);
  poHome->param->RepulsiveCrossSlipEventsCount =
      (unsigned int)floor(dTotalCrossSlipNewStatistics[3] + 0.5);
  poHome->param->TotalBulkCrossSlippedChainsLength =
      dTotalCrossSlipNewStatistics[4];
  poHome->param->TotalSurfaceCrossSlippedChainsLength =
      dTotalCrossSlipNewStatistics[5];
  poHome->param->TotalAttractiveCrossSlippedChainsLength =
      dTotalCrossSlipNewStatistics[6];
  poHome->param->TotalRepulsiveCrossSlippedChainsLength =
      dTotalCrossSlipNewStatistics[7];
  if (poHome->param->HandleCrossSlipPerSlipSystem == 1) {
    poHome->param->SlipSystem1BulkCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[8] + 0.5);
    poHome->param->SlipSystem2BulkCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[9] + 0.5);
    poHome->param->SlipSystem3BulkCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[10] + 0.5);
    poHome->param->SlipSystem4BulkCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[11] + 0.5);
    poHome->param->SlipSystem5BulkCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[12] + 0.5);
    poHome->param->SlipSystem6BulkCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[13] + 0.5);
    poHome->param->SlipSystem7BulkCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[14] + 0.5);
    poHome->param->SlipSystem8BulkCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[15] + 0.5);
    poHome->param->SlipSystem9BulkCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[16] + 0.5);
    poHome->param->SlipSystem10BulkCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[17] + 0.5);
    poHome->param->SlipSystem11BulkCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[18] + 0.5);
    poHome->param->SlipSystem12BulkCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[19] + 0.5);
    poHome->param->OtherSlipSystemBulkCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[20] + 0.5);
    poHome->param->SlipSystem1SurfaceCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[21] + 0.5);
    poHome->param->SlipSystem2SurfaceCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[22] + 0.5);
    poHome->param->SlipSystem3SurfaceCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[23] + 0.5);
    poHome->param->SlipSystem4SurfaceCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[24] + 0.5);
    poHome->param->SlipSystem5SurfaceCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[25] + 0.5);
    poHome->param->SlipSystem6SurfaceCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[26] + 0.5);
    poHome->param->SlipSystem7SurfaceCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[27] + 0.5);
    poHome->param->SlipSystem8SurfaceCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[28] + 0.5);
    poHome->param->SlipSystem9SurfaceCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[29] + 0.5);
    poHome->param->SlipSystem10SurfaceCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[30] + 0.5);
    poHome->param->SlipSystem11SurfaceCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[31] + 0.5);
    poHome->param->SlipSystem12SurfaceCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[32] + 0.5);
    poHome->param->OtherSlipSystemSurfaceCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[33] + 0.5);
    poHome->param->SlipSystem1RepulsiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[34] + 0.5);
    poHome->param->SlipSystem2RepulsiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[35] + 0.5);
    poHome->param->SlipSystem3RepulsiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[36] + 0.5);
    poHome->param->SlipSystem4RepulsiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[37] + 0.5);
    poHome->param->SlipSystem5RepulsiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[38] + 0.5);
    poHome->param->SlipSystem6RepulsiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[39] + 0.5);
    poHome->param->SlipSystem7RepulsiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[40] + 0.5);
    poHome->param->SlipSystem8RepulsiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[41] + 0.5);
    poHome->param->SlipSystem9RepulsiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[42] + 0.5);
    poHome->param->SlipSystem10RepulsiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[43] + 0.5);
    poHome->param->SlipSystem11RepulsiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[44] + 0.5);
    poHome->param->SlipSystem12RepulsiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[45] + 0.5);
    poHome->param->OtherSlipSystemRepulsiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[46] + 0.5);
    poHome->param->SlipSystem1AttractiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[47] + 0.5);
    poHome->param->SlipSystem2AttractiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[48] + 0.5);
    poHome->param->SlipSystem3AttractiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[49] + 0.5);
    poHome->param->SlipSystem4AttractiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[50] + 0.5);
    poHome->param->SlipSystem5AttractiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[51] + 0.5);
    poHome->param->SlipSystem6AttractiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[52] + 0.5);
    poHome->param->SlipSystem7AttractiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[53] + 0.5);
    poHome->param->SlipSystem8AttractiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[54] + 0.5);
    poHome->param->SlipSystem9AttractiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[55] + 0.5);
    poHome->param->SlipSystem10AttractiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[56] + 0.5);
    poHome->param->SlipSystem11AttractiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[57] + 0.5);
    poHome->param->SlipSystem12AttractiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[58] + 0.5);
    poHome->param->OtherSlipSystemAttractiveCrossSlipEventsCount =
        (unsigned int)floor(dTotalCrossSlipNewStatistics[59] + 0.5);
  }
  DSMPI::Barrier();
}
void ParadisCrossSlipServer::Set(Home_t *poParadisHome) {
  Reset();
  // set the nodes
  unsigned int i = 0;
  unsigned int iNodesCount = (unsigned int)poParadisHome->newNodeKeyPtr;
  Node_t *poParadisNode = NULL;
  DislocationNetworkNode *poNetworkNode = NULL;
  DislocationNode *poDislocationNode = NULL;
  vector<DislocationNetworkNode *> vpoNodes;
  vpoNodes.resize(iNodesCount);
  map<unsigned int, unsigned int> miNodeIndices;

  for (i = 0; i < iNodesCount; i++) {
    poParadisNode = poParadisHome->nodeKeys[i];
    // only consider the nodes in this domain
    if ((poParadisNode != NULL) &&
        (poParadisNode->myTag.domainID == poParadisHome->myDomain)) {
      // reset node cross slip prevention flag
      poParadisNode->flags &= ~NODE_NO_CROSS_SLIP;
      poNetworkNode = new DislocationNetworkNode; // ask
      poDislocationNode = poNetworkNode->GetDataPointer();
      poDislocationNode->Set(poParadisNode->x, poParadisNode->y,
                             poParadisNode->z);
      poDislocationNode->SetID(poParadisNode->myTag.index);
      poDislocationNode->SetCategory(poParadisNode->constraint);
      poDislocationNode->SetSurfaceNormal(
          Vector(poParadisNode->dNx, poParadisNode->dNy, poParadisNode->dNz));
      m_oNetwork.AddNode(poNetworkNode);
      vpoNodes[i] = poNetworkNode;
      miNodeIndices.operator[](poParadisNode->myTag.index) = i; // ask
    }
  }
  // set the arms
  unsigned int j = 0;
  unsigned int iArmsCount = 0;
  Node_t *poParadisNeighbour = NULL;
  DislocationNetworkNode *poNetworkNeighbour = NULL;
  DislocationNetworkArm *poNetworkEdge = NULL;
  DislocationSegment *poDislocationSegment;
  for (i = 0; i < iNodesCount; i++) {
    poParadisNode = poParadisHome->nodeKeys[i];
    if ((poParadisNode != NULL) &&
        (poParadisNode->myTag.domainID == poParadisHome->myDomain)) {
      poNetworkNode =
          vpoNodes[miNodeIndices.find(poParadisNode->myTag.index)->second];
      iArmsCount = poParadisNode->numNbrs;
      for (j = 0; j < iArmsCount; j++) {
        poParadisNeighbour = GetNeighborNode(poParadisHome, poParadisNode, j);
        if ((poParadisNeighbour != NULL) &&
            (poParadisNeighbour->myTag.domainID == poParadisHome->myDomain)) {
          poNetworkNeighbour =
              vpoNodes[miNodeIndices.find(poParadisNeighbour->myTag.index)
                           ->second];
          if (poNetworkNode->ConnectTo(poNetworkNeighbour, poNetworkEdge)) {
            poDislocationSegment = poNetworkEdge->GetDataPointer();
            poDislocationSegment->SetSlipPlaneNormal(
                Vector(poParadisNode->nx[j], poParadisNode->ny[j],
                       poParadisNode->nz[j]));
            poDislocationSegment->SetBurgersVector(
                Vector(poParadisNode->burgX[j], poParadisNode->burgY[j],
                       poParadisNode->burgZ[j]));
            m_oNetwork.AddEdge(poNetworkEdge);
          }
        }
      }
    }
  }
  miNodeIndices.clear();
  vpoNodes.clear();
  m_oNetwork.GenerateSequentialNodeIDs();
}
list<DislocationChain *>
ParadisCrossSlipServer::GenerateDislocationChains(Home_t *poHome) {
  // generate chains based on the network topology
  list<GraphChain<DislocationNode, DislocationSegment> *> lpoGraphChains =
      m_oNetwork.GenerateGraphChains();
  list<DislocationChain *> lpoChains;
  list<GraphChain<DislocationNode, DislocationSegment> *>::iterator
      liGraphChains;
  for (liGraphChains = lpoGraphChains.begin();
       liGraphChains != lpoGraphChains.end(); liGraphChains++) {
    lpoChains.push_back(new DislocationChain((*liGraphChains)));
  }
  // split chains based on the node type
  DislocationSystem::NodeCategory eInitialCategory =
      DislocationSystem::UnconstrainedNode;
  DislocationSystem::NodeCategory eCategory =
      DislocationSystem::UnconstrainedNode;
  DislocationChain *poNewChain = NULL;
  unsigned int i = 0;
  unsigned int iChainLength = 0;
  list<DislocationChain *>::iterator liChains;
  DislocationChain *poChain = NULL;
  // split chains based on changes in node types
  for (liChains = lpoChains.begin(); liChains != lpoChains.end(); liChains++) {
    poChain = (*liChains);
    eInitialCategory = (DislocationSystem::NodeCategory)(
        poChain->GetCurrentNode()->GetDataPointer()->GetCategory());
    iChainLength = poChain->GetSize();
    poChain->ResetIterator();
    for (i = 0; i < iChainLength; i++) {
      eCategory = (DislocationSystem::NodeCategory)(
          poChain->GetCurrentNode()->GetDataPointer()->GetCategory());
      if (eCategory != eInitialCategory) {
        if ((eInitialCategory == DislocationSystem::SurfaceNode) ||
            (eCategory == DislocationSystem::SurfaceNode)) {
          poNewChain = poChain->Split();
          if (poNewChain == NULL) {
            poChain->IncrementIterator();
            continue;
          } else {
            lpoChains.push_back(poNewChain);
            break;
          }
        }
      }
      poChain->IncrementIterator();
    }
  }
  // split chains based on the node slip plane
  Vector oInitialVector;
  Vector oCurrentVector;
  DislocationNetworkNode *poCurrentNode = NULL;
  DislocationNetworkNode *poNextNode = NULL;
  DislocationNetworkArm *poEdge = NULL;
  double dTolerance = 1E-6;
  for (liChains = lpoChains.begin(); liChains != lpoChains.end(); liChains++) {
    poChain = (*liChains);
    oInitialVector = poChain->GetSlipPlaneNormal();
    iChainLength = poChain->GetSize();
    poChain->ResetIterator();
    for (i = 0; i < iChainLength; i++) {
      poCurrentNode = poChain->GetCurrentNode();
      poNextNode = poChain->GetNextNode();
      if (poNextNode == NULL) {
        break;
      }
      if (poCurrentNode->IsConnected(poNextNode, poEdge)) {
        oCurrentVector = poEdge->GetDataPointer()->GetSlipPlaneNormal();
      } else {
        break; // this should be unreachable
      }
      // check to see if the normals are the same
      if ((oCurrentVector ^ oInitialVector).Length() > dTolerance) {
        poNewChain = poChain->Split();
        if (poNewChain == NULL) {
          poChain->IncrementIterator();
          continue;
        } else {
          lpoChains.push_back(poNewChain);
          break;
        }
      }
      poChain->IncrementIterator();
    }
  }
  // split chains based on the node Burgers vector
  for (liChains = lpoChains.begin(); liChains != lpoChains.end(); liChains++) {
    poChain = (*liChains);
    oInitialVector = poChain->GetBurgersVector();
    iChainLength = poChain->GetSize();
    poChain->ResetIterator();
    for (i = 0; i < iChainLength; i++) {
      poCurrentNode = poChain->GetCurrentNode();
      poNextNode = poChain->GetNextNode();
      if (poNextNode == NULL) {
        break;
      }
      if (poCurrentNode->IsConnected(poNextNode, poEdge)) {
        oCurrentVector = poEdge->GetDataPointer()->GetBurgersVector();
      } else {
        break; // this should be unreachable
      }
      // check to see if the burgers vectors are the same
      if ((oCurrentVector ^ oInitialVector).Length() > dTolerance) {
        poNewChain = poChain->Split();
        if (poNewChain == NULL) {
          poChain->IncrementIterator();
          continue;
        } else {
          lpoChains.push_back(poNewChain);
          break;
        }
      }
      poChain->IncrementIterator();
    }
  }
  // last step, fix the dislocation lines slip plane normals and burgers vectors
  // slip plane normal
  for (liChains = lpoChains.begin(); liChains != lpoChains.end(); liChains++) {
    oInitialVector = (*liChains)->GetSlipPlaneNormal();
    // fix the initial vector direction
    (*liChains)->ResetIterator();
    poCurrentNode = (*liChains)->GetCurrentNode();
    poNextNode = (*liChains)->GetNextNode();
    if (poCurrentNode->IsConnected(poNextNode, poEdge)) {
      if (poEdge->GetFirst() == poNextNode) {
        oInitialVector.Reverse();
      }
    }
    iChainLength = (*liChains)->GetSize();
    (*liChains)->ResetIterator();
    for (i = 0; i < iChainLength; i++) {
      poCurrentNode = (*liChains)->GetCurrentNode();
      poNextNode = (*liChains)->GetNextNode();
      if (poNextNode == NULL) {
        break;
      }
      if (poCurrentNode->IsConnected(poNextNode, poEdge)) {
        poEdge->GetDataPointer()->SetSlipPlaneNormal(oInitialVector);
      } else {
        break; // this should be unreachable
      }
      (*liChains)->IncrementIterator();
    }
  }
  // burgers vector
  for (liChains = lpoChains.begin(); liChains != lpoChains.end(); liChains++) {
    oInitialVector = (*liChains)->GetBurgersVector();
    // fix the initial vector direction
    (*liChains)->ResetIterator();
    poCurrentNode = (*liChains)->GetCurrentNode();
    poNextNode = (*liChains)->GetNextNode();
    if (poCurrentNode->IsConnected(poNextNode, poEdge)) {
      if (poEdge->GetFirst() == poNextNode) {
        oInitialVector.Reverse();
      }
    }
    iChainLength = (*liChains)->GetSize();
    (*liChains)->ResetIterator();
    for (i = 0; i < iChainLength; i++) {
      poCurrentNode = (*liChains)->GetCurrentNode();
      poNextNode = (*liChains)->GetNextNode();
      if (poNextNode == NULL) {
        break;
      }
      if (poCurrentNode->IsConnected(poNextNode, poEdge)) {
        poEdge->GetDataPointer()->SetBurgersVector(oInitialVector);
      } else {
        break; // this should be unreachable
      }
      (*liChains)->IncrementIterator();
    }
  }
  return lpoChains;
}
void ParadisCrossSlipServer::SeparateScrewChains(
    Param_t *poParam, list<DislocationChain *> *plpoChains,
    const double &dScrewAngleTolerance, bool bIntersectionOnly) {
  list<DislocationChain *>::iterator liChains;
  unsigned int i = 0;
  unsigned int iChainLength = 0;
  Vector oBurgersVector;
  Vector oLineDirection;
  DislocationNetworkNode *poCurrentNode = NULL;
  DislocationNetworkNode *poNextNode = NULL;
  DislocationNetworkArm *poEdge = NULL;
  double dDotProductLowerBound = cos(dScrewAngleTolerance);
  double closetoboundarydScrewAngleTolerance = PI / 180.0 * 1.0;
  double tempdDotProductLowerBound = cos(closetoboundarydScrewAngleTolerance);
  double closetoboundary = 10;
  liChains = plpoChains->begin();
  bool bIsInitiallyScrew = false;
  bool bIsCurrentlyScrew = false;
  unsigned int iNonScrewSegmentsCounter = 0;
  unsigned int iMaxAllowedInternalNonScrewSegments = 3;
  DislocationChain *poNewChain = NULL;
  unsigned int iFirstNodeNeighboursCount = 0;
  unsigned int iLastNodeNeighboursCount = 0;
  ChainType eTempType = InactiveChain;
  DislocationNode *po1end;
  DislocationNode *po2end;
  Point p1end, p2end;
  // start by splitting the chains based on their screw segments
  for (liChains = plpoChains->begin(); liChains != plpoChains->end();
       liChains++) {
    eTempType = (*liChains)->GetType();
    if (eTempType == InactiveChain || eTempType == JunctionChain) {
      continue;
    }
    // skip the chain if we are doing an intersection only check
    if (bIntersectionOnly && (eTempType != IntersectionChain)) {
      continue;
    }
    if (poParam->EnableTwinPlaneCrossSlip == 1) {
      // junjie: when the chain are close to boundaries, tolerance will be
      // strict
      po1end = (*liChains)->GetFirst()->GetDataPointer();
      po2end = (*liChains)->GetLast()->GetDataPointer();
      p1end.Set(po1end->GetX(), po1end->GetY(), po1end->GetZ());
      p2end.Set(po2end->GetX(), po2end->GetY(), po2end->GetZ());
      p1end = GetInLocalCoordinates(poParam, p1end);
      p2end = GetInLocalCoordinates(poParam, p2end);
      // junjie
      if (dDotProductLowerBound != tempdDotProductLowerBound) {
        dDotProductLowerBound = tempdDotProductLowerBound;
      }
      // junjie: when the chain are close to boundaries, tolerance will be
      // strict
      if (fabs(p1end.GetX()) > (poParam->Lx / 2 - closetoboundary) ||
          fabs(p2end.GetX()) > (poParam->Lx / 2 - closetoboundary) ||
          fabs(p1end.GetY()) > (poParam->Ly / 2 - closetoboundary) ||
          fabs(p2end.GetY()) > (poParam->Ly / 2 - closetoboundary) ||
          -(p1end.GetZ()) > (poParam->Lz / 2 - closetoboundary)) {
        tempdDotProductLowerBound = dDotProductLowerBound;
        dDotProductLowerBound = cos(closetoboundarydScrewAngleTolerance);
      }
      // junjie
    }

    oBurgersVector = (*liChains)->GetBurgersVector();
    oBurgersVector.Normalize(); // burgers vectors are NOT normalized by default
    oLineDirection = (*liChains)->GetInitialLineDirection();
    bIsInitiallyScrew =
        (fabs((oLineDirection * oBurgersVector)) > dDotProductLowerBound);
    iChainLength = (*liChains)->GetSize();
    (*liChains)->ResetIterator();
    iNonScrewSegmentsCounter = 0;
    for (i = 0; i < iChainLength; i++) {
      poCurrentNode = (*liChains)->GetCurrentNode();
      poNextNode = (*liChains)->GetNextNode();
      if (poNextNode == NULL) {
        continue;
      }
      oLineDirection.SetByPoints((*poCurrentNode->GetDataPointer()),
                                 (*poNextNode->GetDataPointer()));
      oLineDirection.Normalize();
      bIsCurrentlyScrew =
          (fabs((oLineDirection * oBurgersVector)) > dDotProductLowerBound);
      if (!bIsInitiallyScrew && bIsCurrentlyScrew) {
        poNewChain = (*liChains)->Split();
        if (poNewChain != NULL) {
          plpoChains->push_back(poNewChain);
          break;
        }
      }
      if (bIsInitiallyScrew) {
        if (bIsCurrentlyScrew) {
          iNonScrewSegmentsCounter = 0;
        } else {
          if (iNonScrewSegmentsCounter <= iMaxAllowedInternalNonScrewSegments) {
            iNonScrewSegmentsCounter = iNonScrewSegmentsCounter + 1;
          } else {
            (*liChains)->DecrementIterator(iMaxAllowedInternalNonScrewSegments);
            poNewChain = (*liChains)->Split();
            if (poNewChain != NULL) {
              plpoChains->push_back(poNewChain);
              break;
            }
          }
        }
      }
      (*liChains)->IncrementIterator();
    }
  }
  if (poParam->EnableTwinPlaneCrossSlip == 1) {
    if (dDotProductLowerBound != tempdDotProductLowerBound) {
      dDotProductLowerBound = tempdDotProductLowerBound;
    }
  }

  // set the non screw chains type
  for (liChains = plpoChains->begin(); liChains != plpoChains->end();
       liChains++) {
    oBurgersVector = (*liChains)->GetBurgersVector();
    oBurgersVector.Normalize(); // burgers vectors are NOT normalized by default
    oLineDirection = (*liChains)->GetInitialLineDirection();
    bIsInitiallyScrew =
        (fabs((oLineDirection * oBurgersVector)) > dDotProductLowerBound);
    if (!bIsInitiallyScrew) {
      (*liChains)->OverrideType(NonScrewChain);
    }
  }
  //		// delete all the non screw chains
  // 		liChains = plpoChains->begin();
  // 		while(liChains != plpoChains->end())
  // 		{
  // 			(*liChains)->ResetIterator();
  // 			poCurrentNode = (*liChains)->GetCurrentNode();
  //   		 	poNextNode = (*liChains)->GetNextNode();
  //   		 	if(poNextNode == NULL)
  //   		 	{
  //   		 		continue;
  //   		 	}
  //   		 	if(poCurrentNode->IsConnected(poNextNode,poEdge))
  //   		 	{
  //   		 		oBurgersVector = poEdge->GetDataPointer()->GetBurgersVector();
  //   //
  //   to get the burgers vector of the node
  //   oBurgersVector.Normalize();
  //   		 		oLineDirection.SetByPoints((*poCurrentNode->GetDataPointer()),(*poNextNode->GetDataPointer()));
  //   		 		oLineDirection.Normalize();
  //   		 		bIsInitiallyScrew =
  //   (fabs((oLineDirection*oBurgersVector))
  //   > dDotProductLowerBound);
  //   		 	}
  //   		 	else
  //   		 	{
  //   		 		continue;		// this should be
  //   unreachable
  //   		 	}
  //   		 	if(!bIsInitiallyScrew)
  //   		 	{
  //   		 		delete (*liChains);
  //   		 		liChains = plpoChains->erase(liChains);
  //   		 	}
  //   		 	else
  //   		 	{
  //   		 		liChains++;
  //   		 	}
  // 		}
}
void ParadisCrossSlipServer::RefineScrewChains(
    list<DislocationChain *> *plpoChains, bool bIntersectionOnly) {
  list<DislocationChain *>::iterator liChains;
  DislocationNetworkNode *poFirstNode = NULL;
  DislocationNetworkNode *poLastNode = NULL;
  DislocationChain *poNewChain = NULL;
  bool bIsFirstNodeFree = false;
  bool bIsLastNodeFree = false;
  ChainType eType = InactiveChain;
  unsigned int iFirstNodeNeighboursCount = 0;
  unsigned int iLastNodeNeighboursCount = 0;
  bool bSplitRequired = false;
  for (liChains = plpoChains->begin(); liChains != plpoChains->end();
       liChains++) {
    eType = (*liChains)->GetType();
    if (eType == InactiveChain || eType == NonScrewChain ||
        eType == JunctionChain) {
      continue;
    }
    // skip the chain if we are doing an intersection only check
    if (bIntersectionOnly && (eType != IntersectionChain)) {
      continue;
    }
    poFirstNode = (*liChains)->GetFirst();
    poLastNode = (*liChains)->GetLast();
    iFirstNodeNeighboursCount = poFirstNode->GetAllNeighboursCount();
    iLastNodeNeighboursCount = poLastNode->GetAllNeighboursCount();
    // handle the 13 and 14 chains cases, we will rarely need this because all
    // the 1-nodes shouldn't be free anyway
    bSplitRequired = false;
    if (iFirstNodeNeighboursCount == 1 && iLastNodeNeighboursCount >= 3) {
      bSplitRequired = true;
    } else if (iFirstNodeNeighboursCount >= 3 &&
               iLastNodeNeighboursCount == 1) {
      bSplitRequired = true;
    } else {
      bIsFirstNodeFree =
          (iFirstNodeNeighboursCount <= 2) &&
          (poFirstNode->GetDataPointer()->GetCategory() == UnconstrainedNode);
      bIsLastNodeFree =
          (iLastNodeNeighboursCount <= 2) &&
          (poLastNode->GetDataPointer()->GetCategory() == UnconstrainedNode);
      if (bIsFirstNodeFree || bIsLastNodeFree) {
        continue;
      } else {
        bSplitRequired = true;
      }
    }

    if (bSplitRequired) {
      poNewChain = (*liChains)->SplitEvenly();
      if (poNewChain != NULL) {
        plpoChains->push_back(poNewChain);
      }
    }
  }
}
void ParadisCrossSlipServer::HandleChainCrossSlip(
    Home_t *poParadisHome, list<DislocationChain *> *plpoChains,
    bool bIntersectionOnly) {
  list<DislocationChain *>::iterator liChains;
  ChainType eType = BulkChain;
  char cWrite[500];
  unsigned int iMinimumChainLength = 4;
  for (liChains = plpoChains->begin(); liChains != plpoChains->end();
       liChains++) {
    if ((*liChains)->GetSize() < iMinimumChainLength) {
      continue;
    }
    eType = (*liChains)->GetType();
    // skip the chain if we are doing an intersection only check
    if (bIntersectionOnly && (eType != IntersectionChain)) {
      continue;
    }
    // if the chain cannot cross slip due to a previous cross slip of one of its
    // node, skip it
    if (!CanChainCrossSlip(poParadisHome, (*liChains))) {
      continue;
    }
    if (poParadisHome->param->EnableTwinPlaneCrossSlip == 1) {
      // junjie: if two ends of this chain are on the cs plane, skip this chain
      // for cs determination. This means the activation energy on the cs plane
      // is super large, which is not true. We need MD simulations here to see
      // if the activation energy is different on twin plane
      double A, B, C, D;
      double X_1, Y_1, Z_1;
      double dTol = 1E-6; // set a tolerance for geometrical detection
      DislocationNode *po1end = (*liChains)->GetFirst()->GetDataPointer();
      DislocationNode *po2end = (*liChains)->GetLast()->GetDataPointer();
      Point p1end, p2end;
      A = poParadisHome->param->A;
      B = poParadisHome->param->B;
      C = poParadisHome->param->C;
      X_1 = poParadisHome->param->X_1;
      Y_1 = poParadisHome->param->Y_1;
      Z_1 = poParadisHome->param->Z_1;
      D = -(A * X_1 + B * Y_1 + C * Z_1);
      Plane crossSlipPlane;
      crossSlipPlane.Set(Vector(A, B, C), Point(X_1, Y_1, Z_1));
      p1end.Set(po1end->GetX(), po1end->GetY(), po1end->GetZ());
      p2end.Set(po2end->GetX(), po2end->GetY(), po2end->GetZ());
      if ((crossSlipPlane.GetPointDistance(p1end) < dTol) &&
          (crossSlipPlane.GetPointDistance(p2end) < dTol)) {
        continue;
      }
      // junjie
    }

    if (eType == BulkChain) {
      HandleBulkCrossSlip(poParadisHome, (*liChains));
    } else if (eType == SurfaceChain) {
      HandleSurfaceCrossSlip(poParadisHome, (*liChains));
    } else if (eType == IntersectionChain) {
      HandleIntersectionCrossSlip(poParadisHome, (*liChains), plpoChains);
    }
  }
}
void ParadisCrossSlipServer::HandleBulkCrossSlip(Home_t *poParadisHome,
                                                 DislocationChain *poChain) {
  double dGlidePlaneShearStress = 0.0;
  double dCrossSlipPlaneShearStress = 0.0;
  double dGlideEscaigStress = 0.0;
  double dCrossSlipEscaigStress = 0.0;
  double dTotalChainLength = 0.0;
  GetChainCrossSlipValues(poParadisHome, poChain, dGlidePlaneShearStress,
                          dCrossSlipPlaneShearStress, dGlideEscaigStress,
                          dCrossSlipEscaigStress, dTotalChainLength, false);
  double dBurgersMagnitude = poParadisHome->param->burgMag;
  // now, to make a cross slip decision, we need to make sure that all of the
  // following conditions are satisfied (according to satish)
  // 1. the cross slip shear stress is > 0.5*glide plane shear stress
  // 2. the cross slip shear stress is > Gb/10/L
  // 3. the probability is high enough for cross slip at this instance
  double dCrossSlipShearStressThreshold =
      poParadisHome->param->shearModulus / 2.0 / dTotalChainLength;
  if (fabs(dCrossSlipPlaneShearStress) < dCrossSlipShearStressThreshold) {
    return;
  }
  if (fabs(dCrossSlipPlaneShearStress) < fabs(0.5 * dGlidePlaneShearStress)) {
    return;
  }
  double dActivationEnergy =
      poParadisHome->param->BulkCrossSlipActivationEnergy;
  double dActivationVolume =
      poParadisHome->param->BulkCrossSlipActivationVolumeFactor *
      dBurgersMagnitude * dBurgersMagnitude * dBurgersMagnitude;
  double dBoltzmannConstant = 1.3806503E-23;
  double dTemperature = poParadisHome->param->TempK;
  // modifications requested by satish, february 2013
  double dGamma = poParadisHome->param->BulkStackingFaultEnergy;
  double dDoverBRatio = poParadisHome->param->BulkDoverBRatio;
  double dBe = poParadisHome->param->burgMag / 2.0 / sqrt(3.0);
  double dGammaEffective = dGamma + dCrossSlipEscaigStress * dBe;
  double dTolerance = 1.0E-6;
  if (dGammaEffective < dTolerance) {
    dGammaEffective = dTolerance;
  }
  double dDOverBEffective = dDoverBRatio * dGamma / dGammaEffective;
  if (dDOverBEffective < 1.0) {
    dDOverBEffective = 1.0;
  }
  dActivationEnergy = dActivationEnergy * dGamma / dGammaEffective *
                      sqrt(log(dDOverBEffective) / log(dDoverBRatio));
  double dEscaigStressDifference = dGlideEscaigStress - dCrossSlipEscaigStress;
  double dProbabilityExponent =
      -(dActivationEnergy - dActivationVolume * dEscaigStressDifference) /
      dBoltzmannConstant / dTemperature;
  double dProbabilityExponent_noTemp =
      -(dActivationEnergy - dActivationVolume * dEscaigStressDifference) /
      dBoltzmannConstant;
  double dCrossSlipFrequency = poParadisHome->param->BulkCrossSlipFrequency;
  double dCrossSlipReferenceLength =
      poParadisHome->param->BulkCrossSlipReferenceLength;
  double dCrossSlipActualTime =
      poParadisHome->param->CrossSlipHandlingFrequency *
      poParadisHome->param->deltaTT;
  double dProbabilityBound = dCrossSlipFrequency * dCrossSlipActualTime *
                             dTotalChainLength / dCrossSlipReferenceLength *
                             exp(dProbabilityExponent);
  double dRandomNumber = 0.0;

  if (poParadisHome->param->Ttype <= 0 ||
      poParadisHome->param->stress_timestep <= 0) {
    if (dProbabilityBound < 1.0) {
      dRandomNumber = Randomizer::Random();
      if (dRandomNumber > dProbabilityBound) {
        return;
      }
    }
  } else {
    dRandomNumber = Randomizer::Random();
  }
  if (AdjustChainForCrossSlip(poParadisHome, poChain)) {
    Vector oSlipPlaneNormal = poChain->GetSlipPlaneNormal();
    Vector oBurgersVector = poChain->GetBurgersVector();
    Vector oCrossSlipPlaneNormal =
        GetFCCCrossSlipPlaneNormal(oSlipPlaneNormal, oBurgersVector);
    // Yejun
    if (poParadisHome->param->Ttype > 0 &&
        poParadisHome->param->stress_timestep > 0) {
      CrossSlipChain_ThermalGradient(
          poParadisHome, poChain, oCrossSlipPlaneNormal,
          dCrossSlipFrequency * dCrossSlipActualTime * dTotalChainLength /
              dCrossSlipReferenceLength,
          dProbabilityExponent_noTemp, dRandomNumber);
    } else {
      CrossSlipChain(poParadisHome, poChain, oCrossSlipPlaneNormal);
    }
    poParadisHome->param->BulkCrossSlipEventsCount =
        poParadisHome->param->BulkCrossSlipEventsCount + 1;
    poParadisHome->param->TotalBulkCrossSlippedChainsLength =
        poParadisHome->param->TotalBulkCrossSlippedChainsLength +
        poChain->GetLength();
    // ParadisCrossSlipServer obb; // object to access IdentifyFCCSlipSystem
    // debug purpose... changed to static function, no need for object to access
    // now add the Bulk cross-slip event based on the slip system using
    // IdentifyFCCSlipSystem function
    if (poParadisHome->param->HandleCrossSlipPerSlipSystem ==
        0) // flag to on and off the crossSlip extra statistics
    {
      return;
    }
    if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) == 1) {
      poParadisHome->param->SlipSystem1BulkCrossSlipEventsCount =
          poParadisHome->param->SlipSystem1BulkCrossSlipEventsCount + 1;
      // cout <<"plane 1 accessed"
      // <<"\t"<<poParadisHome->param->SlipSystem1BulkCrossSlipEventsCount
      // <<"\t"
      // <<"Done" <<"\t"<<endl;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               2) {
      poParadisHome->param->SlipSystem2BulkCrossSlipEventsCount =
          poParadisHome->param->SlipSystem2BulkCrossSlipEventsCount + 1;
      // cout <<"plane 2 accessed" <<"\t" <<
      // poParadisHome->param->SlipSystem2BulkCrossSlipEventsCount <<"\t"
      // <<"Done" <<"\t" <<endl;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               3) {
      poParadisHome->param->SlipSystem3BulkCrossSlipEventsCount =
          poParadisHome->param->SlipSystem3BulkCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               4) {
      poParadisHome->param->SlipSystem4BulkCrossSlipEventsCount =
          poParadisHome->param->SlipSystem4BulkCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               5) {
      poParadisHome->param->SlipSystem5BulkCrossSlipEventsCount =
          poParadisHome->param->SlipSystem5BulkCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               6) {
      poParadisHome->param->SlipSystem6BulkCrossSlipEventsCount =
          poParadisHome->param->SlipSystem6BulkCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               7) {
      poParadisHome->param->SlipSystem7BulkCrossSlipEventsCount =
          poParadisHome->param->SlipSystem7BulkCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               8) {
      poParadisHome->param->SlipSystem8BulkCrossSlipEventsCount =
          poParadisHome->param->SlipSystem8BulkCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               9) {
      poParadisHome->param->SlipSystem9BulkCrossSlipEventsCount =
          poParadisHome->param->SlipSystem9BulkCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               10) {
      poParadisHome->param->SlipSystem10BulkCrossSlipEventsCount =
          poParadisHome->param->SlipSystem10BulkCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               11) {
      poParadisHome->param->SlipSystem11BulkCrossSlipEventsCount =
          poParadisHome->param->SlipSystem11BulkCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               12) {
      poParadisHome->param->SlipSystem12BulkCrossSlipEventsCount =
          poParadisHome->param->SlipSystem12BulkCrossSlipEventsCount + 1;
    } else {
      poParadisHome->param->OtherSlipSystemBulkCrossSlipEventsCount =
          poParadisHome->param->OtherSlipSystemBulkCrossSlipEventsCount + 1;
    }
  }
}
void ParadisCrossSlipServer::HandleSurfaceCrossSlip(Home_t *poParadisHome,
                                                    DislocationChain *poChain) {
  double dGlidePlaneShearStress = 0.0;
  double dCrossSlipPlaneShearStress = 0.0;
  double dGlideEscaigStress = 0.0;
  double dCrossSlipEscaigStress = 0.0;
  double dTotalChainLength = 0.0;
  GetChainCrossSlipValues(poParadisHome, poChain, dGlidePlaneShearStress,
                          dCrossSlipPlaneShearStress, dGlideEscaigStress,
                          dCrossSlipEscaigStress, dTotalChainLength, false);
  // now, to make a cross slip decision, we need to make sure that all of the
  // following conditions are satisfied (according to satish)
  // 1. the cross slip shear stress is > 0.9*glide plane shear stress
  // 2. the cross slip shear stress is > Gb/20/L
  // 3. the probability is high enough for cross slip at this instance
  double dBurgersMagnitude = poParadisHome->param->burgMag;
  double dCrossSlipShearStressThreshold =
      poParadisHome->param->shearModulus / 4.0 / dTotalChainLength;
  if (fabs(dCrossSlipPlaneShearStress) < dCrossSlipShearStressThreshold) {
    return;
  }
  if (fabs(dCrossSlipPlaneShearStress) < fabs(0.5 * dGlidePlaneShearStress)) {
    return;
  }
  double dActivationEnergy =
      poParadisHome->param->SurfaceCrossSlipActivationEnergy;
  double dActivationVolume =
      poParadisHome->param->SurfaceCrossSlipActivationVolumeFactor *
      dBurgersMagnitude * dBurgersMagnitude * dBurgersMagnitude;
  double dBoltzmannConstant = 1.3806503E-23;
  double dTemperature = poParadisHome->param->TempK;
  // modifications requested by satish, february 2013
  double dGamma = poParadisHome->param->SurfaceStackingFaultEnergy;
  double dDoverBRatio = poParadisHome->param->SurfaceDoverBRatio;
  double dBe = poParadisHome->param->burgMag / 2.0 / sqrt(3.0);
  double dGammaEffective = dGamma + dCrossSlipEscaigStress * dBe;
  double dTolerance = 1.0E-6;
  if (dGammaEffective < dTolerance) {
    dGammaEffective = dTolerance;
  }
  double dDOverBEffective = dDoverBRatio * dGamma / dGammaEffective;
  if (dDOverBEffective < 1.0) {
    dDOverBEffective = 1.0;
  }
  dActivationEnergy = dActivationEnergy * dGamma / dGammaEffective *
                      sqrt(log(dDOverBEffective) / log(dDoverBRatio));
  double dEscaigStressDifference = dGlideEscaigStress - dCrossSlipEscaigStress;
  double dProbabilityExponent =
      -(dActivationEnergy - dActivationVolume * dEscaigStressDifference) /
      dBoltzmannConstant / dTemperature;
  double dProbabilityExponent_noTemp =
      -(dActivationEnergy - dActivationVolume * dEscaigStressDifference) /
      dBoltzmannConstant;
  double dCrossSlipFrequency = poParadisHome->param->SurfaceCrossSlipFrequency;
  double dCrossSlipReferenceLength =
      poParadisHome->param->SurfaceCrossSlipReferenceLength;
  double dSurfaceCrossSlipLength = poParadisHome->param->SurfaceCrossSlipLength;
  double dCrossSlipActualTime =
      poParadisHome->param->CrossSlipHandlingFrequency *
      poParadisHome->param->deltaTT;
  double dProbabilityBound =
      dCrossSlipFrequency * dCrossSlipActualTime * dSurfaceCrossSlipLength /
      dCrossSlipReferenceLength * exp(dProbabilityExponent);
  double dRandomNumber = 0.0;
  if (poParadisHome->param->Ttype <= 0 ||
      poParadisHome->param->stress_timestep <= 0) {
    if (dProbabilityBound < 1.0) {
      dRandomNumber = Randomizer::Random();
      if (dRandomNumber > dProbabilityBound) {
        return;
      }
    }
  } else {
    dRandomNumber = Randomizer::Random();
  }
  if (AdjustChainForCrossSlip(poParadisHome, poChain)) {
    Vector oSlipPlaneNormal = poChain->GetSlipPlaneNormal();
    Vector oBurgersVector = poChain->GetBurgersVector();
    Vector oCrossSlipPlaneNormal =
        GetFCCCrossSlipPlaneNormal(oSlipPlaneNormal, oBurgersVector);
    // Yejun
    if (poParadisHome->param->Ttype > 0 &&
        poParadisHome->param->stress_timestep > 0) {
      CrossSlipChain_ThermalGradient(
          poParadisHome, poChain, oCrossSlipPlaneNormal,
          dCrossSlipFrequency * dCrossSlipActualTime * dSurfaceCrossSlipLength /
              dCrossSlipReferenceLength,
          dProbabilityExponent_noTemp, dRandomNumber);
    } else {
      CrossSlipChain(poParadisHome, poChain, oCrossSlipPlaneNormal);
    }
    poParadisHome->param->SurfaceCrossSlipEventsCount =
        poParadisHome->param->SurfaceCrossSlipEventsCount + 1;
    poParadisHome->param->TotalSurfaceCrossSlippedChainsLength =
        poParadisHome->param->TotalSurfaceCrossSlippedChainsLength +
        poChain->GetLength();
    // cout <<"surface accessed" <<endl;
    // now add the surface cross-slip event based on the slip system using
    // IdentifyFCCSlipSystem function
    if (poParadisHome->param->HandleCrossSlipPerSlipSystem == 0) {
      return;
    }
    if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) == 1) {
      poParadisHome->param->SlipSystem1SurfaceCrossSlipEventsCount =
          poParadisHome->param->SlipSystem1SurfaceCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               2) {
      poParadisHome->param->SlipSystem2SurfaceCrossSlipEventsCount =
          poParadisHome->param->SlipSystem2SurfaceCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               3) {
      poParadisHome->param->SlipSystem3SurfaceCrossSlipEventsCount =
          poParadisHome->param->SlipSystem3SurfaceCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               4) {
      poParadisHome->param->SlipSystem4SurfaceCrossSlipEventsCount =
          poParadisHome->param->SlipSystem4SurfaceCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               5) {
      poParadisHome->param->SlipSystem5SurfaceCrossSlipEventsCount =
          poParadisHome->param->SlipSystem5SurfaceCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               6) {
      poParadisHome->param->SlipSystem6SurfaceCrossSlipEventsCount =
          poParadisHome->param->SlipSystem6SurfaceCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               7) {
      poParadisHome->param->SlipSystem7SurfaceCrossSlipEventsCount =
          poParadisHome->param->SlipSystem7SurfaceCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               8) {
      poParadisHome->param->SlipSystem8SurfaceCrossSlipEventsCount =
          poParadisHome->param->SlipSystem8SurfaceCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               9) {
      poParadisHome->param->SlipSystem9SurfaceCrossSlipEventsCount =
          poParadisHome->param->SlipSystem9SurfaceCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               10) {
      poParadisHome->param->SlipSystem10SurfaceCrossSlipEventsCount =
          poParadisHome->param->SlipSystem10SurfaceCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               11) {
      poParadisHome->param->SlipSystem11SurfaceCrossSlipEventsCount =
          poParadisHome->param->SlipSystem11SurfaceCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               12) {
      poParadisHome->param->SlipSystem12SurfaceCrossSlipEventsCount =
          poParadisHome->param->SlipSystem12SurfaceCrossSlipEventsCount + 1;
    } else {
      poParadisHome->param->OtherSlipSystemSurfaceCrossSlipEventsCount =
          poParadisHome->param->OtherSlipSystemSurfaceCrossSlipEventsCount + 1;
      // cout <<"others accessed in surface" <<endl;
      // cout <<"Value" <<
      // poParadisHome->param->OtherSlipSystemSurfaceCrossSlipEventsCount
      // <<endl;
    }
  }
}
void ParadisCrossSlipServer::HandleIntersectionCrossSlip(
    Home_t *poParadisHome, DislocationChain *poChain,
    list<DislocationChain *> *plpoAllChains) {
  // check the intersection topology validity, this function handles ONLY 23 and
  // 24 chains
  DislocationNetworkNode *poFirstNode = poChain->GetFirst();
  DislocationNetworkNode *poLastNode = poChain->GetLast();
  unsigned int iFirstNodeNeighboursCount = poFirstNode->GetAllNeighboursCount();
  unsigned int iLastNodeNeighboursCount = poLastNode->GetAllNeighboursCount();
  if (iFirstNodeNeighboursCount != 2 && iLastNodeNeighboursCount != 2) {
    return;
  }
  if (iFirstNodeNeighboursCount == 2 && iLastNodeNeighboursCount == 2) {
    return;
  }
  // now we are sure that we have either a 23 or a 24 chain, check to see if
  // it is topoligically admissible to cross slip it
  bool bIsTopologicallyAdmissibleToCrossSlip = false;
  DislocationChain *poIntersectionChain = NULL;
  Vector oJunctionBurgers;
  Vector oJunctionSlipPlane;
  if (iFirstNodeNeighboursCount == 3 || iLastNodeNeighboursCount == 3) {
    bIsTopologicallyAdmissibleToCrossSlip = Check23IntersectionChainCrossSlip(
        poParadisHome, poChain, plpoAllChains, poIntersectionChain,
        oJunctionBurgers, oJunctionSlipPlane);
  }
  // then, the 24 chains
  if (iFirstNodeNeighboursCount == 4 || iLastNodeNeighboursCount == 4) {
    bIsTopologicallyAdmissibleToCrossSlip = Check24IntersectionChainCrossSlip(
        poParadisHome, poChain, plpoAllChains, poIntersectionChain,
        oJunctionBurgers, oJunctionSlipPlane);
  }
  // proceed only if the cross slip is topologically admissible
  if (!bIsTopologicallyAdmissibleToCrossSlip) {
    return;
  }
  // set the correct burgers vector for the intersecting dislocation
  Vector oIntersectionDirection =
      poIntersectionChain->GetInitialLineDirection();
  Vector oChainDirection = poChain->GetInitialLineDirection();
  double dDotProduct = oIntersectionDirection * oChainDirection;
  if (dDotProduct < 0.0) {
    poIntersectionChain->Reverse();
  }
  Vector oChainBurgersVector = poChain->GetBurgersVector();
  Vector oIntersectionBurgersVector = poIntersectionChain->GetBurgersVector();
  Vector oTemp = oChainBurgersVector + oIntersectionBurgersVector;
  oTemp.Normalize();
  // detect the intersection cross slip type
  if (Is112Vector(oTemp)) {
    HandleRepulsiveCrossSlip(poParadisHome, poChain, poIntersectionChain);
  } else {
    HandleAttractiveCrossSlip(poParadisHome, poChain, poIntersectionChain,
                              oJunctionBurgers, oJunctionSlipPlane);
  }
}
bool ParadisCrossSlipServer::Check23IntersectionChainCrossSlip(
    Home_t *poParadisHome, DislocationChain *poChain,
    list<DislocationChain *> *plpoAllChains,
    DislocationChain *&poIntersectionChain, Vector &oJunctionBurgers,
    Vector &oJunctionSlipPlane) {
  poIntersectionChain = NULL;
  DislocationChain *poJunctionChain = NULL;
  DislocationNetworkNode *poFreeNode = NULL;
  DislocationNetworkNode *poIntersectionNode = NULL;
  if (poChain->GetFirst()->GetAllNeighboursCount() == 2) {
    poFreeNode = poChain->GetFirst();
    poIntersectionNode = poChain->GetLast();
  } else {
    poFreeNode = poChain->GetLast();
    poIntersectionNode = poChain->GetFirst();
  }
  DislocationChain *poTempChain1 = NULL;
  DislocationChain *poTempChain2 = NULL;
  list<DislocationChain *>::iterator liAllChains;

  for (liAllChains = plpoAllChains->begin();
       liAllChains != plpoAllChains->end(); liAllChains++) {
    if ((*liAllChains)->IsEndNode(poIntersectionNode)) {
      if ((*liAllChains) == poChain) {
        continue;
      }
      if (poTempChain1 == NULL) {
        poTempChain1 = (*liAllChains);
      } else {
        poTempChain2 = (*liAllChains);
        break;
      }
    }
  }
  // at this point, we should have got the 2 chains, check if this is not the
  // case
  if (poTempChain1 == NULL || poTempChain2 == NULL) {
    return false;
  }

  if (poTempChain1->GetType() == JunctionChain) {
    if (poTempChain2->GetType() == JunctionChain) {
      return false;
    } else {
      poJunctionChain = poTempChain1;
      poIntersectionChain = poTempChain2;
    }
  } else if (poTempChain2->GetType() == JunctionChain) {
    poJunctionChain = poTempChain2;
    poIntersectionChain = poTempChain1;
  } else {
    return false;
  }

  // now we know for sure which chain is the junction and which is the
  // intersection check to see that the intersection dislocation doesn't lie on
  // the same slip plane as the screw chain, this is why it is an intersection
  // in the 1st place
  Vector oChainBurgers = poChain->GetBurgersVector();
  Vector oChainNormal = poChain->GetSlipPlaneNormal();
  Vector oBurgers1 = poIntersectionChain->GetBurgersVector();
  Vector oNormal1 = poIntersectionChain->GetSlipPlaneNormal();
  double dTolerance = 1E-6;
  if ((oChainNormal ^ oNormal1).Length() < dTolerance) {
    return false;
  }
  // we need to make sure that the screw dislocation is continued at the other
  // side of the junction before anything else. for this check, it only suffices
  // to see if we have a chain ending at the other end of the junction chain
  // that has the same burgers vector and normal as the cross slip chain the
  // same normal requirement takes care of cross slipping only one side of the
  // screw dislocation not both of them
  DislocationNetworkNode *poOtherJunctionNode =
      poJunctionChain->GetOtherEnd(poIntersectionNode);
  if (poOtherJunctionNode == NULL) {
    return false;
  }

  // now we have the other end of the junction, get the 2 chains connected to it
  // we are sure that they should be exactly 2 chains because this is a junction
  // (33 chain)
  for (liAllChains = plpoAllChains->begin();
       liAllChains != plpoAllChains->end(); liAllChains++) {
    if ((*liAllChains)->IsEndNode(poOtherJunctionNode)) {
      if ((*liAllChains) == poJunctionChain) {
        continue;
      }
      if (poTempChain1 == NULL) {
        poTempChain1 = (*liAllChains);
      } else {
        poTempChain2 = (*liAllChains);
        break;
      }
    }
  }

  // we have the two candidate chains, see if any of them is the continuation of
  // the screw dislocation, doesn't matter if both of them are, this is just a
  // continuity check

  oBurgers1 = poTempChain1->GetBurgersVector();
  oNormal1 = poTempChain1->GetSlipPlaneNormal();
  Vector oBurgers2 = poTempChain2->GetBurgersVector();
  Vector oNormal2 = poTempChain2->GetSlipPlaneNormal();

  bool bIsContinuous = false;
  if ((oChainBurgers ^ oBurgers1).Length() < dTolerance) {
    if ((oChainNormal ^ oNormal1).Length() < dTolerance) {
      bIsContinuous = true;
    }
  }

  if ((oChainBurgers ^ oBurgers2).Length() < dTolerance) {
    if ((oChainNormal ^ oNormal2).Length() < dTolerance) {
      bIsContinuous = true;
    }
  }

  if (!bIsContinuous) {
    return false;
  }
  oJunctionBurgers = poJunctionChain->GetBurgersVector();
  oJunctionSlipPlane = poJunctionChain->GetSlipPlaneNormal();
  return true;
}
bool ParadisCrossSlipServer::Check24IntersectionChainCrossSlip(
    Home_t *poParadisHome, DislocationChain *poChain,
    list<DislocationChain *> *plpoAllChains,
    DislocationChain *&poIntersectionChain, Vector &oJunctionBurgers,
    Vector &oJunctionSlipPlane) {
  poIntersectionChain = NULL;
  DislocationNetworkNode *poFreeNode = NULL;
  DislocationNetworkNode *poIntersectionNode = NULL;
  if (poChain->GetFirst()->GetAllNeighboursCount() == 2) {
    poFreeNode = poChain->GetFirst();
    poIntersectionNode = poChain->GetLast();
  } else {
    poFreeNode = poChain->GetLast();
    poIntersectionNode = poChain->GetFirst();
  }
  unsigned int iNeighbouringChainsCount = 3;
  vector<DislocationChain *> vpoChains;
  vpoChains.resize(iNeighbouringChainsCount);
  unsigned int i = 0;
  for (i = 0; i < iNeighbouringChainsCount; i++) {
    vpoChains[i] = NULL;
  }
  list<DislocationChain *>::iterator liAllChains;
  i = 0;
  for (liAllChains = plpoAllChains->begin();
       liAllChains != plpoAllChains->end(); liAllChains++) {
    if ((*liAllChains)->IsEndNode(poIntersectionNode)) {
      if ((*liAllChains) == poChain) {
        continue;
      }
      vpoChains[i] = (*liAllChains);
      i = i + 1;
      if (i == iNeighbouringChainsCount) {
        break;
      }
    }
  }
  // if the 3 chains were not found, return
  if (i != iNeighbouringChainsCount) {
    return false;
  }

  // get the normals and the burgers vectors for all of the chains
  Vector oChainBurgers = poChain->GetBurgersVector();
  Vector oChainNormal = poChain->GetSlipPlaneNormal();
  vector<Vector> voNormals;
  vector<Vector> voBurgers;
  voNormals.resize(iNeighbouringChainsCount);
  voBurgers.resize(iNeighbouringChainsCount);
  for (i = 0; i < iNeighbouringChainsCount; i++) {
    voBurgers[i] = vpoChains[i]->GetBurgersVector();
    voNormals[i] = vpoChains[i]->GetSlipPlaneNormal();
  }
  // look first for the extension of the screw chain
  DislocationChain *poScrewChainExtension = NULL;
  double dTolerance = 1.0E-6;
  bool bFound = false;
  for (i = 0; i < iNeighbouringChainsCount; i++) {
    if ((voNormals[i] ^ oChainNormal).Length() < dTolerance) {
      if ((voBurgers[i] ^ oChainBurgers).Length() < dTolerance) {
        poScrewChainExtension = vpoChains[i];
        bFound = true;
        break;
      }
    }
  }
  // don't continue unless you find an extension
  if (!bFound) {
    return false;
  }

  // now look for an intersection, the other two chains must be on the same
  // plane and they must have the same burgers vector
  DislocationChain *poIntersectionChainExtension = NULL;
  bFound = false;
  bool bSatisfied = false;
  Vector oIntersectionNormal;
  Vector oIntersectionBurgers;
  for (i = 0; i < iNeighbouringChainsCount; i++) {
    if (vpoChains[i] == poScrewChainExtension) {
      continue;
    }
    if ((voNormals[i] ^ oChainNormal).Length() > dTolerance) {
      if (!bFound) {
        poIntersectionChain = vpoChains[i];
        oIntersectionNormal = voNormals[i];
        oIntersectionBurgers = voBurgers[i];
        bFound = true;
      } else {
        if ((voNormals[i] ^ oIntersectionNormal).Length() < dTolerance) {
          if ((voBurgers[i] ^ oIntersectionBurgers).Length() < dTolerance) {
            poIntersectionChainExtension = vpoChains[i];
            bSatisfied = true;
            break;
          }
        }
      }
    }
  }
  if (!bSatisfied) {
    return false;
  }
  oJunctionBurgers =
      poIntersectionChain->GetBurgersVector() + poChain->GetBurgersVector();
  Vector oJunctionDirection =
      poIntersectionChain->GetSlipPlaneNormal() ^ oChainNormal;
  oJunctionDirection.Normalize();
  oJunctionSlipPlane = oJunctionDirection ^ oJunctionBurgers;
  oJunctionSlipPlane.Normalize();
  return true;
}
void ParadisCrossSlipServer::HandleRepulsiveCrossSlip(
    Home_t *poParadisHome, DislocationChain *poChain,
    DislocationChain *poIntersectionChain) {
  Vector oSlipPlaneNormal = poChain->GetSlipPlaneNormal();
  Vector oBurgersVector = poChain->GetBurgersVector();
  Vector oCrossSlipPlaneNormal =
      GetFCCCrossSlipPlaneNormal(oSlipPlaneNormal, oBurgersVector);
  Vector oIntersectionNormal = poIntersectionChain->GetSlipPlaneNormal();
  Vector oSlipIntersectionCommon = oSlipPlaneNormal ^ oIntersectionNormal;
  Vector oCrossSlipIntersectionCommon =
      oCrossSlipPlaneNormal ^ oIntersectionNormal;
  double dTolerance = 1.0E-6;
  if (oSlipIntersectionCommon.Length() < dTolerance) {
    return;
  }
  if (oCrossSlipIntersectionCommon.Length() < dTolerance) {
    return;
  }
  oSlipIntersectionCommon.Normalize();
  oCrossSlipIntersectionCommon.Normalize();
  Vector oIntersectionDirection =
      poIntersectionChain->GetInitialLineDirection();
  if (fabs(oSlipIntersectionCommon * oIntersectionDirection) <
      fabs(oCrossSlipIntersectionCommon * oIntersectionDirection)) {
    // spontaneous cross slip, this is NOT thermally activated, so no need to
    // do any probability calculations
    if (AdjustChainForCrossSlip(poParadisHome, poChain)) {
      Vector oSlipPlaneNormal = poChain->GetSlipPlaneNormal();
      Vector oBurgersVector = poChain->GetBurgersVector();
      Vector oCrossSlipPlaneNormal =
          GetFCCCrossSlipPlaneNormal(oSlipPlaneNormal, oBurgersVector);
      CrossSlipChain(poParadisHome, poChain, oCrossSlipPlaneNormal);
      poParadisHome->param->RepulsiveCrossSlipEventsCount =
          poParadisHome->param->RepulsiveCrossSlipEventsCount + 1;
      poParadisHome->param->TotalRepulsiveCrossSlippedChainsLength =
          poParadisHome->param->TotalRepulsiveCrossSlippedChainsLength +
          poChain->GetLength();
      // cout <<"repulsive accessed" <<endl;
      // now add the Repulsive cross-slip event based on the slip system using
      // IdentifyFCCSlipSystem function
      if (poParadisHome->param->HandleCrossSlipPerSlipSystem == 0) {
        return;
      }
      if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) == 1) {
        poParadisHome->param->SlipSystem1RepulsiveCrossSlipEventsCount =
            poParadisHome->param->SlipSystem1RepulsiveCrossSlipEventsCount + 1;
      } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
                 2) {
        poParadisHome->param->SlipSystem2RepulsiveCrossSlipEventsCount =
            poParadisHome->param->SlipSystem2RepulsiveCrossSlipEventsCount + 1;
      } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
                 3) {
        poParadisHome->param->SlipSystem3RepulsiveCrossSlipEventsCount =
            poParadisHome->param->SlipSystem3RepulsiveCrossSlipEventsCount + 1;
      } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
                 4) {
        poParadisHome->param->SlipSystem4RepulsiveCrossSlipEventsCount =
            poParadisHome->param->SlipSystem4RepulsiveCrossSlipEventsCount + 1;
      } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
                 5) {
        poParadisHome->param->SlipSystem5RepulsiveCrossSlipEventsCount =
            poParadisHome->param->SlipSystem5RepulsiveCrossSlipEventsCount + 1;
      } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
                 6) {
        poParadisHome->param->SlipSystem6RepulsiveCrossSlipEventsCount =
            poParadisHome->param->SlipSystem6RepulsiveCrossSlipEventsCount + 1;
      } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
                 7) {
        poParadisHome->param->SlipSystem7RepulsiveCrossSlipEventsCount =
            poParadisHome->param->SlipSystem7RepulsiveCrossSlipEventsCount + 1;
      } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
                 8) {
        poParadisHome->param->SlipSystem8RepulsiveCrossSlipEventsCount =
            poParadisHome->param->SlipSystem8RepulsiveCrossSlipEventsCount + 1;
      } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
                 9) {
        poParadisHome->param->SlipSystem9RepulsiveCrossSlipEventsCount =
            poParadisHome->param->SlipSystem9RepulsiveCrossSlipEventsCount + 1;
      } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
                 10) {
        poParadisHome->param->SlipSystem10RepulsiveCrossSlipEventsCount =
            poParadisHome->param->SlipSystem10RepulsiveCrossSlipEventsCount + 1;
      } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
                 11) {
        poParadisHome->param->SlipSystem11RepulsiveCrossSlipEventsCount =
            poParadisHome->param->SlipSystem11RepulsiveCrossSlipEventsCount + 1;
      } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
                 12) {
        poParadisHome->param->SlipSystem12RepulsiveCrossSlipEventsCount =
            poParadisHome->param->SlipSystem12RepulsiveCrossSlipEventsCount + 1;
      } else {
        poParadisHome->param->OtherSlipSystemRepulsiveCrossSlipEventsCount =
            poParadisHome->param->OtherSlipSystemRepulsiveCrossSlipEventsCount +
            1;
        // cout <<"others accessed in repulsive" <<endl;
        // cout <<"Value" <<
        // poParadisHome->param->OtherSlipSystemRepulsiveCrossSlipEventsCount
        // <<endl;
      }
    }
  }
}
void ParadisCrossSlipServer::HandleAttractiveCrossSlip(
    Home_t *poParadisHome, DislocationChain *poChain,
    DislocationChain *poIntersectionChain, Vector &oJunctionBurgers,
    Vector &oJunctionSlipPlane) {
  Vector oChainBurgersVector = poChain->GetBurgersVector();
  Vector oIntersectionBurgersVector = poIntersectionChain->GetBurgersVector();
  Vector oTemp = oChainBurgersVector + oIntersectionBurgersVector;
  oTemp.Normalize();

  double dBurgersMagnitude = poParadisHome->param->burgMag;

  double dActivationEnergy = 0.0;
  double dActivationVolume = 0.0;
  double dCrossSlipLength = 0.0;
  double dCrossSlipReferenceLength = 0.0;
  double dCrossSlipFrequency = 0.0;
  double dGamma = 0.0;
  double dDoverBRatio = 0.0;

  if (Is100Vector(oTemp)) {
    // hirth lock
    dActivationEnergy = poParadisHome->param->HirthCrossSlipActivationEnergy;
    dActivationVolume =
        poParadisHome->param->HirthCrossSlipActivationVolumeFactor *
        dBurgersMagnitude * dBurgersMagnitude * dBurgersMagnitude;
    dCrossSlipLength = poParadisHome->param->HirthCrossSlipLength;
    dCrossSlipReferenceLength =
        poParadisHome->param->HirthCrossSlipReferenceLength;
    dCrossSlipFrequency = poParadisHome->param->HirthCrossSlipFrequency;
    dGamma = poParadisHome->param->HirthStackingFaultEnergy;
    dDoverBRatio = poParadisHome->param->HirthDoverBRatio;
  } else if (Is110Vector(oTemp)) {
    // attractive cross slip - lomer cottrel or glide lock
    if (Is111Vector(oJunctionSlipPlane)) {
      // glide lock
      dActivationEnergy =
          poParadisHome->param->GlideLockCrossSlipActivationEnergy;
      dActivationVolume =
          poParadisHome->param->GlideLockCrossSlipActivationVolumeFactor *
          dBurgersMagnitude * dBurgersMagnitude * dBurgersMagnitude;
      dCrossSlipLength = poParadisHome->param->GlideLockCrossSlipLength;
      dCrossSlipReferenceLength =
          poParadisHome->param->GlideLockCrossSlipReferenceLength;
      dCrossSlipFrequency = poParadisHome->param->GlideLockCrossSlipFrequency;
      dGamma = poParadisHome->param->GlideLockStackingFaultEnergy;
      dDoverBRatio = poParadisHome->param->GlideLockDoverBRatio;
    } else if (Is100Vector(oJunctionSlipPlane) ||
               Is110Vector(oJunctionSlipPlane)) {
      // LC lock
      dActivationEnergy = poParadisHome->param->LCLockCrossSlipActivationEnergy;
      dActivationVolume =
          poParadisHome->param->LCLockCrossSlipActivationVolumeFactor *
          dBurgersMagnitude * dBurgersMagnitude * dBurgersMagnitude;
      dCrossSlipLength = poParadisHome->param->LCLockCrossSlipLength;
      dCrossSlipReferenceLength =
          poParadisHome->param->LCLockCrossSlipReferenceLength;
      dCrossSlipFrequency = poParadisHome->param->LCLockCrossSlipFrequency;
      dGamma = poParadisHome->param->LCLockStackingFaultEnergy;
      dDoverBRatio = poParadisHome->param->LCLockDoverBRatio;
    } else {
      return;
    }
  }

  double dGlidePlaneShearStress = 0.0;
  double dCrossSlipPlaneShearStress = 0.0;
  double dGlideEscaigStress = 0.0;
  double dCrossSlipEscaigStress = 0.0;
  double dTotalChainLength = 0.0;
  GetChainCrossSlipValues(poParadisHome, poChain, dGlidePlaneShearStress,
                          dCrossSlipPlaneShearStress, dGlideEscaigStress,
                          dCrossSlipEscaigStress, dTotalChainLength, true);
  // now, to make a cross slip decision, we need to make sure that all of the
  // following conditions are satisfied (according to satish)
  // 1. the cross slip shear stress is < 0.5*glide plane shear stress
  // 2. the cross slip shear stress is > Gb/10/L
  // 3. the probability is high enough for cross slip at this instance

  double dCrossSlipShearStressThreshold =
      poParadisHome->param->shearModulus / 2.0 / dTotalChainLength;
  if (fabs(dCrossSlipPlaneShearStress) < dCrossSlipShearStressThreshold) {
    return;
  }
  if (fabs(dCrossSlipPlaneShearStress) < fabs(0.5 * dGlidePlaneShearStress)) {
    return;
  }

  double dBoltzmannConstant = 1.3806503E-23;
  double dTemperature = poParadisHome->param->TempK;
  // modifications requested by satish, february 2013
  double dBe = poParadisHome->param->burgMag / 2.0 / sqrt(3.0);
  double dGammaEffective = dGamma + dCrossSlipEscaigStress * dBe;
  double dTolerance = 1.0E-6;
  if (dGammaEffective < dTolerance) {
    dGammaEffective = dTolerance;
  }
  double dDOverBEffective = dDoverBRatio * dGamma / dGammaEffective;
  if (dDOverBEffective < 1.0) {
    dDOverBEffective = 1.0;
  }
  dActivationEnergy = dActivationEnergy * dGamma / dGammaEffective *
                      sqrt(log(dDOverBEffective) / log(dDoverBRatio));
  double dEscaigStressDifference = dGlideEscaigStress - dCrossSlipEscaigStress;

  double dProbabilityExponent =
      -(dActivationEnergy - dActivationVolume * dEscaigStressDifference) /
      dBoltzmannConstant / dTemperature;
  double dProbabilityExponent_noTemp =
      -(dActivationEnergy - dActivationVolume * dEscaigStressDifference) /
      dBoltzmannConstant;
  // intersection cross slip is checked every time step regardless of the actual
  // cross slip handling frequency
  // double dCrossSlipActualTime = poParadisHome->param->deltaTT;
  double dCrossSlipActualTime =
      poParadisHome->param->CrossSlipHandlingFrequency *
      poParadisHome->param->deltaTT;
  double dProbabilityBound = dCrossSlipFrequency * dCrossSlipActualTime *
                             dCrossSlipLength / dCrossSlipReferenceLength *
                             exp(dProbabilityExponent);
  double dRandomNumber = 0.0;

  if (poParadisHome->param->Ttype <= 0 ||
      poParadisHome->param->stress_timestep <= 0) {
    if (dProbabilityBound < 1.0) {
      dRandomNumber = Randomizer::Random();
      if (dRandomNumber > dProbabilityBound) {
        return;
      }
    }
  } else {
    dRandomNumber = Randomizer::Random();
  }

  if (AdjustChainForCrossSlip(poParadisHome, poChain)) {
    Vector oSlipPlaneNormal = poChain->GetSlipPlaneNormal();
    Vector oBurgersVector = poChain->GetBurgersVector();
    Vector oCrossSlipPlaneNormal =
        GetFCCCrossSlipPlaneNormal(oSlipPlaneNormal, oBurgersVector);
    // Yejun
    if (poParadisHome->param->Ttype > 0 &&
        poParadisHome->param->stress_timestep > 0) {
      CrossSlipChain_ThermalGradient(
          poParadisHome, poChain, oCrossSlipPlaneNormal,
          dCrossSlipFrequency * dCrossSlipActualTime * dCrossSlipLength /
              dCrossSlipReferenceLength,
          dProbabilityExponent_noTemp, dRandomNumber);
    } else {
      CrossSlipChain(poParadisHome, poChain, oCrossSlipPlaneNormal);
    }
    poParadisHome->param->AttractiveCrossSlipEventsCount =
        poParadisHome->param->AttractiveCrossSlipEventsCount + 1;
    poParadisHome->param->TotalAttractiveCrossSlippedChainsLength =
        poParadisHome->param->TotalAttractiveCrossSlippedChainsLength +
        poChain->GetLength();
    // cout <<"attractive accessed" <<endl;
    // now add the Attractive cross-slip event based on the slip system using
    // IdentifyFCCSlipSystem function
    if (poParadisHome->param->HandleCrossSlipPerSlipSystem == 0) {
      return;
    }
    if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) == 1) {
      poParadisHome->param->SlipSystem1AttractiveCrossSlipEventsCount =
          poParadisHome->param->SlipSystem1AttractiveCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               2) {
      poParadisHome->param->SlipSystem2AttractiveCrossSlipEventsCount =
          poParadisHome->param->SlipSystem2AttractiveCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               3) {
      poParadisHome->param->SlipSystem3AttractiveCrossSlipEventsCount =
          poParadisHome->param->SlipSystem3AttractiveCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               4) {
      poParadisHome->param->SlipSystem4AttractiveCrossSlipEventsCount =
          poParadisHome->param->SlipSystem4AttractiveCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               5) {
      poParadisHome->param->SlipSystem5AttractiveCrossSlipEventsCount =
          poParadisHome->param->SlipSystem5AttractiveCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               6) {
      poParadisHome->param->SlipSystem6AttractiveCrossSlipEventsCount =
          poParadisHome->param->SlipSystem6AttractiveCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               7) {
      poParadisHome->param->SlipSystem7AttractiveCrossSlipEventsCount =
          poParadisHome->param->SlipSystem7AttractiveCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               8) {
      poParadisHome->param->SlipSystem8AttractiveCrossSlipEventsCount =
          poParadisHome->param->SlipSystem8AttractiveCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               9) {
      poParadisHome->param->SlipSystem9AttractiveCrossSlipEventsCount =
          poParadisHome->param->SlipSystem9AttractiveCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               10) {
      poParadisHome->param->SlipSystem10AttractiveCrossSlipEventsCount =
          poParadisHome->param->SlipSystem10AttractiveCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               11) {
      poParadisHome->param->SlipSystem11AttractiveCrossSlipEventsCount =
          poParadisHome->param->SlipSystem11AttractiveCrossSlipEventsCount + 1;
    } else if (IdentifyFCCSlipSystem(oCrossSlipPlaneNormal, oBurgersVector) ==
               12) {
      poParadisHome->param->SlipSystem12AttractiveCrossSlipEventsCount =
          poParadisHome->param->SlipSystem12AttractiveCrossSlipEventsCount + 1;
    } else {
      poParadisHome->param->OtherSlipSystemAttractiveCrossSlipEventsCount =
          poParadisHome->param->OtherSlipSystemAttractiveCrossSlipEventsCount +
          1;
      // cout <<"others accessed in attractive" <<endl;
      // cout <<"Value" <<
      // poParadisHome->param->OtherSlipSystemAttractiveCrossSlipEventsCount
      // <<endl;
    }
  }
}
void ParadisCrossSlipServer::GetChainCrossSlipValues(
    Home_t *poParadisHome, DislocationChain *poChain,
    double &dAverageGlidePlaneShearStress,
    double &dAverageCrossSlipPlaneShearStress, double &dGlideEscaigStress,
    double &dCrossSlipEscaigStress, double &dTotalLength,
    bool bExcludeLocalEscaigStress) {
  Vector oSlipPlaneNormal = poChain->GetSlipPlaneNormal();
  // the normal has to follow the convention
  if (oSlipPlaneNormal.GetX() * oSlipPlaneNormal.GetY() *
          oSlipPlaneNormal.GetZ() <
      0.0) {
    oSlipPlaneNormal = oSlipPlaneNormal * (-1.0);
  }
  Vector oBurgersVector = poChain->GetBurgersVector();
  Vector oCrossSlipPlaneNormal =
      GetFCCCrossSlipPlaneNormal(oSlipPlaneNormal, oBurgersVector);
  if (oCrossSlipPlaneNormal.GetX() * oCrossSlipPlaneNormal.GetY() *
          oCrossSlipPlaneNormal.GetZ() <
      0.0) {
    oCrossSlipPlaneNormal = oCrossSlipPlaneNormal * (-1.0);
  }
  // follow the convention of the coordinate system proposed by Satish
  Vector oGPXAxis = oBurgersVector;
  Vector oGPTemp = oSlipPlaneNormal ^ oBurgersVector;
  Vector oCSXAxis = oBurgersVector;
  Vector oCSTemp = oCrossSlipPlaneNormal ^ oBurgersVector;
  bool bReverseSign = false;
  double dTolerance = 1.0E-6;
  // fix glide plane x axis
  if ((fabs(oGPXAxis.GetX()) < dTolerance) &&
      (oGPTemp.GetX() * oSlipPlaneNormal.GetX() < 0.0)) {
    bReverseSign = true;
  }
  if ((fabs(oGPXAxis.GetY()) < dTolerance) &&
      (oGPTemp.GetY() * oSlipPlaneNormal.GetY() < 0.0)) {
    bReverseSign = true;
  }
  if ((fabs(oGPXAxis.GetZ()) < dTolerance) &&
      (oGPTemp.GetZ() * oSlipPlaneNormal.GetZ() < 0.0)) {
    bReverseSign = true;
  }
  if (bReverseSign) {
    oGPXAxis = oGPXAxis * (-1.0);
  }
  // same for the cross slip plane x axis
  bReverseSign = false;
  if ((fabs(oCSXAxis.GetX()) < dTolerance) &&
      (oCSTemp.GetX() * oCrossSlipPlaneNormal.GetX() < 0.0)) {
    bReverseSign = true;
  }
  if ((fabs(oCSXAxis.GetY()) < dTolerance) &&
      (oCSTemp.GetY() * oCrossSlipPlaneNormal.GetY() < 0.0)) {
    bReverseSign = true;
  }
  if ((fabs(oCSXAxis.GetZ()) < dTolerance) &&
      (oCSTemp.GetZ() * oCrossSlipPlaneNormal.GetZ() < 0.0)) {
    bReverseSign = true;
  }
  if (bReverseSign) {
    oCSXAxis = oCSXAxis * (-1.0);
  }
  oGPXAxis.Normalize();
  oCSXAxis.Normalize();

  DislocationNetworkNode *poCurrentChainNode = NULL;
  DislocationNetworkNode *poNextChainNode = NULL;

  unsigned int i = 0;
  unsigned int iChainLength = poChain->GetSize();
  DislocationNode *poCurrentDislocationNode = NULL;
  DislocationNode *poNextDislocationNode = NULL;
  Vector oTotalForce(0.0, 0.0, 0.0);
  Node_t *poNode = NULL;
  Node_t *poNeighbour = NULL;
  int iArmID = 0;
  double dCurrentSegmentLength = 0.0;
  double dTotalChainLength = 0.0;
  poChain->ResetIterator();
  for (i = 0; i < iChainLength; i++) {
    poCurrentChainNode = poChain->GetCurrentNode();
    poNextChainNode = poChain->GetNextNode();
    if (poNextChainNode == NULL) {
      continue;
    }
    poCurrentDislocationNode = poCurrentChainNode->GetDataPointer();
    poNextDislocationNode = poNextChainNode->GetDataPointer();

    poNode = GetNodeFromIndex(poParadisHome, poParadisHome->myDomain,
                              poCurrentDislocationNode->GetID());
    poNeighbour = GetNodeFromIndex(poParadisHome, poParadisHome->myDomain,
                                   poNextDislocationNode->GetID());

    if ((poNode == NULL) || (poNeighbour == NULL)) {
      continue;
    }
    iArmID = GetArmID(poNode, poNeighbour);
    if (iArmID == -1) {
      continue;
    }
    dCurrentSegmentLength =
        poCurrentDislocationNode->Distance(*poNextDislocationNode);
    dTotalChainLength = dTotalChainLength + dCurrentSegmentLength;
    oTotalForce =
        oTotalForce +
        Vector(poNode->fX + poNeighbour->fX, poNode->fY + poNeighbour->fY,
               poNode->fZ + poNeighbour->fZ) *
            0.5;
    poChain->IncrementIterator();
  }
  // according to Satish, we need the force component normal to the
  // dislocation line on both the glide and cross slip planes (the edge
  // dislocation burgers vector direction)
  Vector oGlidePlaneDislocationNormal = oBurgersVector ^ oSlipPlaneNormal;
  Vector oCrossSlipPlaneDislocationNormal =
      oBurgersVector ^ oCrossSlipPlaneNormal;
  oGlidePlaneDislocationNormal.Normalize();
  oCrossSlipPlaneDislocationNormal.Normalize();
  double dBurgersMagnitude = poParadisHome->param->burgMag;
  dAverageGlidePlaneShearStress =
      (oTotalForce * oGlidePlaneDislocationNormal) / dTotalChainLength;
  dAverageCrossSlipPlaneShearStress =
      (oTotalForce * oCrossSlipPlaneDislocationNormal) / dTotalChainLength;
  // get the escaig stress in the glide plane
  dGlideEscaigStress = GetEscaigStress(oTotalForce, oGPXAxis, oSlipPlaneNormal,
                                       dBurgersMagnitude);
  // get the escaig stress in the cross slip plane
  dCrossSlipEscaigStress = GetEscaigStress(
      oTotalForce, oCSXAxis, oCrossSlipPlaneNormal, dBurgersMagnitude);
  dTotalLength = dTotalChainLength;
}
double ParadisCrossSlipServer::GetEscaigStress(
    const Vector &oForce, const Vector &oBurgersDirection,
    const Vector &oNormal,
    const double &dBurgersMagnitude) // escaig stress calculation
{
  double daB[3];
  daB[0] = oBurgersDirection.GetX();
  daB[1] = oBurgersDirection.GetY();
  daB[2] = oBurgersDirection.GetZ();

  Matrix oS(3, 3);
  unsigned int i = 0;
  unsigned int j = 0;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      oS.Set(i + 1, j + 1, daB[i] * daB[j]);
    }
  }
  oS = oS * dBurgersMagnitude;
  oS = oS.GetInverse();

  double daS[3][3];
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      daS[i][j] = oS.Get(i + 1, j + 1);
    }
  }

  double daN[3];
  daN[0] = oNormal.GetX();
  daN[1] = oNormal.GetY();
  daN[2] = oNormal.GetZ();
  double daR[3][3][3];
  unsigned int k = 0;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        daR[i][j][k] = daN[i] * daB[j] * daN[k];
      }
    }
  }

  double daFC[3];
  daFC[0] = 0.0;
  daFC[1] = 0.0;
  daFC[2] = 0.0;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      daFC[0] = daFC[0] + daS[i][j] * daR[i][j][0];
      daFC[1] = daFC[1] + daS[i][j] * daR[i][j][1];
      daFC[2] = daFC[2] + daS[i][j] * daR[i][j][2];
    }
  }
  double dEscaigStress = daFC[0] * oForce.GetX() + daFC[1] * oForce.GetY() +
                         daFC[2] * oForce.GetZ();
  return dEscaigStress;
}
Vector ParadisCrossSlipServer::GetFCCCrossSlipPlaneNormal(
    const Vector &oSlipPlaneNormal, const Vector &oBurgersVector) {
  Vector oNormal(0.0, 0.0, 0.0);
  double dTolerance = 1E-6;
  Vector oNormalDirection = oSlipPlaneNormal.GetDirection();
  Vector oBurgersDirection = oBurgersVector.GetDirection();
  if (fabs(oNormalDirection * oBurgersDirection) < dTolerance) {
    if (fabs(oBurgersDirection.GetX()) < dTolerance) {
      oNormal.Set(-oNormalDirection.GetX(), oNormalDirection.GetY(),
                  oNormalDirection.GetZ());
    } else if (fabs(oBurgersDirection.GetY()) < dTolerance) {
      oNormal.Set(oNormalDirection.GetX(), -oNormalDirection.GetY(),
                  oNormalDirection.GetZ());
    } else if (fabs(oBurgersDirection.GetZ()) < dTolerance) {
      oNormal.Set(oNormalDirection.GetX(), oNormalDirection.GetY(),
                  -oNormalDirection.GetZ());
    }
  }
  return oNormal;
}
bool ParadisCrossSlipServer::GetCrossSlipFitLine(DislocationChain *poChain,
                                                 Line &oFitLine) {
  DislocationNetworkNode *poNode1 = poChain->GetFirst();
  DislocationNetworkNode *poNode2 = poChain->GetLast();
  unsigned int iConstraint1 = 0;
  unsigned int iConstraint2 = 0;
  Vector oTemp;
  GetNodeDynamicConstraint(poNode1, iConstraint1, oTemp);
  GetNodeDynamicConstraint(poNode2, iConstraint2, oTemp);
  bool bNode1Moveable =
      (poNode1->GetDataPointer()->GetCategory() == UnconstrainedNode) &&
      (iConstraint1 == 1);
  bool bNode2Moveable =
      (poNode2->GetDataPointer()->GetCategory() == UnconstrainedNode) &&
      (iConstraint2 == 1);
  Point oPivotPoint;

  if (bNode1Moveable) {
    if (bNode2Moveable) {
      unsigned int iChainSize = poChain->GetSize();
      unsigned int i = 0;
      poChain->ResetIterator();
      oPivotPoint.Set(0.0, 0.0, 0.0);
      for (i = 0; i < iChainSize; i++) {
        poNode1 = poChain->GetCurrentNode();
        oPivotPoint = oPivotPoint + *(poNode1->GetDataPointer());
        poChain->IncrementIterator();
      }
      oPivotPoint = oPivotPoint * (1.0 / (double)iChainSize);
    } else {
      oPivotPoint = *(poNode2->GetDataPointer());
    }
  } else {
    if (bNode2Moveable) {
      oPivotPoint = *(poNode1->GetDataPointer());
    } else {
      return false;
    }
  }

  Vector oBurgersVector = poChain->GetBurgersVector();
  oFitLine = Line(oBurgersVector, oPivotPoint);
  return true;
}
bool ParadisCrossSlipServer::AdjustChainForCrossSlip(
    Home_t *poParadisHome, DislocationChain *poChain) {
  Line oBestFitLine;
  if (!GetCrossSlipFitLine(poChain, oBestFitLine))
    return false;
  // make sure that all the dislocation chain nodes are NOT inside any
  // precipitates
  unsigned int i = 0;
  unsigned int iChainLength = poChain->GetSize();
  poChain->ResetIterator();
  Polyhedron *poPrecipitate = NULL;
  for (i = 0; i < iChainLength; i++) {
    if (poParadisHome->poPrecipitateServer->IsPointInsideAPrecipitate(
            poChain->GetCurrentNode()->GetDataPointer(), poPrecipitate)) {
      return false;
    }
    poChain->IncrementIterator();
  }

  DislocationNetworkNode *poChainNode = NULL;
  DislocationNode *poDislocationNode = NULL;
  Point oProjection;
  poChain->ResetIterator();
  for (i = 0; i < iChainLength; i++) {
    poChainNode = poChain->GetCurrentNode();
    poDislocationNode = poChainNode->GetDataPointer();
    oProjection = oBestFitLine.GetPointProjection(*poDislocationNode);
    poDislocationNode->Set(oProjection.GetX(), oProjection.GetY(),
                           oProjection.GetZ());
    poChain->IncrementIterator();
  }
  return true;
}
void ParadisCrossSlipServer::CrossSlipChain(
    Home_t *poParadisHome, DislocationChain *poChain,
    const Vector &oCrossSlipPlaneNormal) // needed function
{
  DislocationNetworkNode *poCurrentChainNode = NULL;
  DislocationNetworkNode *poNextChainNode = NULL;
  DislocationNode *poCurrentDislocationNode = NULL;
  DislocationNode *poNextDislocationNode = NULL;
  Node_t *poNode = NULL;
  Node_t *poNextNode = NULL;
  Vector oNormal;
  Vector oBurgers;
  // now loop over all of the chain nodes, reposition them and switch their
  // glide planes
  unsigned int i = 0;
  unsigned int iChainLength = poChain->GetSize();
  double daTemp[3] = {0.0, 0.0, 0.0};
  poChain->ResetIterator();
  printf("@ %d : cross slipping chain of length %d\n", poParadisHome->myDomain,
         iChainLength);
  for (i = 0; i < iChainLength; i++) {
    poCurrentChainNode = poChain->GetCurrentNode();
    poNextChainNode = poChain->GetNextNode();
    if (poNextChainNode == NULL) {
      continue;
    }
    poCurrentDislocationNode = poCurrentChainNode->GetDataPointer();
    poNextDislocationNode = poNextChainNode->GetDataPointer();

    poNode = GetNodeFromIndex(poParadisHome, poParadisHome->myDomain,
                              poCurrentDislocationNode->GetID());
    poNextNode = GetNodeFromIndex(poParadisHome, poParadisHome->myDomain,
                                  poNextDislocationNode->GetID());

    if ((poNode == NULL) || (poNextNode == NULL)) {
      continue;
    }

    daTemp[0] =
        poCurrentDislocationNode->GetX(); // get coardinates of the current node
    daTemp[1] = poCurrentDislocationNode->GetY();
    daTemp[2] = poCurrentDislocationNode->GetZ();
    RepositionNode(poParadisHome, daTemp, &(poNode->myTag), 1);
    poNode->flags |= NODE_NO_CROSS_SLIP;

    daTemp[0] =
        poNextDislocationNode->GetX(); // get coardinates of the next node
    daTemp[1] = poNextDislocationNode->GetY();
    daTemp[2] = poNextDislocationNode->GetZ();
    RepositionNode(poParadisHome, daTemp, &(poNextNode->myTag), 1);
    poNextNode->flags |= NODE_NO_CROSS_SLIP;

    daTemp[0] = oCrossSlipPlaneNormal.GetX(); // get crossSlip plane normal
    daTemp[1] = oCrossSlipPlaneNormal.GetY();
    daTemp[2] = oCrossSlipPlaneNormal.GetZ();
    ResetGlidePlane(poParadisHome, daTemp, &(poNode->myTag),
                    &(poNextNode->myTag), 1);
    // oNormal = poChain->GetSlipPlaneNormal();
    // oBurgers = poChain->GetBurgersVector();
    // IdentifyFCCSlipSystem(oNormal,oBurgers);

    poChain->IncrementIterator();
  }
  printf("@ %d : done cross slipping chain\n", poParadisHome->myDomain);
}

// Yejun
void ParadisCrossSlipServer::CrossSlipChain_ThermalGradient(
    Home_t *poParadisHome, DislocationChain *poChain,
    const Vector &oCrossSlipPlaneNormal, double preProb, double Exponent_noTemp,
    double randProb) // needed function
{
  DislocationNetworkNode *poCurrentChainNode = NULL;
  DislocationNetworkNode *poNextChainNode = NULL;
  DislocationNode *poCurrentDislocationNode = NULL;
  DislocationNode *poNextDislocationNode = NULL;
  Node_t *poNode = NULL;
  Node_t *poNextNode = NULL;
  Vector oNormal;
  Vector oBurgers;
  // now loop over all of the chain nodes, reposition them and switch their
  // glide planes
  unsigned int i = 0;
  unsigned int iChainLength = poChain->GetSize();
  double daTemp[3] = {0.0, 0.0, 0.0};
  double daTemp1[3] = {0.0, 0.0, 0.0};
  double dTemperature = 0.0;
  double pos_mid[3];
  int n_x, n_y, n_z;
  int TempCounter = 0;
  int s_timestep;
  double xd, yd, zd;
  double pos_stress0[3], pos_stress1[3];
  double c000, c001, c010, c011, c100, c101, c110, c111;

  poChain->ResetIterator();

  for (i = 0; i < iChainLength; i++) {
    poCurrentChainNode = poChain->GetCurrentNode();
    poNextChainNode = poChain->GetNextNode();
    if (poNextChainNode == NULL) {
      continue;
    }
    poCurrentDislocationNode = poCurrentChainNode->GetDataPointer();
    poNextDislocationNode = poNextChainNode->GetDataPointer();

    poNode = GetNodeFromIndex(poParadisHome, poParadisHome->myDomain,
                              poCurrentDislocationNode->GetID());
    poNextNode = GetNodeFromIndex(poParadisHome, poParadisHome->myDomain,
                                  poNextDislocationNode->GetID());

    if ((poNode == NULL) || (poNextNode == NULL)) {
      continue;
    }

    daTemp[0] =
        poCurrentDislocationNode->GetX(); // get coardinates of the current node
    daTemp[1] = poCurrentDislocationNode->GetY();
    daTemp[2] = poCurrentDislocationNode->GetZ();

    daTemp1[0] =
        poNextDislocationNode->GetX(); // get coardinates of the next node
    daTemp1[1] = poNextDislocationNode->GetY();
    daTemp1[2] = poNextDislocationNode->GetZ();
    pos_mid[0] = 0.5 * daTemp[0] + 0.5 * daTemp1[0];
    pos_mid[1] = 0.5 * daTemp[1] + 0.5 * daTemp1[1];
    pos_mid[2] = 0.5 * daTemp[2] + 0.5 * daTemp1[2];

    GetInLocalCoordinates(poParadisHome->param, pos_mid[0], pos_mid[1],
                          pos_mid[2]);
    if (fabs(pos_mid[0]) > 0.5 * poParadisHome->param->stress_Dim[0] ||
        fabs(pos_mid[1]) > 0.5 * poParadisHome->param->stress_Dim[1] ||
        fabs(pos_mid[2]) > 0.5 * poParadisHome->param->stress_Dim[2]) {
      dTemperature += poParadisHome->param->TempK;
    } else {
      n_x = ceil((pos_mid[0] + 0.5 * poParadisHome->param->stress_Dim[0]) /
                 poParadisHome->param->stress_Dim[0] *
                 (double)(poParadisHome->param->GP_x - 1));
      n_y = ceil((pos_mid[1] + 0.5 * poParadisHome->param->stress_Dim[1]) /
                 poParadisHome->param->stress_Dim[1] *
                 (double)(poParadisHome->param->GP_y - 1));
      n_z = ceil((pos_mid[2] + 0.5 * poParadisHome->param->stress_Dim[2]) /
                 poParadisHome->param->stress_Dim[2] *
                 (double)(poParadisHome->param->GP_z - 1));

      n_x = fmax(n_x, 1);
      n_y = fmax(n_y, 1);
      n_z = fmax(n_z, 1);

      s_timestep = poParadisHome->param->stress_timestep - 1;
      pos_stress0[0] =
          poParadisHome->stress[n_x - 1][n_y - 1][n_z - 1][s_timestep].x;
      pos_stress0[1] =
          poParadisHome->stress[n_x - 1][n_y - 1][n_z - 1][s_timestep].y;
      pos_stress0[2] =
          poParadisHome->stress[n_x - 1][n_y - 1][n_z - 1][s_timestep].z;
      pos_stress1[0] = poParadisHome->stress[n_x][n_y][n_z][s_timestep].x;
      pos_stress1[1] = poParadisHome->stress[n_x][n_y][n_z][s_timestep].y;
      pos_stress1[2] = poParadisHome->stress[n_x][n_y][n_z][s_timestep].z;

      // trilinear interpolation
      xd = (pos_mid[0] - pos_stress0[0]) / (pos_stress1[0] - pos_stress0[0]);
      yd = (pos_mid[1] - pos_stress0[1]) / (pos_stress1[1] - pos_stress0[1]);
      zd = (pos_mid[2] - pos_stress0[2]) / (pos_stress1[2] - pos_stress0[2]);
      c000 = poParadisHome->stress[n_x - 1][n_y - 1][n_z - 1][s_timestep].Temp;
      c001 = poParadisHome->stress[n_x - 1][n_y - 1][n_z][s_timestep].Temp;
      c010 = poParadisHome->stress[n_x - 1][n_y][n_z - 1][s_timestep].Temp;
      c011 = poParadisHome->stress[n_x - 1][n_y][n_z][s_timestep].Temp;
      c100 = poParadisHome->stress[n_x][n_y - 1][n_z - 1][s_timestep].Temp;
      c101 = poParadisHome->stress[n_x][n_y - 1][n_z][s_timestep].Temp;
      c110 = poParadisHome->stress[n_x][n_y][n_z - 1][s_timestep].Temp;
      c111 = poParadisHome->stress[n_x][n_y][n_z][s_timestep].Temp;

      dTemperature += ((c000 * (1 - xd) + c100 * xd) * (1 - yd) +
                       (c010 * (1 - xd) + c110 * xd) * yd) *
                          (1 - zd) +
                      ((c001 * (1 - xd) + c101 * xd) * (1 - yd) +
                       (c011 * (1 - xd) + c111 * xd) * yd) *
                          zd;
      // printf("%d %lf\n", i,
      // ((c000*(1-xd)+c100*xd)*(1-yd)+(c010*(1-xd)+c110*xd)*yd)*(1-zd)+((c001*(1-xd)+c101*xd)*(1-yd)+(c011*(1-xd)+c111*xd)*yd)*zd);
    }
    TempCounter++;
    poChain->IncrementIterator();
  }

  dTemperature = dTemperature / (double)TempCounter;
  double dProbabilityBound = preProb * exp(Exponent_noTemp / dTemperature);
  if (dProbabilityBound >= 1.0 ||
      (dProbabilityBound < 1.0 && randProb <= dProbabilityBound)) {
    poChain->ResetIterator();
    printf("@ %d : cross slipping chain of length %d at T= %lf K\n",
           poParadisHome->myDomain, iChainLength, dTemperature);

    /*if (poParadisHome->param->BoundaryType != RIGID_BOUNDARY)
    {
            for (i = 0; i < iChainLength; i++)
            {
                    poCurrentChainNode = poChain->GetCurrentNode();
                    poNextChainNode = poChain->GetNextNode();
                    if (poNextChainNode == NULL)
                    {
                            continue;
                    }
                    poCurrentDislocationNode =
    poCurrentChainNode->GetDataPointer(); poNextDislocationNode =
    poNextChainNode->GetDataPointer();

                    poNode = GetNodeFromIndex(poParadisHome,
    poParadisHome->myDomain, poCurrentDislocationNode->GetID()); poNextNode =
    GetNodeFromIndex(poParadisHome, poParadisHome->myDomain,
    poNextDislocationNode->GetID());

                    if ((poNode == NULL) || (poNextNode == NULL))
                    {
                            continue;
                    }

                    daTemp[0] = poCurrentDislocationNode->GetX(); // get
    coardinates of the current node daTemp[1] =
    poCurrentDislocationNode->GetY(); daTemp[2] =
    poCurrentDislocationNode->GetZ(); RepositionNode(poParadisHome, daTemp,
    &(poNode->myTag), 1); poNode->flags |= NODE_NO_CROSS_SLIP;

                    daTemp[0] = poNextDislocationNode->GetX(); // get
    coardinates of the next node daTemp[1] = poNextDislocationNode->GetY();
                    daTemp[2] = poNextDislocationNode->GetZ();
                    RepositionNode(poParadisHome, daTemp, &(poNextNode->myTag),
    1); poNextNode->flags |= NODE_NO_CROSS_SLIP;

                    daTemp[0] = oCrossSlipPlaneNormal.GetX(); // get crossSlip
    plane normal daTemp[1] = oCrossSlipPlaneNormal.GetY(); daTemp[2] =
    oCrossSlipPlaneNormal.GetZ(); ResetGlidePlane(poParadisHome, daTemp,
    &(poNode->myTag), &(poNextNode->myTag), 1);
                    //oNormal = poChain->GetSlipPlaneNormal();
                    //oBurgers = poChain->GetBurgersVector();
                    //IdentifyFCCSlipSystem(oNormal,oBurgers);

                    poChain->IncrementIterator();
            }
    }

    else
    {
            for (i = 0; i < iChainLength; i++)
            {
                    poCurrentChainNode = poChain->GetCurrentNode();
                    poNextChainNode = poChain->GetNextNode();
                    if (poNextChainNode == NULL)
                    {
                            continue;
                    }
                    poCurrentDislocationNode =
    poCurrentChainNode->GetDataPointer(); poNextDislocationNode =
    poNextChainNode->GetDataPointer();

                    poNode = GetNodeFromIndex(poParadisHome,
    poParadisHome->myDomain, poCurrentDislocationNode->GetID()); poNextNode =
    GetNodeFromIndex(poParadisHome, poParadisHome->myDomain,
    poNextDislocationNode->GetID());

                    if ((poNode == NULL) || (poNextNode == NULL))
                    {
                            continue;
                    }

                    daTemp[0] = poCurrentDislocationNode->GetX(); // get
    coardinates of the current node daTemp[1] =
    poCurrentDislocationNode->GetY(); daTemp[2] =
    poCurrentDislocationNode->GetZ(); daTemp1[0] =
    poNextDislocationNode->GetX(); // get coardinates of the next node
                    daTemp1[1] = poNextDislocationNode->GetY();
                    daTemp1[2] = poNextDislocationNode->GetZ();

                    GetInLocalCoordinates(poParadisHome->param, daTemp[0],
    daTemp[1], daTemp[2]); GetInLocalCoordinates(poParadisHome->param,
    daTemp1[0], daTemp1[1], daTemp1[2]);

                    if (fabs(daTemp[0]) <= 0.5 *
    poParadisHome->param->Dimensions[0] && fabs(daTemp[1]) <= 0.5 *
    poParadisHome->param->Dimensions[1] && fabs(daTemp[2]) <= 0.5 *
    poParadisHome->param->Dimensions[2] && fabs(daTemp1[0]) <= 0.5 *
    poParadisHome->param->Dimensions[0] && fabs(daTemp1[1]) <= 0.5 *
    poParadisHome->param->Dimensions[1] && fabs(daTemp1[2]) <= 0.5 *
    poParadisHome->param->Dimensions[2])
                    {
                            GetInGlobalCoordinates(poParadisHome->param,
    daTemp[0], daTemp[1], daTemp[2]);
                            GetInGlobalCoordinates(poParadisHome->param,
    daTemp1[0], daTemp1[1], daTemp1[2]);

                            RepositionNode(poParadisHome, daTemp,
    &(poNode->myTag), 1); RepositionNode(poParadisHome, daTemp1,
    &(poNextNode->myTag), 1); daTemp[0] = oCrossSlipPlaneNormal.GetX(); // get
    crossSlip plane normal daTemp[1] = oCrossSlipPlaneNormal.GetY(); daTemp[2] =
    oCrossSlipPlaneNormal.GetZ(); ResetGlidePlane(poParadisHome, daTemp,
    &(poNode->myTag), &(poNextNode->myTag), 1);
                            //oNormal = poChain->GetSlipPlaneNormal();
                            //oBurgers = poChain->GetBurgersVector();
                            //IdentifyFCCSlipSystem(oNormal,oBurgers);
                    }
                    poNode->flags |= NODE_NO_CROSS_SLIP;
                    poNextNode->flags |= NODE_NO_CROSS_SLIP;

                    poChain->IncrementIterator();
            }
    }*/
    for (i = 0; i < iChainLength; i++) {
      poCurrentChainNode = poChain->GetCurrentNode();
      poNextChainNode = poChain->GetNextNode();
      if (poNextChainNode == NULL) {
        continue;
      }
      poCurrentDislocationNode = poCurrentChainNode->GetDataPointer();
      poNextDislocationNode = poNextChainNode->GetDataPointer();

      poNode = GetNodeFromIndex(poParadisHome, poParadisHome->myDomain,
                                poCurrentDislocationNode->GetID());
      poNextNode = GetNodeFromIndex(poParadisHome, poParadisHome->myDomain,
                                    poNextDislocationNode->GetID());

      if ((poNode == NULL) || (poNextNode == NULL)) {
        continue;
      }

      daTemp[0] = poCurrentDislocationNode
                      ->GetX(); // get coardinates of the current node
      daTemp[1] = poCurrentDislocationNode->GetY();
      daTemp[2] = poCurrentDislocationNode->GetZ();
      RepositionNode(poParadisHome, daTemp, &(poNode->myTag), 1);
      poNode->flags |= NODE_NO_CROSS_SLIP;

      daTemp[0] =
          poNextDislocationNode->GetX(); // get coardinates of the next node
      daTemp[1] = poNextDislocationNode->GetY();
      daTemp[2] = poNextDislocationNode->GetZ();
      RepositionNode(poParadisHome, daTemp, &(poNextNode->myTag), 1);
      poNextNode->flags |= NODE_NO_CROSS_SLIP;

      daTemp[0] = oCrossSlipPlaneNormal.GetX(); // get crossSlip plane normal
      daTemp[1] = oCrossSlipPlaneNormal.GetY();
      daTemp[2] = oCrossSlipPlaneNormal.GetZ();
      ResetGlidePlane(poParadisHome, daTemp, &(poNode->myTag),
                      &(poNextNode->myTag), 1);
      // oNormal = poChain->GetSlipPlaneNormal();
      // oBurgers = poChain->GetBurgersVector();
      // IdentifyFCCSlipSystem(oNormal,oBurgers);

      poChain->IncrementIterator();
    }

    printf("@ %d : done cross slipping chain at T=%lf K\n",
           poParadisHome->myDomain, dTemperature);
  }
}

unsigned int ParadisCrossSlipServer::IdentifyFCCSlipSystem(
    Vector oNormalVector,
    Vector oBurgersVector) // defined as static in the .h file
{
  oNormalVector.Normalize();
  oBurgersVector.Normalize();
  double dTolerance = 1E-3;
  Vector oReferenceNormal(0.0, 0.0, 0.0);
  Vector oReferenceBurgers(0.0, 0.0, 0.0);

  // slip plane 1,1,1
  oReferenceNormal.Set(1.0, 1.0, 1.0); // check how it is done
  oReferenceNormal.Normalize();        // remember to normalize

  oReferenceBurgers.Set(1.0, -1.0, 0.0);
  oReferenceBurgers.Normalize();
  if (oReferenceNormal.IsSame(oNormalVector, dTolerance) ||
      oReferenceNormal.IsOpposite(oNormalVector, dTolerance)) {
    if (oReferenceBurgers.IsSame(oBurgersVector, dTolerance) ||
        oReferenceBurgers.IsOpposite(
            oBurgersVector,
            dTolerance)) // to check on tolerance and on being opoosite
    {
      return 1;
    }
  }

  oReferenceBurgers.Set(1.0, 0.0, -1.0);
  oReferenceBurgers.Normalize();
  if (oReferenceNormal.IsSame(oNormalVector, dTolerance) ||
      oReferenceNormal.IsOpposite(oNormalVector, dTolerance)) {
    if (oReferenceBurgers.IsSame(oBurgersVector, dTolerance) ||
        oReferenceBurgers.IsOpposite(oBurgersVector, dTolerance)) {
      return 2;
    }
  }

  oReferenceBurgers.Set(0.0, 1.0, -1.0);
  oReferenceBurgers.Normalize();
  if (oReferenceNormal.IsSame(oNormalVector, dTolerance) ||
      oReferenceNormal.IsOpposite(oNormalVector, dTolerance)) {
    if (oReferenceBurgers.IsSame(oBurgersVector, dTolerance) ||
        oReferenceBurgers.IsOpposite(oBurgersVector, dTolerance)) {
      return 3;
    }
  }

  // slip plane -1,1,1
  oReferenceNormal.Set(-1.0, 1.0, 1.0);
  oReferenceNormal.Normalize();

  oReferenceBurgers.Set(1.0, 1.0, 0.0);
  oReferenceBurgers.Normalize();
  if (oReferenceNormal.IsSame(oNormalVector, dTolerance) ||
      oReferenceNormal.IsOpposite(oNormalVector, dTolerance)) {
    if (oReferenceBurgers.IsSame(oBurgersVector, dTolerance) ||
        oReferenceBurgers.IsOpposite(oBurgersVector, dTolerance)) {
      return 4;
    }
  }

  oReferenceBurgers.Set(1.0, 0.0, 1.0);
  oReferenceBurgers.Normalize();
  if (oReferenceNormal.IsSame(oNormalVector, dTolerance) ||
      oReferenceNormal.IsOpposite(oNormalVector, dTolerance)) {
    if (oReferenceBurgers.IsSame(oBurgersVector, dTolerance) ||
        oReferenceBurgers.IsOpposite(oBurgersVector, dTolerance)) {
      return 5;
    }
  }

  oReferenceBurgers.Set(0.0, 1.0, -1.0);
  oReferenceBurgers.Normalize();
  if (oReferenceNormal.IsSame(oNormalVector, dTolerance) ||
      oReferenceNormal.IsOpposite(oNormalVector, dTolerance)) {
    if (oReferenceBurgers.IsSame(oBurgersVector, dTolerance) ||
        oReferenceBurgers.IsOpposite(oBurgersVector, dTolerance)) {
      return 6;
    }
  }

  // slip plane 1,-1,1
  oReferenceNormal.Set(1.0, -1.0, 1.0);
  oReferenceNormal.Normalize();

  oReferenceBurgers.Set(1.0, 1.0, 0.0);
  oReferenceBurgers.Normalize();
  if (oReferenceNormal.IsSame(oNormalVector, dTolerance) ||
      oReferenceNormal.IsOpposite(oNormalVector, dTolerance)) {
    if (oReferenceBurgers.IsSame(oBurgersVector, dTolerance) ||
        oReferenceBurgers.IsOpposite(oBurgersVector, dTolerance)) {
      return 7;
    }
  }

  oReferenceBurgers.Set(1.0, 0.0, -1.0);
  oReferenceBurgers.Normalize();
  if (oReferenceNormal.IsSame(oNormalVector, dTolerance) ||
      oReferenceNormal.IsOpposite(oNormalVector, dTolerance)) {
    if (oReferenceBurgers.IsSame(oBurgersVector, dTolerance) ||
        oReferenceBurgers.IsOpposite(oBurgersVector, dTolerance)) {
      return 8;
    }
  }

  oReferenceBurgers.Set(0.0, 1.0, 1.0);
  oReferenceBurgers.Normalize();
  if (oReferenceNormal.IsSame(oNormalVector, dTolerance) ||
      oReferenceNormal.IsOpposite(oNormalVector, dTolerance)) {
    if (oReferenceBurgers.IsSame(oBurgersVector, dTolerance) ||
        oReferenceBurgers.IsOpposite(oBurgersVector, dTolerance)) {
      return 9;
    }
  }

  // slip plane 1,1,-1
  oReferenceNormal.Set(1.0, 1.0, -1.0);
  oReferenceNormal.Normalize();

  oReferenceBurgers.Set(1.0, -1.0, 0.0);
  oReferenceBurgers.Normalize();
  if (oReferenceNormal.IsSame(oNormalVector, dTolerance) ||
      oReferenceNormal.IsOpposite(oNormalVector, dTolerance)) {
    if (oReferenceBurgers.IsSame(oBurgersVector, dTolerance) ||
        oReferenceBurgers.IsOpposite(oBurgersVector, dTolerance)) {
      return 10;
    }
  }

  oReferenceBurgers.Set(1.0, 0.0, 1.0);
  oReferenceBurgers.Normalize();
  if (oReferenceNormal.IsSame(oNormalVector, dTolerance) ||
      oReferenceNormal.IsOpposite(oNormalVector, dTolerance)) {
    if (oReferenceBurgers.IsSame(oBurgersVector, dTolerance) ||
        oReferenceBurgers.IsOpposite(oBurgersVector, dTolerance)) {
      return 11;
    }
  }

  oReferenceBurgers.Set(0.0, 1.0, 1.0);
  oReferenceBurgers.Normalize();
  if (oReferenceNormal.IsSame(oNormalVector, dTolerance) ||
      oReferenceNormal.IsOpposite(oNormalVector, dTolerance)) {
    if (oReferenceBurgers.IsSame(oBurgersVector, dTolerance) ||
        oReferenceBurgers.IsOpposite(oBurgersVector, dTolerance)) {
      return 12;
    }
  }

  return 13;
}
void ParadisCrossSlipServer::DisposeChains(
    list<DislocationChain *> *plpoChains) {
  list<DislocationChain *>::iterator liChains;
  for (liChains = plpoChains->begin(); liChains != plpoChains->end();
       liChains++) {
    (*liChains)->Collapse();
    delete (*liChains);
  }
  plpoChains->clear();
}
bool ParadisCrossSlipServer::Is111Vector(Vector oVector) {
  oVector.Normalize();
  double dX = oVector.GetX();
  double dY = oVector.GetY();
  double dZ = oVector.GetZ();
  double dTolerance = 1.0E-6;
  if (fabs(fabs(dX) - fabs(dY)) < dTolerance) {
    if (fabs(fabs(dY) - fabs(dZ)) < dTolerance) {
      return true;
    }
  }
  return false;
}
bool ParadisCrossSlipServer::Is110Vector(Vector oVector) {
  oVector.Normalize();
  double dX = oVector.GetX();
  double dY = oVector.GetY();
  double dZ = oVector.GetZ();
  double dTolerance = 1.0E-6;
  if (fabs(fabs(dX) - fabs(dY)) < dTolerance) {
    if (fabs(dZ) < dTolerance) {
      return true;
    }
  }

  if (fabs(fabs(dY) - fabs(dZ)) < dTolerance) {
    if (fabs(dX) < dTolerance) {
      return true;
    }
  }

  if (fabs(fabs(dZ) - fabs(dX)) < dTolerance) {
    if (fabs(dY) < dTolerance) {
      return true;
    }
  }
  return false;
}
bool ParadisCrossSlipServer::Is112Vector(Vector oVector) {
  oVector.Normalize();
  double dX = oVector.GetX();
  double dY = oVector.GetY();
  double dZ = oVector.GetZ();
  double dTolerance = 1.0E-6;
  double dOne = 1.0 / sqrt(6.0);
  double dTwo = 2.0 * dOne;
  if (fabs(fabs(dX) - dTwo) < dTolerance) {
    if (fabs(fabs(dY) - dOne) < dTolerance) {
      if (fabs(fabs(dZ) - dOne) < dTolerance) {
        return true;
      }
    }
  }

  if (fabs(fabs(dY) - dTwo) < dTolerance) {
    if (fabs(fabs(dZ) - dOne) < dTolerance) {
      if (fabs(fabs(dX) - dOne) < dTolerance) {
        return true;
      }
    }
  }

  if (fabs(fabs(dZ) - dTwo) < dTolerance) {
    if (fabs(fabs(dX) - dOne) < dTolerance) {
      if (fabs(fabs(dY) - dOne) < dTolerance) {
        return true;
      }
    }
  }
  return false;
}
bool ParadisCrossSlipServer::Is100Vector(Vector oVector) {
  oVector.Normalize();
  double dX = oVector.GetX();
  double dY = oVector.GetY();
  double dZ = oVector.GetZ();
  double dTolerance = 1.0E-6;
  if (fabs(fabs(dX) - 1.0) < dTolerance) {
    if (fabs(dY) < dTolerance) {
      if (fabs(dZ) < dTolerance) {
        return true;
      }
    }
  }

  if (fabs(fabs(dY) - 1.0) < dTolerance) {
    if (fabs(dZ) < dTolerance) {
      if (fabs(dX) < dTolerance) {
        return true;
      }
    }
  }

  if (fabs(fabs(dZ) - 1.0) < dTolerance) {
    if (fabs(dX) < dTolerance) {
      if (fabs(dY) < dTolerance) {
        return true;
      }
    }
  }
  return false;
}
void ParadisCrossSlipServer::GetNodeDynamicConstraint(
    DislocationNetworkNode *poNode, unsigned int &iConstraintType,
    Vector &oConstraintVector) {
  oConstraintVector.Set(0.0, 0.0, 0.0);
  iConstraintType = 0;
  double dTolerance = 1.0E-6;
  list<DislocationNetworkArm *> *plpoEdges = poNode->GetEdges();
  list<DislocationNetworkArm *>::iterator liEdges;
  unsigned int iNeighboursCount = plpoEdges->size();
  // get the node dynamic constraint
  if (iNeighboursCount == 0) {
    oConstraintVector.Set(0.0, 0.0, 0.0);
    iConstraintType = 0;
  } else if (iNeighboursCount == 1) {
    iConstraintType = 1;
    liEdges = plpoEdges->begin();
    oConstraintVector = (*liEdges)->GetDataPointer()->GetSlipPlaneNormal();
  } else if (iNeighboursCount == 2) {
    liEdges = plpoEdges->begin();
    Vector oNormal1 = (*liEdges)->GetDataPointer()->GetSlipPlaneNormal();
    liEdges++;
    Vector oNormal2 = (*liEdges)->GetDataPointer()->GetSlipPlaneNormal();
    oNormal1.Normalize();
    oNormal2.Normalize();
    Vector oCrossProduct = oNormal1 ^ oNormal2;
    if (oCrossProduct.Length() < dTolerance) {
      iConstraintType = 1;
      oConstraintVector = oNormal1;
    } else {
      iConstraintType = 2;
      oConstraintVector = oCrossProduct;
    }
  } else {
    list<Vector> loNormals;
    list<Vector>::iterator liNormals;
    Vector oTempNormal;
    Vector oCrossProduct;
    bool bAddNormal = false;
    for (liEdges = plpoEdges->begin(); liEdges != plpoEdges->end(); liEdges++) {
      oTempNormal = (*liEdges)->GetDataPointer()->GetSlipPlaneNormal();
      oTempNormal.Normalize();
      bAddNormal = false;
      if (loNormals.empty()) {
        loNormals.push_back(oTempNormal);
      } else {
        for (liNormals = loNormals.begin(); liNormals != loNormals.end();
             liNormals++) // junjie: ?you might adding two same normals
        {
          oCrossProduct = oTempNormal ^ (*liNormals);
          if (oCrossProduct.Length() > dTolerance) {
            bAddNormal = true;
            break;
          }
        }
        if (bAddNormal) {
          loNormals.push_back(oTempNormal);
        }
      }
    }

    unsigned int iDynamicConstraintsCount = (unsigned int)loNormals.size();
    if (iDynamicConstraintsCount == 1) {
      iConstraintType = 1;
      oConstraintVector = loNormals.front();
    } else if (iDynamicConstraintsCount == 2) {
      iConstraintType = 2;
      oConstraintVector = loNormals.front() ^ loNormals.back();
    } else {
      iConstraintType = 3;
      oConstraintVector.Set(0.0, 0.0, 0.0);
    }
  }
  oConstraintVector.Normalize();
}
bool ParadisCrossSlipServer::CanChainCrossSlip(Home_t *poParadisHome,
                                               DislocationChain *poChain) {
  DislocationNode *poCurrentDislocationNode = NULL;
  Node_t *poNode = NULL;
  unsigned int i = 0;
  unsigned int iChainLength = poChain->GetSize();
  poChain->ResetIterator();
  for (i = 0; i < iChainLength; i++) {
    poCurrentDislocationNode = poChain->GetCurrentNode()->GetDataPointer();
    poNode = GetNodeFromIndex(poParadisHome, poParadisHome->myDomain,
                              poCurrentDislocationNode->GetID());
    if (poNode == NULL) {
      continue;
    }
    if (poNode->flags & NODE_NO_CROSS_SLIP) {
      return false;
    }
    poChain->IncrementIterator();
  }
  return true;
}
