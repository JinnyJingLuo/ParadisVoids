// Ahmed M. Hussein

#ifndef PARADISCROSSSLIPSERVER_H_
#define PARADISCROSSSLIPSERVER_H_

#include "Home.h"
#include "DislocationChain.h"
#include "Line.h"
#include "Graph.h"
#include "list"

using namespace std;
using namespace DislocationSystem;
using namespace GeometrySystem;
using namespace GraphSystem;

class ParadisCrossSlipServer {
public:
  static ParadisCrossSlipServer *GetInstance();
  static void Free();
  ~ParadisCrossSlipServer();
  void Set(Home_t *poParadisHome);
  void Reset();
  void HandleCrossSlip(Home_t *poHome);
  static bool Is111Vector(Vector oVector);
  static bool Is110Vector(Vector oVector);
  static bool Is112Vector(Vector oVector);
  static bool Is100Vector(Vector oVector);
  static unsigned int IdentifyFCCSlipSystem(Vector oNormalVector,
                                            Vector oBurgersVector);

private:
protected:
  static ParadisCrossSlipServer *m_poServer;
  ParadisCrossSlipServer();
  void Initialize();
  void GetChainCrossSlipValues(Home_t *poParadisHome, DislocationChain *poChain,
                               double &dAverageGlidePlaneShearStress,
                               double &dAverageCrossSlipPlaneShearStress,
                               double &dGlideEscaigStress,
                               double &dCrossSlipEscaigStress,
                               double &dTotalLength,
                               bool bExcludeLocalEscaigStress);
  double GetEscaigStress(const Vector &oForce, const Vector &oBurgersDirection,
                         const Vector &oNormal,
                         const double &dBurgersMagnitude);
  list<DislocationChain *> GenerateDislocationChains(Home_t *poHome);
  void SeparateScrewChains(Param_t *poParam,
                           list<DislocationChain *> *plpoChains,
                           const double &dScrewAngleTolerance,
                           bool bIntersectionOnly = false);
  void RefineScrewChains(list<DislocationChain *> *plpoChains,
                         bool bIntersectionOnly = false);
  void HandleChainCrossSlip(Home_t *poParadisHome,
                            list<DislocationChain *> *plpoChains,
                            bool bIntersectionOnly = false);
  void HandleBulkCrossSlip(Home_t *poParadisHome, DislocationChain *poChain);
  void HandleSurfaceCrossSlip(Home_t *poParadisHome, DislocationChain *poChain);
  void HandleIntersectionCrossSlip(Home_t *poParadisHome,
                                   DislocationChain *poChain,
                                   list<DislocationChain *> *plpoAllChains);
  bool Check23IntersectionChainCrossSlip(
      Home_t *poParadisHome, DislocationChain *poChain,
      list<DislocationChain *> *plpoAllChains,
      DislocationChain *&poIntersectionChain, Vector &oJunctionBurgers,
      Vector &oJunctionSlipPlane);
  bool Check24IntersectionChainCrossSlip(
      Home_t *poParadisHome, DislocationChain *poChain,
      list<DislocationChain *> *plpoAllChains,
      DislocationChain *&poIntersectionChain, Vector &oJunctionBurgers,
      Vector &oJunctionSlipPlane);
  void HandleRepulsiveCrossSlip(Home_t *poParadisHome,
                                DislocationChain *poChain,
                                DislocationChain *poIntersectionChain);
  void HandleAttractiveCrossSlip(Home_t *poParadisHome,
                                 DislocationChain *poChain,
                                 DislocationChain *poIntersectionChain,
                                 Vector &oJunctionBurgers,
                                 Vector &oJunctionSlipPlane);
  Vector GetFCCCrossSlipPlaneNormal(const Vector &oSlipPlaneNormal,
                                    const Vector &oBurgersVector);
  bool AdjustChainForCrossSlip(Home_t *poParadisHome,
                               DislocationChain *poChain);
  void CrossSlipChain(Home_t *poParadisHome, DislocationChain *poChain,
                      const Vector &oCrossSlipPlaneNormal);
  void CrossSlipChain_ThermalGradient(Home_t *poParadisHome,
                                      DislocationChain *poChain,
                                      const Vector &oCrossSlipPlaneNormal,
                                      double preProb, double Exponent_noTemp,
                                      double randProb);
  bool GetCrossSlipFitLine(DislocationChain *poChain, Line &oFitLine);
  bool CanChainCrossSlip(Home_t *poParadisHome, DislocationChain *poChain);
  void DisposeChains(list<DislocationChain *> *plpoChains);
  void GetNodeDynamicConstraint(DislocationNetworkNode *poNode,
                                unsigned int &iConstraintType,
                                Vector &oConstraintVector);
  Graph<DislocationNode, DislocationSegment> m_oNetwork;
};

#endif
