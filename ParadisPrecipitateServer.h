#ifndef PARADISPRECIPITATESERVER_H_
#define PARADISPRECIPITATESERVER_H_

#include "Typedefs.h"
#include "PrecipitateStructure.h"
#include "APBEventCommunicator.h"

using namespace std;
using namespace EZ;

class APBPoint {
public:
  APBPoint();
  ~APBPoint();
  APBPoint(const APBPoint &oPoint);
  APBPoint &operator=(const APBPoint &oPoint);
  bool Shear(Vector oBurgersVector);
  unsigned int GetOwningDomain() const;
  void SetOwningDomain(const unsigned int &iDomain);
  bool IsSheared() const;
  void SetShearable();
  double m_dX;
  double m_dY;
  unsigned int m_iX;
  unsigned int m_iY;
  unsigned int m_iSingletonCycle;
  void Write(FILE *fpFile);
  void SetShear(const int &iXShear, const int &iYShear, const int &iZShear);
  int GetXShear() const;
  int GetYShear() const;
  int GetZShear() const;

private:
protected:
  void Initialize();
  void CheckSuperDislocated();
  unsigned int m_iOwningDomain;
  bool m_bIsShearable;
  int m_iXShear;
  int m_iYShear;
  int m_iZShear;
};

class APB {
public:
  APB();
  ~APB();
  APB(const APB &oAPB);
  APB &operator=(const APB &oAPB);
  void Reset();
  void SetSlipPlaneNormal(const Vector &oVector);
  void SetResolution(const unsigned int &iResolution);
  void SetXMin(const double &dValue);
  void SetXMax(const double &dValue);
  void SetYMin(const double &dValue);
  void SetYMax(const double &dValue);
  void SetPlaneSpacing(const double &dValue);
  void GeneratePoints(PrecipitateStructure *poPrecipitateStructure);
  Vector GetSlipPlaneNormal();
  double GetPlaneSpacing() const;
  unsigned int GetResolution() const;
  double GetXMin() const;
  double GetXMax() const;
  double GetYMin() const;
  double GetYMax() const;
  void Update(Node_t *poNode1, Node_t *poNode2,
              APBEventCommunicator *poAPBEventCommunicator);
  void Write(FILE *fpFile);
  void Read(FILE *fpFile, const unsigned int &iResolution,
            PrecipitateStructure *poPrecipitateStructure,
            unsigned int &iSlipPlaneID);
  void UpdateTriangle(const Point &oPoint1, const Point &oPoint2,
                      const Point &oPoint3, const Vector &oBurgersVector,
                      APBEventCommunicator *poAPBEventCommunicator);
  void UpdateAPBCellOwnership(Home_t *poHome);
  void ShearCell(const unsigned int &iCellX, const unsigned int &iCellY,
                 const int &iXShear, const int &iYShear, const int &iZShear);
  void Report(const unsigned int &iProcessorID) const;
  APBPoint *GetContainingCell(const Point &oPoint);
  Point GetAPBPoint(APBPoint *poPoint);
  unsigned int GetSlipPlaneID();
  static unsigned int GetFCCSlipPlaneID(const Vector &oSlipPlaneNormal);
  static Vector GetFCCSlipPlaneNormalFromID(const unsigned int &iSlipPlaneID);
  void CleanSingletons(const unsigned int &iCycle,
                       const unsigned int &iProcessorID,
                       APBEventCommunicator *poAPBEventCommunicator);

private:
protected:
  void Initialize();
  bool IsAPBPointInSweepQuad(APBPoint *poPoint, Node_t *poNode1,
                             Node_t *poNode2);
  bool IsAPBPointInSweepTri(APBPoint *poAPBPoint, const Point &oPoint1,
                            const Point &oPoint2, const Point &oPoint3);
  void SetFCCSlipPlaneNormalFromID(const unsigned int &iSlipPlaneID);
  bool IsSingleton(const unsigned int &iX, const unsigned int &iY);
  void ClearPoints();
  double m_dXMin;
  double m_dXMax;
  double m_dYMin;
  double m_dYMax;
  double m_dXSpacing;
  double m_dYSpacing;
  unsigned int m_iResolution;
  APBPoint **m_ppoPoints;
  Vector m_oSlipPlaneNormal;
  double m_dPlaneSpacing;
};

class ParadisPrecipitateServer {
public:
  static ParadisPrecipitateServer *CreateInstance();
  ~ParadisPrecipitateServer();
  void Reset();
  bool Set(Home_t *poHome);
  void CheckNodes(Home_t *poHome);
  void ComputeAPBForces(Home_t *poHome, Node_t *poNode1, Node_t *poNode2,
                        double pdForce1[3], double pdForce2[3]);
  void AddAPBShearing(Home_t *poHome);
  void UpdateAPB(Home_t *poHome);
  static unsigned int IndentifyFCCSlipSystem(Vector oNormalVector,
                                             Vector oBurgersVector);
  static Vector GetFCCSlipPlaneNormal(const unsigned int &iPlaneID);
  static void GetFCCSlipSystemVectors(const unsigned int &iSlipSystem,
                                      Vector &oNormalVector,
                                      Vector &oBurgersVector);
  void WriteAPB(const string &sFileName) const;
  bool IsNodeInsideAPrecipitate(Node_t *poNode,
                                Polyhedron *&poContainingPrecipitate) const;
  bool IsPointInsideAPrecipitate(Point *poPoint,
                                 Polyhedron *&poContainingPrecipitate) const;
  void UpdateTriangleSweep(Home_t *poHome, const Point &oPoint1,
                           const Point &oPoint2, const Point &oPoint3,
                           const Vector &oBurgersVector);
  void ReadAPBs(const string &sAPBFileName);
  void UpdateAPBCellOwnership(Home_t *poHome);
  APB *GetAPBByPlane(const unsigned int &iPlaneID,
                     const double &dPlaneSpacing) const;
  APB *GetAPB(const unsigned int &iSlipSystemID,
              const double &dPlaneSpacing) const;
  APB *GetAPB(Node_t *poNode1, Node_t *poNode2) const;
  APB *GetAPB(Node_t *poNode, const Vector &oSlipPlaneNormal) const;
  APB *AddAPBByPlane(const unsigned int &iPlaneID, const double &dPlaneSpacing);
  APB *AddAPB(const unsigned int &iSlipSystemID, const double &dPlaneSpacing);
  void Report(const unsigned int &iProcessorID) const;
  bool IsSegmentShearCapable(Node_t *poNode1, Node_t *poNode2,
                             Polyhedron *&poPrecipitate) const;
  unsigned int GetPrecipitatesCount() const;
  unsigned int GetTotalNucleatedCount() const;
  bool WritePrecipitateStructureSnapshot(const string &sFileName) const;

private:
protected:
  static ParadisPrecipitateServer *m_poParadisPrecipitateServerInstance;
  static double m_dCellSize;
  ParadisPrecipitateServer();
  void Initialize();
  static bool IsSuperDislocation(Node_t *poNode,
                                 const unsigned int &iNeighbourIndex);
  void RetractNode(Home_t *poHome, Node_t *poNode,
                   Polyhedron *poContainingPrecipitate = NULL) const;
  bool IsNodeShearCapable(Node_t *poNode, Polyhedron *&poPrecipitate) const;
  static unsigned int GetFCCSlipPlaneID(const unsigned int &iSlipSystemID);
  void NucleatePrecipitates(Home_t *poHome);
  unsigned int GetNucleationCount(const double &dExpectedCount) const;
  Point GetRandomNucleationCenter(const Home_t *poHome) const;

  PrecipitateStructure m_oStructure;
  vector<list<APB *>> m_vlpAPB;
  APBEventCommunicator *m_poAPBEventCommunicator;
  unsigned int m_iDomainID;
  double m_dXMin;
  double m_dXMax;
  double m_dYMin;
  double m_dYMax;
  static unsigned int m_iAPBPointResolution;
  bool m_bShearablePrecipitates;
  bool m_bPairShearing;
  bool m_bUseNucleationMode;
  double m_dNucleationRate;
  double m_dNucleationRadius;
  int m_iNucleationMaxCount;
  unsigned int m_iTotalNucleatedCount;
};

#endif
