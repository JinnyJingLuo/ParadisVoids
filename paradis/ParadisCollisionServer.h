// Ahmed M. Hussein

#ifndef PARADISCOLLISIONSERVER_H_
#define PARADISCOLLISIONSERVER_H_

#include "Typedefs.h"
#include "Matrix.h"
#include "Vector.h"
#include "ParadisSurface.h"
#include "stdio.h"

using namespace std;
using namespace EZ;
using namespace GeometrySystem;

class _home;
class _node;

class AnnihilationSegment {
public:
  AnnihilationSegment();
  AnnihilationSegment(const AnnihilationSegment &oSegment);
  ~AnnihilationSegment();
  AnnihilationSegment &operator=(const AnnihilationSegment &oSegment);
  void Reset();
  void SetPoints(const Point &oStart, const Point &oEnd);
  void SetBurgersVector(const Vector &oBurgersVector);
  void SetTimeStep(const unsigned int &iTimeStep);
  void Write(FILE *fpFile) const;

private:
protected:
  void Initialize();
  Point m_oStartPoint;
  Point m_oEndPoint;
  Vector m_oBurgersVector;
  unsigned int m_iTimeStep;
};

class ParadisCollisionServer {
public:
  static ParadisCollisionServer *CreateInstance();
  ~ParadisCollisionServer();
  void Reset();
  void Set(Home_t *poHome);
  void HandleCollisions(Home_t *poHome);

private:
protected:
  static ParadisCollisionServer *m_poParadisCollisionServerInstance;
  ParadisCollisionServer();
  void Initialize();
  bool IsNodeCollisionValid(Node_t *poNode) const;
  bool AreSegmentsCloseAndApproaching(Home_t *poHome, Node_t *poNode1,
                                      Node_t *poNode2, Node_t *poNode3,
                                      Node_t *poNode4, Point &oNearPoint1,
                                      Point &oNearPoint2) const;
  bool AreNodesCloseAndApproaching(Home_t *poHome, Node_t *poNode1,
                                   Node_t *poNode2) const;
  bool IsTriangleCollapsing(Home_t *poHome, Node_t *poNode1, Node_t *poNode2,
                            Node_t *poNode3) const;
  bool GetPlanePlaneCollisionPoint(Point *poNode1, const Vector &oNormal1,
                                   Point *poNode2, const Vector &oNormal2,
                                   Point &oCollisionPoint) const;
  bool GetLinePlaneCollisionPoint(Point *poNode1, const Vector &oDirection,
                                  Point *poNode2, const Vector &oNormal,
                                  Point &oCollisionPoint) const;
  bool GetLineLineCollisionPoint(Point *poNode1, const Vector &oDirection1,
                                 Point *poNode2, const Vector &oDirection2,
                                 Point &oCollisionPoint) const;
  bool DoesExactCollisionPointExist(Point *poNode1,
                                    const unsigned int &iConstraintType1,
                                    const Vector &oConstraintVector1,
                                    Point *poNode2,
                                    const unsigned int &iConstraintType2,
                                    const Vector &oConstraintVector2,
                                    Point &oExactCollisionPoint) const;
  bool GetExactCollisionPoint(Home_t *poHome, Node_t *poNode1, Node_t *poNode2,
                              Point &oCollisionPoint) const;
  void HandleAnnihilations(Home_t *poHome);
  void HandleSegmentSegmentCollisions(Home_t *poHome) const;
  void HandleHingeJointCollisions(Home_t *poHome) const;
  void HandleTriangularLoops(Home_t *poHome) const;
  void WriteAnnihilatedSegments();
  list<AnnihilationSegment *> m_lpoAnnihilatedSegments;
  string m_sAnnihilatedSegmentsFileName;
  FILE *m_fpAnnihilatedSegments;
};

#endif
