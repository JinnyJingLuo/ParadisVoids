// Ahmed M. Hussein

#ifndef PARADISSURFACE_H_
#define PARADISSURFACE_H_

#include "Typedefs.h"
#include "Vector.h"
#include "TriPatch.h"
#include "list"
#include "FEMElement.h"
#include "SurfaceSegment.h"

using namespace std;
using namespace EZ;
using namespace GeometrySystem;
using namespace FEMSystem;

class _home;

class ParadisSurface {
public:
  ParadisSurface();
  ParadisSurface(const ParadisSurface &oSurface);
  ~ParadisSurface();
  ParadisSurface &operator=(const ParadisSurface &oSurface);
  void Reset();
  virtual void Set(Home_t *poHome) = 0;
  virtual void CheckNodes(Home_t *poHome);
  void StoreSurfaceArms(Home_t *poHome);
  void StoreSurfaceNodesMotion(Home_t *poHome);
  void WriteSurfaceSegments(Home_t *poHome);
  static void GetNodeDynamicConstraint(Node_t *poNode,
                                       unsigned int &iConstraintType,
                                       Vector &oConstraintVector);

private:
protected:
  void Initialize();
  TriPatch *GetNearestTriangle(const Point &oPoint, Point &oNearestPoint) const;
  TriPatch *GetNearestTriangleOnPlane(const Plane &oPlane, const Point &oPoint,
                                      Point &oNearestPoint) const;
  TriPatch *GetNearestTriangleOnLine(const Line &oLine, const Point &oPoint,
                                     Point &oNearestPoint) const;
  virtual bool IsPointInside(const Point &oPoint) const;
  virtual bool IsPointOnSurface(const Point &oPoint) const;
  virtual double GetLeastSurfaceDistance(const Point &oPoint) const;
  virtual void GenerateTriangulations() = 0;
  void ClearTriangulations();
  void GenerateSurfaceTriangulationFromFEMNodes(
      const vector<FEMNode *> &vpoNodes,
      const vector<FEMElement *> &vpoElement);
  void ReduceSurfaceSemgnts();
  void HandleRigidBoundary(Home_t *poHome, Node_t *poNode);

  list<GenericNode *> m_lpoPoints;
  list<TriPatch *> m_lpoTriangles;
  list<SurfaceSegment *> m_lpoSurfaceSegments;
  double m_dXMin;
  double m_dXMax;
  double m_dYMin;
  double m_dYMax;
  double m_dZMin;
  double m_dZMax;
  double m_dTolerance;
  double m_dVolume;
};

#endif
