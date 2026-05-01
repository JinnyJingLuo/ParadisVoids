#ifndef SURFACESEGMENT_H_
#define SURFACESEGMENT_H_

#include "DislocationNode.h"
#include "DislocationSegment.h"
#include "vector"

using namespace std;
using namespace DislocationSystem;

class SurfaceSegment {
public:
  SurfaceSegment();
  SurfaceSegment(const SurfaceSegment &oSegment);
  ~SurfaceSegment();
  SurfaceSegment &operator=(const SurfaceSegment &oSegment);
  void Reset();
  static void SetVirtualNodeSpacing(const double &dSpacing);
  static void SetPoissonsRatio(const double &dRatio);
  static void SetGaussPointsCount(const unsigned int &iCount);
  void Set(const DislocationNode &oNode1, const DislocationNode &oNode2,
           const DislocationSegment &oSegment);
  Vector GetPointDisplacement(const Point &oPoint) const;
  void Print() const;
  bool IsPostMergeable(SurfaceSegment *poSegment) const;
  void PostMerge(SurfaceSegment *poSegment);
  void Write(FILE *fpFile) const;
  double GetLength() const;
  const Vector *GetBurgersVector() const;
  const Vector *GetSlipPlaneNormal() const;
  const DislocationNode *GetFirstNode() const;
  const DislocationNode *GetSecondNode() const;

private:
protected:
  void Initialize();
  static Vector GetSegmentDisplacement(const Point &oPoint,
                                       const Point &oEndPoint1,
                                       const Point &oEndPoint2,
                                       const Vector &oBurgersVector);
  bool IsPostTouching(SurfaceSegment *poSegment,
                      const double &dTolerance = 1.0E-2) const;

  static unsigned int GaussPointsCount;
  static double VirtualNodeSpacing;
  static double PoissonsRatio;
  static vector<double> GaussPointsLocations;
  static vector<double> GaussPointsWeights;

  Point m_oOrigin;
  DislocationNode m_oFirstOriginalNode;
  DislocationNode m_oSecondOriginalNode;
  Point m_oFirstVirtualNode;
  Point m_oSecondVirtualNode;
  Vector m_oBurgersVector;
  Vector m_oSlipPlaneNormal;
};

#endif
