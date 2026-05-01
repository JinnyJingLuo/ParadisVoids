// Paradis Processor Project
// Ahmed M. Hussein (mailto : am.hussin@gmail.com)
// June 2012

#ifndef DISLOCATION_H_
#define DISLOCATION_H_

#include "DislocationSegment.h"
#include "DislocationNode.h"
#include "Matrix.h"
#include "Graph.h"

using namespace GraphSystem;

namespace DislocationSystem {
class Dislocation {
public:
  static Vector GetPointDisplacement(
      GraphEdge<DislocationNode, DislocationSegment> *poDislocationSegment,
      Point *poPoint, bool bIsDebug = false);
  static Matrix GetSegmentPointStress(const Point &oSegmentStart,
                                      const Point &oSegmentEnd,
                                      const Vector &oBurgersVector,
                                      const double &dShearModulus,
                                      const double &dPoissonRatio,
                                      const Point &oFieldPoint);
  static Matrix GetSegmentPointStrain(const Point &oSegmentStart,
                                      const Point &oSegmentEnd,
                                      const Vector &oBurgersVector,
                                      const double &dPoissonRatio,
                                      const Point &oFieldPoint);

private:
protected:
};
} // namespace DislocationSystem

#endif
