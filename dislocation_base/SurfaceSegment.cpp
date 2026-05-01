#include "SurfaceSegment.h"
#include "MathServices.h"

unsigned int SurfaceSegment::GaussPointsCount = 0;
double SurfaceSegment::VirtualNodeSpacing = 0.0;
double SurfaceSegment::PoissonsRatio = 0.0;
vector<double> SurfaceSegment::GaussPointsLocations;
vector<double> SurfaceSegment::GaussPointsWeights;

SurfaceSegment::SurfaceSegment() { Initialize(); }
SurfaceSegment::SurfaceSegment(const SurfaceSegment &oSegment) {
  *this = oSegment;
}
SurfaceSegment::~SurfaceSegment() { Reset(); }
SurfaceSegment &SurfaceSegment::operator=(const SurfaceSegment &oSegment) {
  m_oOrigin = oSegment.m_oOrigin;
  m_oFirstOriginalNode = oSegment.m_oFirstOriginalNode;
  m_oSecondOriginalNode = oSegment.m_oSecondOriginalNode;
  m_oFirstVirtualNode = oSegment.m_oFirstVirtualNode;
  m_oSecondVirtualNode = oSegment.m_oSecondVirtualNode;
  m_oBurgersVector = oSegment.m_oBurgersVector;
  m_oSlipPlaneNormal = oSegment.m_oSlipPlaneNormal;
  return *this;
}
void SurfaceSegment::Reset() {}
void SurfaceSegment::SetVirtualNodeSpacing(const double &dSpacing) {
  VirtualNodeSpacing = dSpacing;
}
void SurfaceSegment::SetPoissonsRatio(const double &dRatio) {
  PoissonsRatio = dRatio;
}
void SurfaceSegment::SetGaussPointsCount(const unsigned int &iCount) {
  GaussPointsLocations.clear();
  GaussPointsWeights.clear();
  GaussPointsCount = iCount;
  MathServices::GenerateGaussPoints(GaussPointsLocations, GaussPointsWeights,
                                    GaussPointsCount);
}
void SurfaceSegment::Set(const DislocationNode &oNode1,
                         const DislocationNode &oNode2,
                         const DislocationSegment &oSegment) {
  // the origin is always zero
  m_oOrigin.Set(0.0, 0.0, 0.0);
  // the original nodes are the dislocation nodes
  m_oFirstOriginalNode = oNode1;
  m_oSecondOriginalNode = oNode2;

  // the virtual nodes must be projected in the surface normal direction but
  // they have to stay on their slip plane
  Vector oSlipNormal = oSegment.GetSlipPlaneNormal();
  m_oBurgersVector = oSegment.GetBurgersVector();
  m_oSlipPlaneNormal = oSegment.GetSlipPlaneNormal();

  Vector oDisplacement = oNode1.GetSurfaceNormal() * VirtualNodeSpacing;
  double dDropHeight = oDisplacement * oSlipNormal;
  oDisplacement = oDisplacement - oSlipNormal * dDropHeight;
  m_oFirstVirtualNode = m_oFirstOriginalNode + oDisplacement;

  oDisplacement = oNode2.GetSurfaceNormal() * VirtualNodeSpacing;
  dDropHeight = oDisplacement * oSlipNormal;
  oDisplacement = oDisplacement - oSlipNormal * dDropHeight;
  m_oSecondVirtualNode = m_oSecondOriginalNode + oDisplacement;
}
Vector SurfaceSegment::GetPointDisplacement(const Point &oPoint) const {
  Vector oDisplacement = GetSegmentDisplacement(
      oPoint, m_oOrigin, m_oFirstOriginalNode, m_oBurgersVector);
  oDisplacement = oDisplacement +
                  GetSegmentDisplacement(oPoint, m_oFirstOriginalNode,
                                         m_oFirstVirtualNode, m_oBurgersVector);
  oDisplacement = oDisplacement + GetSegmentDisplacement(
                                      oPoint, m_oFirstVirtualNode,
                                      m_oSecondVirtualNode, m_oBurgersVector);
  oDisplacement = oDisplacement + GetSegmentDisplacement(
                                      oPoint, m_oSecondVirtualNode,
                                      m_oSecondOriginalNode, m_oBurgersVector);
  oDisplacement =
      oDisplacement + GetSegmentDisplacement(oPoint, m_oSecondOriginalNode,
                                             m_oOrigin, m_oBurgersVector);
  return oDisplacement;
}
Vector SurfaceSegment::GetSegmentDisplacement(const Point &oPoint,
                                              const Point &oEndPoint1,
                                              const Point &oEndPoint2,
                                              const Vector &oBurgersVector) {
  Vector oDisplacement(0.0, 0.0, 0.0);
  unsigned int i = 0;
  double dEta = 0.0;
  double dWeight = 0.0;
  Point oGaussPoint;
  Vector oTangent;
  double dR = 0.0;
  double dTolerance = 1.0E-6;
  // double dMinR = 2.0*dTolerance;
  Vector oRVector;
  double dX = 0.0;
  double dY = 0.0;
  double dZ = 0.0;
  double dBX = oBurgersVector.GetX();
  double dBY = oBurgersVector.GetY();
  double dBZ = oBurgersVector.GetZ();
  double dTemp1 = 0.0;
  double dTemp2 = 0.0;
  double dTemp3 = 0.0;
  double dTemp4 = 0.0;
  double dUX = 0.0;
  double dUY = 0.0;
  double dUZ = 0.0;
  double dK = 0.5 / (1.0 - PoissonsRatio);
  oTangent.SetByPoints(oEndPoint1, oEndPoint2);
  // solution according to the expression in Peach paper
  for (i = 0; i < GaussPointsCount; i++) {
    dEta = 0.5 * (GaussPointsLocations[i] + 1.0);
    dWeight = GaussPointsWeights[i];
    oGaussPoint = oEndPoint1 * (1.0 - dEta) + oEndPoint2 * dEta;
    oRVector = Vector(oGaussPoint, oPoint);

    dX = oRVector.GetX();
    dY = oRVector.GetY();
    dZ = oRVector.GetZ();
    dR = oRVector.Length();
    if (dR < dTolerance) {
      continue;
    }
    // displacement x component

    dTemp1 = dY * dZ / 3.0 / dR *
             (1 / (dX * dX + dZ * dZ) - 1 / (dX * dX + dY * dY)) * dBX;
    dTemp2 = dX * dZ / 3.0 / dR *
             (1 / (dX * dX + dY * dY) - 1 / (dZ * dZ + dY * dY)) * dBX;
    dTemp3 = dY * dX / 3.0 / dR *
             (1 / (dY * dY + dZ * dZ) - 1 / (dX * dX + dZ * dZ)) * dBX;
    // check for NaN
    if (dTemp1 != dTemp1) {
      dTemp1 = 0.0;
    }
    if (dTemp2 != dTemp2) {
      dTemp2 = 0.0;
    }
    if (dTemp3 != dTemp3) {
      dTemp3 = 0.0;
    }

    dTemp1 = dTemp1 + dX * dK / dR / dR / dR * (-dBY * dZ + dBZ * dY);
    dTemp2 = dTemp2 + dX * dK / dR / dR / dR * (-dBZ * dX + dBX * dZ) -
             (1.0 - dK) / dR * dBZ;
    dTemp3 = dTemp3 + dX * dK / dR / dR / dR * (-dBX * dY + dBY * dX) +
             (1.0 - dK) / dR * dBY;

    dTemp1 = dTemp1 * oTangent.GetX();
    dTemp2 = dTemp2 * oTangent.GetY();
    dTemp3 = dTemp3 * oTangent.GetZ();

    dTemp4 = dTemp1 + dTemp2 + dTemp3;
    dTemp4 = dTemp4 * dWeight;
    dUX = dUX + dTemp4;

    // displacement y component
    dTemp1 = dZ * dX / 3.0 / dR *
             (1 / (dY * dY + dX * dX) - 1 / (dY * dY + dZ * dZ)) * dBY;
    dTemp2 = dY * dX / 3.0 / dR *
             (1 / (dY * dY + dZ * dZ) - 1 / (dX * dX + dZ * dZ)) * dBY;
    dTemp3 = dZ * dY / 3.0 / dR *
             (1 / (dZ * dZ + dX * dX) - 1 / (dY * dY + dX * dX)) * dBY;
    // check for NaN
    if (dTemp1 != dTemp1) {
      dTemp1 = 0.0;
    }
    if (dTemp2 != dTemp2) {
      dTemp2 = 0.0;
    }
    if (dTemp3 != dTemp3) {
      dTemp3 = 0.0;
    }
    dTemp1 = dTemp1 + dY * dK / dR / dR / dR * (-dBZ * dX + dBX * dZ);
    dTemp2 = dTemp2 + dY * dK / dR / dR / dR * (-dBX * dY + dBY * dX) -
             (1.0 - dK) / dR * dBX;
    dTemp3 = dTemp3 + dY * dK / dR / dR / dR * (-dBY * dZ + dBZ * dY) +
             (1.0 - dK) / dR * dBZ;

    dTemp1 = dTemp1 * oTangent.GetY();
    dTemp2 = dTemp2 * oTangent.GetZ();
    dTemp3 = dTemp3 * oTangent.GetX();
    dTemp4 = dTemp1 + dTemp2 + dTemp3;
    dTemp4 = dTemp4 * dWeight;
    dUY = dUY + dTemp4;

    // displacement z component
    dTemp1 = dX * dY / 3.0 / dR *
             (1 / (dZ * dZ + dY * dY) - 1 / (dZ * dZ + dX * dX)) * dBZ;
    dTemp2 = dZ * dY / 3.0 / dR *
             (1 / (dZ * dZ + dX * dX) - 1 / (dY * dY + dX * dX)) * dBZ;
    dTemp3 = dX * dZ / 3.0 / dR *
             (1 / (dX * dX + dY * dY) - 1 / (dZ * dZ + dY * dY)) * dBZ;
    // check for NaN
    if (dTemp1 != dTemp1) {
      dTemp1 = 0.0;
    }
    if (dTemp2 != dTemp2) {
      dTemp2 = 0.0;
    }
    if (dTemp3 != dTemp3) {
      dTemp3 = 0.0;
    }
    dTemp1 = dTemp1 + dZ * dK / dR / dR / dR * (-dBX * dY + dBY * dX);
    dTemp2 = dTemp2 + dZ * dK / dR / dR / dR * (-dBY * dZ + dBZ * dY) -
             (1.0 - dK) / dR * dBY;
    dTemp3 = dTemp3 + dZ * dK / dR / dR / dR * (-dBZ * dX + dBX * dZ) +
             (1.0 - dK) / dR * dBX;

    dTemp1 = dTemp1 * oTangent.GetZ();
    dTemp2 = dTemp2 * oTangent.GetX();
    dTemp3 = dTemp3 * oTangent.GetY();
    dTemp4 = dTemp1 + dTemp2 + dTemp3;
    dTemp4 = dTemp4 * dWeight;
    dUZ = dUZ + dTemp4;
  }
  double dLength = oTangent.Length();
  oDisplacement.Set(dUX / 4.0 / PI / dLength, dUY / 4.0 / PI / dLength,
                    dUZ / 4.0 / PI / dLength);
  return oDisplacement;
}
void SurfaceSegment::Print() const {
  printf("(%lf,%lf,%lf) -> (%lf,%lf,%lf) -> (%lf,%lf,%lf) -> (%lf,%lf,%lf) -> "
         "(%lf,%lf,%lf) @ (%lf,%lf,%lf) / (%lf,%lf,%lf)\n",
         m_oOrigin.GetX(), m_oOrigin.GetY(), m_oOrigin.GetZ(),
         m_oFirstOriginalNode.GetX(), m_oFirstOriginalNode.GetY(),
         m_oFirstOriginalNode.GetZ(), m_oFirstVirtualNode.GetX(),
         m_oFirstVirtualNode.GetY(), m_oFirstVirtualNode.GetZ(),
         m_oSecondVirtualNode.GetX(), m_oSecondVirtualNode.GetY(),
         m_oSecondVirtualNode.GetZ(), m_oSecondOriginalNode.GetX(),
         m_oSecondOriginalNode.GetY(), m_oSecondOriginalNode.GetZ(),
         m_oBurgersVector.GetX(), m_oBurgersVector.GetY(),
         m_oBurgersVector.GetZ(), m_oSlipPlaneNormal.GetX(),
         m_oSlipPlaneNormal.GetY(), m_oSlipPlaneNormal.GetZ());
}
bool SurfaceSegment::IsPostMergeable(SurfaceSegment *poSegment) const {
  if (!IsPostTouching(poSegment)) {
    return false;
  }

  if (!m_oBurgersVector.IsSameDirection(poSegment->m_oBurgersVector)) {
    return false;
  }
  if (!m_oSlipPlaneNormal.IsSameDirection(poSegment->m_oSlipPlaneNormal)) {
    return false;
  }

  Vector oNormal = m_oFirstOriginalNode.GetSurfaceNormal();
  if (oNormal.IsSameDirection(m_oSecondOriginalNode.GetSurfaceNormal())) {
    if (oNormal.IsSameDirection(
            poSegment->m_oFirstOriginalNode.GetSurfaceNormal())) {
      if (oNormal.IsSameDirection(
              poSegment->m_oSecondOriginalNode.GetSurfaceNormal())) {
        return true;
      }
    }
  }
  return false;
}
void SurfaceSegment::PostMerge(SurfaceSegment *poSegment) {
  m_oSecondOriginalNode = poSegment->m_oSecondOriginalNode;
  m_oSecondVirtualNode = poSegment->m_oSecondVirtualNode;
}
void SurfaceSegment::Write(FILE *fpFile) const {
  fprintf(fpFile, "0,0 : (0,0) -> (0,0) : (%lf,%lf,%lf) / (%lf,%lf,%lf)\n",
          m_oSlipPlaneNormal.GetX(), m_oSlipPlaneNormal.GetY(),
          m_oSlipPlaneNormal.GetZ(), m_oBurgersVector.GetX(),
          m_oBurgersVector.GetY(), m_oBurgersVector.GetZ());
  Vector oNormal1 = m_oFirstOriginalNode.GetSurfaceNormal();
  Vector oNormal2 = m_oSecondOriginalNode.GetSurfaceNormal();
  fprintf(fpFile,
          "(%lf,%lf,%lf) -> (%lf,%lf,%lf) : (%lf,%lf,%lf) -> (%lf,%lf,%lf)\n",
          m_oFirstOriginalNode.GetX(), m_oFirstOriginalNode.GetY(),
          m_oFirstOriginalNode.GetZ(), m_oSecondOriginalNode.GetX(),
          m_oSecondOriginalNode.GetY(), m_oSecondOriginalNode.GetZ(),
          oNormal1.GetX(), oNormal1.GetY(), oNormal1.GetZ(), oNormal2.GetX(),
          oNormal2.GetY(), oNormal2.GetZ());
}
void SurfaceSegment::Initialize() {
  m_oOrigin.Set(0.0, 0.0, 0.0);
  m_oFirstOriginalNode.Set(0.0, 0.0, 0.0);
  m_oSecondOriginalNode.Set(0.0, 0.0, 0.0);
  m_oFirstVirtualNode.Set(0.0, 0.0, 0.0);
  m_oSecondVirtualNode.Set(0.0, 0.0, 0.0);
  m_oBurgersVector.Set(0.0, 0.0, 0.0);
  m_oSlipPlaneNormal.Set(0.0, 0.0, 0.0);
}
bool SurfaceSegment::IsPostTouching(SurfaceSegment *poSegment,
                                    const double &dTolerance) const {
  if (m_oSecondOriginalNode.GetDistanceSquared(
          poSegment->m_oFirstOriginalNode) < dTolerance) {
    return true;
  }
  return false;
}
double SurfaceSegment::GetLength() const {
  return m_oFirstOriginalNode.Distance(m_oSecondOriginalNode);
}
const Vector *SurfaceSegment::GetBurgersVector() const {
  return &m_oBurgersVector;
}
const Vector *SurfaceSegment::GetSlipPlaneNormal() const {
  return &m_oSlipPlaneNormal;
}
const DislocationNode *SurfaceSegment::GetFirstNode() const {
  return &m_oFirstOriginalNode;
}
const DislocationNode *SurfaceSegment::GetSecondNode() const {
  return &m_oSecondOriginalNode;
}
