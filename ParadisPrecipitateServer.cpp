#include "ParadisPrecipitateServer.h"
#include "Tools.h"
#include "DSMPI.h"
#include "Comm.h"
#include "ParadisSurface.h"
#include "ParadisCrossSlipServer.h"
#include "Randomizer.h"
#include "TriPatch.h"
#include "Decomp.h"

using namespace EZ;
using namespace SupportSystem;

// problem orientation rules
// 1. the dislocations are input in the global coordinates
// 2. the precipitates are input in the global coordinates
// 3. the APBs are generated in the global coordinates (which means that the
// spacing and the plane normal are both in global coordinates)
// 4. the APBs are written and read in global coordinates
// 5. when shearing, all Burgers vector are in global coordinates
// 6. problem boundaries are in local coordinates
// 7. the precipitate structure is all global
// 8. domain decomposition is in local coordinates

APBPoint::APBPoint() { Initialize(); }
APBPoint::~APBPoint() {}
APBPoint::APBPoint(const APBPoint &oPoint) { *this = oPoint; }
APBPoint &APBPoint::operator=(const APBPoint &oPoint) {
  m_dX = oPoint.m_dX;
  m_dY = oPoint.m_dY;
  m_iX = oPoint.m_iX;
  m_iY = oPoint.m_iY;
  m_iOwningDomain = oPoint.m_iOwningDomain;
  m_bIsShearable = oPoint.m_bIsShearable;
  m_iXShear = oPoint.m_iXShear;
  m_iYShear = oPoint.m_iYShear;
  m_iZShear = oPoint.m_iZShear;
  m_iSingletonCycle = oPoint.m_iSingletonCycle;
  return *this;
}
bool APBPoint::Shear(Vector oBurgersVector) {
  if (!m_bIsShearable)
    return false;
  double dTolerance = 1.0E-3;
  // add the new shear to the current shear
  double dTemp = oBurgersVector.GetX();
  if (dTemp > dTolerance)
    m_iXShear++;
  else if (dTemp < -dTolerance)
    m_iXShear--;

  dTemp = oBurgersVector.GetY();
  if (dTemp > dTolerance)
    m_iYShear++;
  else if (dTemp < -dTolerance)
    m_iYShear--;

  dTemp = oBurgersVector.GetZ();
  if (dTemp > dTolerance)
    m_iZShear++;
  else if (dTemp < -dTolerance)
    m_iZShear--;
  // make sure super dislocation shears are properly accounted for
  CheckSuperDislocated();
  return true;
}
unsigned int APBPoint::GetOwningDomain() const { return m_iOwningDomain; }
void APBPoint::SetOwningDomain(const unsigned int &iDomain) {
  m_iOwningDomain = iDomain;
}
bool APBPoint::IsSheared() const {
  if ((m_iXShear == 0) && (m_iYShear == 0) && (m_iZShear == 0))
    return false;
  return true;
}
void APBPoint::SetShearable() { m_bIsShearable = true; }
void APBPoint::Write(FILE *fpFile) {
  fwrite(&m_iX, sizeof(int), 1, fpFile);
  fwrite(&m_iY, sizeof(int), 1, fpFile);
  fwrite(&m_iXShear, sizeof(int), 1, fpFile);
  fwrite(&m_iYShear, sizeof(int), 1, fpFile);
  fwrite(&m_iZShear, sizeof(int), 1, fpFile);
}
void APBPoint::SetShear(const int &iXShear, const int &iYShear,
                        const int &iZShear) {
  m_iXShear = iXShear;
  m_iYShear = iYShear;
  m_iZShear = iZShear;
}
void APBPoint::Initialize() {
  m_dX = 0.0;
  m_dY = 0.0;
  m_iX = 0;
  m_iY = 0;
  m_iOwningDomain = 0;
  m_bIsShearable = false;
  m_iXShear = 0;
  m_iYShear = 0;
  m_iZShear = 0;
  m_iSingletonCycle = 0;
}
void APBPoint::CheckSuperDislocated() {
  int iX = m_iXShear;
  int iY = m_iYShear;
  int iZ = m_iZShear;
  int iV = 0;

  // reduce the Burgers vector by any combination of super dislocations' Burgers
  // vector
  iV = min(abs(iY), abs(iZ)) / 2;
  if ((iY > 0) && (iZ > 0))
    iV = -2 * iV;
  else if ((iY < 0) && (iZ < 0))
    iV = 2 * iV;
  else
    iV = 0;
  iY += iV;
  iZ += iV;

  iV = min(abs(iX), abs(iZ)) / 2;
  if ((iX > 0) && (iZ > 0))
    iV = -2 * iV;
  else if ((iX < 0) && (iZ < 0))
    iV = 2 * iV;
  else
    iV = 0;
  iX += iV;
  iZ += iV;

  iV = min(abs(iX), abs(iY)) / 2;
  if ((iX > 0) && (iY > 0))
    iV = -2 * iV;
  else if ((iX < 0) && (iY < 0))
    iV = 2 * iV;
  else
    iV = 0;
  iX += iV;
  iY += iV;

  iV = min(abs(iY), abs(iZ)) / 2;
  if ((iY > 0) && (iZ < 0))
    iV = -2 * iV;
  else if ((iY < 0) && (iZ > 0))
    iV = 2 * iV;
  else
    iV = 0;
  iY += iV;
  iZ -= iV;

  iV = min(abs(iX), abs(iZ)) / 2;
  if ((iX > 0) && (iZ < 0))
    iV = -2 * iV;
  else if ((iX < 0) && (iZ > 0))
    iV = 2 * iV;
  else
    iV = 0;
  iX += iV;
  iZ -= iV;

  iV = min(abs(iX), abs(iY)) / 2;
  if ((iX > 0) && (iY < 0))
    iV = -2 * iV;
  else if ((iX < 0) && (iY > 0))
    iV = 2 * iV;
  else
    iV = 0;
  iX += iV;
  iY -= iV;

  // if we end up with a zero dislocation or a <100> dislocation, then there is
  // no shearing, otherwise, set the shearing vector to the reduced vector
  bool bReset = false;
  if (((iX == 0) && (iY == 0)) || ((iX == 0) && (iZ == 0)) ||
      ((iY == 0) && (iZ == 0))) {
    m_iXShear = 0;
    m_iYShear = 0;
    m_iZShear = 0;
  } else {
    m_iXShear = iX;
    m_iYShear = iY;
    m_iZShear = iZ;
  }
}
int APBPoint::GetXShear() const { return m_iXShear; }
int APBPoint::GetYShear() const { return m_iYShear; }
int APBPoint::GetZShear() const { return m_iZShear; }

APB::APB() { Initialize(); }
APB::~APB() { Reset(); }
APB::APB(const APB &oAPB) { *this = oAPB; }
APB &APB::operator=(const APB &oAPB) {
  Reset();
  m_dXMin = oAPB.m_dXMin;
  m_dXMax = oAPB.m_dXMax;
  m_dYMin = oAPB.m_dYMin;
  m_dYMax = oAPB.m_dYMax;
  m_iResolution = oAPB.m_iResolution;
  m_oSlipPlaneNormal = oAPB.m_oSlipPlaneNormal;
  m_dPlaneSpacing = oAPB.m_dPlaneSpacing;
  m_dXSpacing = oAPB.m_dXSpacing;
  m_dYSpacing = oAPB.m_dYSpacing;
  m_ppoPoints = new APBPoint *[m_iResolution];
  unsigned int i = 0;
  unsigned int j = 0;
  for (i = 0; i < m_iResolution; i++) {
    m_ppoPoints[i] = new APBPoint[m_iResolution];
    for (j = 0; j < m_iResolution; j++) {
      m_ppoPoints[i][j] = oAPB.m_ppoPoints[i][j];
    }
  }
  return *this;
}
void APB::Reset() {
  ClearPoints();
  Initialize();
}
void APB::ClearPoints() {
  unsigned int i = 0;
  for (i = 0; i < m_iResolution; i++) {
    delete[] m_ppoPoints[i];
  }
  delete[] m_ppoPoints;
}
void APB::Initialize() {
  m_dXMin = 0.0;
  m_dXMax = 0.0;
  m_dYMin = 0.0;
  m_dYMax = 0.0;
  m_iResolution = 0;
  m_ppoPoints = NULL;
  m_oSlipPlaneNormal.Reset();
  m_dPlaneSpacing = 0.0;
  m_dXSpacing = 0.0;
  m_dYSpacing = 0.0;
}
Vector APB::GetSlipPlaneNormal() { return m_oSlipPlaneNormal; }
double APB::GetPlaneSpacing() const { return m_dPlaneSpacing; }
unsigned int APB::GetResolution() const { return m_iResolution; }
double APB::GetXMin() const { return m_dXMin; }
double APB::GetXMax() const { return m_dXMax; }
double APB::GetYMin() const { return m_dYMin; }
double APB::GetYMax() const { return m_dYMax; }
void APB::SetSlipPlaneNormal(const Vector &oVector) {
  m_oSlipPlaneNormal = oVector;
  m_oSlipPlaneNormal.Normalize();
}
void APB::SetResolution(const unsigned int &iResolution) {
  m_iResolution = iResolution;
}
void APB::SetXMin(const double &dValue) { m_dXMin = dValue; }
void APB::SetXMax(const double &dValue) { m_dXMax = dValue; }
void APB::SetYMin(const double &dValue) { m_dYMin = dValue; }
void APB::SetYMax(const double &dValue) { m_dYMax = dValue; }
void APB::SetPlaneSpacing(const double &dValue) { m_dPlaneSpacing = dValue; }
void APB::GeneratePoints(PrecipitateStructure *poPrecipitateStructure) {
  // generate a uniform grid of x y points between the limits
  // initialize all the points to non shearable
  // this function HAS to be called right after the creation of the APB object
  m_ppoPoints = new APBPoint *[m_iResolution];
  unsigned int i = 0;
  unsigned int j = 0;
  m_dXSpacing = (m_dXMax - m_dXMin) / (double)(m_iResolution - 1);
  m_dYSpacing = (m_dYMax - m_dYMin) / (double)(m_iResolution - 1);
  for (i = 0; i < m_iResolution; i++) {
    m_ppoPoints[i] = new APBPoint[m_iResolution];
    for (j = 0; j < m_iResolution; j++) {
      m_ppoPoints[i][j].m_dX = i * m_dXSpacing + m_dXMin;
      m_ppoPoints[i][j].m_dY = j * m_dYSpacing + m_dYMin;
      // see if the point is inside the precipitate
      if (poPrecipitateStructure->IsPointInsidePrecipitate(
              GetAPBPoint(&m_ppoPoints[i][j]))) {
        m_ppoPoints[i][j].SetShearable();
      }
      m_ppoPoints[i][j].m_iX = i;
      m_ppoPoints[i][j].m_iY = j;
    }
  }
}
void APB::Update(Node_t *poNode1, Node_t *poNode2,
                 APBEventCommunicator *poAPBEventCommunicator) {
  // get the sweeping quad bounding box
  double dXMin =
      min(min(poNode1->x, poNode1->oldx), min(poNode2->x, poNode2->oldx));
  double dXMax =
      max(max(poNode1->x, poNode1->oldx), max(poNode2->x, poNode2->oldx));
  double dYMin =
      min(min(poNode1->y, poNode1->oldy), min(poNode2->y, poNode2->oldy));
  double dYMax =
      max(max(poNode1->y, poNode1->oldy), max(poNode2->y, poNode2->oldy));
  if (dXMin < m_dXMin)
    dXMin = m_dXMin;
  if (dXMax < m_dXMin)
    dXMax = m_dXMin;
  if (dYMin < m_dYMin)
    dYMin = m_dYMin;
  if (dYMax < m_dYMin)
    dYMax = m_dYMin;
  if (dXMin > m_dXMax)
    dXMin = m_dXMax;
  if (dXMax > m_dXMax)
    dXMax = m_dXMax;
  if (dYMin > m_dYMax)
    dYMin = m_dYMax;
  if (dYMax > m_dYMax)
    dYMax = m_dYMax;
  // get the array region covered by this bounding box and test all the points
  // inside
  int iMin = (unsigned int)floor((dXMin - m_dXMin) / m_dXSpacing) - 1;
  int iMax = (unsigned int)ceil((dXMax - m_dXMin) / m_dXSpacing) + 1;
  int jMin = (unsigned int)floor((dYMin - m_dYMin) / m_dYSpacing) - 1;
  int jMax = (unsigned int)ceil((dYMax - m_dYMin) / m_dYSpacing) + 1;
  if (iMin < 0)
    iMin = 0;
  if (jMin < 0)
    jMin = 0;
  if (iMax >= m_iResolution)
    iMax = m_iResolution - 1;
  if (jMax >= m_iResolution)
    jMax = m_iResolution - 1;
  unsigned int i = 0;
  unsigned int j = 0;
  APBCellShearingEvent *poAPBShearingEvent = NULL;
  unsigned int iNeighbourID = GetArmID(poNode1, poNode2);
  for (i = iMin; i <= iMax; i++) {
    for (j = jMin; j <= jMax; j++) {
      if (IsAPBPointInSweepQuad(&m_ppoPoints[i][j], poNode1, poNode2)) {
        if (m_ppoPoints[i][j].Shear(Vector(poNode1->burgX[iNeighbourID],
                                           poNode1->burgY[iNeighbourID],
                                           poNode1->burgZ[iNeighbourID]))) {
          // create an APB flipping event and add it to the list of events
          poAPBShearingEvent = new APBCellShearingEvent;
          poAPBShearingEvent->SetAPB(this);
          poAPBShearingEvent->SetCellPoint(&m_ppoPoints[i][j]);
          poAPBEventCommunicator->AddEvent(poAPBShearingEvent);
        }
      }
    }
  }
}
void APB::Write(FILE *fpFile) {
  // write the slip plane id and spacing first
  unsigned int iSlipPlaneID = GetSlipPlaneID();
  fwrite(&iSlipPlaneID, sizeof(unsigned int), 1, fpFile);
  fwrite(&m_dPlaneSpacing, sizeof(double), 1, fpFile);
  unsigned int i = 0;
  unsigned int j = 0;
  unsigned int iCount = 0;
  for (i = 0; i < m_iResolution; i++) {
    for (j = 0; j < m_iResolution; j++) {
      if (m_ppoPoints[i][j].IsSheared())
        iCount++;
    }
  }

  fwrite(&iCount, sizeof(unsigned int), 1, fpFile);
  for (i = 0; i < m_iResolution; i++) {
    for (j = 0; j < m_iResolution; j++) {
      if (m_ppoPoints[i][j].IsSheared())
        m_ppoPoints[i][j].Write(fpFile);
    }
  }
  fflush(fpFile);
}
void APB::Read(FILE *fpFile, const unsigned int &iResolution,
               PrecipitateStructure *poPrecipitateStructure,
               unsigned int &iSlipPlaneID) {
  ClearPoints();
  SetResolution(iResolution);
  unsigned int iCount = 0;
  fread(&iSlipPlaneID, sizeof(unsigned int), 1, fpFile);
  SetFCCSlipPlaneNormalFromID(iSlipPlaneID);
  fread(&m_dPlaneSpacing, sizeof(double), 1, fpFile);
  fread(&iCount, sizeof(unsigned int), 1, fpFile);
  GeneratePoints(poPrecipitateStructure);
  // read the points
  int i = 0;
  int j = 0;
  unsigned int k = 0;
  int iX = 0;
  int iY = 0;
  int iZ = 0;
  for (k = 0; k < iCount; k++) {
    fread(&i, sizeof(int), 1, fpFile);
    fread(&j, sizeof(int), 1, fpFile);
    fread(&iX, sizeof(int), 1, fpFile);
    fread(&iY, sizeof(int), 1, fpFile);
    fread(&iZ, sizeof(int), 1, fpFile);
    m_ppoPoints[i][j].SetShear(iX, iY, iZ);
  }
}
bool APB::IsAPBPointInSweepQuad(APBPoint *poPoint, Node_t *poNode1,
                                Node_t *poNode2) {
  // 	GenericNode oP1(poNode1->oldx,poNode1->oldy,poNode1->oldz);
  // 	GenericNode oP2(poNode2->oldx,poNode2->oldy,poNode2->oldz);
  // 	GenericNode oP3(poNode1->x,poNode1->y,poNode1->z);
  // 	GenericNode oP4(poNode2->x,poNode2->y,poNode2->z);
  //
  // 	TriPatch oTri1(&oP1,&oP2,&oP3);
  // 	TriPatch oTri2(&oP4,&oP2,&oP3);
  // 	Point oAPBPoint = GetAPBPoint(poPoint);

  // Sometimes, the segments are off plane due to the APB selection process
  // which returns the nearest APB within a certain tolerance. Hence, it is
  // better to check whether the segment's projection sweeps the APB point's
  // projection instead of checking the actual 3D sweep. Plus, it helps avoiding
  // some roundoff errors since some numbers are strictly zero in the projection
  // check case
  GenericNode oP1(poNode1->oldx, poNode1->oldy, 0.0);
  GenericNode oP2(poNode2->oldx, poNode2->oldy, 0.0);
  GenericNode oP3(poNode1->x, poNode1->y, 0.0);
  GenericNode oP4(poNode2->x, poNode2->y, 0.0);
  TriPatch oTri1(&oP1, &oP2, &oP3);
  TriPatch oTri2(&oP4, &oP2, &oP3);

  Point oAPBPoint = GetAPBPoint(poPoint);
  oAPBPoint.SetZ(0.0);

  bool bIsInTri1 = oTri1.IsPointInTriangle(oAPBPoint);
  bool bIsInTri2 = oTri2.IsPointInTriangle(oAPBPoint);
  if ((bIsInTri1 && !bIsInTri2) || (!bIsInTri1 && bIsInTri2))
    return true;
  return false;
}
bool APB::IsAPBPointInSweepTri(APBPoint *poAPBPoint, const Point &oPoint1,
                               const Point &oPoint2, const Point &oPoint3) {
  GenericNode oP1 = oPoint1;
  GenericNode oP2 = oPoint2;
  GenericNode oP3 = oPoint3;
  oP1.SetZ(0.0);
  oP2.SetZ(0.0);
  oP3.SetZ(0.0);
  TriPatch oTri(&oP1, &oP2, &oP3);

  Point oAPBPoint = GetAPBPoint(poAPBPoint);
  oAPBPoint.SetZ(0.0);
  return oTri.IsPointInTriangle(oAPBPoint);
}
Point APB::GetAPBPoint(APBPoint *poPoint) {
  Point oPoint;
  oPoint.SetX(poPoint->m_dX);
  oPoint.SetY(poPoint->m_dY);
  double dZ = m_dPlaneSpacing - m_oSlipPlaneNormal.GetX() * poPoint->m_dX -
              m_oSlipPlaneNormal.GetY() * poPoint->m_dY;
  dZ = dZ / m_oSlipPlaneNormal.GetZ();
  oPoint.SetZ(dZ);
  return oPoint;
}
void APB::UpdateTriangle(const Point &oPoint1, const Point &oPoint2,
                         const Point &oPoint3, const Vector &oBurgersVector,
                         APBEventCommunicator *poAPBEventCommunicator) {
  // get the sweeping quad bounding box
  double dXMin = min(min(oPoint1.GetX(), oPoint2.GetX()), oPoint3.GetX());
  double dXMax = max(max(oPoint1.GetX(), oPoint2.GetX()), oPoint3.GetX());
  double dYMin = min(min(oPoint1.GetY(), oPoint2.GetY()), oPoint3.GetY());
  double dYMax = max(max(oPoint1.GetY(), oPoint2.GetY()), oPoint3.GetY());
  if (dXMin < m_dXMin)
    dXMin = m_dXMin;
  if (dXMax < m_dXMin)
    dXMax = m_dXMin;
  if (dYMin < m_dYMin)
    dYMin = m_dYMin;
  if (dYMax < m_dYMin)
    dYMax = m_dYMin;
  if (dXMin > m_dXMax)
    dXMin = m_dXMax;
  if (dXMax > m_dXMax)
    dXMax = m_dXMax;
  if (dYMin > m_dYMax)
    dYMin = m_dYMax;
  if (dYMax > m_dYMax)
    dYMax = m_dYMax;
  // get the array region covered by this bounding box and test all the points
  // inside
  int iMin = (unsigned int)floor((dXMin - m_dXMin) / m_dXSpacing) - 1;
  int iMax = (unsigned int)ceil((dXMax - m_dXMin) / m_dXSpacing) + 1;
  int jMin = (unsigned int)floor((dYMin - m_dYMin) / m_dYSpacing) - 1;
  int jMax = (unsigned int)ceil((dYMax - m_dYMin) / m_dYSpacing) + 1;
  if (iMin < 0)
    iMin = 0;
  if (jMin < 0)
    jMin = 0;
  if (iMax >= m_iResolution)
    iMax = m_iResolution - 1;
  if (jMax >= m_iResolution)
    jMax = m_iResolution - 1;
  unsigned int i = 0;
  unsigned int j = 0;
  APBCellShearingEvent *poAPBShearingEvent = NULL;

  for (i = iMin; i <= iMax; i++) {
    for (j = jMin; j <= jMax; j++) {
      if (IsAPBPointInSweepTri(&m_ppoPoints[i][j], oPoint1, oPoint2, oPoint3)) {
        if (m_ppoPoints[i][j].Shear(oBurgersVector)) {
          // create an APB flipping event and add it to the list of events
          poAPBShearingEvent = new APBCellShearingEvent;
          poAPBShearingEvent->SetAPB(this);
          poAPBShearingEvent->SetCellPoint(&m_ppoPoints[i][j]);
          poAPBEventCommunicator->AddEvent(poAPBShearingEvent);
        }
      }
    }
  }
}
void APB::UpdateAPBCellOwnership(Home_t *poHome) {
  unsigned int i = 0;
  unsigned int j = 0;
  Point oCellPoint;
  unsigned int iDomain = 0;

  for (i = 0; i < m_iResolution; i++) {
    for (j = 0; j < m_iResolution; j++) {
      oCellPoint = GetAPBPoint(&m_ppoPoints[i][j]);
      iDomain = FindCoordDomain(poHome, &oCellPoint);
      m_ppoPoints[i][j].SetOwningDomain(iDomain);
    }
  }
}
void APB::ShearCell(const unsigned int &iCellX, const unsigned int &iCellY,
                    const int &iXShear, const int &iYShear,
                    const int &iZShear) {
  m_ppoPoints[iCellX][iCellY].SetShear(iXShear, iYShear, iZShear);
}
void APB::Report(const unsigned int &iProcessorID) const {
  unsigned int i = 0;
  unsigned int j = 0;
  unsigned int iCount = 0;
  for (i = 0; i < m_iResolution; i++) {
    for (j = 0; j < m_iResolution; j++) {
      if (m_ppoPoints[i][j].IsSheared())
        iCount++;
    }
  }
  printf("@ %d : n = (%e,%e,%e), sheared : %d\n", iProcessorID,
         m_oSlipPlaneNormal.GetX(), m_oSlipPlaneNormal.GetY(),
         m_oSlipPlaneNormal.GetZ(), iCount);
  fflush(NULL);
}
APBPoint *APB::GetContainingCell(const Point &oPoint) {
  int i = (unsigned int)floor((oPoint.GetX() - m_dXMin) / m_dXSpacing + 0.5);
  int j = (unsigned int)floor((oPoint.GetY() - m_dYMin) / m_dYSpacing + 0.5);
  return &m_ppoPoints[i][j];
}
unsigned int APB::GetSlipPlaneID() {
  return GetFCCSlipPlaneID(m_oSlipPlaneNormal);
}
unsigned int APB::GetFCCSlipPlaneID(const Vector &oSlipPlaneNormal) {
  double dTolerance = 1E-3;
  Vector oReferenceNormal(0.0, 0.0, 0.0);

  oReferenceNormal.Set(1.0, 1.0, 1.0);
  oReferenceNormal.Normalize();
  if (oReferenceNormal.IsSame(oSlipPlaneNormal, dTolerance) ||
      oReferenceNormal.IsOpposite(oSlipPlaneNormal, dTolerance))
    return 1;

  oReferenceNormal.Set(-1.0, 1.0, 1.0);
  oReferenceNormal.Normalize();
  if (oReferenceNormal.IsSame(oSlipPlaneNormal, dTolerance) ||
      oReferenceNormal.IsOpposite(oSlipPlaneNormal, dTolerance))
    return 2;

  oReferenceNormal.Set(1.0, -1.0, 1.0);
  oReferenceNormal.Normalize();
  if (oReferenceNormal.IsSame(oSlipPlaneNormal, dTolerance) ||
      oReferenceNormal.IsOpposite(oSlipPlaneNormal, dTolerance))
    return 3;

  oReferenceNormal.Set(1.0, 1.0, -1.0);
  oReferenceNormal.Normalize();
  if (oReferenceNormal.IsSame(oSlipPlaneNormal, dTolerance) ||
      oReferenceNormal.IsOpposite(oSlipPlaneNormal, dTolerance))
    return 4;

  return 0;
}
Vector APB::GetFCCSlipPlaneNormalFromID(const unsigned int &iSlipPlaneID) {
  Vector oSlipPlaneNormal(0.0, 0.0, 0.0);
  if (iSlipPlaneID == 1)
    oSlipPlaneNormal.Set(1.0, 1.0, 1.0);
  if (iSlipPlaneID == 2)
    oSlipPlaneNormal.Set(-1.0, 1.0, 1.0);
  if (iSlipPlaneID == 3)
    oSlipPlaneNormal.Set(1.0, -1.0, 1.0);
  if (iSlipPlaneID == 4)
    oSlipPlaneNormal.Set(1.0, 1.0, -1.0);
  oSlipPlaneNormal.Normalize();
  return oSlipPlaneNormal;
}
void APB::SetFCCSlipPlaneNormalFromID(const unsigned int &iSlipPlaneID) {
  m_oSlipPlaneNormal = GetFCCSlipPlaneNormalFromID(iSlipPlaneID);
}
void APB::CleanSingletons(const unsigned int &iCycle,
                          const unsigned int &iProcessorID,
                          APBEventCommunicator *poAPBEventCommunicator) {
  unsigned int i = 0;
  unsigned int j = 0;
  // start by detecting all singletons and setting their cycles
  for (i = 1; i < m_iResolution - 1; i++) {
    for (j = 1; j < m_iResolution - 1; j++) {
      // see if the point i,j is a singleton, if not, reset its singleton cycle,
      // if it is, see if its singleton cycle is 0 or not, if it is 0, set it to
      // the current cycle so that the system can start counting the number of
      // cycles the cell has been a singleton and cleans it at the right time,
      // if it is not 0, this means that it has been set before, leave it as it
      // is
      if (IsSingleton(i, j)) {
        if (m_ppoPoints[i][j].m_iSingletonCycle == 0)
          m_ppoPoints[i][j].m_iSingletonCycle = iCycle;
      } else
        m_ppoPoints[i][j].m_iSingletonCycle = 0;
    }
  }
  // in a different loop, see if the point has been a singleton for too long, if
  // this is the case, clean it. Do NOT clean the cells that are not owned by
  // this processor
  unsigned int iSingletonLifeTime = 500;
  APBCellShearingEvent *poAPBShearingEvent = NULL;
  for (i = 1; i < m_iResolution - 1; i++) {
    for (j = 1; j < m_iResolution - 1; j++) {
      if (m_ppoPoints[i][j].m_iSingletonCycle == 0)
        continue;
      if (m_ppoPoints[i][j].GetOwningDomain() != iProcessorID)
        continue;
      if ((iCycle - m_ppoPoints[i][j].m_iSingletonCycle) > iSingletonLifeTime) {
        // this cell has been a singleton for too long, unshear it
        m_ppoPoints[i][j].SetShear(0, 0, 0);
        // inform all the other domains about it
        poAPBShearingEvent = new APBCellShearingEvent;
        poAPBShearingEvent->SetAPB(this);
        poAPBShearingEvent->SetCellPoint(&m_ppoPoints[i][j]);
        poAPBEventCommunicator->AddEvent(poAPBShearingEvent);
      }
    }
  }
}
bool APB::IsSingleton(const unsigned int &iX, const unsigned int &iY) {
  // a singleton is a sheared cell surrounded by all unsheared cells
  if (m_ppoPoints[iX][iY].IsSheared()) {
    if (!m_ppoPoints[iX - 1][iY - 1].IsSheared() &&
        !m_ppoPoints[iX - 1][iY].IsSheared() &&
        !m_ppoPoints[iX - 1][iY + 1].IsSheared() &&
        !m_ppoPoints[iX][iY - 1].IsSheared() &&
        !m_ppoPoints[iX][iY + 1].IsSheared() &&
        !m_ppoPoints[iX + 1][iY - 1].IsSheared() &&
        !m_ppoPoints[iX + 1][iY].IsSheared() &&
        !m_ppoPoints[iX + 1][iY + 1].IsSheared()) {
      return true;
    }
  }
  return false;
}

ParadisPrecipitateServer
    *ParadisPrecipitateServer::m_poParadisPrecipitateServerInstance = NULL;
unsigned int ParadisPrecipitateServer::m_iAPBPointResolution = 150;
ParadisPrecipitateServer *ParadisPrecipitateServer::CreateInstance() {
  if (m_poParadisPrecipitateServerInstance == NULL) {
    m_poParadisPrecipitateServerInstance = new ParadisPrecipitateServer;
  }
  return m_poParadisPrecipitateServerInstance;
}
ParadisPrecipitateServer::~ParadisPrecipitateServer() { Reset(); }
void ParadisPrecipitateServer::Reset() {
  m_oStructure.Reset();
  unsigned int i = 0;
  list<APB *>::iterator liAPB;
  unsigned int iSlipPlanesCount = m_vlpAPB.size();
  for (i = 0; i < iSlipPlanesCount; i++) {
    for (liAPB = m_vlpAPB[i].begin(); liAPB != m_vlpAPB[i].end(); liAPB++) {
      if ((*liAPB) != NULL) {
        delete (*liAPB);
      }
    }
    m_vlpAPB[i].clear();
  }
  m_vlpAPB.clear();
  m_poAPBEventCommunicator->Reset();
}
bool ParadisPrecipitateServer::Set(Home_t *poHome) {
  Reset();

  m_bShearablePrecipitates = (poHome->param->ShearablePrecipitates != 0);
  m_bPairShearing = (poHome->param->PairPrecipitateShearing != 0);
  m_dNucleationRate = poHome->param->PrecipitateNucleationRate;
  m_dNucleationRadius = poHome->param->PrecipitateNucleationRadius;
  m_iNucleationMaxCount = poHome->param->PrecipitateNucleationMaxCount;
  if (m_dNucleationRate < 0.0)
    m_dNucleationRate = 0.0;
  if (m_dNucleationRadius < 0.0)
    m_dNucleationRadius = 0.0;
  m_bUseNucleationMode =
      (!m_bShearablePrecipitates) && (m_dNucleationRate > 0.0) &&
      (m_dNucleationRadius > 0.0);

  m_iDomainID = poHome->myDomain;
  m_poAPBEventCommunicator->SetDomainID(m_iDomainID);
  m_poAPBEventCommunicator->SetDomainsCount(poHome->numDomains);

  m_dXMin = -0.5 * poHome->param->Dimensions[0];
  m_dXMax = 0.5 * poHome->param->Dimensions[0];
  m_dYMin = -0.5 * poHome->param->Dimensions[1];
  m_dYMax = 0.5 * poHome->param->Dimensions[1];

  string sPrecipitatesFileName = string(poHome->param->PrecipitatesFileName);
  if (sPrecipitatesFileName.find_first_not_of(" \t\r\n") == string::npos) {
    sPrecipitatesFileName = "";
  }
  bool bStructureLoaded = false;
  if (!sPrecipitatesFileName.empty()) {
    // try to read the precipitates file
    if (m_oStructure.Set(sPrecipitatesFileName, poHome->param->Dimensions[0],
                         poHome->param->Dimensions[1],
                         poHome->param->Dimensions[2])) {
      bStructureLoaded = true;
    } else {
      char cFilePath[256];
      sprintf(cFilePath, "../%s", sPrecipitatesFileName.c_str());
      if (m_oStructure.Set(string(cFilePath), poHome->param->Dimensions[0],
                           poHome->param->Dimensions[1],
                           poHome->param->Dimensions[2])) {
        bStructureLoaded = true;
      }
    }
  }

  if (!bStructureLoaded && !m_bUseNucleationMode) {
    return false;
  }

  if (bStructureLoaded) {
    // processor 0 writes the precipitate structure
    if (poHome->myDomain == 0) {
      m_oStructure.WritePrecipitateStructureVTK();
    }
    DSMPI::Barrier();
  }

  // finally, allocate the APB array if required
  if (m_bShearablePrecipitates) {
    unsigned int iFCCSlipPlanesCount = 4;
    m_vlpAPB.resize(iFCCSlipPlanesCount);
    unsigned int i = 0;
    for (i = 0; i < iFCCSlipPlanesCount; i++) {
      m_vlpAPB[i].clear();
    }
  }

  ReadAPBs(poHome->param->APBFileName);
  if (!m_oStructure.IsEmpty()) {
    AxisAlignedBoundingBox oBox = m_oStructure.GetBoundingBox();
    m_dXMin = oBox.GetXMin();
    m_dXMax = oBox.GetXMax();
    m_dYMin = oBox.GetYMin();
    m_dYMax = oBox.GetYMax();
  }
  m_iTotalNucleatedCount = m_oStructure.GetPrecipitatesCount();
  return true;
}
unsigned int
ParadisPrecipitateServer::GetNucleationCount(const double &dExpectedCount) const {
  if (dExpectedCount <= 0.0)
    return 0;

  // Use exact Poisson sampling for small means.
  if (dExpectedCount < 30.0) {
    const double dLimit = exp(-dExpectedCount);
    double dProduct = 1.0;
    unsigned int iCount = 0;
    do {
      iCount++;
      dProduct *= Randomizer::Random();
    } while (dProduct > dLimit);
    return iCount - 1;
  }

  // Normal approximation for large means.
  double dCount = Randomizer::RandomNormal(dExpectedCount, sqrt(dExpectedCount));
  if (dCount <= 0.0)
    return 0;
  return (unsigned int)floor(dCount + 0.5);
}
Point ParadisPrecipitateServer::GetRandomNucleationCenter(
    const Home_t *poHome) const {
  const double dHalfX = 0.5 * poHome->param->Dimensions[0];
  const double dHalfY = 0.5 * poHome->param->Dimensions[1];
  const double dHalfZ = 0.5 * poHome->param->Dimensions[2];
  return Point(Randomizer::Random(-dHalfX, dHalfX),
               Randomizer::Random(-dHalfY, dHalfY),
               Randomizer::Random(-dHalfZ, dHalfZ));
}
void ParadisPrecipitateServer::NucleatePrecipitates(Home_t *poHome) {
  if (!m_bUseNucleationMode)
    return;

  const double dExpectedCount = m_dNucleationRate * poHome->param->deltaTT;
  if (dExpectedCount <= 0.0)
    return;

  int iNucleationCount = 0;
  if (m_iDomainID == 0) {
    iNucleationCount = (int)GetNucleationCount(dExpectedCount);
    if (m_iNucleationMaxCount >= 0) {
      int iRemaining =
          m_iNucleationMaxCount - (int)m_oStructure.GetPrecipitatesCount();
      if (iRemaining <= 0) {
        iNucleationCount = 0;
      } else if (iNucleationCount > iRemaining) {
        iNucleationCount = iRemaining;
      }
    }
  }
  DSMPI::Broadcast(&iNucleationCount, 1, MPI_INT, 0);

  unsigned int iAddedCount = 0;
  for (int i = 0; i < iNucleationCount; i++) {
    double daCenter[3] = {0.0, 0.0, 0.0};
    if (m_iDomainID == 0) {
      Point oCenter = GetRandomNucleationCenter(poHome);
      daCenter[0] = oCenter.GetX();
      daCenter[1] = oCenter.GetY();
      daCenter[2] = oCenter.GetZ();
    }
    DSMPI::Broadcast(daCenter, 3, MPI_DOUBLE, 0);
    if (m_oStructure.AddSphericalPrecipitate(
            Point(daCenter[0], daCenter[1], daCenter[2]), m_dNucleationRadius)) {
      iAddedCount++;
    }
  }
  m_iTotalNucleatedCount += iAddedCount;
}
void ParadisPrecipitateServer::CheckNodes(Home_t *poHome) {
  NucleatePrecipitates(poHome);
  if (m_oStructure.IsEmpty()) {
    return;
  }
  ClearOpList(poHome);
  unsigned int iSize = (unsigned int)poHome->newNodeKeyPtr;
  unsigned int i = 0;
  Node_t *poNode = NULL;
  for (i = 0; i < iSize; i++) {
    poNode = poHome->nodeKeys[i];
    if (poNode == NULL)
      continue;
    if ((poNode->constraint == PINNED_NODE) ||
        (poNode->constraint >=SURFACE_NODE))
      continue;
    // if the precipitates are not shearable, retract any nodes that are inside
    // the precipitates to the surface
    RetractNode(poHome, poNode);
  }
  CommSendRemesh(poHome);
  FixRemesh(poHome);
  DSMPI::Barrier();
}
void ParadisPrecipitateServer::ComputeAPBForces(Home_t *poHome, Node_t *poNode1,
                                                Node_t *poNode2,
                                                double pdForce1[3],
                                                double pdForce2[3]) {
  pdForce1[0] = 0.0;
  pdForce1[1] = 0.0;
  pdForce1[2] = 0.0;
  pdForce2[0] = 0.0;
  pdForce2[1] = 0.0;
  pdForce2[2] = 0.0;

  if (m_oStructure.IsEmpty())
    return;

  // if the segment is a super dislocation, the APB forces are zero
  int iID = GetArmID(poNode1, poNode2);
  if (IsSuperDislocation(poNode1, iID))
    return;

  // if the segment's midpoint is not inside a precipitate, the APB forces are
  // zero
  Point oMidPoint(0.5 * (poNode1->x + poNode2->x),
                  0.5 * (poNode1->y + poNode2->y),
                  0.5 * (poNode1->z + poNode2->z));
  Polyhedron *poPrecipitate = m_oStructure.GetContainingPrecipitate(&oMidPoint);
  if (poPrecipitate == NULL)
    return;

  // get the APB
  APB *poAPB = GetAPB(poNode1, poNode2);
  if (poAPB == NULL)
    return;

  // get the center of the APB cell containing segment's midpoint
  APBPoint *poAPBCell = poAPB->GetContainingCell(oMidPoint);
  Point oCenterPoint = poAPB->GetAPBPoint(poAPBCell);

  // the force vector direction points to the side where the cell center is
  // (but it doesn't point towards the cell center)
  Vector oSegment(poNode2->x - poNode1->x, poNode2->y - poNode1->y,
                  poNode2->z - poNode1->z);
  double dSegmentLength = oSegment.Length();
  oSegment.Normalize();
  Vector oCenter(oMidPoint, oCenterPoint);
  Vector oForceDirection = oCenter - oSegment * (oCenter * oSegment);
  oForceDirection.Normalize();

  // if the cell is sheared, then the force is towards the cell's center and
  // vice versa because the APB force works to minimize the sheared region
  if (!poAPBCell->IsSheared())
    oForceDirection.Reverse();

  // the force direction is guaranteed to be on plane, scale it by the magnitude
  // of the APB force
  Vector oForce = oForceDirection *
                  (0.5 * poHome->param->APBEnrgy / poHome->param->burgMag);
  // we are looking for the force per unit length, multiply the obtained force
  // by the segment's length
  oForce = oForce * dSegmentLength;
  pdForce1[0] = oForce.GetX();
  pdForce1[1] = oForce.GetY();
  pdForce1[2] = oForce.GetZ();
  pdForce2[0] = pdForce1[0];
  pdForce2[1] = pdForce1[1];
  pdForce2[2] = pdForce1[2];
}
void ParadisPrecipitateServer::AddAPBShearing(Home_t *poHome) {
  if (m_oStructure.IsEmpty())
    return;
  if (poHome->param->APBEnrgy <= 0.0)
    return;
  double dAPBShearStress = poHome->param->APBEnrgy / poHome->param->burgMag;
  unsigned int i = 0;
  unsigned int j = 0;
  Node_t *poNode = NULL;
  Node_t *poNeighbour = NULL;
  Polyhedron *poPrecipitate = NULL;
  unsigned int iConstraintType = 0;
  Vector oSlipPlaneNormal;
  bool bIsSuperDislocation = false;
  double dLength = 0.0;
  double dTemp = 0.0;
  APB *poAPB = NULL;
  APBPoint *poAPBCell = NULL;
  Point oCenterPoint;
  Vector oNodeForce;
  Vector oNodeForceDirection;
  double dAPBForceMagnitude = 0.0;
  Vector oNodePoint;
  Vector oCenter;
  bool bIsAPBIncreasing = false;
  Vector oAPBForce;
  // loop over all the arms and subtract the APB resistance forces from the
  // total node force
  for (i = 0; i < poHome->newNodeKeyPtr; i++) {
    poNode = poHome->nodeKeys[i];
    if (poNode == NULL)
      continue;
    // if the node is not inside a precipitate, don't do anything to it
    if (!IsNodeInsideAPrecipitate(poNode, poPrecipitate))
      continue;
    // only nodes moving on planes will be resisted as they shear through the
    // precipitate
    ParadisSurface::GetNodeDynamicConstraint(poNode, iConstraintType,
                                             oSlipPlaneNormal);
    if (iConstraintType != 1)
      continue;
    // check to see whether the node belongs to a super dislocation, also get
    // the node's effective length
    oNodePoint.Set(poNode->x, poNode->y, poNode->z);
    bIsSuperDislocation = false;
    dLength = 0.0;
    for (j = 0; j < poNode->numNbrs; j++) {
      if (IsSuperDislocation(poNode, j)) {
        bIsSuperDislocation = true;
        break;
      }
      poNeighbour = GetNeighborNode(poHome, poNode, j);
      dTemp = (poNeighbour->x - poNode->x) * (poNeighbour->x - poNode->x);
      dTemp += (poNeighbour->y - poNode->y) * (poNeighbour->y - poNode->y);
      dTemp += (poNeighbour->z - poNode->z) * (poNeighbour->z - poNode->z);
      dLength += 0.5 * sqrt(dTemp);
    }
    // if the node belongs to a super dislocation, the APB forces are zero
    if (bIsSuperDislocation)
      continue;
    // get the APB on which this node is moving
    poAPB = GetAPB(poNode, oSlipPlaneNormal);
    if (poAPB == NULL)
      continue;
    // get the center of the APB cell containing the node
    poAPBCell = poAPB->GetContainingCell(oNodePoint);
    oCenterPoint = poAPB->GetAPBPoint(poAPBCell);
    oNodeForce.Set(poNode->fX, poNode->fY, poNode->fZ);
    // get the component of this force on the slip plane normal, this is the
    // effective force component that will be resisted by the APB force
    oNodeForce =
        oNodeForce - oSlipPlaneNormal * (oNodeForce * oSlipPlaneNormal);
    oNodeForceDirection = oNodeForce;
    oNodeForceDirection.Normalize();
    // The magnitude of the APB friction force is f = gamma*L/b
    dAPBForceMagnitude = fabs(dLength * dAPBShearStress);
    // the total force can work on either increasing or decreasing the total
    // sheared area in the case of decreasing it, then add the APB force to the
    // total force, in case of increasing it, if the APB force is greater than
    // the total force, zero out the force, if it is less, subtract it from the
    // total force get the vector pointing to the cell's center from the node
    // and normalize it
    oCenter.SetByPoints(oNodePoint, oCenterPoint);
    oCenter.Normalize();
    // assume unsheared cells first, if the node is moving towards the center of
    // the cell, then its motion increases the APB, otherwise, it decreases it
    if (oNodeForceDirection * oCenter > 0.0)
      bIsAPBIncreasing = true;
    else
      bIsAPBIncreasing = false;
    // flip the case if the cell is sheared
    if (poAPBCell->IsSheared())
      bIsAPBIncreasing = !bIsAPBIncreasing;

    // finally, if the APB area is increasing, reduce the total force,
    // otherwise, add the APB force to the total force
    if (bIsAPBIncreasing) {
      if (poHome->param->ShearablePrecipitates == 0) {
        // the dislocation wants to shear through an unshearable precipitate,
        // stop it
        poNode->fX = 0.0;
        poNode->fY = 0.0;
        poNode->fZ = 0.0;
        continue;
      } else {
        // precipitates can be sheared
        if (dAPBForceMagnitude > oNodeForce.Length()) {
          poNode->fX = 0.0;
          poNode->fY = 0.0;
          poNode->fZ = 0.0;
          continue;
        }
        // get the APB force vector and subtract it from the total force
        oAPBForce = oNodeForceDirection * dAPBForceMagnitude;
        poNode->fX -= oAPBForce.GetX();
        poNode->fY -= oAPBForce.GetY();
        poNode->fZ -= oAPBForce.GetZ();
        // printf("shear happens !!!\n");
      }
    } else {
      if (poHome->param->ShearablePrecipitates == 0) {
        // precipitates cannot be sheared, however, the dislocation is
        // unshearing them, so let it be and don't change the arm forces
      } else {
        // get the APB force vector and add it to the total force
        oAPBForce = oNodeForceDirection * dAPBForceMagnitude;
        poNode->fX += oAPBForce.GetX();
        poNode->fY += oAPBForce.GetY();
        poNode->fZ += oAPBForce.GetZ();
      }
    }
  }
}
/*void ParadisPrecipitateServer::AddAPBShearing(Home_t* poHome)
{
        if(m_oStructure.IsEmpty())				return;
        if(poHome->param->APBEnrgy <= 0.0)		return;
        double dAPBShearStress = poHome->param->APBEnrgy/poHome->param->burgMag;
        // if there is no friction in the system, return
        unsigned int i = 0;
        unsigned int j = 0;
        Node_t* poNode = NULL;
        Node_t* poNeighbour = NULL;
        Point oMidPoint;
        APB* poAPB = NULL;
        APBPoint* poAPBCell = NULL;
        Point oCenterPoint;
        Vector oSlipPlaneNormal;
        Vector oLineDirection;
        Vector oBurgersVector;
        double dLength = 0.0;
        Vector oAPBForce;
        Vector oArmForce;
        Vector oArmForceDirection;
        double dAPBForceMagnitude = 0.0;
        Polyhedron* poPrecipitate = NULL;
        Vector oCenter;
        bool bIsAPBIncreasing = false;
        unsigned int iArmID = 0;
        // loop over all the arms and subtraction the friction forces from the
total arm force for(i = 0 ; i < poHome->newNodeKeyPtr ; i++)
        {
                poNode = poHome->nodeKeys[i];
                if(poNode == NULL)				continue;
                if(!IsNodeInsideAPrecipitate(poNode,poPrecipitate))
continue; for(j = 0 ; j < poNode->numNbrs ; j++)
                {
                        // if the segment is a super dislocation, the APB forces
are zero if(IsSuperDislocation(poNode,j))
continue; poNeighbour = GetNeighborNode(poHome,poNode,j); if(poNeighbour ==
NULL)				continue;
                        // if the segment's midpoint is not inside a
precipitate, the APB forces are zero oMidPoint.Set(0.5*(poNode->x +
poNeighbour->x),0.5*(poNode->y + poNeighbour->y),0.5*(poNode->z +
poNeighbour->z)); poPrecipitate =
m_oStructure.GetContainingPrecipitate(&oMidPoint); if(poPrecipitate == NULL)
continue;
                        // get the APB
                        if(!m_bPairShearing)
                        {
                                poAPB = GetAPB(poNode,poNeighbour);
                                if(poAPB == NULL)
continue;
                                // get the center of the APB cell containing
segment's midpoint poAPBCell = poAPB->GetContainingCell(oMidPoint); oCenterPoint
= poAPB->GetAPBPoint(poAPBCell);
                        }

                        // get the slip plane normal and burgers vector
                        oSlipPlaneNormal.Set(poNode->nx[j],poNode->ny[j],poNode->nz[j]);
                        oSlipPlaneNormal.Normalize();
                        oBurgersVector.Set(poNode->burgX[j],poNode->burgY[j],poNode->burgZ[j]);

                        // get the dislocation segment line direction
                        oLineDirection.SetByPoints(Point(poNode->x,poNode->y,poNode->z),Point(poNeighbour->x,poNeighbour->y,poNeighbour->z));
                        dLength = oLineDirection.Length();
                        oLineDirection.Normalize();

                        // get the total force on this arm (excluding the APB
force), the force should come
                        // from BOTH nodes
                        iArmID = GetArmID(poNeighbour,poNode);
                        oArmForce.Set(poNode->armfx[j],poNode->armfy[j],poNode->armfz[j]);
                        printf("single : %e\n",oArmForce.Length());
                        oArmForce.Set(poNode->armfx[j] +
poNeighbour->armfx[iArmID],poNode->armfy[j] +
poNeighbour->armfy[iArmID],poNode->armfz[j] + poNeighbour->armfz[iArmID]);

                        //printf("xforce ws:
%e,%e\n",poNode->armfx[j],poNeighbour->armfx[iArmID]);
                        //printf("yforce ws:
%e,%e\n",poNode->armfy[j],poNeighbour->armfy[iArmID]);
                        //printf("zforce ws:
%e,%e\n",poNode->armfz[j],poNeighbour->armfz[iArmID]);

                        //printf("xforce : %e,%e\n",poNode->sigbRem[3*j +
0],poNeighbour->sigbRem[3*iArmID + 0]);
                        //printf("yforce : %e,%e\n",poNode->sigbRem[3*j +
1],poNeighbour->sigbRem[3*iArmID + 1]);
                        //printf("zforce : %e,%e\n",poNode->sigbRem[3*j +
2],poNeighbour->sigbRem[3*iArmID + 2]);
                        // get the component of this force on the slip plane
normal, this is the effective
                        // force component that will be resisted by the APB
force oArmForce = oArmForce - oSlipPlaneNormal*(oArmForce*oSlipPlaneNormal);
                        // get the component of the resolved force in the
direction normal to the segment's
                        // direction, this is the effective	force component
that will be resisted by the APB force oArmForce = oArmForce -
oLineDirection*(oArmForce*oLineDirection);
                        // get the projected arm force direction
                        oArmForceDirection = oArmForce;
                        oArmForceDirection.Normalize();

                        // The magnitude of the APB friction force (positive or
negative) is f_i = tau b_j l s_i s_j
                        // where tau is the APB stress (gamma/b), b_j are the
Burgers vector components, l is
                        // the segment length and s_i are the components of the
segment displacement normal
                        // to the dislocation segment.
                        //dAPBForceMagnitude =
fabs(dLength*dAPBShearStress*(oBurgersVector*oArmForceDirection));
                        dAPBForceMagnitude = fabs(dLength*dAPBShearStress);

            printf("check : %e : %e :
%e\n",dAPBShearStress,dAPBForceMagnitude,oArmForce.Length());
                        if(m_bPairShearing)
                        {
                                // if the chain ID is odd, then it is a leading
dislocation, otherwise, it is a
                                // trailing one, the assumption here is that the
leading dislocation will always
                                // move far away from the trailing, shearing
more areas in the process, the trailing
                                // will always tend to move further from the
leading, but the APB force will keep it
                                // local. In other words, the leading edge
always works so as to increase the APB region
                                // and the trailing edge works so as to always
decrease it if((poNode->piChainID[j]%2) == 0)
                                {
                                        bIsAPBIncreasing = false;
                                }
                                else
                                {
                                        bIsAPBIncreasing = true;
                                }
                        }
                        else
                        {
                                // the total force can work on either increasing
or decreasing the total sheared area
                                // in the case of decreasing it, then add the
APB force to the total force, in case
                                // of decreasing it, if the APB force is greater
than the total force, zero out the
                                // force, if it is less, subtract it from the
total force
                                // get the vector pointing to the cell's center
from the segment's midpoint and normalize it
                                oCenter.SetByPoints(oMidPoint,oCenterPoint);
                                oCenter.Normalize();
                                // see if the total force is in the direction of
the center or opposite to it, assume
                                // unsheared cell first
                                if(oArmForceDirection*oCenter > 0.0)
bIsAPBIncreasing = true; else
bIsAPBIncreasing = false;
                                // flip the case if the cell is sheared
                                if(poAPBCell->IsSheared())
bIsAPBIncreasing = !bIsAPBIncreasing;
                        }

                        // if the APB area is increasing, reduce the total
force, otherwise, add the APB to the total force if(bIsAPBIncreasing)
                        {
                                if(poHome->param->ShearablePrecipitates == 0)
                                {
                                        // the dislocation wants to shear
through an unshearable precipitate, stop it poNode->armfx[j] = 0.0;
                                        poNode->armfy[j] = 0.0;
                                        poNode->armfz[j] = 0.0;
                                        continue;
                                }
                                else
                                {
                                        // precipitates can be sheared
                                        if(dAPBForceMagnitude >
oArmForce.Length())
                                        {
                                                poNode->armfx[j] = 0.0;
                                                poNode->armfy[j] = 0.0;
                                                poNode->armfz[j] = 0.0;
                                                continue;
                                        }
                                        // get the APB force vector and subtract
it from the total force oAPBForce = oArmForceDirection*dAPBForceMagnitude;
                                        poNode->armfx[j] -= oAPBForce.GetX();
                                        poNode->armfy[j] -= oAPBForce.GetY();
                                        poNode->armfz[j] -= oAPBForce.GetZ();
                                }
                        }
                        else
                        {
                                if(poHome->param->ShearablePrecipitates == 0)
                                {
                                        // precipitates cannot be sheared,
however, the dislocation is unshearing them, so
                                        // let it be and don't change the arm
forces
                                }
                                else
                                {
                                        // get the APB force vector and add it
to the total force oAPBForce = oArmForceDirection*dAPBForceMagnitude;
                                        poNode->armfx[j] += oAPBForce.GetX();
                                        poNode->armfy[j] += oAPBForce.GetY();
                                        poNode->armfz[j] += oAPBForce.GetZ();
                                }
                        }
                }
        }
        // update the nodal forces
        for(i = 0 ; i < poHome->newNodeKeyPtr ; i++)
        {
                poNode = poHome->nodeKeys[i];
                if(poNode == NULL)				continue;
                poNode->fX = 0.0;
                poNode->fY = 0.0;
                poNode->fZ = 0.0;
                for(j = 0 ; j < poNode->numNbrs ; j++)
                {
                        poNode->fX += poNode->armfx[j];
                        poNode->fY += poNode->armfy[j];
                        poNode->fZ += poNode->armfz[j];
                }
        }
}*/
void ParadisPrecipitateServer::UpdateAPB(Home_t *poHome) {
  if (!m_bShearablePrecipitates)
    return;
  // this method updates the APB, hence, it scans all the dislocation motion and
  // does 2 things
  // 1. if the dislocation moving has an active APB, it updates it
  // 2. if the corresponding APB is not active, it adds it
  // it does so by driving lower level methods
  if (m_oStructure.IsEmpty())
    return;
  if (m_bPairShearing)
    return;

  unsigned int iSize = (unsigned int)poHome->newNodeKeyPtr;
  unsigned int i = 0;
  unsigned int j = 0;
  Node_t *poNode = NULL;
  Node_t *poNeighbour = NULL;
  Polyhedron *poPrecipitate = NULL;
  unsigned int iSlipSystem = 0;
  Vector oSlipPlaneNormal;
  Vector oBurgersVector;
  APB *poAPB = NULL;
  double dPlaneSpacing = 0.0;
  double dTolerance = 1.0E-3;
  APBCreationEvent *poAPBCreationEvent = NULL;
  for (i = 0; i < iSize; i++) {
    poNode = poHome->nodeKeys[i];
    if (poNode == NULL) {
      continue;
    }
    // if the node has just crossed the PBC boundary, skip the check
    if (DidNodeCrossPBC(poHome->param, poNode))
      continue;
    // look for all of its neighbours
    for (j = 0; j < poNode->numNbrs; j++) {
      // avoid handling the same segment twice
      if (fabs(poNode->burgX[j]) < dTolerance) {
        if (fabs(poNode->burgY[j]) < dTolerance) {
          if (poNode->burgZ[j] <= dTolerance) {
            continue;
          }
        } else if (poNode->burgY[j] < 0.0) {
          continue;
        }
      } else if (poNode->burgX[j] < 0.0) {
        continue;
      }

      // non <110> Burgers dislocations don't shear APBs
      oBurgersVector.Set(poNode->burgX[j], poNode->burgY[j], poNode->burgZ[j]);
      if (!ParadisCrossSlipServer::Is110Vector(oBurgersVector))
        continue;
      // non <111> plane dislocations don't shear APBs
      oSlipPlaneNormal.Set(poNode->nx[j], poNode->ny[j], poNode->nz[j]);
      oSlipPlaneNormal.Normalize();
      if (!ParadisCrossSlipServer::Is111Vector(oSlipPlaneNormal))
        continue;
      // super dislocations don't shear APBs
      if (IsSuperDislocation(poNode, j))
        continue;

      // now we have a dislocation that might shear a precipitate
      poNeighbour = GetNeighborNode(poHome, poNode, j);
      // if the nodes spanned, or span the PBC boundary, skip the check
      if (DidNodesSpanPBC(poHome->param, poNode, poNeighbour))
        continue;
      if (DoNodesSpanPBC(poHome->param, poNode, poNeighbour))
        continue;

      if (!IsSegmentShearCapable(poNode, poNeighbour, poPrecipitate))
        continue;
      // now we have a plane sweep that is at least partially inside one
      // precipitate
      iSlipSystem = IndentifyFCCSlipSystem(oSlipPlaneNormal, oBurgersVector);
      // if the slip system is not an FCC slip system, skip it
      if ((iSlipSystem > 12) || (iSlipSystem < 1))
        continue;
      // this step is redundant for the most part, however, in the case of
      // inverted slip plane normal, we need to get the standard normal vector
      // so as to properly compute the plane spacing which will be used later in
      // the comparison with the already existing apbs
      GetFCCSlipSystemVectors(iSlipSystem, oSlipPlaneNormal, oBurgersVector);
      dPlaneSpacing =
          oSlipPlaneNormal * Vector(poNode->x, poNode->y, poNode->z);
      // loop over the APB list and see if the current dislocation belongs to
      // any of their planes
      poAPB = GetAPB(iSlipSystem, dPlaneSpacing);
      if (poAPB == NULL) {
        // no active APBs were found, create a new one
        poAPB = AddAPB(iSlipSystem, dPlaneSpacing);
        // if the APB cannot be added, skip this node
        if (poAPB == NULL)
          continue;
        // create an APB creation event and add it to the list of events
        poAPBCreationEvent = new APBCreationEvent;
        poAPBCreationEvent->SetAPB(poAPB);
        m_poAPBEventCommunicator->AddEvent(poAPBCreationEvent);
      }
      // in either case, update the APB
      poAPB->Update(poNode, poNeighbour, m_poAPBEventCommunicator);
    }
  }
  // when done, communicate all updates in the APBs, notice that this updates
  // all the triangle sweeps from the last step as well
  m_poAPBEventCommunicator->Communicate(this);
  // when done updating and communicating everything, clean all the persisting
  // singletons and communicate the cleaning events
  // unsigned int iFCCSlipPlanesCount = (unsigned int)m_vlpAPB.size();
  // list< APB* >::const_iterator liAPB;
  // for(i = 0 ; i < iFCCSlipPlanesCount ; i++)
  //{
  // for(liAPB = m_vlpAPB[i].begin() ; liAPB != m_vlpAPB[i].end() ; liAPB++)
  //{
  //(*liAPB)->CleanSingletons(poHome->cycle,poHome->myDomain,m_poAPBEventCommunicator);
  //}
  //}
  // m_poAPBEventCommunicator->Communicate(this);
}
ParadisPrecipitateServer::ParadisPrecipitateServer() { Initialize(); }
void ParadisPrecipitateServer::Initialize() {
  m_oStructure.Reset();
  m_vlpAPB.clear();
  m_poAPBEventCommunicator = APBEventCommunicator::GetCommunicator();
  m_iDomainID = 0;
  m_dXMin = 0.0;
  m_dXMax = 0.0;
  m_dYMin = 0.0;
  m_dYMax = 0.0;
  m_bShearablePrecipitates = true;
  m_bPairShearing = false;
  m_bUseNucleationMode = false;
  m_dNucleationRate = 0.0;
  m_dNucleationRadius = 0.0;
  m_iNucleationMaxCount = -1;
  m_iTotalNucleatedCount = 0;
}
void ParadisPrecipitateServer::RetractNode(
    Home_t *poHome, Node_t *poNode, Polyhedron *poContainingPrecipitate) const {
  if (m_oStructure.IsEmpty())
    return;
  Point oNodePoint(poNode->x, poNode->y, poNode->z);
  // get the containing precipitates, if any
  Polyhedron *poPrecipitate = NULL;
  if (poContainingPrecipitate == NULL) {
    poPrecipitate = m_oStructure.GetContainingPrecipitate(&oNodePoint);
  } else {
    poPrecipitate = poContainingPrecipitate;
  }

  if (poPrecipitate == NULL) {
    return;
  }
  // push the nodes inside the precipitate back outside
  unsigned int iConstraintType = 0;
  Vector oConstraintVector;
  ParadisSurface::GetNodeDynamicConstraint(poNode, iConstraintType,
                                           oConstraintVector);
  Point oNearestPoint = oNodePoint;
  double dDistance = 0.0;
  if (iConstraintType == 0) {
    oNearestPoint = poPrecipitate->GetNearestPoint(oNodePoint, dDistance);
  } else if (iConstraintType == 1) {
    Plane oPlane(oConstraintVector, oNodePoint);
    poPrecipitate->GetNearestPointOnPlane(oNodePoint, oPlane, oNearestPoint,
                                          dDistance);
  } else if (iConstraintType == 2) {
    Line oLine(oConstraintVector, oNodePoint);
    poPrecipitate->GetNearestPointOnLine(oNodePoint, oLine, oNearestPoint,
                                         dDistance);
  }

  double daPos[3] = {0.0, 0.0, 0.0};
  daPos[0] = oNearestPoint.GetX();
  daPos[1] = oNearestPoint.GetY();
  daPos[2] = oNearestPoint.GetZ();
  RepositionNode(poHome, daPos, &(poNode->myTag), 1);
}
unsigned int
ParadisPrecipitateServer::IndentifyFCCSlipSystem(Vector oNormalVector,
                                                 Vector oBurgersVector) {
  oNormalVector.Normalize();
  oBurgersVector.Normalize();
  double dTolerance = 1E-3;
  Vector oReferenceNormal(0.0, 0.0, 0.0);
  Vector oReferenceBurgers(0.0, 0.0, 0.0);

  // slip plane 1,1,1
  oReferenceNormal.Set(1.0, 1.0, 1.0);
  oReferenceNormal.Normalize();

  oReferenceBurgers.Set(1.0, -1.0, 0.0);
  oReferenceBurgers.Normalize();
  if (oReferenceNormal.IsSame(oNormalVector, dTolerance) ||
      oReferenceNormal.IsOpposite(oNormalVector, dTolerance)) {
    if (oReferenceBurgers.IsSame(oBurgersVector, dTolerance) ||
        oReferenceBurgers.IsOpposite(oBurgersVector, dTolerance)) {
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
Vector
ParadisPrecipitateServer::GetFCCSlipPlaneNormal(const unsigned int &iPlaneID) {
  Vector oNormalVector;
  if ((iPlaneID < 1) || (iPlaneID > 4))
    return oNormalVector;
  oNormalVector.Set(1.0, 1.0, 1.0);

  if (iPlaneID == 2)
    oNormalVector.SetX(-1.0);
  else if (iPlaneID == 3)
    oNormalVector.SetY(-1.0);
  else if (iPlaneID == 4)
    oNormalVector.SetZ(-1.0);
  oNormalVector.Normalize();
  return oNormalVector;
}
void ParadisPrecipitateServer::GetFCCSlipSystemVectors(
    const unsigned int &iSlipSystem, Vector &oNormalVector,
    Vector &oBurgersVector) {
  if ((iSlipSystem < 1) || (iSlipSystem > 12)) {
    return;
  }
  if ((iSlipSystem == 1) || (iSlipSystem == 2) || (iSlipSystem == 3)) {
    oNormalVector.Set(1.0, 1.0, 1.0);
    if (iSlipSystem == 1)
      oBurgersVector.Set(1.0, -1.0, 0.0);
    if (iSlipSystem == 2)
      oBurgersVector.Set(1.0, 0.0, -1.0);
    if (iSlipSystem == 3)
      oBurgersVector.Set(0.0, 1.0, -1.0);
  }
  if ((iSlipSystem == 4) || (iSlipSystem == 5) || (iSlipSystem == 6)) {
    oNormalVector.Set(-1.0, 1.0, 1.0);
    if (iSlipSystem == 4)
      oBurgersVector.Set(1.0, 1.0, 0.0);
    if (iSlipSystem == 5)
      oBurgersVector.Set(1.0, 0.0, 1.0);
    if (iSlipSystem == 6)
      oBurgersVector.Set(0.0, 1.0, -1.0);
  }
  if ((iSlipSystem == 7) || (iSlipSystem == 8) || (iSlipSystem == 9)) {
    oNormalVector.Set(1.0, -1.0, 1.0);
    if (iSlipSystem == 7)
      oBurgersVector.Set(1.0, 1.0, 0.0);
    if (iSlipSystem == 8)
      oBurgersVector.Set(1.0, 0.0, -1.0);
    if (iSlipSystem == 9)
      oBurgersVector.Set(0.0, 1.0, 1.0);
  }
  if ((iSlipSystem == 10) || (iSlipSystem == 11) || (iSlipSystem == 12)) {
    oNormalVector.Set(1.0, 1.0, -1.0);
    if (iSlipSystem == 10)
      oBurgersVector.Set(1.0, -1.0, 0.0);
    if (iSlipSystem == 11)
      oBurgersVector.Set(1.0, 0.0, 1.0);
    if (iSlipSystem == 12)
      oBurgersVector.Set(0.0, 1.0, 1.0);
  }
  oNormalVector.Normalize();
  oBurgersVector.Normalize();
}
bool ParadisPrecipitateServer::IsNodeShearCapable(
    Node_t *poNode, Polyhedron *&poPrecipitate) const {
  if (m_oStructure.IsEmpty())
    return false;
  poPrecipitate = NULL;
  if (poNode == NULL) {
    return false;
  }
  // if the node is not inside any precipitate at the moment and hasn't been
  // inside any precipitates in the last time step, skip it
  Point oPosition(poNode->x, poNode->y, poNode->z);
  Point oOldPosition(poNode->oldx, poNode->oldy, poNode->oldz);
  poPrecipitate = m_oStructure.GetContainingPrecipitate(&oPosition);
  Polyhedron *poOldPrecipitate =
      m_oStructure.GetContainingPrecipitate(&oOldPosition);
  // if the node was out and is still out, skip it
  if (poPrecipitate == NULL) {
    if (poOldPrecipitate == NULL) {
      return false;
    } else {
      poPrecipitate = poOldPrecipitate;
    }
  } else {
    if (poOldPrecipitate == NULL) {
      poOldPrecipitate = poPrecipitate;
    } else if (poPrecipitate != poOldPrecipitate) {
      // if the node switched precipitates between steps (an extremely rare
      // situation), it is too difficult to handle it, just skip it for now
      return false;
    }
  }
  return true;
}
bool ParadisPrecipitateServer::IsSegmentShearCapable(
    Node_t *poNode1, Node_t *poNode2, Polyhedron *&poPrecipitate) const {
  if (m_oStructure.IsEmpty())
    return false;
  poPrecipitate = NULL;
  bool bNode1Capable = IsNodeShearCapable(poNode1, poPrecipitate);
  Polyhedron *poTempPrecipitate = NULL;
  bool bNode2Capable = IsNodeShearCapable(poNode2, poTempPrecipitate);
  if (bNode1Capable || bNode2Capable) {
    // at least one of them is shear capable
    // if only one of them is shear capable, use its precipitate
    if (!bNode2Capable) {
      return true;
    }
    if (!bNode1Capable) {
      poPrecipitate = poTempPrecipitate;
      return true;
    }
    // if the two of them are shear capable, make sure that they are in the same
    // precipitate
    if (poPrecipitate == poTempPrecipitate) {
      // the two of them are shear capable and they are in the same precipitate
      return true;
    }
  }
  return false;
}
void ParadisPrecipitateServer::WriteAPB(const string &sFileName) const {
  if (!m_bShearablePrecipitates)
    return;
  if (m_oStructure.IsEmpty())
    return;
  if (m_bPairShearing)
    return;
  // since the cell shearing is being constantly communicated, only processor 0
  // writes the apb data
  MPI_Barrier(MPI_COMM_WORLD);
  if (m_iDomainID == 0) {
    FILE *fpFile = fopen(sFileName.c_str(), "w");
    // first, write the total number of APB planes in the system
    unsigned int iAPBCount = 0;
    unsigned int i = 0;
    unsigned int iFCCSlipPlanesCount = (unsigned int)m_vlpAPB.size();
    for (i = 0; i < iFCCSlipPlanesCount; i++) {
      iAPBCount = iAPBCount + (unsigned int)m_vlpAPB[i].size();
    }
    fwrite(&iAPBCount, sizeof(unsigned int), 1, fpFile);
    // then write the resolution used in the APB planes
    fwrite(&m_iAPBPointResolution, sizeof(unsigned int), 1, fpFile);
    // then write the x-y limits of the APBs
    fwrite(&m_dXMin, sizeof(double), 1, fpFile);
    fwrite(&m_dXMax, sizeof(double), 1, fpFile);
    fwrite(&m_dYMin, sizeof(double), 1, fpFile);
    fwrite(&m_dYMax, sizeof(double), 1, fpFile);
    // finally, loop over the APB planes and write their cell shearing data
    list<APB *>::const_iterator liAPB;
    for (i = 0; i < iFCCSlipPlanesCount; i++) {
      for (liAPB = m_vlpAPB[i].begin(); liAPB != m_vlpAPB[i].end(); liAPB++) {
        (*liAPB)->Write(fpFile);
      }
    }
    fclose(fpFile);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}
bool ParadisPrecipitateServer::IsNodeInsideAPrecipitate(
    Node_t *poNode, Polyhedron *&poContainingPrecipitate) const {
  Point oPosition(poNode->x, poNode->y, poNode->z);
  return IsPointInsideAPrecipitate(&oPosition, poContainingPrecipitate);
}
bool ParadisPrecipitateServer::IsPointInsideAPrecipitate(
    Point *poPoint, Polyhedron *&poContainingPrecipitate) const {
  if (m_oStructure.IsEmpty())
    return false;
  poContainingPrecipitate = m_oStructure.GetContainingPrecipitate(poPoint);
  if (poContainingPrecipitate == NULL)
    return false;
  return true;
}
void ParadisPrecipitateServer::UpdateTriangleSweep(
    Home_t *poHome, const Point &oPoint1, const Point &oPoint2,
    const Point &oPoint3, const Vector &oBurgersVector) {
  if (!m_bShearablePrecipitates)
    return;
  if (m_oStructure.IsEmpty())
    return;
  if (m_bPairShearing)
    return;
  // this function takes in a triangle and a precipitate where the points of the
  // triangle lie inside, gets the proper APB plane where this triangle is on
  // and checks to see whether the triangle sweeps any of the APB plane points,
  // if so, it flips them
  Vector oNormalVector = Vector(oPoint1, oPoint2) ^ Vector(oPoint1, oPoint3);
  oNormalVector.Normalize();
  unsigned int iSlipSystem =
      IndentifyFCCSlipSystem(oNormalVector, oBurgersVector);
  // if the slip system is not an FCC slip system, skip it
  if ((iSlipSystem > 12) || (iSlipSystem < 1))
    return;
  Vector oTempBurgersVector;
  GetFCCSlipSystemVectors(iSlipSystem, oNormalVector, oTempBurgersVector);
  double dPlaneSpacing = oNormalVector * ((Vector)oPoint1);
  APB *poAPB = GetAPB(iSlipSystem, dPlaneSpacing);
  if (poAPB == NULL)
    return;
  poAPB->UpdateTriangle(oPoint1, oPoint2, oPoint3, oBurgersVector,
                        m_poAPBEventCommunicator);
}
void ParadisPrecipitateServer::ReadAPBs(const string &sAPBFileName) {
  // if the precipitates are not shearable, return
  if (!m_bShearablePrecipitates)
    return;
  if (m_bPairShearing)
    return;
  // if the file name is empty, return
  if (sAPBFileName.empty())
    return;
  FILE *fpFile = fopen(sAPBFileName.c_str(), "r");
  unsigned int iAPBCount = 0;
  // read the APB count
  fread(&iAPBCount, sizeof(unsigned int), 1, fpFile);
  // read the resolution used in the APB planes
  fread(&m_iAPBPointResolution, sizeof(unsigned int), 1, fpFile);
  // then read the x-y limits of the APBs
  fread(&m_dXMin, sizeof(double), 1, fpFile);
  fread(&m_dXMax, sizeof(double), 1, fpFile);
  fread(&m_dYMin, sizeof(double), 1, fpFile);
  fread(&m_dYMax, sizeof(double), 1, fpFile);
  // finally, loop over the APBs and read them
  unsigned int i = 0;
  APB *poAPB = NULL;
  unsigned int iSlipPlaneID = 0;
  for (i = 0; i < iAPBCount; i++) {
    poAPB = new APB;
    poAPB->SetXMin(m_dXMin);
    poAPB->SetXMax(m_dXMax);
    poAPB->SetYMin(m_dYMin);
    poAPB->SetYMax(m_dYMax);
    poAPB->Read(fpFile, m_iAPBPointResolution, &m_oStructure, iSlipPlaneID);
    m_vlpAPB[iSlipPlaneID - 1].push_back(poAPB);
  }
  fclose(fpFile);
}
bool ParadisPrecipitateServer::IsSuperDislocation(
    Node_t *poNode, const unsigned int &iNeighbourIndex) {
  Vector oBurgersVector(poNode->burgX[iNeighbourIndex],
                        poNode->burgY[iNeighbourIndex],
                        poNode->burgZ[iNeighbourIndex]);
  if (ParadisCrossSlipServer::Is100Vector(oBurgersVector))
    return true;
  if (ParadisCrossSlipServer::Is110Vector(oBurgersVector)) {
    double dMagnitude = oBurgersVector.Length();
    // if the magnitude is even, then it is a super partial
    if (((int)floor(dMagnitude + 0.5) % 2) == 0)
      return true;
  }
  return false;
}
void ParadisPrecipitateServer::UpdateAPBCellOwnership(Home_t *poHome) {
  if (m_oStructure.IsEmpty())
    return;
  unsigned int iSlipSystemsCount = m_vlpAPB.size();
  unsigned int i = 0;
  list<APB *>::iterator liAPB;

  for (i = 0; i < iSlipSystemsCount; i++) {
    for (liAPB = m_vlpAPB[i].begin(); liAPB != m_vlpAPB[i].end(); liAPB++) {
      (*liAPB)->UpdateAPBCellOwnership(poHome);
    }
  }
}
APB *ParadisPrecipitateServer::GetAPBByPlane(
    const unsigned int &iPlaneID, const double &dPlaneSpacing) const {
  APB *poAPB = NULL;
  if (m_vlpAPB.size() > 0) {
    const list<APB *> *plpAPBs = &m_vlpAPB[iPlaneID - 1];
    list<APB *>::const_iterator liAPB;
    double dTolerance = 1.0;
    for (liAPB = plpAPBs->begin(); liAPB != plpAPBs->end(); liAPB++) {
      if (fabs((*liAPB)->GetPlaneSpacing() - dPlaneSpacing) < dTolerance) {
        poAPB = (*liAPB);
        break;
      }
    }
  }
  return poAPB;
}
APB *ParadisPrecipitateServer::GetAPB(const unsigned int &iSlipSystemID,
                                      const double &dPlaneSpacing) const {
  return GetAPBByPlane(GetFCCSlipPlaneID(iSlipSystemID), dPlaneSpacing);
}
APB *ParadisPrecipitateServer::GetAPB(Node_t *poNode,
                                      const Vector &oSlipPlaneNormal) const {
  Vector oWorkingNormal = oSlipPlaneNormal;
  oWorkingNormal.Normalize();
  unsigned int iSlipPlaneID = APB::GetFCCSlipPlaneID(oWorkingNormal);
  if ((iSlipPlaneID < 1) || (iSlipPlaneID > 4))
    return NULL;
  // this step is redundant for the most part, however, in the case of inverted
  // slip plane normal, we need to get the standard normal vector so as to
  // properly compute the plane spacing which will be used later in the
  // comparison with the already existing apbs
  oWorkingNormal = APB::GetFCCSlipPlaneNormalFromID(iSlipPlaneID);
  double dPlaneSpacing =
      oWorkingNormal * Vector(poNode->x, poNode->y, poNode->z);
  return GetAPBByPlane(iSlipPlaneID, dPlaneSpacing);
}
APB *ParadisPrecipitateServer::GetAPB(Node_t *poNode1, Node_t *poNode2) const {
  int iArmID = GetArmID(poNode1, poNode2);
  Vector oSlipPlaneNormal(poNode1->nx[iArmID], poNode1->ny[iArmID],
                          poNode1->nz[iArmID]);
  oSlipPlaneNormal.Normalize();
  Vector oBurgersVector(poNode1->burgX[iArmID], poNode1->burgY[iArmID],
                        poNode1->burgZ[iArmID]);
  unsigned int iSlipSystem =
      IndentifyFCCSlipSystem(oSlipPlaneNormal, oBurgersVector);
  // if the slip system is not an FCC slip system, skip it
  if ((iSlipSystem > 12) || (iSlipSystem < 1))
    return NULL;

  // this step is redundant for the most part, however, in the case of inverted
  // slip plane normal, we need to get the standard normal vector so as to
  // properly compute the plane spacing which will be used later in the
  // comparison with the already existing apbs
  GetFCCSlipSystemVectors(iSlipSystem, oSlipPlaneNormal, oBurgersVector);
  double dPlaneSpacing =
      oSlipPlaneNormal * Vector(poNode1->x, poNode1->y, poNode1->z);
  return GetAPB(iSlipSystem, dPlaneSpacing);
}
APB *ParadisPrecipitateServer::AddAPBByPlane(const unsigned int &iPlaneID,
                                             const double &dPlaneSpacing) {
  // check plane id
  if ((iPlaneID < 1) || (iPlaneID > 4))
    return NULL;
  // APB doesn't exist, create and add it
  APB *poAPB = new APB;
  poAPB->SetSlipPlaneNormal(GetFCCSlipPlaneNormal(iPlaneID));
  poAPB->SetPlaneSpacing(dPlaneSpacing);
  poAPB->SetResolution(m_iAPBPointResolution);

  poAPB->SetXMin(m_dXMin);
  poAPB->SetXMax(m_dXMax);
  poAPB->SetYMin(m_dYMin);
  poAPB->SetYMax(m_dYMax);

  poAPB->GeneratePoints(&m_oStructure);
  m_vlpAPB[iPlaneID - 1].push_back(poAPB);
  return poAPB;
}
APB *ParadisPrecipitateServer::AddAPB(const unsigned int &iSlipSystemID,
                                      const double &dPlaneSpacing) {
  return AddAPBByPlane(GetFCCSlipPlaneID(iSlipSystemID), dPlaneSpacing);
}
void ParadisPrecipitateServer::Report(const unsigned int &iProcessorID) const {
  unsigned int iSlipPlanesCount = m_vlpAPB.size();
  unsigned int i = 0;
  list<APB *>::const_iterator liAPB;
  for (i = 0; i < iSlipPlanesCount; i++) {
    for (liAPB = m_vlpAPB[i].begin(); liAPB != m_vlpAPB[i].end(); liAPB++) {
      (*liAPB)->Report(iProcessorID);
    }
  }
}
unsigned int ParadisPrecipitateServer::GetPrecipitatesCount() const {
  return m_oStructure.GetPrecipitatesCount();
}
unsigned int ParadisPrecipitateServer::GetTotalNucleatedCount() const {
  return m_iTotalNucleatedCount;
}
bool ParadisPrecipitateServer::WritePrecipitateStructureSnapshot(
    const string &sFileName) const {
  if (m_iDomainID != 0) {
    return true;
  }
  m_oStructure.WritePrecipitateStructure(sFileName);
  return true;
}
unsigned int
ParadisPrecipitateServer::GetFCCSlipPlaneID(const unsigned int &iSlipSystemID) {
  unsigned int iSlipPlaneID = 0;
  if ((iSlipSystemID == 1) || (iSlipSystemID == 2) || (iSlipSystemID == 3))
    iSlipPlaneID = 1;
  if ((iSlipSystemID == 4) || (iSlipSystemID == 5) || (iSlipSystemID == 6))
    iSlipPlaneID = 2;
  if ((iSlipSystemID == 7) || (iSlipSystemID == 8) || (iSlipSystemID == 9))
    iSlipPlaneID = 3;
  if ((iSlipSystemID == 10) || (iSlipSystemID == 11) || (iSlipSystemID == 12))
    iSlipPlaneID = 4;
  return iSlipPlaneID;
}
