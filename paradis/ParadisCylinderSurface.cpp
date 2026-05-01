#include "ParadisCylinderSurface.h"
#include "FEMMesh.h"
#include "Home.h"

using namespace std;

ParadisCylinderSurface::ParadisCylinderSurface() { Initialize(); }
ParadisCylinderSurface::ParadisCylinderSurface(
    const ParadisCylinderSurface &oSurface) {
  *this = oSurface;
}
ParadisCylinderSurface::~ParadisCylinderSurface() { Reset(); }
ParadisCylinderSurface &ParadisCylinderSurface::
operator=(const ParadisCylinderSurface &oSurface) {
  Reset();
  ParadisSurface::operator=(oSurface);
  m_poCylinder = (Cylinder *)(oSurface.m_poCylinder->Clone());
  return *this;
}
void ParadisCylinderSurface::Reset() {
  ParadisSurface::Reset();
  if (m_poCylinder != NULL) {
    delete m_poCylinder;
  }
  m_poCylinder = NULL;
}
void ParadisCylinderSurface::Set(Home_t *poHome) {
  Reset();
  m_dXMin = -0.5 * poHome->param->Dimensions[0];
  m_dXMax = 0.5 * poHome->param->Dimensions[0];
  m_dYMin = -0.5 * poHome->param->Dimensions[1];
  m_dYMax = 0.5 * poHome->param->Dimensions[1];
  m_dZMin = -0.5 * poHome->param->Dimensions[2];
  m_dZMax = 0.5 * poHome->param->Dimensions[2];
  double dToleranceFactor = 1.0E-6;
  double dMinimumDimension = m_dXMax - m_dXMin;
  if (dMinimumDimension > (m_dYMax - m_dYMin)) {
    dMinimumDimension = m_dYMax - m_dZMin;
  }
  if (dMinimumDimension > (m_dZMax - m_dYMin)) {
    dMinimumDimension = m_dZMax - m_dZMin;
  }
  m_dTolerance = dToleranceFactor * dMinimumDimension;
  m_dXMin = m_dXMin - m_dTolerance;
  m_dXMax = m_dXMax + m_dTolerance;
  m_dYMin = m_dYMin - m_dTolerance;
  m_dYMax = m_dYMax + m_dTolerance;
  m_dZMin = m_dZMin - m_dTolerance;
  m_dZMax = m_dZMax + m_dTolerance;
  double dRadius = min((m_dXMax - m_dXMin), (m_dYMax - m_dYMin));
  double dLength = m_dZMax - m_dZMin;
  m_poCylinder = new Cylinder();
  m_poCylinder->SetRadius(dRadius);
  m_poCylinder->SetLength(dLength);
  m_poCylinder->SetResolution(5, 18, 3);
  m_dVolume = m_poCylinder->GetVolume();
  GenerateTriangulations();
  poHome->param->simVol = m_dVolume;
  poHome->param->burgVolFactor =
      1.0 /
      (poHome->param->burgMag * poHome->param->burgMag * poHome->param->simVol);
}
void ParadisCylinderSurface::Initialize() {
  ParadisSurface::Initialize();
  m_poCylinder = NULL;
}
bool ParadisCylinderSurface::IsPointInside(const Point &oPoint) const {
  return m_poCylinder->IsPointInside(oPoint);
}
bool ParadisCylinderSurface::IsPointOnSurface(const Point &oPoint) const {
  return m_poCylinder->IsOnSurface(oPoint);
}
void ParadisCylinderSurface::GenerateTriangulations() {
  // generate an FEM mesh for the cylinder
  vector<FEMNode *> vpoNodes;
  vector<FEMElement *> vpoElements;
  FEMMesh::GenerateMeshFromCylinder(m_poCylinder, vpoNodes, vpoElements,
                                    SolidFEMNode, SolidFEMElement);
  GenerateSurfaceTriangulationFromFEMNodes(vpoNodes, vpoElements);
}
