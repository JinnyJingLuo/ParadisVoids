#include "ParadisBlockSurface.h"
#include "Comm.h"
#include "Home.h"

ParadisBlockSurface::ParadisBlockSurface() { Initialize(); }
ParadisBlockSurface::ParadisBlockSurface(const ParadisBlockSurface &oSurface) {
  *this = oSurface;
}
ParadisBlockSurface::~ParadisBlockSurface() { Reset(); }
ParadisBlockSurface &ParadisBlockSurface::
operator=(const ParadisBlockSurface &oSurface) {
  Reset();
  ParadisSurface::operator=(oSurface);
  return *this;
}
void ParadisBlockSurface::Reset() {
  ParadisSurface::Reset();
  Initialize();
}
void ParadisBlockSurface::Set(Home_t *poHome) {
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
  m_dVolume = (m_dXMax - m_dXMin) * (m_dYMax - m_dYMin) * (m_dZMax - m_dZMin);
  GenerateTriangulations();
  poHome->param->simVol = m_dVolume;
  poHome->param->burgVolFactor =
      1.0 /
      (poHome->param->burgMag * poHome->param->burgMag * poHome->param->simVol);
}
void ParadisBlockSurface::Initialize() { ParadisSurface::Initialize(); }
bool ParadisBlockSurface::IsPointInside(const Point &oPoint) const {
  double dX = oPoint.GetX();
  if ((dX >= m_dXMin) && (dX <= m_dXMax)) {
    double dY = oPoint.GetY();
    if ((dY >= m_dYMin) && (dY <= m_dYMax)) {
      double dZ = oPoint.GetZ();
      if ((dZ >= m_dZMin) && (dZ <= m_dZMax)) {
        return true;
      }
    }
  }
  return false;
}
bool ParadisBlockSurface::IsPointOnSurface(const Point &oPoint) const {
  double dX = oPoint.GetX();
  double dY = oPoint.GetY();
  double dZ = oPoint.GetZ();
  if ((dX >= m_dXMin) && (dX <= m_dXMax)) {
    if ((dY >= m_dYMin) && (dY <= m_dYMax)) {
      if ((fabs(dZ - m_dZMin) < m_dTolerance) ||
          (fabs(dZ - m_dZMax) < m_dTolerance)) {
        return true;
      }
    }
  }

  if ((dX >= m_dXMin) && (dX <= m_dXMax)) {
    if ((dZ >= m_dZMin) && (dZ <= m_dZMax)) {
      if ((fabs(dY - m_dYMin) < m_dTolerance) ||
          (fabs(dY - m_dYMax) < m_dTolerance)) {
        return true;
      }
    }
  }

  if ((dY >= m_dYMin) && (dY <= m_dYMax)) {
    if ((dZ >= m_dZMin) && (dZ <= m_dZMax)) {
      if ((fabs(dX - m_dXMin) < m_dTolerance) ||
          (fabs(dX - m_dXMax) < m_dTolerance)) {
        return true;
      }
    }
  }
  return false;
}
void ParadisBlockSurface::GenerateTriangulations() {
  GenericNode *poPointNNN = new GenericNode(m_dXMin, m_dYMin, m_dZMin);
  GenericNode *poPointPNN = new GenericNode(m_dXMax, m_dYMin, m_dZMin);
  GenericNode *poPointPPN = new GenericNode(m_dXMax, m_dYMax, m_dZMin);
  GenericNode *poPointNPN = new GenericNode(m_dXMin, m_dYMax, m_dZMin);
  GenericNode *poPointNNP = new GenericNode(m_dXMin, m_dYMin, m_dZMax);
  GenericNode *poPointPNP = new GenericNode(m_dXMax, m_dYMin, m_dZMax);
  GenericNode *poPointPPP = new GenericNode(m_dXMax, m_dYMax, m_dZMax);
  GenericNode *poPointNPP = new GenericNode(m_dXMin, m_dYMax, m_dZMax);
  m_lpoPoints.push_back(poPointNNN);
  m_lpoPoints.push_back(poPointPNN);
  m_lpoPoints.push_back(poPointPPN);
  m_lpoPoints.push_back(poPointNPN);
  m_lpoPoints.push_back(poPointNNP);
  m_lpoPoints.push_back(poPointPNP);
  m_lpoPoints.push_back(poPointPPP);
  m_lpoPoints.push_back(poPointNPP);

  m_lpoTriangles.push_back(new TriPatch(poPointNNN, poPointPNN, poPointNNP));
  m_lpoTriangles.push_back(new TriPatch(poPointPNN, poPointPNP, poPointNNP));

  m_lpoTriangles.push_back(new TriPatch(poPointNPN, poPointNPP, poPointPPP));
  m_lpoTriangles.push_back(new TriPatch(poPointNPN, poPointPPP, poPointPPN));

  m_lpoTriangles.push_back(new TriPatch(poPointNPN, poPointNNN, poPointNNP));
  m_lpoTriangles.push_back(new TriPatch(poPointNPN, poPointNNP, poPointNPP));

  m_lpoTriangles.push_back(new TriPatch(poPointPNN, poPointPPN, poPointPPP));
  m_lpoTriangles.push_back(new TriPatch(poPointPNN, poPointPPP, poPointPNP));

  m_lpoTriangles.push_back(new TriPatch(poPointNNN, poPointPPN, poPointPNN));
  m_lpoTriangles.push_back(new TriPatch(poPointNNN, poPointNPN, poPointPPN));

  m_lpoTriangles.push_back(new TriPatch(poPointPNP, poPointNPP, poPointNNP));
  m_lpoTriangles.push_back(new TriPatch(poPointPNP, poPointPPP, poPointNPP));
}
