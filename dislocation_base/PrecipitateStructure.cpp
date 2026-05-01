#include "PrecipitateStructure.h"
#include "Tools.h"
#include "Dodecahedron.h"
#include "math.h"
#include "string"
#include "list"
#include "float.h"

using namespace std;
using namespace SupportSystem;

double PrecipitateStructure::m_dCellSize = 1000.0;
PrecipitateStructure::PrecipitateStructure() { Initialize(); }
PrecipitateStructure::PrecipitateStructure(
    const PrecipitateStructure &oStructure)
{
  *this = oStructure;
}
PrecipitateStructure::~PrecipitateStructure() { Reset(); }
PrecipitateStructure &PrecipitateStructure::
operator=(const PrecipitateStructure &oStructure)
{
  return *this;
}
void PrecipitateStructure::Reset()
{
  unsigned int iSize = (unsigned int)m_vpoPrecipitates.size();
  unsigned int i = 0;
  for (i = 0; i < iSize; i++)
  {
    if (m_vpoPrecipitates[i] != NULL)
    {
      delete m_vpoPrecipitates[i];
    }
  }
  Initialize();
}
bool PrecipitateStructure::Set(const string &sFileName)
{
  Reset();
  if (!ReadPrecipitatesFile(sFileName))
  {
    return false;
  }
  return true;
}
bool PrecipitateStructure::Set(const string &sFileName,
                               const double &dXDimension,
                               const double &dYDimension,
                               const double &dZDimension)
{
  if (!Set(sFileName))
    return false;
  // estimate the number of cells, the cells are preferably cubic with a side
  // length of m_dCellSize
  m_iXCellsCount = (unsigned int)ceil(dXDimension / m_dCellSize);
  m_iYCellsCount = (unsigned int)ceil(dYDimension / m_dCellSize);
  m_iZCellsCount = (unsigned int)ceil(dZDimension / m_dCellSize);
  m_dXShift = -0.5 * dXDimension;
  m_dYShift = -0.5 * dYDimension;
  m_dZShift = -0.5 * dZDimension;
  // compute the intersection of each cell with the precipitate
  unsigned int iCellsCount = m_iXCellsCount * m_iYCellsCount * m_iZCellsCount;
  m_vsIntersectingPrecipitatesIDs.resize(iCellsCount);
  m_viIntersectingPrecipitatesCounts.resize(iCellsCount);
  unsigned int i = 0;
  for (i = 0; i < iCellsCount; i++)
  {
    m_vsIntersectingPrecipitatesIDs[i] = "";
    m_viIntersectingPrecipitatesCounts[i] = 0;
  }
  // ComputeCellPrecipitateIntersection();
  return true;
}
void PrecipitateStructure::Initialize()
{
  m_iXCellsCount = 0;
  m_iYCellsCount = 0;
  m_iZCellsCount = 0;
  m_vpoPrecipitates.clear();
  m_vsIntersectingPrecipitatesIDs.clear();
  m_viIntersectingPrecipitatesCounts.clear();
  m_dXShift = 0.0;
  m_dYShift = 0.0;
  m_dZShift = 0.0;
}
bool PrecipitateStructure::ReadPrecipitatesFile(const string &sFileName)
{
  Reset();
  unsigned int iPrecipitatesCount = 0;
  FILE *fpFile = fopen(sFileName.c_str(), "r");
  if (fpFile == NULL)
  {
    return false;
  }
  string sRead = GetRealString(256, fpFile);
  sscanf(sRead.c_str(), "%d\n", &iPrecipitatesCount);
  unsigned int iPointsCount = 0;
  unsigned int i = 0;
  unsigned int j = 0;
  double dX = 0.0;
  double dY = 0.0;
  double dZ = 0.0;
  list<Point> loPoints;
  Polyhedron *poPrecipitate = NULL;
  m_vpoPrecipitates.resize(iPrecipitatesCount);
  for (i = 0; i < iPrecipitatesCount; i++)
  {
    iPointsCount = 0;
    sRead = GetRealString(256, fpFile);
    sscanf(sRead.c_str(), "%d\n", &iPointsCount);
    loPoints.clear();
    // read the precipitate points
    for (j = 0; j < iPointsCount; j++)
    {
      sRead = GetRealString(256, fpFile);
      sscanf(sRead.c_str(), "%lf\t%lf\t%lf\n", &dX, &dY, &dZ);
      loPoints.push_back(Point(dX, dY, dZ));
    }
    // create the precipitate
    poPrecipitate = new Polyhedron;
    if (!poPrecipitate->CreateAsHullFromPoints(&loPoints))
    {
      loPoints.clear();
      printf("error: precipitate generation failed\n");
      fclose(fpFile);
      return false;
    }
    loPoints.clear();
    poPrecipitate->SetID(i + 1);
    m_vpoPrecipitates[i] = poPrecipitate;
  }
  fclose(fpFile);
  return true;
}
unsigned int
PrecipitateStructure::GetPointCellIndex(const Point *poPoint) const
{
  // the indices are zero based
  unsigned int iXCellIndex =
      (unsigned int)floor((poPoint->GetX() - m_dXShift) / m_dCellSize);
  unsigned int iYCellIndex =
      (unsigned int)floor((poPoint->GetY() - m_dYShift) / m_dCellSize);
  unsigned int iZCellIndex =
      (unsigned int)floor((poPoint->GetZ() - m_dZShift) / m_dCellSize);
  if (iXCellIndex >= m_iXCellsCount)
  {
    iXCellIndex = m_iXCellsCount - 1;
  }
  if (iYCellIndex >= m_iYCellsCount)
  {
    iYCellIndex = m_iYCellsCount - 1;
  }
  if (iZCellIndex >= m_iZCellsCount)
  {
    iZCellIndex = m_iZCellsCount - 1;
  }
  unsigned int iCellIndex = iZCellIndex * (m_iYCellsCount * m_iXCellsCount) +
                            iYCellIndex * m_iXCellsCount + iXCellIndex;
  return iCellIndex;
}
void PrecipitateStructure::ComputeCellPrecipitateIntersection()
{
  unsigned int i = 0;
  unsigned int j = 0;
  unsigned int k = 0;
  unsigned int l = 0;
  unsigned int iPrecipitatesCount = (unsigned int)m_vpoPrecipitates.size();
  double dX = 0.0;
  double dY = 0.0;
  double dZ = 0.0;
  unsigned int iCellIndex = 0;
  unsigned int iCornersCount = 0;
  char cIndexString[16] = "";
  string sIndexString;
  Polyhedron oCellPolyhedron;
  Point oCellNNN;
  Point oCellPNN;
  Point oCellPPN;
  Point oCellNPN;
  Point oCellNNP;
  Point oCellPNP;
  Point oCellPPP;
  Point oCellNPP;
  list<Point> loCellPoints;
  for (l = 0; l < iPrecipitatesCount; l++)
  {
    printf("in precipitate %d of %d\n", l + 1, iPrecipitatesCount);
    sprintf(cIndexString, "%d:", m_vpoPrecipitates[l]->GetID());
    sIndexString = string(cIndexString);
    for (k = 0; k < m_iZCellsCount; k++)
    {
      for (j = 0; j < m_iYCellsCount; j++)
      {
        for (i = 0; i < m_iXCellsCount; i++)
        {
          iCellIndex =
              k * (m_iYCellsCount * m_iXCellsCount) + j * m_iXCellsCount + i;

          loCellPoints.clear();
          loCellPoints.push_back(Point(i * m_dCellSize + m_dXShift,
                                       j * m_dCellSize + m_dYShift,
                                       k * m_dCellSize + m_dZShift));
          loCellPoints.push_back(Point((i + 1) * m_dCellSize + m_dXShift,
                                       j * m_dCellSize + m_dYShift,
                                       k * m_dCellSize + m_dZShift));
          loCellPoints.push_back(Point((i + 1) * m_dCellSize + m_dXShift,
                                       (j + 1) * m_dCellSize + m_dYShift,
                                       k * m_dCellSize + m_dZShift));
          loCellPoints.push_back(Point(i * m_dCellSize + m_dXShift,
                                       (j + 1) * m_dCellSize + m_dYShift,
                                       k * m_dCellSize + m_dZShift));
          loCellPoints.push_back(Point(i * m_dCellSize + m_dXShift,
                                       j * m_dCellSize + m_dYShift,
                                       (k + 1) * m_dCellSize + m_dZShift));
          loCellPoints.push_back(Point((i + 1) * m_dCellSize + m_dXShift,
                                       j * m_dCellSize + m_dYShift,
                                       (k + 1) * m_dCellSize + m_dZShift));
          loCellPoints.push_back(Point((i + 1) * m_dCellSize + m_dXShift,
                                       (j + 1) * m_dCellSize + m_dYShift,
                                       (k + 1) * m_dCellSize + m_dZShift));
          loCellPoints.push_back(Point(i * m_dCellSize + m_dXShift,
                                       (j + 1) * m_dCellSize + m_dYShift,
                                       (k + 1) * m_dCellSize + m_dZShift));
          oCellPolyhedron.CreateAsHullFromPoints(&loCellPoints);

          if (oCellPolyhedron.IsIntersecting(*m_vpoPrecipitates[l]))
          {
            m_vsIntersectingPrecipitatesIDs[iCellIndex] =
                m_vsIntersectingPrecipitatesIDs[iCellIndex] + sIndexString;
            printf("current string @ %d : %s\n", iCellIndex,
                   m_vsIntersectingPrecipitatesIDs[iCellIndex].c_str());
            m_viIntersectingPrecipitatesCounts[iCellIndex]++;
          }
        }
      }
    }
  }
}
Polyhedron *
PrecipitateStructure::GetContainingPrecipitate(const Point *poPoint) const
{
  // 	unsigned int iCellIndex = GetPointCellIndex(poPoint);
  // 	// see if this cell has no precipitates, return NULL
  // 	if(m_viIntersectingPrecipitatesCounts[iCellIndex] == 0)
  // 	{
  // 		return NULL;
  // 	}
  // 	unsigned int iIndex = 0;
  // 	// if the cell is completely contained in a precipitate, return it
  // 	if(m_viIntersectingPrecipitatesCounts[iCellIndex] == 1)
  // 	{
  // 		sscanf(m_vsIntersectingPrecipitatesIDs[iCellIndex].c_str(),"%d:",&iIndex);
  // 		return m_vpoPrecipitates[iIndex - 1];
  // 	}
  // 	// more than 1 precipitate intersect this cell, get them all
  // 	string sRemainingString = m_vsIntersectingPrecipitatesIDs[iCellIndex];
  // 	unsigned int iLocation = 0;
  // 	string sIndex = "";
  // 	list< Polyhedron* > lpoIntersectingPrecipitates;
  // 	lpoIntersectingPrecipitates.clear();
  // 	while(sRemainingString.size() > 1)
  // 	{
  // 		iLocation = sRemainingString.find_first_of(":");
  // 		sIndex = sRemainingString.substr(0,iLocation);
  // 		sRemainingString = sRemainingString.substr(iLocation + 1);
  // 		sscanf(sIndex.c_str(),"%d",&iIndex);
  // 		lpoIntersectingPrecipitates.push_back(m_vpoPrecipitates[iIndex -
  // 1]);
  // 	}
  //
  // 	// loop over the precipitates and see which one contains the point in
  // question, if any 	list< Polyhedron* >::iterator liPrecipitates;
  // Polyhedron*
  // poContainingPrecipitate = NULL; 	for(liPrecipitates =
  // lpoIntersectingPrecipitates.begin() ; liPrecipitates !=
  // lpoIntersectingPrecipitates.end() ; liPrecipitates++)
  // 	{
  // 		if((*liPrecipitates)->IsPointInside(*poPoint))
  // 		{
  // 			poContainingPrecipitate = (*liPrecipitates);
  // 			break;
  // 		}
  // 	}
  // 	lpoIntersectingPrecipitates.clear();
  // 	return poContainingPrecipitate;

  // for the naive way of doing it, just use the following
  unsigned int i = 0;
  for (i = 0; i < m_vpoPrecipitates.size(); i++)
  {
    if (m_vpoPrecipitates[i]->IsPointInside(*poPoint))
    {
      return m_vpoPrecipitates[i];
    }
  }
  return NULL;
}
void PrecipitateStructure::WritePrecipitateStructure(
    const string &sFileName) const
{
  unsigned int i = 0;
  unsigned int iSize = (unsigned int)m_vpoPrecipitates.size();
  FILE *fpFile = fopen(sFileName.c_str(), "w");
  fprintf(fpFile, "%d\n", iSize);
  for (i = 0; i < iSize; i++)
  {
    m_vpoPrecipitates[i]->WriteNodes(fpFile);
  }
  fclose(fpFile);
}
void PrecipitateStructure::WritePrecipitateStructureVTK(
    const string &sFileName) const
{
  unsigned int i = 0;
  unsigned int iSize = (unsigned int)m_vpoPrecipitates.size();
  char cFileName[256];
  for (i = 0; i < iSize; i++)
  {
    if (sFileName.empty())
    {
      sprintf(cFileName, "precipitate_%04d.vtk", i + 1);
    }
    else
    {
      sprintf(cFileName, "%s_%04d.vtk", sFileName.c_str(), i + 1);
    }
    m_vpoPrecipitates[i]->WriteParaview(cFileName);
  }
}
bool PrecipitateStructure::IsEmpty() const
{
  if (m_vpoPrecipitates.size() == 0)
  {
    return true;
  }
  return false;
}
const vector<Polyhedron *> *PrecipitateStructure::GetPrecipitates() const
{
  return &m_vpoPrecipitates;
}
unsigned int PrecipitateStructure::GetPrecipitatesCount() const
{
  return (unsigned int)m_vpoPrecipitates.size();
}
bool PrecipitateStructure::AddSphericalPrecipitate(const Point &oCenter,
                                                   const double &dRadius)
{
  if (dRadius <= 0.0)
  {
    return false;
  }

  Dodecahedron oSphereApproximation;
  oSphereApproximation.Set(2.0 * dRadius, false);
  oSphereApproximation.SetOrigin(oCenter);

  list<GenericNode *> *ploNodes = oSphereApproximation.GetTriangulationPoints();
  list<Point> loPoints;
  list<GenericNode *>::iterator liNodes;
  unsigned int iNodeID = 0;
  unsigned int iAttempt = 0;
  const unsigned int iMaxAttempts = 3;
  Polyhedron *poPrecipitate = NULL;
  bool bSuccess = false;

  for (iAttempt = 0; iAttempt < iMaxAttempts; iAttempt++)
  {
    loPoints.clear();
    iNodeID = 0;
    double dEpsilon = 0.0;
    if (iAttempt > 0)
    {
      dEpsilon = dRadius * 1.0e-8 * (double)iAttempt;
    }
    for (liNodes = ploNodes->begin(); liNodes != ploNodes->end(); liNodes++)
    {
      Point oPoint((*liNodes)->GetX(), (*liNodes)->GetY(), (*liNodes)->GetZ());
      if (dEpsilon > 0.0)
      {
        int iSignX = (iNodeID % 2 == 0) ? 1 : -1;
        int iSignY = (iNodeID % 3 == 0) ? 1 : -1;
        int iSignZ = (iNodeID % 5 == 0) ? 1 : -1;
        oPoint.Shift(iSignX * dEpsilon, iSignY * dEpsilon, iSignZ * dEpsilon);
      }
      loPoints.push_back(oPoint);
      iNodeID++;
    }

    poPrecipitate = new Polyhedron;
    if (poPrecipitate->CreateAsHullFromPoints(&loPoints))
    {
      bSuccess = true;
      break;
    }
    delete poPrecipitate;
    poPrecipitate = NULL;
  }

  if (!bSuccess || (poPrecipitate == NULL))
  {
    return false;
  }

  poPrecipitate->SetID((unsigned int)m_vpoPrecipitates.size() + 1);
  m_vpoPrecipitates.push_back(poPrecipitate);
  return true;
}
bool PrecipitateStructure::IsPointInsidePrecipitate(const Point &oPoint) const
{
  const Point *poPoint = &oPoint;
  if (GetContainingPrecipitate(poPoint) == NULL)
  {
    return false;
  }
  return true;
}
AxisAlignedBoundingBox PrecipitateStructure::GetBoundingBox() const
{
  double dXMin = DBL_MAX;
  double dXMax = -DBL_MAX;
  double dYMin = DBL_MAX;
  double dYMax = -DBL_MAX;
  double dZMin = DBL_MAX;
  double dZMax = -DBL_MAX;

  unsigned int i = 0;
  const AxisAlignedBoundingBox *poBox = NULL;
  for (i = 0; i < m_vpoPrecipitates.size(); i++)
  {
    poBox = m_vpoPrecipitates[i]->GetBox();
    if (poBox->GetXMin() < dXMin)
      dXMin = poBox->GetXMin();
    if (poBox->GetXMax() > dXMax)
      dXMax = poBox->GetXMax();
    if (poBox->GetYMin() < dYMin)
      dYMin = poBox->GetYMin();
    if (poBox->GetYMax() > dYMax)
      dYMax = poBox->GetYMax();
    if (poBox->GetZMin() < dZMin)
      dZMin = poBox->GetZMin();
    if (poBox->GetZMax() > dZMax)
      dZMax = poBox->GetZMax();
  }
  AxisAlignedBoundingBox oBox;
  oBox.SetXMin(dXMin);
  oBox.SetXMax(dXMax);
  oBox.SetYMin(dYMin);
  oBox.SetYMax(dYMax);
  oBox.SetZMin(dZMin);
  oBox.SetZMax(dZMax);
  return oBox;
}
void PrecipitateStructure::RemoveNonIntersectingPrecipitates(
    const Plane &oPlane)
{
  // see which precipitates intersect the plane
  list<Polyhedron *> lpoSurvivors;
  unsigned int i = 0;
  for (i = 0; i < m_vpoPrecipitates.size(); i++)
  {
    if (!m_vpoPrecipitates[i]->IsIntersecting(oPlane))
    {
      // delete those that do not
      delete m_vpoPrecipitates[i];
    }
    else
    {
      // keep those that do
      lpoSurvivors.push_back(m_vpoPrecipitates[i]);
    }
  }
  m_vpoPrecipitates.clear();
  m_vpoPrecipitates.resize(lpoSurvivors.size());
  // copy whatever survives to the new precipitates vector
  list<Polyhedron *>::iterator liSurvivors;
  i = 0;
  for (liSurvivors = lpoSurvivors.begin(); liSurvivors != lpoSurvivors.end();
       liSurvivors++)
  {
    m_vpoPrecipitates[i] = (*liSurvivors);
    m_vpoPrecipitates[i]->SetID(i + 1);
    i++;
  }
  lpoSurvivors.clear();
}
void PrecipitateStructure::RemoveExternalPrecipitates(
    Polyhedron *poVolumePolyhedron)
{
  // see which precipitates are internal
  list<Polyhedron *> lpoSurvivors;
  unsigned int i = 0;
  for (i = 0; i < m_vpoPrecipitates.size(); i++)
  {
    if (!poVolumePolyhedron->Contains(*m_vpoPrecipitates[i]))
    {
      // delete those that do not
      delete m_vpoPrecipitates[i];
    }
    else
    {
      // keep those that do
      lpoSurvivors.push_back(m_vpoPrecipitates[i]);
    }
  }
  m_vpoPrecipitates.clear();
  m_vpoPrecipitates.resize(lpoSurvivors.size());
  // copy whatever survives to the new precipitates vector
  list<Polyhedron *>::iterator liSurvivors;
  i = 0;
  for (liSurvivors = lpoSurvivors.begin(); liSurvivors != lpoSurvivors.end();
       liSurvivors++)
  {
    m_vpoPrecipitates[i] = (*liSurvivors);
    m_vpoPrecipitates[i]->SetID(i + 1);
    i++;
  }
  lpoSurvivors.clear();
}
void PrecipitateStructure::PlaneCut(const Plane &oPlane, const bool &bPerturb,
                                    const double &dPerturbationDistance)
{
  // see which precipitates intersect the plane or lie below it
  list<Polyhedron *> lpoSurvivors;
  unsigned int i = 0;
  Polyhedron *poNewPrecipitate = NULL;
  int iClassification = 0;
  for (i = 0; i < m_vpoPrecipitates.size(); i++)
  {
    iClassification = m_vpoPrecipitates[i]->ClassifyPlane(oPlane);
    if (iClassification == 1)
    {
      // the precipitate is above the plane, get rid of it
      delete m_vpoPrecipitates[i];
      m_vpoPrecipitates[i] = NULL;
    }
    else if (iClassification == -1)
    {
      // the precipitate is below the plane, keep it
      lpoSurvivors.push_back(m_vpoPrecipitates[i]);
    }
    else
    {
      // the precipitate intersects the plane, cut it
      poNewPrecipitate = m_vpoPrecipitates[i]->PlaneCut(oPlane, bPerturb,
                                                        dPerturbationDistance);
      lpoSurvivors.push_back(poNewPrecipitate);
      delete m_vpoPrecipitates[i];
    }
  }
  // copy whatever survives to the new precipitates vector
  m_vpoPrecipitates.clear();
  m_vpoPrecipitates.resize(lpoSurvivors.size());
  list<Polyhedron *>::iterator liSurvivors;
  i = 0;
  for (liSurvivors = lpoSurvivors.begin(); liSurvivors != lpoSurvivors.end();
       liSurvivors++)
  {
    m_vpoPrecipitates[i] = (*liSurvivors);
    m_vpoPrecipitates[i]->SetID(i + 1);
    i++;
  }
  lpoSurvivors.clear();
}
Polyhedron *PrecipitateStructure::GeneratePrecipitateInBox(
    AxisAlignedBoundingBox *poBox, const double &dVolumeFraction,
    const unsigned int &iPointsCount) const
{
  AxisAlignedBoundingBox oBox = poBox->GetSubVolume(dVolumeFraction);
  unsigned int i = 0;
  list<Point> loPoints;
  for (i = 0; i < iPointsCount; i++)
  {
    loPoints.push_back(oBox.GenerateRandomPoint());
  }
  Polyhedron *poPrecipitate = new Polyhedron;
  if (!poPrecipitate->CreateAsHullFromPoints(&loPoints))
  {
    delete poPrecipitate;
    poPrecipitate = NULL;
  }
  loPoints.clear();
  return poPrecipitate;
}
void PrecipitateStructure::GenerateInBox(AxisAlignedBoundingBox *poBox,
                                         const unsigned int &iXCount,
                                         const unsigned int &iYCount,
                                         const unsigned int &iZCount,
                                         const double &dVolumeFraction)
{
  Reset();
  list<AxisAlignedBoundingBox *> lpoBoxes =
      poBox->NonUniformPartition(iXCount, iYCount, iZCount);
  list<AxisAlignedBoundingBox *>::iterator liBoxes;
  Polyhedron *poPrecipitate = NULL;
  unsigned int iPointsCount = 500;
  unsigned int iCount = 0;
  m_vpoPrecipitates.resize(lpoBoxes.size());
  for (liBoxes = lpoBoxes.begin(); liBoxes != lpoBoxes.end(); liBoxes++)
  {
    poPrecipitate =
        GeneratePrecipitateInBox((*liBoxes), dVolumeFraction, iPointsCount);
    if (poPrecipitate == NULL)
      continue;
    poPrecipitate->SetID(iCount + 1);
    m_vpoPrecipitates[iCount] = poPrecipitate;
    iCount++;
    delete (*liBoxes);
    printf("in %d\n", iCount);
  }
}
