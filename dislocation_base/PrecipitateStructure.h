#ifndef PRECIPITATESTRUCTURE_H_
#define PRECIPITATESTRUCTURE_H_

#include "Polyhedron.h"
#include "vector"

using namespace std;
using namespace GeometrySystem;

class PrecipitateStructure {
public:
  PrecipitateStructure();
  PrecipitateStructure(const PrecipitateStructure &oStructure);
  ~PrecipitateStructure();
  PrecipitateStructure &operator=(const PrecipitateStructure &oStructure);
  void Reset();
  bool Set(const string &sFileName);
  bool Set(const string &sFileName, const double &dXDimension,
           const double &dYDimension, const double &dZDimension);
  void WritePrecipitateStructure(const string &sFileName) const;
  void WritePrecipitateStructureVTK(const string &sFileName = "") const;
  bool IsEmpty() const;
  Polyhedron *GetContainingPrecipitate(const Point *poPoint) const;
  const vector<Polyhedron *> *GetPrecipitates() const;
  unsigned int GetPrecipitatesCount() const;
  void WriteAPB() const;
  bool IsPointInsidePrecipitate(const Point &oPoint) const;
  AxisAlignedBoundingBox GetBoundingBox() const;
  bool AddSphericalPrecipitate(const Point &oCenter, const double &dRadius);
  void RemoveNonIntersectingPrecipitates(const Plane &oPlane);
  void RemoveExternalPrecipitates(Polyhedron *poVolumePolyhedron);
  void PlaneCut(const Plane &oPlane, const bool &bPerturb = false,
                const double &dPerturbationDistance = 0.0);
  void GenerateInBox(AxisAlignedBoundingBox *poBox, const unsigned int &iXCount,
                     const unsigned int &iYCount, const unsigned int &iZCount,
                     const double &dVolumeFraction);

private:
protected:
  static double m_dCellSize;
  void Initialize();
  bool ReadPrecipitatesFile(const string &sFileName);
  unsigned int GetPointCellIndex(const Point *poPoint) const;
  void ComputeCellPrecipitateIntersection();
  Polyhedron *
  GeneratePrecipitateInBox(AxisAlignedBoundingBox *poBox,
                           const double &dVolumeFraction,
                           const unsigned int &iPointsCount = 1000) const;

  unsigned int m_iXCellsCount;
  unsigned int m_iYCellsCount;
  unsigned int m_iZCellsCount;
  vector<Polyhedron *> m_vpoPrecipitates;
  vector<string> m_vsIntersectingPrecipitatesIDs;
  vector<unsigned int> m_viIntersectingPrecipitatesCounts;
  double m_dXShift;
  double m_dYShift;
  double m_dZShift;
};

#endif
