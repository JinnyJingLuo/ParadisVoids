#ifndef PARADISEXTERNALLOADSERVER_H_
#define PARADISEXTERNALLOADSERVER_H_

#include "FEMSolver.h"
#include "Typedefs.h"

using namespace std;

enum ParadisLoadType {
  NullLoadType = 0,
  ConstantStrainRateLoad = 1,
  CyclicStrainLoad = 2,
  CyclicStrainLoadLinearTemperature = 3
};

class _home;

class ParadisExternalLoadServer {
public:
  static ParadisExternalLoadServer *CreateInstance();
  ~ParadisExternalLoadServer();
  void Reset();
  bool Set(Home_t *poHome);
  void RunFEM(const double &dTime, Home_t *poHome);
  FEMSolver *GetFEMSolver() const;
  double GetTimeStep() const;
  double GetCurrentTime() const;
  void GetSegmentForce(Home_t *poHome, double bx, double by, double bz,
                       double x1, double y1, double z1, double x2, double y2,
                       double z2, double f1[3], double f2[3]);
  static Matrix
  ComputeSegmentStressAtPoint(const Point &oStart, const Point &oEnd,
                              const Vector &oBurgersVector, const Point &oPoint,
                              const double &dMu, const double &dNu);

private:
protected:
  static ParadisExternalLoadServer *m_poParadisExternalLoadServerInstance;
  ParadisExternalLoadServer();
  void Initialize();
  bool InitializeFEM(const string &sFileName);

  bool GenerateHeader(FILE *fpFile, Home_t *poHome) const;
  bool GenerateLoads(FILE *fpFile, Home_t *poHome) const;
  bool GenerateMaterials(FILE *fpFile, Home_t *poHome) const;

  bool GenerateGeometry(FILE *fpFile, Home_t *poHome) const;
  bool GenerateSolidCylindricalGeometry(FILE *fpFile, Home_t *poHome) const;
  bool GenerateThermoMechanicalCylindricalGeometry(FILE *fpFile,
                                                   Home_t *poHome) const;
  bool GenerateSolidBlockGeometry(FILE *fpFile, Home_t *poHome) const;
  bool GenerateThermoMechanicalBlockGeometry(FILE *fpFile,
                                             Home_t *poHome) const;

  void ComputeSurfaceStresses(Home_t *poHome);

  MainDataStructure *m_poData;
  FEMSolver *m_poFEMSolver;

  ParadisLoadType m_eLoadType;
  double m_dStrainRate;
  double m_dLoadAmplitude;
  double m_dLoadFrequency;
  double m_dInitialTemperature;
  double m_dHeatingRate;
  double m_dTimeStep;
  unsigned int m_iProcessID;
  bool m_bUseFEM;
  list<FEMBoundaryElementFace *> *m_plpoBoundaryFaces;
};

#endif
