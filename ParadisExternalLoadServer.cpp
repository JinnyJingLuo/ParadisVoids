#include "ParadisExternalLoadServer.h"
#include "Tools.h"
#include "Block.h"
#include "DSMPI.h"
#include "FEMMesh.h"
#include "FEMAxialTimeLinearTotalStrainLoad.h"
#include "FEMAxialTimeSinusoidalTotalStrainLoad.h"
#include "FEMImplicitDynamicsSolidSolver.h"
#include "FEMBoundaryElementFace.h"
#include "Home.h"
#include "MathServices.h"

using namespace SupportSystem;
using namespace FEMSystem;
using namespace EZ;

ParadisExternalLoadServer
    *ParadisExternalLoadServer::m_poParadisExternalLoadServerInstance = NULL;
ParadisExternalLoadServer *ParadisExternalLoadServer::CreateInstance() {
  if (m_poParadisExternalLoadServerInstance == NULL) {
    m_poParadisExternalLoadServerInstance = new ParadisExternalLoadServer;
  }
  return m_poParadisExternalLoadServerInstance;
}
ParadisExternalLoadServer::ParadisExternalLoadServer() { Initialize(); }
ParadisExternalLoadServer::~ParadisExternalLoadServer() { Reset(); }
void ParadisExternalLoadServer::Reset() {
  if (m_poData != NULL) {
    m_poData->Reset();
    delete m_poData;
  }
  if (m_poFEMSolver != NULL) {
    m_poFEMSolver->Reset();
    delete m_poFEMSolver;
  }
  m_plpoBoundaryFaces = NULL;
}
bool ParadisExternalLoadServer::Set(Home_t *poHome) {
  // read relevant DD parameters
  m_iProcessID = poHome->myDomain;
  if (poHome->param->EnableFEM != 0) {
    m_bUseFEM = true;
  } else {
    m_bUseFEM = false;
  }

  if (poHome->param->loadType == 1) {
    m_eLoadType = ConstantStrainRateLoad;
  } else if (poHome->param->loadType == 7) {
    m_eLoadType = CyclicStrainLoad;
  } else if (poHome->param->loadType == 9) {
    m_eLoadType = CyclicStrainLoadLinearTemperature;
  } else {
    if (m_bUseFEM) {
      printf("error: unknwon FEM load type\n");
      return false;
    }
  }

  m_dStrainRate = poHome->param->eRate;
  m_dLoadAmplitude = poHome->param->CyclicLoadAmplitude;
  m_dLoadFrequency = poHome->param->CyclicLoadFrequency;
  m_dInitialTemperature = poHome->param->InitialTemperature;
  m_dHeatingRate = poHome->param->HeatingRate;

  // now this run is either a fresh run or a restart run
  // for fresh runs, generate the fem problem file and use it, for restart runs,
  // do NOT run unless there is an FEM restart file available. the name of the
  // FEM restart file should be given in paradis control file parameter
  // FEMRestartFileName
  string sFEMInputFile = poHome->param->FEMRestartFileName;
  bool bIsValidModel = true;
  if (m_bUseFEM) {
    if (sFEMInputFile.empty()) {
      // fresh run, generate the FEM problem and write it in a file
      // only process 0 writes the input file
      sFEMInputFile = "fem_model.txt";
      if (m_iProcessID == 0) {
        FILE *fpFEMInput = fopen(sFEMInputFile.c_str(), "w");
        // file organization
        // 1. header
        GenerateHeader(fpFEMInput, poHome);
        // 2. loads
        bIsValidModel = GenerateLoads(fpFEMInput, poHome);
        // 3. materials
        bIsValidModel = bIsValidModel && GenerateMaterials(fpFEMInput, poHome);
        // 4. nodes (type, position, surface, initial conditions, boundary
        // conditions)
        // 5. elements (connectivity, boundary conditions)
        bIsValidModel = bIsValidModel && GenerateGeometry(fpFEMInput, poHome);
        fclose(fpFEMInput);
      }
      DSMPI::Barrier();
    } else {
      // restart run, use the supplied FEM restart file, if this file doesn't
      // exist, then the model is not valid
      sFEMInputFile = poHome->param->FEMRestartFileName;
      FILE *fpTestFile = fopen(sFEMInputFile.c_str(), "r");
      if (fpTestFile == NULL) {
        printf("error: couldn't find FEM restart file %s\n",
               sFEMInputFile.c_str());
        return false;
      }
      fclose(fpTestFile);
    }
  }

  int iIsValidModel = 0;
  if (bIsValidModel) {
    iIsValidModel = 1;
  }
  DSMPI::Broadcast(&iIsValidModel, 1, MPI_INT, 0);
  if (iIsValidModel == 0) {
    return false;
  }
  // read the created file and build the FEM model
  if (m_bUseFEM) {
    InitializeFEM(sFEMInputFile);
  }
  // finally, get the surface Gauss points and store them
  unsigned int iGaussPointsCount = 0;
  return true;
}
void ParadisExternalLoadServer::RunFEM(const double &dTime, Home_t *poHome) {
  // apply the plastic strain, the total strain load is always load number 2
  // (index 1) in the FEM load vector
  double al = poHome->param->edotdir[0];
  double am = poHome->param->edotdir[1];
  double an = poHome->param->edotdir[2];
  double amag = sqrt(al * al + am * am + an * an);
  al /= amag;
  am /= amag;
  an /= amag;

  double dAxialPlasticStrain = poHome->param->totpStn[0] * al * al +
                               poHome->param->totpStn[1] * am * am +
                               poHome->param->totpStn[2] * an * an +
                               2.0 * (poHome->param->totpStn[3] * am * an +
                                      poHome->param->totpStn[4] * an * al +
                                      poHome->param->totpStn[5] * al * am);
  if (m_eLoadType == ConstantStrainRateLoad) {
    ((FEMAxialTimeLinearTotalStrainLoad *)(m_poData->GetLoads()->at(1)))
        ->SetCurrentPlasticStrain(dAxialPlasticStrain);
  } else if ((m_eLoadType == CyclicStrainLoad) ||
             (m_eLoadType == CyclicStrainLoadLinearTemperature)) {
    ((FEMAxialTimeSinusoidalTotalStrainLoad *)(m_poData->GetLoads()->at(1)))
        ->SetCurrentPlasticStrain(dAxialPlasticStrain);
  }

  // update the FEM main data structure
  m_poData->UpdateTime(dTime);
  m_poData->IncrementCurrentOutputCount();

  // run the FEM solver for 1 step
  // first, calculate image forces
  ComputeSurfaceStresses(poHome);
  // apply the loads
  m_poData->ApplyLoads(dTime);
  // apply the image forces
  m_poData->ApplyImageForces();
  // solve the system
  m_poFEMSolver->Solve(dTime);
  // update the output count and write the output
  unsigned int iOutputCount = m_poData->GetCurrentOutputCount();
  // only process 0 writes the output
  if (m_iProcessID == 0) {
    m_poFEMSolver->WriteFEMSolution(iOutputCount);
    m_poFEMSolver->WriteFEMSolutionToParaview(iOutputCount);
  }
  sprintf(poHome->param->FEMRestartFileName, "%s",
          (m_poFEMSolver->GetOutputFileName(iOutputCount)).c_str());
  DSMPI::Barrier();
}
FEMSolver *ParadisExternalLoadServer::GetFEMSolver() const {
  return m_poFEMSolver;
}
double ParadisExternalLoadServer::GetTimeStep() const { return m_dTimeStep; }
double ParadisExternalLoadServer::GetCurrentTime() const {
  if (m_bUseFEM) {
    return m_poData->GetCurrentTime();
  }
  return 0.0;
}
void ParadisExternalLoadServer::GetSegmentForce(Home_t *poHome, double bx,
                                                double by, double bz, double x1,
                                                double y1, double z1, double x2,
                                                double y2, double z2,
                                                double f1[3], double f2[3]) {
  // if the FEM is enabled, get the stress from the FEM server, otherwise, use
  // the stress passed to this function
  double strb[3];

  if (poHome->param->EnableFEM == 1) {

    // get FEM stress at the midpoint of the segment, convert coordinates to
    // microns
    double dFactor = 0.5; //*poHome->param->burgMag;
    Point oMidPoint(dFactor * (x1 + x2), dFactor * (y1 + y2),
                    dFactor * (z1 + z2));
    unsigned int iStatus;
    Matrix oStress = ((FEMImplicitDynamicsSolidSolver *)m_poFEMSolver)
                         ->GetStress(&oMidPoint, iStatus);
    // the stress coming from the FEM is in kg /micron sec^2, multiply by 1.0E6
    // to convert it to Pa
    oStress = oStress * 1.0E6;
    strb[0] = oStress.Get(1, 1) * bx + oStress.Get(1, 2) * by +
              oStress.Get(1, 3) * bz;
    strb[1] = oStress.Get(2, 1) * bx + oStress.Get(2, 2) * by +
              oStress.Get(2, 3) * bz;
    strb[2] = oStress.Get(3, 1) * bx + oStress.Get(3, 2) * by +
              oStress.Get(3, 3) * bz;
  } else if (poHome->param->EnableFEM == 2) {
    // get FEM stress at the midpoint of the segment, convert coordinates to
    // microns
    //		double dFactor = 5.0E5*poHome->param->burgMag;
    // update KS: Now the coordinates are in units of b. So no need to scale
    // here as this way
    // it is conform with the ParaDIS input.
    Point oMidPoint(0.5 * (x1 + x2), 0.5 * (y1 + y2), 0.5 * (z1 + z2));
    unsigned int iStatus;
    Matrix oStress = ((FEMImplicitDynamicsSolidSolver *)m_poFEMSolver)
                         ->GetStress(&oMidPoint, iStatus);
    // the stress coming from the FEM is in kg /micron sec^2, multiply by 1.0E6
    // to convert it to Pa
    oStress = oStress * 1e6;
    strb[0] = oStress.Get(1, 1) * bx + oStress.Get(1, 2) * by +
              oStress.Get(1, 3) * bz;
    strb[1] = oStress.Get(2, 1) * bx + oStress.Get(2, 2) * by +
              oStress.Get(2, 3) * bz;
    strb[2] = oStress.Get(3, 1) * bx + oStress.Get(3, 2) * by +
              oStress.Get(3, 3) * bz;

    strb[0] += poHome->param->appliedStress[0] * bx +
               poHome->param->appliedStress[5] * by +
               poHome->param->appliedStress[4] * bz;
    strb[1] += poHome->param->appliedStress[5] * bx +
               poHome->param->appliedStress[1] * by +
               poHome->param->appliedStress[3] * bz;
    strb[2] += poHome->param->appliedStress[4] * bx +
               poHome->param->appliedStress[3] * by +
               poHome->param->appliedStress[2] * bz;

  } else {
    strb[0] = poHome->param->appliedStress[0] * bx +
              poHome->param->appliedStress[5] * by +
              poHome->param->appliedStress[4] * bz;
    strb[1] = poHome->param->appliedStress[5] * bx +
              poHome->param->appliedStress[1] * by +
              poHome->param->appliedStress[3] * bz;
    strb[2] = poHome->param->appliedStress[4] * bx +
              poHome->param->appliedStress[3] * by +
              poHome->param->appliedStress[2] * bz;
    if (poHome->param->SType > 0 && poHome->param->stress_timestep > 0) {
      real8 pos_mid[3];
      real8 extStress[6][5];
      pos_mid[X] = (x1 + x2) * 0.5;
      pos_mid[Y] = (y1 + y2) * 0.5;
      pos_mid[Z] = (z1 + z2) * 0.5;
      // printf("%lf %lf", pos_mid[X], pos_mid[Y]);
      GetInLocalCoordinates(poHome->param, pos_mid[X], pos_mid[Y], pos_mid[Z]);
      real8 posm[3];
      // printf("%lf %lf\n", pos_mid[X], pos_mid[Y]);
      int n_x, n_y, n_z;
      n_x = ceil((pos_mid[X] + 0.5 * poHome->param->stress_Dim[0]) /
                 poHome->param->stress_Dim[0] *
                 (double)(poHome->param->GP_x - 1));
      n_y = ceil((pos_mid[Y] + 0.5 * poHome->param->stress_Dim[1]) /
                 poHome->param->stress_Dim[1] *
                 (double)(poHome->param->GP_y - 1));
      n_z = ceil((pos_mid[Z] + 0.5 * poHome->param->stress_Dim[2]) /
                 poHome->param->stress_Dim[2] *
                 (double)(poHome->param->GP_z - 1));

      if (n_x >= 0 && n_x <= poHome->param->GP_x - 1 && n_y >= 0 &&
          n_y <= poHome->param->GP_y - 1 && n_z >= 0 &&
          n_z <= poHome->param->GP_z - 1) {
        n_x = fmax(n_x, 1);
        n_y = fmax(n_y, 1);
        n_z = fmax(n_z, 1);

        int s_timestep = poHome->param->stress_timestep - 1;
        real8 pos_stress0[3], pos_stress1[3];
        pos_stress0[0] =
            poHome->stress[n_x - 1][n_y - 1][n_z - 1][s_timestep].x;
        pos_stress0[1] =
            poHome->stress[n_x - 1][n_y - 1][n_z - 1][s_timestep].y;
        pos_stress0[2] =
            poHome->stress[n_x - 1][n_y - 1][n_z - 1][s_timestep].z;
        pos_stress1[0] = poHome->stress[n_x][n_y][n_z][s_timestep].x;
        pos_stress1[1] = poHome->stress[n_x][n_y][n_z][s_timestep].y;
        pos_stress1[2] = poHome->stress[n_x][n_y][n_z][s_timestep].z;
        // GetInLocalCoordinates(poHome->param,pos_stress1[0], pos_stress1[1],
        // pos_stress1[2]);
        real8 xd, yd, zd;
        // trilinear interpolation
        xd = (pos_mid[X] - pos_stress0[0]) / (pos_stress1[0] - pos_stress0[0]);
        yd = (pos_mid[Y] - pos_stress0[1]) / (pos_stress1[1] - pos_stress0[1]);
        zd = (pos_mid[Z] - pos_stress0[2]) / (pos_stress1[2] - pos_stress0[2]);

        real8 extStress0, extStress1, extStress2, extStress3, extStress4,
            extStress5;
        real8 c000, c001, c010, c011, c100, c101, c110, c111;
        /* Vxyz =	V000 (1 - x) (1 - y) (1 - z) +
                V100 x (1 - y) (1 - z) +
                V010 (1 - x) y (1 - z) +
                V001 (1 - x) (1 - y) z +
                V101 x (1 - y) z +
                V011 (1 - x) y z +
                V110 x y (1 - z) +
                V111 x y z
        */
        c000 = poHome->stress[n_x - 1][n_y - 1][n_z - 1][s_timestep].s_xx;
        c001 = poHome->stress[n_x - 1][n_y - 1][n_z][s_timestep].s_xx;
        c010 = poHome->stress[n_x - 1][n_y][n_z - 1][s_timestep].s_xx;
        c011 = poHome->stress[n_x - 1][n_y][n_z][s_timestep].s_xx;
        c100 = poHome->stress[n_x][n_y - 1][n_z - 1][s_timestep].s_xx;
        c101 = poHome->stress[n_x][n_y - 1][n_z][s_timestep].s_xx;
        c110 = poHome->stress[n_x][n_y][n_z - 1][s_timestep].s_xx;
        c111 = poHome->stress[n_x][n_y][n_z][s_timestep].s_xx;
        extStress0 = ((c000 * (1 - xd) + c100 * xd) * (1 - yd) +
                      (c010 * (1 - xd) + c110 * xd) * yd) *
                         (1 - zd) +
                     ((c001 * (1 - xd) + c101 * xd) * (1 - yd) +
                      (c011 * (1 - xd) + c111 * xd) * yd) *
                         zd;

        c000 = poHome->stress[n_x - 1][n_y - 1][n_z - 1][s_timestep].s_yy;
        c001 = poHome->stress[n_x - 1][n_y - 1][n_z][s_timestep].s_yy;
        c010 = poHome->stress[n_x - 1][n_y][n_z - 1][s_timestep].s_yy;
        c011 = poHome->stress[n_x - 1][n_y][n_z][s_timestep].s_yy;
        c100 = poHome->stress[n_x][n_y - 1][n_z - 1][s_timestep].s_yy;
        c101 = poHome->stress[n_x][n_y - 1][n_z][s_timestep].s_yy;
        c110 = poHome->stress[n_x][n_y][n_z - 1][s_timestep].s_yy;
        c111 = poHome->stress[n_x][n_y][n_z][s_timestep].s_yy;
        extStress1 = ((c000 * (1 - xd) + c100 * xd) * (1 - yd) +
                      (c010 * (1 - xd) + c110 * xd) * yd) *
                         (1 - zd) +
                     ((c001 * (1 - xd) + c101 * xd) * (1 - yd) +
                      (c011 * (1 - xd) + c111 * xd) * yd) *
                         zd;

        c000 = poHome->stress[n_x - 1][n_y - 1][n_z - 1][s_timestep].s_zz;
        c001 = poHome->stress[n_x - 1][n_y - 1][n_z][s_timestep].s_zz;
        c010 = poHome->stress[n_x - 1][n_y][n_z - 1][s_timestep].s_zz;
        c011 = poHome->stress[n_x - 1][n_y][n_z][s_timestep].s_zz;
        c100 = poHome->stress[n_x][n_y - 1][n_z - 1][s_timestep].s_zz;
        c101 = poHome->stress[n_x][n_y - 1][n_z][s_timestep].s_zz;
        c110 = poHome->stress[n_x][n_y][n_z - 1][s_timestep].s_zz;
        c111 = poHome->stress[n_x][n_y][n_z][s_timestep].s_zz;
        extStress2 = ((c000 * (1 - xd) + c100 * xd) * (1 - yd) +
                      (c010 * (1 - xd) + c110 * xd) * yd) *
                         (1 - zd) +
                     ((c001 * (1 - xd) + c101 * xd) * (1 - yd) +
                      (c011 * (1 - xd) + c111 * xd) * yd) *
                         zd;

        c000 = poHome->stress[n_x - 1][n_y - 1][n_z - 1][s_timestep].s_yz;
        c001 = poHome->stress[n_x - 1][n_y - 1][n_z][s_timestep].s_yz;
        c010 = poHome->stress[n_x - 1][n_y][n_z - 1][s_timestep].s_yz;
        c011 = poHome->stress[n_x - 1][n_y][n_z][s_timestep].s_yz;
        c100 = poHome->stress[n_x][n_y - 1][n_z - 1][s_timestep].s_yz;
        c101 = poHome->stress[n_x][n_y - 1][n_z][s_timestep].s_yz;
        c110 = poHome->stress[n_x][n_y][n_z - 1][s_timestep].s_yz;
        c111 = poHome->stress[n_x][n_y][n_z][s_timestep].s_yz;
        extStress3 = ((c000 * (1 - xd) + c100 * xd) * (1 - yd) +
                      (c010 * (1 - xd) + c110 * xd) * yd) *
                         (1 - zd) +
                     ((c001 * (1 - xd) + c101 * xd) * (1 - yd) +
                      (c011 * (1 - xd) + c111 * xd) * yd) *
                         zd;

        c000 = poHome->stress[n_x - 1][n_y - 1][n_z - 1][s_timestep].s_xz;
        c001 = poHome->stress[n_x - 1][n_y - 1][n_z][s_timestep].s_xz;
        c010 = poHome->stress[n_x - 1][n_y][n_z - 1][s_timestep].s_xz;
        c011 = poHome->stress[n_x - 1][n_y][n_z][s_timestep].s_xz;
        c100 = poHome->stress[n_x][n_y - 1][n_z - 1][s_timestep].s_xz;
        c101 = poHome->stress[n_x][n_y - 1][n_z][s_timestep].s_xz;
        c110 = poHome->stress[n_x][n_y][n_z - 1][s_timestep].s_xz;
        c111 = poHome->stress[n_x][n_y][n_z][s_timestep].s_xz;
        extStress4 = ((c000 * (1 - xd) + c100 * xd) * (1 - yd) +
                      (c010 * (1 - xd) + c110 * xd) * yd) *
                         (1 - zd) +
                     ((c001 * (1 - xd) + c101 * xd) * (1 - yd) +
                      (c011 * (1 - xd) + c111 * xd) * yd) *
                         zd;

        c000 = poHome->stress[n_x - 1][n_y - 1][n_z - 1][s_timestep].s_xy;
        c001 = poHome->stress[n_x - 1][n_y - 1][n_z][s_timestep].s_xy;
        c010 = poHome->stress[n_x - 1][n_y][n_z - 1][s_timestep].s_xy;
        c011 = poHome->stress[n_x - 1][n_y][n_z][s_timestep].s_xy;
        c100 = poHome->stress[n_x][n_y - 1][n_z - 1][s_timestep].s_xy;
        c101 = poHome->stress[n_x][n_y - 1][n_z][s_timestep].s_xy;
        c110 = poHome->stress[n_x][n_y][n_z - 1][s_timestep].s_xy;
        c111 = poHome->stress[n_x][n_y][n_z][s_timestep].s_xy;
        extStress5 = ((c000 * (1 - xd) + c100 * xd) * (1 - yd) +
                      (c010 * (1 - xd) + c110 * xd) * yd) *
                         (1 - zd) +
                     ((c001 * (1 - xd) + c101 * xd) * (1 - yd) +
                      (c011 * (1 - xd) + c111 * xd) * yd) *
                         zd;
        strb[0] += extStress0 * bx + extStress5 * by + extStress4 * bz;
        strb[1] += extStress5 * bx + extStress1 * by + extStress3 * bz;
        strb[2] += extStress4 * bx + extStress3 * by + extStress2 * bz;
      }
    }
  }

  double xi[3];
  xi[0] = x2 - x1;
  xi[1] = y2 - y1;
  xi[2] = z2 - z1;

  double ft[3];
  cross(strb, xi, ft);

  int i;
  for (i = 0; i < 3; i++) {
    f1[i] = ft[i] * 0.5;
    f2[i] = ft[i] * 0.5;
  }
}
void ParadisExternalLoadServer::Initialize() {
  m_poData = NULL;
  m_poFEMSolver = NULL;
  m_eLoadType = NullLoadType;
  m_dStrainRate = 0.0;
  m_dLoadAmplitude = 0.0;
  m_dLoadFrequency = 0.0;
  m_dInitialTemperature = 0.0;
  m_dHeatingRate = 0.0;
  m_dTimeStep = 0.0;
  m_iProcessID = 0;
  m_bUseFEM = false;
  m_plpoBoundaryFaces = NULL;
}
bool ParadisExternalLoadServer::InitializeFEM(const string &sFileName) {
  m_poData = MainDataStructure::CreateInstance();
  m_poData->SetInputFileName(sFileName);
  // read FEM problem
  printf("Generating FEM model\n");
  printf("Reading FEM input from file %s\n", sFileName.c_str());
  m_poData->ReadInput();
  m_poData->GenerateBoundaryFaces();
  m_plpoBoundaryFaces = m_poData->GetBoundaryFaces();
  printf("Done generating FEM model\n");
  m_poFEMSolver = FEMSolver::CreateSolverByPhysicsIndex(
      m_poData->GetProblemPhysics(), m_poData->GetProblemType());
  if (m_poFEMSolver == NULL) {
    printf("error: unknown problem physics %d\n",
           m_poData->GetProblemPhysics());
    return false;
  }
  m_poFEMSolver->SetDataStructure(m_poData);
  // initialize the FEM solver
  // initialize the matrices using the time step given in the input file, during
  // initialization, the time step will be readjusted and stored in the solver
  // data
  m_poFEMSolver->InitializeMatrices(m_poData->GetTimeStep());
  m_dTimeStep = m_poFEMSolver->GetTimeStep();
  printf("FEM time step : %25.20f\n", m_dTimeStep);
  // now the time step has been adjusted, store it since it will be used to tell
  // the FEM solver when to run
  return true;
}
bool ParadisExternalLoadServer::GenerateHeader(FILE *fpFile,
                                               Home_t *poHome) const {
  fprintf(fpFile, "* output base file name\n");
  fprintf(fpFile, "fem_out\n");
  fprintf(fpFile, "* problem type, physics, current time, time step, target "
                  "time, this output count\n");
  // the FEM target time is not important at all for this use case, can be
  // anything
  if (m_eLoadType == CyclicStrainLoadLinearTemperature) {
    fprintf(fpFile, "3,3,0.0,%e,1.0,0\n", poHome->param->FEMTimeStep);
  } else {
    fprintf(fpFile, "3,2,0.0,%e,1.0,0\n", poHome->param->FEMTimeStep);
  }
  return true;
}
bool ParadisExternalLoadServer::GenerateLoads(FILE *fpFile,
                                              Home_t *poHome) const {
  // loads count (either 2 or 3)
  fprintf(fpFile, "* Loads : \n");
  unsigned int iLoadsCount = 2;
  if (m_eLoadType == CyclicStrainLoadLinearTemperature) {
    iLoadsCount = 3;
  }
  fprintf(fpFile, "%d\n", iLoadsCount);

  // first load (constant with zero)
  fprintf(fpFile, "1\n");   // 1 for constant load
  fprintf(fpFile, "0.0\n"); // 0.0 for zero load

  // second load
  double dBurgersVectorMagnitude = poHome->param->burgMag;
  // write everything in microns
  double dLength =
      1.0E6 * dBurgersVectorMagnitude * poHome->param->Dimensions[2];
  if (m_eLoadType == ConstantStrainRateLoad) {
    fprintf(fpFile, "6\n"); // 6 for axial time linear total strain load
    // the axial direction is always the z direction
    fprintf(fpFile, "%20.15f\t%20.15f\n", dLength, m_dStrainRate);
  } else if (m_eLoadType == CyclicStrainLoad) {
    fprintf(fpFile, "7\n"); // 7 for axial time sinusoidal total strain load
    fprintf(fpFile, "%20.15f\t%20.15f\t%20.15f\n", dLength, m_dLoadAmplitude,
            m_dLoadFrequency);
  } else if (m_eLoadType == CyclicStrainLoadLinearTemperature) {
    fprintf(fpFile, "7\n"); // 7 for axial time sinusoidal total strain load
    fprintf(fpFile, "%20.15f\t%20.15f\t%20.15f\n", dLength, m_dLoadAmplitude,
            m_dLoadFrequency);
    fprintf(fpFile, "5\n"); // 5 for linear temperature load
    fprintf(fpFile, "%20.15f\t%20.15f\n", m_dInitialTemperature,
            m_dHeatingRate);
  } else {
    printf("error: unrecognized FEM load : %d\n", m_eLoadType);
    return false;
  }
  return true;
}
bool ParadisExternalLoadServer::GenerateMaterials(FILE *fpFile,
                                                  Home_t *poHome) const {
  fprintf(fpFile, "* Materials : \n");
  unsigned int iMaterailsCount = 1;
  unsigned int iLinearIsotropicMaterialType = 1;
  fprintf(fpFile, "%d\n", iMaterailsCount);
  // the one and only material
  fprintf(fpFile, "%d\n", iLinearIsotropicMaterialType);
  // set the stress unit to kg/micron sec^2
  double dYoungsModulus = poHome->param->YoungsModulus * 1.0E-6;
  // set the density unit to kg/micron^3
  double dDensity = poHome->param->MassDensity * 1.0E-18;
  // set the thermal conductivity unit to kg micron/s^3/K
  double dThermalConductivity = poHome->param->ThermalConductivity * 1.0E6;
  // set the specific heat capacity unit to micron^2/sec^2/K
  double dHeatCapacity = poHome->param->SpecificHeatCapacity * 1.0E12;

  // properties E nu rho k alpha c To
  fprintf(fpFile, "%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\n", dYoungsModulus,
          poHome->param->PoissonsRatio, dDensity, dThermalConductivity,
          poHome->param->ThermalExpansionCoefficient, dHeatCapacity,
          poHome->param->ReferenceTemperature);
  return true;
}
bool ParadisExternalLoadServer::GenerateGeometry(FILE *fpFile,
                                                 Home_t *poHome) const {
  if (poHome->param->GeometryType == 1) {
    if (m_eLoadType == CyclicStrainLoadLinearTemperature) {
      return GenerateThermoMechanicalBlockGeometry(fpFile, poHome);
    } else {
      return GenerateSolidBlockGeometry(fpFile, poHome);
    }
  } else if (poHome->param->GeometryType == 4) {
    if (m_eLoadType == CyclicStrainLoadLinearTemperature) {
      return GenerateThermoMechanicalCylindricalGeometry(fpFile, poHome);
    } else {
      return GenerateSolidCylindricalGeometry(fpFile, poHome);
    }
  } else {
    printf("error: unrecognized FEM geometry : %d\n",
           poHome->param->GeometryType);
    return false;
  }
}
bool ParadisExternalLoadServer::GenerateSolidBlockGeometry(
    FILE *fpFile, Home_t *poHome) const {
  Block *poBlock = new Block;
  double dXLength = poHome->param->Dimensions[0];
  double dYLength = poHome->param->Dimensions[1];
  double dZLength = poHome->param->Dimensions[2];
  double dMinDimension = min(min(dXLength, dYLength), dZLength);
  // make sure the shortest dimension has at least 5 elements and that the
  // elements aspect ratio is close to 1
  double dElementLength = 0.2 * dMinDimension;
  unsigned int iXResolution = (unsigned int)ceil(dXLength / dElementLength);
  unsigned int iYResolution = (unsigned int)ceil(dYLength / dElementLength);
  unsigned int iZResolution = (unsigned int)ceil(dZLength / dElementLength);
  // write everything in microns
  double dFactor = 1;
  1.0E6 * poHome->param->burgMag;
  dXLength = dXLength * dFactor; // scaling factor
  dYLength = dYLength * dFactor;
  dZLength = dZLength * dFactor;
  poBlock->SetXLength(dXLength);
  poBlock->SetYLength(dYLength);
  poBlock->SetZLength(dZLength);

  poBlock->SetResolution(iXResolution, iYResolution, iZResolution);

  vector<FEMNode *> vpoNodes;
  vector<FEMElement *> vpoElements;
  FEMMesh::GenerateMeshFromBlock(poBlock, vpoNodes, vpoElements, SolidFEMNode,
                                 SolidFEMElement);
  delete poBlock;
  unsigned int iNodesCount = (unsigned int)vpoNodes.size();
  unsigned int iElementsCount = (unsigned int)vpoElements.size();
  unsigned int i = 0;
  unsigned int iTemp = 0;
  char cTempString[1024];
  double dTolerance = 1.0E-6;

  // nodes section
  {
    fprintf(fpFile, "* Nodes : \n");
    fprintf(fpFile, "%d\n", iNodesCount);
    for (i = 0; i < iNodesCount; i++) {
      iTemp = 0;
      if (vpoNodes[i]->IsOnSurface()) {
        iTemp = 1;
      }
      fprintf(fpFile, "2\n");
      if (fabs(vpoNodes[i]->GetZ() + dZLength / 2.0) < dTolerance) {
        fprintf(fpFile, "%e,%e,%e,1,1,1,1,1,1,%d\n", vpoNodes[i]->GetX(),
                vpoNodes[i]->GetY(), vpoNodes[i]->GetZ(), iTemp);
      } else if (fabs(vpoNodes[i]->GetZ() - dZLength / 2.0) < dTolerance) {
        fprintf(fpFile, "%e,%e,%e,1,1,1,1,1,2,%d\n", vpoNodes[i]->GetX(),
                vpoNodes[i]->GetY(), vpoNodes[i]->GetZ(), iTemp);
      } else {
        fprintf(fpFile, "%e,%e,%e,0,1,0,1,0,1,%d\n", vpoNodes[i]->GetX(),
                vpoNodes[i]->GetY(), vpoNodes[i]->GetZ(), iTemp);
      }
      // initial displacements, velocities and accelerations
      fprintf(fpFile, "0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0\n");
      // forces and stresses
      fprintf(fpFile, "0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0\n");
    }
  }

  // elements section
  {
    unsigned int j = 0;
    vector<FEMNode *> *pvpoNodes = NULL;
    fprintf(fpFile, "* Elements : \n");
    fprintf(fpFile, "%d\n", iElementsCount);
    unsigned int iFacesCount = 6;
    Point oCenter;
    for (i = 0; i < iElementsCount; i++) {
      pvpoNodes = vpoElements[i]->GetNodes();
      iTemp = (unsigned int)pvpoNodes->size();
      // element geometry type, element physics type, element material type
      fprintf(fpFile, "2,2,1\n");
      fprintf(fpFile, "%d", pvpoNodes->at(0)->GetID());
      for (j = 1; j < iTemp; j++) {
        fprintf(fpFile, ",%d", pvpoNodes->at(j)->GetID());
      }
      fprintf(fpFile, "\n");

      sprintf(cTempString, "");
      for (j = 1; j <= iFacesCount; j++) {
        oCenter = vpoElements[i]->GetGeometry()->GetFaceCenter(j);
        if (fabs(oCenter.GetZ() - dZLength / 2.0) < dTolerance) {
          sprintf(cTempString, "%s1,1,1,", cTempString);
        } else {
          sprintf(cTempString, "%s1,1,1,", cTempString);
        }
      }
      sprintf(cTempString, "%s1,1,1", cTempString); // body forces
      fprintf(fpFile, "%s\n", cTempString);
      vpoElements[i]->DeleteGeometry();
      delete vpoElements[i];
    }
  }

  for (i = 0; i < iNodesCount; i++) {
    delete vpoNodes[i];
  }
  return true;
}
bool ParadisExternalLoadServer::GenerateThermoMechanicalBlockGeometry(
    FILE *fpFile, Home_t *poHome) const {
  Block *poBlock = new Block;
  double dXLength = poHome->param->Dimensions[0];
  double dYLength = poHome->param->Dimensions[1];
  double dZLength = poHome->param->Dimensions[2];
  double dMinDimension = min(min(dXLength, dYLength), dZLength);
  // make sure the shortest dimension has at least 5 elements and that the
  // elements aspect ratio is close to 1
  double dElementLength = 0.2 * dMinDimension;
  unsigned int iXResolution = (unsigned int)ceil(dXLength / dElementLength);
  unsigned int iYResolution = (unsigned int)ceil(dYLength / dElementLength);
  unsigned int iZResolution = (unsigned int)ceil(dZLength / dElementLength);
  // write everything in microns
  double dFactor = 1.0E6 * poHome->param->burgMag;
  dXLength = dXLength * dFactor;
  dYLength = dYLength * dFactor;
  dZLength = dZLength * dFactor;
  poBlock->SetXLength(dXLength);
  poBlock->SetYLength(dYLength);
  poBlock->SetZLength(dZLength);

  poBlock->SetResolution(iXResolution, iYResolution, iZResolution);

  vector<FEMNode *> vpoNodes;
  vector<FEMElement *> vpoElements;
  FEMMesh::GenerateMeshFromBlock(poBlock, vpoNodes, vpoElements,
                                 ThermoMechanicalFEMNode,
                                 ThermoMechanicalFEMElement);
  delete poBlock;
  unsigned int iNodesCount = (unsigned int)vpoNodes.size();
  unsigned int iElementsCount = (unsigned int)vpoElements.size();
  unsigned int i = 0;
  unsigned int iTemp = 0;
  char cTempString[1024];
  double dTolerance = 1.0E-6;

  // nodes section
  {
    fprintf(fpFile, "* Nodes : \n");
    fprintf(fpFile, "%d\n", iNodesCount);
    for (i = 0; i < iNodesCount; i++) {
      iTemp = 0;
      if (vpoNodes[i]->IsOnSurface()) {
        iTemp = 1;
      }
      fprintf(fpFile, "3\n");
      if (fabs(vpoNodes[i]->GetZ() + dZLength / 2.0) < dTolerance) {
        fprintf(fpFile, "%e,%e,%e,1,1,1,1,1,1,1,1,%d\n", vpoNodes[i]->GetX(),
                vpoNodes[i]->GetY(), vpoNodes[i]->GetZ(), iTemp);
      } else if (fabs(vpoNodes[i]->GetZ() - dZLength / 2.0) < dTolerance) {
        fprintf(fpFile, "%e,%e,%e,1,1,1,1,1,2,1,3,%d\n", vpoNodes[i]->GetX(),
                vpoNodes[i]->GetY(), vpoNodes[i]->GetZ(), iTemp);
      } else {
        fprintf(fpFile, "%e,%e,%e,0,1,0,1,0,1,0,1,%d\n", vpoNodes[i]->GetX(),
                vpoNodes[i]->GetY(), vpoNodes[i]->GetZ(), iTemp);
      }
      // initial displacements, velocities, accelerations, temperatures, heating
      // rates, ...
      fprintf(fpFile, "0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0\n");
      // forces, fluxes and stresses
      fprintf(fpFile, "0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0\n");
    }
  }

  // elements section
  {
    unsigned int j = 0;
    vector<FEMNode *> *pvpoNodes = NULL;
    fprintf(fpFile, "* Elements : \n");
    fprintf(fpFile, "%d\n", iElementsCount);
    unsigned int iFacesCount = 6;
    Point oCenter;
    for (i = 0; i < iElementsCount; i++) {
      pvpoNodes = vpoElements[i]->GetNodes();
      iTemp = (unsigned int)pvpoNodes->size();
      // element geometry type, element physics type, element material type
      fprintf(fpFile, "2,3,1\n");
      fprintf(fpFile, "%d", pvpoNodes->at(0)->GetID());
      for (j = 1; j < iTemp; j++) {
        fprintf(fpFile, ",%d", pvpoNodes->at(j)->GetID());
      }
      fprintf(fpFile, "\n");

      sprintf(cTempString, "");
      for (j = 1; j <= iFacesCount; j++) {
        oCenter = vpoElements[i]->GetGeometry()->GetFaceCenter(j);
        if (fabs(oCenter.GetZ() - dZLength / 2.0) < dTolerance) {
          sprintf(cTempString, "%s1,1,1,1,", cTempString);
        } else {
          sprintf(cTempString, "%s1,1,1,1,", cTempString);
        }
      }
      sprintf(cTempString, "%s1,1,1,1", cTempString); // body forces
      fprintf(fpFile, "%s\n", cTempString);
      vpoElements[i]->DeleteGeometry();
      delete vpoElements[i];
    }
  }

  for (i = 0; i < iNodesCount; i++) {
    delete vpoNodes[i];
  }
  return true;
}
bool ParadisExternalLoadServer::GenerateSolidCylindricalGeometry(
    FILE *fpFile, Home_t *poHome) const {

  return true;
}
bool ParadisExternalLoadServer::GenerateThermoMechanicalCylindricalGeometry(
    FILE *fpFile, Home_t *poHome) const {

  return true;
}
void ParadisExternalLoadServer::ComputeSurfaceStresses(Home_t *poHome) {
  printf("computing surface stresses\n");
  // in this function, we need to compute the stresses at all of the Gauss
  // points on the surface from all the dislocation network, which is
  // distributed over several processors here is how we will do it
  // 1. count the number of Gauss points in the system using a loop
  unsigned int iTotalGaussPointsCount = 0;
  list<FEMBoundaryElementFace *>::iterator liFaces;
  for (liFaces = m_plpoBoundaryFaces->begin();
       liFaces != m_plpoBoundaryFaces->end(); liFaces++) {
    iTotalGaussPointsCount =
        iTotalGaussPointsCount + (*liFaces)->GetGaussPointsCount();
  }
  // 2. allocate an array of size 6*gauss points count to hold the stresses
  double *pdLocalStresses = new double[6 * iTotalGaussPointsCount]();
  // 3. compute the stresses from the local dislocation network for all of the
  // gauss point and store it in the local array
  unsigned int i = 0;
  unsigned int j = 0;
  unsigned int k = 0;
  Node_t *poNode = NULL;
  Node_t *poNeighbour = NULL;
  unsigned int iGaussPointsCount = 0;
  unsigned int iNodesCount = poHome->newNodeKeyPtr;
  vector<Point> *pvpoGaussPoints = NULL;
  Matrix oStress(3, 3);
  double dMu = poHome->param->shearModulus;
  double dNu = poHome->param->pois;
  Point oStart;
  Point oEnd;
  Vector oBurgersVector;
  Point oPoint;
  unsigned int iCount = 0;
  // the coordinates of the points in FEM are in microns, multiply them by a
  // factor to get them in b units
  double dCoordinatesFactor = 1.0E-6 / poHome->param->burgMag;
  // loop over the boundary elements
  for (liFaces = m_plpoBoundaryFaces->begin();
       liFaces != m_plpoBoundaryFaces->end(); liFaces++) {
    pvpoGaussPoints = (*liFaces)->GetGaussPoints();
    iGaussPointsCount = pvpoGaussPoints->size();
    // loop over the Gauss points
    for (i = 0; i < iGaussPointsCount; i++) {
      oStress.ZeroValues();
      oPoint = pvpoGaussPoints->at(i);
      oPoint = oPoint * dCoordinatesFactor;
      // loop over the dislocation nodes
      for (j = 0; j < iNodesCount; j++) {
        poNode = poHome->nodeKeys[j];
        if (poNode == NULL) {
          continue;
        }
        oStart.Set(poNode->x, poNode->y, poNode->z);
        // loop over the dislocation segments
        for (k = 0; k < poNode->numNbrs; k++) {
          // make sure that the segment is owned by this domain
          if (poNode->myTag.domainID < poNode->nbrTag[k].domainID) {
            continue;
          }
          if (poNode->myTag.domainID == poNode->nbrTag[k].domainID) {
            if (poNode->myTag.index < poNode->nbrTag[k].index) {
              continue;
            }
          }
          poNeighbour = GetNeighborNode(poHome, poNode, k);
          oEnd.Set(poNeighbour->x, poNeighbour->y, poNeighbour->z);
          oBurgersVector.Set(poNode->burgX[k], poNode->burgY[k],
                             poNode->burgZ[k]);
          // compute the stress field at from this segment
          oStress = oStress + ComputeSegmentStressAtPoint(oStart, oEnd,
                                                          oBurgersVector,
                                                          oPoint, dMu, dNu);
          // the stress in paradis is in Pa, the stress in FEM is in the units
          // of kg/micron sec^2, divide this stress by 10^6 before sending it to
          // FEM
          oStress = oStress * (1.0E-6);
        }
      }
      // now we have the total stress at this Gauss point, pass it to the FEM
      // boundary element face object
      pdLocalStresses[iCount++] = oStress.Get(1, 1);
      pdLocalStresses[iCount++] = oStress.Get(1, 2);
      pdLocalStresses[iCount++] = oStress.Get(1, 3);
      pdLocalStresses[iCount++] = oStress.Get(2, 2);
      pdLocalStresses[iCount++] = oStress.Get(2, 3);
      pdLocalStresses[iCount++] = oStress.Get(3, 3);
    }
  }
  // 4. mpi_all_reduce sum the stresses
  double *pdTotalStresses = new double[6 * iTotalGaussPointsCount]();
  DSMPI::AllReduce(pdLocalStresses, pdTotalStresses, 6 * iTotalGaussPointsCount,
                   MPI_DOUBLE, MPI_SUM);
  // 5. set these stresses to the fem boundary element face objects
  iCount = 0;
  for (liFaces = m_plpoBoundaryFaces->begin();
       liFaces != m_plpoBoundaryFaces->end(); liFaces++) {
    pvpoGaussPoints = (*liFaces)->GetGaussPoints();
    iGaussPointsCount = pvpoGaussPoints->size();
    // loop over the Gauss points
    for (i = 0; i < iGaussPointsCount; i++) {
      oStress.Set(1, 1, pdTotalStresses[iCount++]);
      oStress.Set(1, 2, pdTotalStresses[iCount++]);
      oStress.Set(1, 3, pdTotalStresses[iCount++]);
      oStress.Set(2, 2, pdTotalStresses[iCount++]);
      oStress.Set(2, 3, pdTotalStresses[iCount++]);
      oStress.Set(3, 3, pdTotalStresses[iCount++]);
      oStress.Set(2, 1, oStress.Get(1, 2));
      oStress.Set(3, 1, oStress.Get(1, 3));
      oStress.Set(3, 2, oStress.Get(3, 2));
      (*liFaces)->SetGaussPointStress(i, oStress);
    }
  }
  delete[] pdLocalStresses;
  delete[] pdTotalStresses;
  DSMPI::Barrier();
}
Matrix ParadisExternalLoadServer::ComputeSegmentStressAtPoint(
    const Point &oStart, const Point &oEnd, const Vector &oBurgersVector,
    const Point &oPoint, const double &dMu, const double &dNu) {
  // compute the stress field from a dislocation segment at a point
  // the general 3D stress field for a mixed dislocation was used in this
  // function the segment is always straight, which simplifies the evaluation of
  // the integrals the integration is exact. the expression of the integrals
  // were obtained from Wolfram integrator
  Matrix oStress(3, 3);
  double dTolerance = 1.0E-6;
  double dToleranceSquared = 1.0E-12;

  unsigned int i = 0;

  double dY[3] = {0.0, 0.0, 0.0};
  double dQ[3] = {0.0, 0.0, 0.0};
  double dP[3] = {0.0, 0.0, 0.0};
  double dM[3] = {0.0, 0.0, 0.0};

  dY[0] = oPoint.GetX();
  dY[1] = oPoint.GetY();
  dY[2] = oPoint.GetZ();
  dQ[0] = oStart.GetX();
  dQ[1] = oStart.GetY();
  dQ[2] = oStart.GetZ();
  dP[0] = oEnd.GetX();
  dP[1] = oEnd.GetY();
  dP[2] = oEnd.GetZ();

  double dC = 0.0;
  for (i = 0; i < 3; i++) {
    dM[i] = dP[i] - dQ[i];
    dC = dC + dM[i] * dM[i];
  }
  if (dC < dToleranceSquared) {
    // zero length segment, return
    return oStress;
  }
  // get the actual length
  double dL = sqrt(dC);
  // if the point lies on the dislocation line or its extension, move it a
  // little bit outwards (to a normal distance of 1b) to avoid the singularity
  // in the stress calculation. This is the only approximation used in this
  // function, everything else is exact
  double dProjection = 0.0;
  for (i = 0; i < 3; i++) {
    dM[i] = dM[i] / dL;
    dProjection = dProjection + (dY[i] - dQ[i]) * dM[i];
  }
  double dTest = 0.0;
  double dNormalDirection[3] = {0.0, 0.0, 0.0};
  for (i = 0; i < 3; i++) {
    dNormalDirection[i] = dY[i] - dQ[i] - dProjection * dM[i];
    dTest = dTest + dNormalDirection[i] * dNormalDirection[i];
  }
  dTest = sqrt(dTest);
  double dTubeRadius = 0.1;
  // dTest is the the normal distance between the point and the line, if it is
  // zero, return the point location so that it lands on the 1b radius cylinder
  // with the dislocation being the axis
  Vector oLineDirection(dM[0], dM[1], dM[2]);
  Vector oRandomDirection;
  if (fabs(dTest) < dTolerance) {
    while (true) {
      oRandomDirection.Set(Randomizer::Random(-1, 1), Randomizer::Random(-1, 1),
                           Randomizer::Random(-1, 1));
      oRandomDirection.Normalize();
      oRandomDirection = oRandomDirection ^ oLineDirection;
      if (oRandomDirection.Length() > dTolerance) {
        oRandomDirection.Normalize();
        break;
      }
    }
    dY[0] = dY[0] + (dTubeRadius - dTest) * oRandomDirection.GetX();
    dY[1] = dY[1] + (dTubeRadius - dTest) * oRandomDirection.GetY();
    dY[2] = dY[2] + (dTubeRadius - dTest) * oRandomDirection.GetZ();
  }
  // if it is a small nonzero number, move the point outwards to 1b away
  else if (dTest < dTubeRadius) {
    dY[0] = dY[0] + (dTubeRadius - dTest) * dNormalDirection[0] / dTest;
    dY[1] = dY[1] + (dTubeRadius - dTest) * dNormalDirection[1] / dTest;
    dY[2] = dY[2] + (dTubeRadius - dTest) * dNormalDirection[2] / dTest;
  }

  double dA = 0.0;
  double dB = 0.0;
  double dG[3] = {0.0, 0.0, 0.0};
  double dH[3] = {0.0, 0.0, 0.0};
  for (i = 0; i < 3; i++) {
    dG[i] = dY[i] - dQ[i];
    dA = dA + dG[i] * dG[i];
    dB = dB - dM[i] * dG[i];
    dH[i] = -dM[i] * dL;
  }
  dB = 2.0 * dL * dB;

  // first, evaluate the expressions needed for the r,ppi  and second term of
  // the r,ijk integration
  double dD[3] = {0.0, 0.0, 0.0};
  double dE[3] = {0.0, 0.0, 0.0};
  for (i = 0; i < 3; i++) {
    dD[i] = 2.0 * (2.0 * dA * dH[i] - dB * dG[i]);
    dE[i] = 2.0 * (dB * dH[i] - 2.0 * dC * dG[i]);
  }

  double dI[3] = {0.0, 0.0, 0.0};
  double dT = dB * dB - 4.0 * dA * dC;
  double dR = sqrt(dA + dB + dC);
  double dS = sqrt(dA);
  for (i = 0; i < 3; i++) {
    dI[i] = (dD[i] + dE[i]) / dR - dD[i] / dS;
    dI[i] = dI[i] / dT;
  }

  // second, evaluate the expressions needed for the first term of the r,ijk
  // integration
  unsigned int iIndexI[10] = {1, 1, 2, 3, 1, 1, 2, 2, 3, 3};
  unsigned int iIndexJ[10] = {2, 1, 2, 3, 1, 1, 2, 2, 3, 3};
  unsigned int iIndexK[10] = {3, 1, 2, 3, 2, 3, 1, 3, 1, 2};
  unsigned int iReverseIndex[3][3][3] = {{{1, 4, 5}, {4, 6, 0}, {5, 0, 8}},
                                         {{4, 6, 0}, {6, 2, 7}, {0, 7, 9}},
                                         {{5, 0, 8}, {0, 7, 9}, {8, 9, 3}}};

  double dF[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double dN[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double dU[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double dV[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  unsigned int j = 0;
  unsigned int k = 0;
  unsigned int iIndex = 0;
  for (iIndex = 0; iIndex < 10; iIndex++) {
    i = iIndexI[iIndex] - 1;
    j = iIndexJ[iIndex] - 1;
    k = iIndexK[iIndex] - 1;
    dF[iIndex] = dG[i] * dG[j] * dG[k];
    dN[iIndex] =
        dG[i] * dG[j] * dH[k] + dG[i] * dG[k] * dH[j] + dG[k] * dG[j] * dH[i];
    dU[iIndex] =
        dG[i] * dH[j] * dH[k] + dG[j] * dH[k] * dH[i] + dG[k] * dH[j] * dH[i];
    dV[iIndex] = dH[i] * dH[j] * dH[k];
  }

  double dK = 0.0;
  double dO = 0.0;
  double dJ[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double dKDen = 2.0 / (3.0 * dT * dT * dR * dR * dR);
  double dODen = -2.0 / (3.0 * dT * dT * dS * dS * dS);
  double dTemp1 = 0.0;
  double dTemp2 = 0.0;
  double dTemp3 = 0.0;
  double dTemp4 = 0.0;
  for (i = 0; i < 10; i++) {
    dTemp1 = dB * dB * dB * (-dF[i] - 3.0 * dN[i] + 3.0 * dU[i] + dV[i]);
    dTemp2 = -8.0 * (2.0 * dA * dA * dA * dV[i] - 2.0 * dC * dC * dC * dF[i] -
                     dA * dC * dC * (3.0 * dF[i] + dU[i]) +
                     dA * dA * dC * (dN[i] + 3.0 * dV[i]));
    dTemp3 = -2.0 * dB * dB * (-(dC * (3.0 * dF[i] - 6.0 * dN[i] + dU[i])) +
                               dA * (dN[i] - 6.0 * dU[i] + 3.0 * dV[i]));
    dTemp4 = 4.0 * dB * (-2.0 * dC * dC * (-3.0 * dF[i] + dN[i]) +
                         2.0 * dA * dA * (dU[i] - 3.0 * dV[i]) +
                         3.0 * dA * dC * (dF[i] - dN[i] + dU[i] - dV[i]));
    dK = (dTemp1 + dTemp2 + dTemp3 + dTemp4) * dKDen;
    dTemp1 = dB * dB * dB * dF[i];
    dTemp2 = 8.0 * (2.0 * dA * dA * dA * dV[i] + dA * dA * dC * dN[i]);
    dTemp3 = 2.0 * dB * dB * dA * dN[i];
    dTemp4 = -4.0 * dB * (2.0 * dA * dA * dU[i] + 3.0 * dA * dC * dF[i]);
    dO = (dTemp1 + dTemp2 + dTemp3 + dTemp4) * dODen;
    dJ[i] = dK - dO;
  }

  // now calculate the stress
  double dBurgers[3] = {0.0, 0.0, 0.0};
  dBurgers[0] = oBurgersVector.GetX();
  dBurgers[1] = oBurgersVector.GetY();
  dBurgers[2] = oBurgersVector.GetZ();
  unsigned int iAlpha = 0;
  unsigned int iBeta = 0;
  double dNuFactor = 1.0 / (1.0 - dNu);
  for (iAlpha = 1; iAlpha <= 3; iAlpha++) {
    for (iBeta = iAlpha; iBeta <= 3; iBeta++) {
      // compute first term
      dTemp1 = 0.0;
      for (i = 0; i < 3; i++) {
        dTemp2 = 0.0;
        for (j = 0; j < 3; j++) {
          dTemp2 = dTemp2 -
                   dBurgers[j] *
                       (MathServices::PermutationSymbol(i + 1, j + 1, iAlpha) *
                            dM[iBeta - 1] +
                        MathServices::PermutationSymbol(i + 1, j + 1, iBeta) *
                            dM[iAlpha - 1]);
        }
        dTemp1 = dTemp1 + dTemp2 * dI[i];
      }
      // compute second term
      dTemp4 = 0.0;
      for (i = 0; i < 3; i++) {
        dTemp2 = 0.0;
        for (j = 0; j < 3; j++) {
          dTemp3 = 0.0;
          for (k = 0; k < 3; k++) {
            dTemp3 =
                dTemp3 +
                MathServices::PermutationSymbol(i + 1, j + 1, k + 1) * dM[k];
          }
          dTemp2 = dTemp2 + dBurgers[j] * dTemp3;
        }
        dTemp4 =
            dTemp4 +
            dTemp2 *
                (MathServices::KroneckerDelta(iAlpha, iBeta) * dI[i] -
                 MathServices::KroneckerDelta(iAlpha, i + 1) * dI[iBeta - 1] -
                 MathServices::KroneckerDelta(iBeta, i + 1) * dI[iAlpha - 1] +
                 3.0 * dJ[iReverseIndex[i][iAlpha - 1][iBeta - 1]]);
      }
      oStress.Set(iAlpha, iBeta, dTemp1 + dNuFactor * dTemp4);
    }
  }
  oStress = oStress * (dMu * dL / 4.0 / PI);
  oStress.Set(2, 1, oStress.Get(1, 2));
  oStress.Set(3, 1, oStress.Get(1, 3));
  oStress.Set(3, 2, oStress.Get(2, 3));
  return oStress;
}
