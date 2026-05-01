// Ahmed M. Hussein
#include "FEMImplicitDynamicsThermoMechanicalSolver.h"
#include "FEMThermoMechanicalElement.h"
#include "FEMThermoMechanicalNode.h"
#include "Tools.h"
#include "cmath"


namespace FEMSystem
{
	FEMImplicitDynamicsThermoMechanicalSolver::FEMImplicitDynamicsThermoMechanicalSolver()
	{
		Initialize();
	}
	FEMImplicitDynamicsThermoMechanicalSolver::FEMImplicitDynamicsThermoMechanicalSolver(const FEMImplicitDynamicsThermoMechanicalSolver& oSolver)
	{
		*this = oSolver;
	}
	FEMImplicitDynamicsThermoMechanicalSolver::~FEMImplicitDynamicsThermoMechanicalSolver()
	{
		Reset();
	}
	FEMImplicitDynamicsThermoMechanicalSolver& FEMImplicitDynamicsThermoMechanicalSolver::operator=(const FEMImplicitDynamicsThermoMechanicalSolver& oSolver)
	{
		FEMSolver::operator=(oSolver);
		m_oForces1 = oSolver.m_oForces1;
		m_oForces2 = oSolver.m_oForces2;
		m_oCurrentDisplacements1 = oSolver.m_oCurrentDisplacements1;
		m_oCurrentDisplacements2 = oSolver.m_oCurrentDisplacements2;
		m_oPreviousDisplacements1 = oSolver.m_oPreviousDisplacements1;
		m_oPreviousDisplacements2 = oSolver.m_oPreviousDisplacements2;
		m_oCurrentVelocities1 = oSolver.m_oCurrentVelocities1;
		m_oCurrentVelocities2 = oSolver.m_oCurrentVelocities2;
		m_oPreviousVelocities1 = oSolver.m_oPreviousVelocities1;
		m_oPreviousVelocities2 = oSolver.m_oPreviousVelocities2;
		m_oCurrentAccelerations1 = oSolver.m_oCurrentAccelerations1;
		m_oCurrentAccelerations2 = oSolver.m_oCurrentAccelerations2;
		m_oPreviousAccelerations1 = oSolver.m_oPreviousAccelerations1;
		m_oPreviousAccelerations2 = oSolver.m_oPreviousAccelerations2;
		
		m_viUnknownDisplacementsIndices = oSolver.m_viUnknownDisplacementsIndices;
		m_viKnownDisplacementsIndices = oSolver.m_viKnownDisplacementsIndices;
		m_viUnknownDisplacementsIndicesReverseMap = oSolver.m_viUnknownDisplacementsIndicesReverseMap;
		m_viKnownDisplacementsIndicesReverseMap = oSolver.m_viKnownDisplacementsIndicesReverseMap;
		
		m_oDynamicMatrix11 = oSolver.m_oDynamicMatrix11;
		m_oDynamicMatrix12 = oSolver.m_oDynamicMatrix12;
		m_oDynamicMatrix21 = oSolver.m_oDynamicMatrix21;
		m_oDynamicMatrix22 = oSolver.m_oDynamicMatrix22;
		
		m_oMassMatrix11 = oSolver.m_oMassMatrix11;
		m_oMassMatrix12 = oSolver.m_oMassMatrix12;
		m_oMassMatrix21 = oSolver.m_oMassMatrix21;
		m_oMassMatrix22 = oSolver.m_oMassMatrix22;
		
		m_oDampingMatrix11 = oSolver.m_oDampingMatrix11;
		m_oDampingMatrix12 = oSolver.m_oDampingMatrix12;
		m_oDampingMatrix21 = oSolver.m_oDampingMatrix21;
		m_oDampingMatrix22 = oSolver.m_oDampingMatrix22;
		return *this;
	}
	void FEMImplicitDynamicsThermoMechanicalSolver::Reset()
	{
		FEMSolver::Reset();
		Initialize();
	}
	void FEMImplicitDynamicsThermoMechanicalSolver::Initialize()
	{
		FEMSolver::Initialize();
		m_oForces1.Reset();
		m_oForces2.Reset();
		m_oCurrentDisplacements1.Reset();
		m_oCurrentDisplacements2.Reset();
		m_oPreviousDisplacements1.Reset();
		m_oPreviousDisplacements2.Reset();
		m_oCurrentVelocities1.Reset();
		m_oCurrentVelocities2.Reset();
		m_oPreviousVelocities1.Reset();
		m_oPreviousVelocities2.Reset();
		m_oCurrentAccelerations1.Reset();
		m_oCurrentAccelerations2.Reset();
		m_oPreviousAccelerations1.Reset();
		m_oPreviousAccelerations2.Reset();
		
		m_viUnknownDisplacementsIndices.clear();
		m_viKnownDisplacementsIndices.clear();
		m_viUnknownDisplacementsIndicesReverseMap.clear();
		m_viKnownDisplacementsIndicesReverseMap.clear();
		
		m_oDynamicMatrix11.Reset();
		m_oDynamicMatrix12.Reset();
		m_oDynamicMatrix21.Reset();
		m_oDynamicMatrix22.Reset();
		
		m_oMassMatrix11.Reset();
		m_oMassMatrix12.Reset();
		m_oMassMatrix21.Reset();
		m_oMassMatrix22.Reset();
		
		m_oDampingMatrix11.Reset();
		m_oDampingMatrix12.Reset();
		m_oDampingMatrix21.Reset();
		m_oDampingMatrix22.Reset();
	}		
	void FEMImplicitDynamicsThermoMechanicalSolver::InitializeMatrices(const double& dTimeStep)
	{
		PrintOnScreen("Initializing FEM problem");
		vector<FEMNode*>* pvpoNodes = m_poDataStructure->GetNodes();
 		unsigned int iSize = pvpoNodes->size();
 		unsigned int i = 0;
 		// get the exact node count
 		unsigned int iDOFCount = 0;
 		for(i = 0; i < iSize ; i++)
 		{
 			iDOFCount = iDOFCount + pvpoNodes->at(i)->GetDOFCount();
 		}
 		m_viUnknownDisplacementsIndices.reserve(iDOFCount);
 		m_viKnownDisplacementsIndices.reserve(iDOFCount);
 		// the maximum DOF index cannot be higher than the number of DOFs
 		m_viUnknownDisplacementsIndicesReverseMap.resize(iDOFCount);
 		m_viKnownDisplacementsIndicesReverseMap.resize(iDOFCount);
 		FEMThermoMechanicalNode* poNode = NULL;
 		unsigned int iIndex = 0;
 		unsigned int iKnownIndicesCount = 0;
 		unsigned int iUnknownIndicesCount = 0;
 		
 		for(i = 0; i < iSize ; i++)
 		{
 			poNode = (FEMThermoMechanicalNode*)pvpoNodes->at(i);
 			
 			iIndex = poNode->GetXDOF()->GetIndex();
 			if(poNode->GetXDOF()->GetCondition())
 			{
 				iKnownIndicesCount = iKnownIndicesCount + 1;
 				m_viKnownDisplacementsIndices.push_back(iIndex);
 				m_viKnownDisplacementsIndicesReverseMap[iIndex - 1] = iKnownIndicesCount;
 				m_viUnknownDisplacementsIndicesReverseMap[iIndex - 1] = 0;
 			}
 			else
 			{
 				iUnknownIndicesCount = iUnknownIndicesCount + 1;
 				m_viUnknownDisplacementsIndices.push_back(iIndex);
 				m_viUnknownDisplacementsIndicesReverseMap[iIndex - 1] = iUnknownIndicesCount;
 				m_viKnownDisplacementsIndicesReverseMap[iIndex - 1] = 0;
 			}

			iIndex = poNode->GetYDOF()->GetIndex();
 			if(poNode->GetYDOF()->GetCondition())
 			{
 				iKnownIndicesCount = iKnownIndicesCount + 1;
 				m_viKnownDisplacementsIndices.push_back(iIndex);
 				m_viKnownDisplacementsIndicesReverseMap[iIndex - 1] = iKnownIndicesCount;
 				m_viUnknownDisplacementsIndicesReverseMap[iIndex - 1] = 0;
 			}
 			else
 			{
 				iUnknownIndicesCount = iUnknownIndicesCount + 1;
 				m_viUnknownDisplacementsIndices.push_back(iIndex);
 				m_viUnknownDisplacementsIndicesReverseMap[iIndex - 1] = iUnknownIndicesCount;
 				m_viKnownDisplacementsIndicesReverseMap[iIndex - 1] = 0;
 			}
 			
 			iIndex = poNode->GetZDOF()->GetIndex();
 			if(poNode->GetZDOF()->GetCondition())
 			{
 				iKnownIndicesCount = iKnownIndicesCount + 1;
 				m_viKnownDisplacementsIndices.push_back(iIndex);
 				m_viKnownDisplacementsIndicesReverseMap[iIndex - 1] = iKnownIndicesCount;
 				m_viUnknownDisplacementsIndicesReverseMap[iIndex - 1] = 0;
 			}
 			else
 			{
 				iUnknownIndicesCount = iUnknownIndicesCount + 1;
 				m_viUnknownDisplacementsIndices.push_back(iIndex);
 				m_viUnknownDisplacementsIndicesReverseMap[iIndex - 1] = iUnknownIndicesCount;
 				m_viKnownDisplacementsIndicesReverseMap[iIndex - 1] = 0;
 			}
 			
 			iIndex = poNode->GetTDOF()->GetIndex();
 			if(poNode->GetTDOF()->GetCondition())
 			{
 				iKnownIndicesCount = iKnownIndicesCount + 1;
 				m_viKnownDisplacementsIndices.push_back(iIndex);
 				m_viKnownDisplacementsIndicesReverseMap[iIndex - 1] = iKnownIndicesCount;
 				m_viUnknownDisplacementsIndicesReverseMap[iIndex - 1] = 0;
 			}
 			else
 			{
 				iUnknownIndicesCount = iUnknownIndicesCount + 1;
 				m_viUnknownDisplacementsIndices.push_back(iIndex);
 				m_viUnknownDisplacementsIndicesReverseMap[iIndex - 1] = iUnknownIndicesCount;
 				m_viKnownDisplacementsIndicesReverseMap[iIndex - 1] = 0;
 			}
 		}
 	 	
		// set the Newmark method parameters
		double dDelta = 0.5;
		double dAlpha = 0.25;
		m_dTimeStep = dTimeStep;
		m_dA0 = 1.0/dAlpha/dTimeStep/dTimeStep;
		m_dA1 = dDelta/dAlpha/dTimeStep;
		m_dA2 = 1.0/dAlpha/dTimeStep;
		m_dA3 = 0.5/dAlpha - 1.0;
		m_dA4 = dDelta/dAlpha - 1.0;
		m_dA5 = 0.5*dTimeStep*(m_dA4 - 1.0);
		m_dA6 = dTimeStep*(1.0 - dDelta);
		m_dA7 = dDelta*dTimeStep;
		GenerateMatrices(iDOFCount);
		// now resize the force, displacement, velocity and acceleration arrays
		unsigned int iSize1 = (unsigned int)m_viUnknownDisplacementsIndices.size();
		unsigned int iSize2 = (unsigned int)m_viKnownDisplacementsIndices.size();
		m_oForces1.SetSize(iSize1,1);
		m_oForces2.SetSize(iSize2,1);
		m_oCurrentDisplacements1.SetSize(iSize1,1);
		m_oCurrentDisplacements2.SetSize(iSize2,1);
		m_oPreviousDisplacements1.SetSize(iSize1,1);
		m_oPreviousDisplacements2.SetSize(iSize2,1);
		m_oCurrentVelocities1.SetSize(iSize1,1);
		m_oCurrentVelocities2.SetSize(iSize2,1);
		m_oPreviousVelocities1.SetSize(iSize1,1);
		m_oPreviousVelocities2.SetSize(iSize2,1);
		m_oCurrentAccelerations1.SetSize(iSize1,1);
		m_oCurrentAccelerations2.SetSize(iSize2,1);
		m_oPreviousAccelerations1.SetSize(iSize1,1);
		m_oPreviousAccelerations2.SetSize(iSize2,1);


		// initialize the displacements, velocities, acceleratios and forces
		FEMDegreeOfFreedom* poDOF = NULL;
		Vector oDisplacement;
		Vector oVelocity;
		Vector oAcceleration;
		Vector oForce;
		for(i = 0; i < iSize ; i++)
		{
			poNode = (FEMThermoMechanicalNode*)pvpoNodes->at(i);
			oDisplacement = poNode->GetDisplacement();
			oVelocity = poNode->GetVelocity();
			oAcceleration = poNode->GetAcceleration();
			oForce = poNode->GetForce();
			
			poDOF = poNode->GetXDOF();
			iIndex = poDOF->GetIndex() - 1;
			if(poDOF->GetCondition())
			{
				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
				m_oCurrentDisplacements2.Set(iIndex,1,oDisplacement.GetX());
				m_oCurrentVelocities2.Set(iIndex,1,oVelocity.GetX());
				m_oCurrentAccelerations2.Set(iIndex,1,oAcceleration.GetX());
				m_oForces2.Set(iIndex,1,oForce.GetX());
			}
			else
			{
				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
				m_oCurrentDisplacements1.Set(iIndex,1,oDisplacement.GetX());
				m_oCurrentVelocities1.Set(iIndex,1,oVelocity.GetX());
				m_oCurrentAccelerations1.Set(iIndex,1,oAcceleration.GetX());
				m_oForces1.Set(iIndex,1,oForce.GetX());
			}
	
			poDOF = poNode->GetYDOF();
			iIndex = poDOF->GetIndex() - 1;
			if(poDOF->GetCondition())
			{
				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
				m_oCurrentDisplacements2.Set(iIndex,1,oDisplacement.GetY());
				m_oCurrentVelocities2.Set(iIndex,1,oVelocity.GetY());
				m_oCurrentAccelerations2.Set(iIndex,1,oAcceleration.GetY());
				m_oForces2.Set(iIndex,1,oForce.GetY());
			}
			else
			{
				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
				m_oCurrentDisplacements1.Set(iIndex,1,oDisplacement.GetY());
				m_oCurrentVelocities1.Set(iIndex,1,oVelocity.GetY());
				m_oCurrentAccelerations1.Set(iIndex,1,oAcceleration.GetY());
				m_oForces1.Set(iIndex,1,oForce.GetY());
			}
	
			poDOF = poNode->GetZDOF();
			iIndex = poDOF->GetIndex() - 1;
			if(poDOF->GetCondition())
			{
				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
				m_oCurrentDisplacements2.Set(iIndex,1,oDisplacement.GetZ());
				m_oCurrentVelocities2.Set(iIndex,1,oVelocity.GetZ());
				m_oCurrentAccelerations2.Set(iIndex,1,oAcceleration.GetZ());
				m_oForces2.Set(iIndex,1,oForce.GetZ());
			}
			else
			{
				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
				m_oCurrentDisplacements1.Set(iIndex,1,oDisplacement.GetZ());
				m_oCurrentVelocities1.Set(iIndex,1,oVelocity.GetZ());
				m_oCurrentAccelerations1.Set(iIndex,1,oAcceleration.GetZ());
				m_oForces1.Set(iIndex,1,oForce.GetZ());
			}
			
			poDOF = poNode->GetTDOF();
			iIndex = poDOF->GetIndex() - 1;
			if(poDOF->GetCondition())
			{
				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
				m_oCurrentDisplacements2.Set(iIndex,1,poDOF->GetPrimaryValue());
				m_oCurrentVelocities2.Set(iIndex,1,poNode->GetHeatingRate());
				m_oCurrentAccelerations2.Set(iIndex,1,poNode->GetHeatingAcceleration());
				m_oForces2.Set(iIndex,1,poDOF->GetSecondaryValue());
			}
			else
			{
				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
				m_oCurrentDisplacements1.Set(iIndex,1,poDOF->GetPrimaryValue());
				m_oCurrentVelocities1.Set(iIndex,1,poNode->GetHeatingRate());
				m_oCurrentAccelerations1.Set(iIndex,1,poNode->GetHeatingAcceleration());
				m_oForces1.Set(iIndex,1,poDOF->GetSecondaryValue());
			}
		}
		PrintOnScreen("Done initializing");
	}
	void FEMImplicitDynamicsThermoMechanicalSolver::Solve(const double& dTime)
	{
		PrintOnScreen("Generating nodal forces and displacements");
		GenerateNodalStateVectors();
		PrintOnScreen("Done");
		PrintOnScreen("Solving the system");
		// get the unknown displacements
		Matrix oMRHS1 = m_oPreviousDisplacements1*m_dA0 + m_oPreviousVelocities1*m_dA2 + m_oPreviousAccelerations1*m_dA3;
		Matrix oMRHS2 = m_oPreviousDisplacements2*m_dA0 + m_oPreviousVelocities2*m_dA2 + m_oPreviousAccelerations2*m_dA3;
		Matrix oCRHS1 = m_oPreviousDisplacements1*m_dA1 + m_oPreviousVelocities1*m_dA4 + m_oPreviousAccelerations1*m_dA5;
		Matrix oCRHS2 = m_oPreviousDisplacements2*m_dA1 + m_oPreviousVelocities2*m_dA4 + m_oPreviousAccelerations2*m_dA5;
		Matrix oRHS = m_oForces1 - m_oDynamicMatrix12*m_oCurrentDisplacements2 + m_oMassMatrix11*oMRHS1 + m_oMassMatrix12*oMRHS2 + m_oDampingMatrix11*oCRHS1 + m_oDampingMatrix12*oCRHS2;
		m_oCurrentDisplacements1 = m_oDynamicMatrix11.SolveGMRES(oRHS,&m_oCurrentDisplacements1);
		// get the unknown forces
		m_oForces2 = m_oDynamicMatrix21*m_oCurrentDisplacements1 + m_oDynamicMatrix22*m_oCurrentDisplacements2;
		m_oForces2 = m_oForces2 - m_oMassMatrix21*oMRHS1 - m_oMassMatrix22*oMRHS2 - m_oDampingMatrix21*oCRHS1 - m_oDampingMatrix22*oCRHS2;
		// get the accelerations
		m_oCurrentAccelerations1 = (m_oCurrentDisplacements1 - m_oPreviousDisplacements1)*m_dA0;
		m_oCurrentAccelerations1 = m_oCurrentAccelerations1 - m_oPreviousVelocities1*m_dA2;
		m_oCurrentAccelerations1 = m_oCurrentAccelerations1 - m_oPreviousAccelerations1*m_dA3;
		m_oCurrentAccelerations2 = (m_oCurrentDisplacements2 - m_oPreviousDisplacements2)*m_dA0;
		m_oCurrentAccelerations2 = m_oCurrentAccelerations2 - m_oPreviousVelocities2*m_dA2;
		m_oCurrentAccelerations2 = m_oCurrentAccelerations2 - m_oPreviousAccelerations2*m_dA3;
		// get the velocities
		m_oCurrentVelocities1 = m_oPreviousVelocities1 + m_oPreviousAccelerations1*m_dA6 + m_oCurrentAccelerations1*m_dA7;
		m_oCurrentVelocities2 = m_oPreviousVelocities2 + m_oPreviousAccelerations2*m_dA6 + m_oCurrentAccelerations2*m_dA7;
		UpdateNodalStateVectors();
		PrintOnScreen("Done");
		PrintOnScreen("Calculating nodal stresses and fluxes");
		GenerateNodalStressesMatrix();
		GenerateNodalFluxesMatrix();
		PrintOnScreen("Done");
	}
	Matrix FEMImplicitDynamicsThermoMechanicalSolver::GetStress(Point* poPoint,unsigned int& iStatus)
	{
		Matrix oStress(3,3);
		vector<double> vdCoordinates;
		FEMElement* poElement = NULL;
		iStatus = DeterminePointLocation(poPoint,poElement,vdCoordinates);
		Matrix oNaturalCoordinates(3,1);
		if(iStatus != 0)
		{
			oNaturalCoordinates.Set(1,1,vdCoordinates[0]);
			oNaturalCoordinates.Set(2,1,vdCoordinates[1]);
			oNaturalCoordinates.Set(3,1,vdCoordinates[2]);
			oStress = ((FEMThermoMechanicalElement*)poElement)->GetStress(oNaturalCoordinates);
		}
		return oStress;
	}
	Vector FEMImplicitDynamicsThermoMechanicalSolver::GetDisplacement(Point* poPoint,unsigned int& iStatus)
	{
		Vector oDisplacement(0.0,0.0,0.0);
		vector<double> vdCoordinates;
		FEMElement* poElement = NULL;
		iStatus = DeterminePointLocation(poPoint,poElement,vdCoordinates);
		Matrix oNaturalCoordinates(3,1);
		if(iStatus != 0)
		{
			oNaturalCoordinates.Set(1,1,vdCoordinates[0]);
			oNaturalCoordinates.Set(2,1,vdCoordinates[1]);
			oNaturalCoordinates.Set(3,1,vdCoordinates[2]);
			oDisplacement = ((FEMThermoMechanicalElement*)poElement)->GetDisplacement(oNaturalCoordinates);
		}
		return oDisplacement;
	}
	double FEMImplicitDynamicsThermoMechanicalSolver::GetTemperature(Point* poPoint,unsigned int& iStatus)
	{
		double dTemperature = 0.0;
		vector<double> vdCoordinates;
		FEMElement* poElement = NULL;
		iStatus = DeterminePointLocation(poPoint,poElement,vdCoordinates);
		Matrix oNaturalCoordinates(3,1);
		if(iStatus != 0)
		{
			oNaturalCoordinates.Set(1,1,vdCoordinates[0]);
			oNaturalCoordinates.Set(2,1,vdCoordinates[1]);
			oNaturalCoordinates.Set(3,1,vdCoordinates[2]);
			dTemperature = ((FEMThermoMechanicalElement*)poElement)->GetTemperature(oNaturalCoordinates);
		}
		return dTemperature;
	}
	Vector FEMImplicitDynamicsThermoMechanicalSolver::GetFlux(Point* poPoint,unsigned int& iStatus)
	{
		Vector oFlux(0.0,0.0,0.0);
		vector<double> vdCoordinates;
		FEMElement* poElement = NULL;
		iStatus = DeterminePointLocation(poPoint,poElement,vdCoordinates);
		Matrix oNaturalCoordinates(3,1);
		if(iStatus != 0)
		{
			oNaturalCoordinates.Set(1,1,vdCoordinates[0]);
			oNaturalCoordinates.Set(2,1,vdCoordinates[1]);
			oNaturalCoordinates.Set(3,1,vdCoordinates[2]);
			oFlux = ((FEMThermoMechanicalElement*)poElement)->GetFlux(oNaturalCoordinates);
		}
		return oFlux;
	}
	void FEMImplicitDynamicsThermoMechanicalSolver::WriteFEMSolutionToParaview(const unsigned int& iStep) const
	{
		char cWrite[512];
		string sOutputBaseFileName = m_poDataStructure->GetOutputBaseFileName();
		sprintf(cWrite,"%s_%d.vtk",sOutputBaseFileName.c_str(),iStep);
		string sOutputFileName = cWrite;
		FILE* fpOutput = fopen(sOutputFileName.c_str(),"w");
		fprintf(fpOutput,"# vtk DataFile Version 3.0\n");
		fprintf(fpOutput,"output from step %d\n",iStep);
		fprintf(fpOutput,"ASCII\n");
		fprintf(fpOutput,"DATASET UNSTRUCTURED_GRID\n");
		
		unsigned int i = 0;
		unsigned int iNodesCount = m_poDataStructure->GetNodesCount();
		fprintf(fpOutput,"POINTS %d double\n",iNodesCount);
		// write the nodes
		FEMSolidNode* poTempNode = NULL;
		for(i = 0; i < iNodesCount ; i++)
		{
			poTempNode = (FEMSolidNode*)(m_poDataStructure->GetNode(i));
			fprintf(fpOutput,"%e\t\t%e\t\t%e\n",poTempNode->GetX(),poTempNode->GetY(),poTempNode->GetZ());
		}
		// write the elements connectivity
		unsigned int iNodesPerElement = 20;
		unsigned int iElementsCount = m_poDataStructure->GetElementsCount();
		unsigned int iTotalDataSize = 0;
		for(i = 0; i < iElementsCount ; i++)
		{
			iTotalDataSize = iTotalDataSize + ((FEMSolidElement*)(m_poDataStructure->GetElement(i)))->GetGeometry()->GetNodesCount();
		}
		iTotalDataSize = iTotalDataSize + iElementsCount;
		fprintf(fpOutput,"CELLS %d %d\n",iElementsCount,iTotalDataSize);
		FEMSolidElement* poTempElement = NULL;
		FEMElementGeometry* poGeometry = NULL;
		unsigned int iElementNodesCount = 0;
		unsigned int j = 0;
		for(i = 0; i < iElementsCount ; i++)
		{
			poTempElement = (FEMSolidElement*)(m_poDataStructure->GetElement(i));
			poGeometry = poTempElement->GetGeometry();
			iElementNodesCount = poGeometry->GetNodesCount();
			fprintf(fpOutput,"%d\t\t",iElementNodesCount);
			for(j = 0; j < iElementNodesCount ; j++)
			{
				fprintf(fpOutput,"%d\t\t",poGeometry->GetNode(j)->GetID() - 1);
			}
			fprintf(fpOutput,"\n");
		}
		// write the elements types
		fprintf(fpOutput,"CELL_TYPES %d\n",iElementsCount);
		for(i = 0 ; i < iElementsCount ; i++)
		{
			fprintf(fpOutput,"25\n");
		}

		// now write the nodal displacements and forces
		fprintf(fpOutput,"POINT_DATA %d\n",iNodesCount);
		fprintf(fpOutput,"SCALARS node_id double 1\n");
		fprintf(fpOutput,"LOOKUP_TABLE default\n");
		for(i = 0 ; i < iNodesCount ; i++)
		{
			fprintf(fpOutput,"%f\n",(float)i);
		}
		
		fprintf(fpOutput,"VECTORS displacement double\n");
		Vector oTemp;
		for(i = 0 ; i < iNodesCount ; i++)
		{
			oTemp = ((FEMThermoMechanicalNode*)(m_poDataStructure->GetNode(i)))->GetDisplacement();
			fprintf(fpOutput,"%e\t\t%e\t\t%e\n",oTemp.GetX(),oTemp.GetY(),oTemp.GetZ());
		}

		fprintf(fpOutput,"VECTORS force double\n");
		for(i = 0 ; i < iNodesCount ; i++)
		{
			oTemp = ((FEMThermoMechanicalNode*)(m_poDataStructure->GetNode(i)))->GetForce();
			fprintf(fpOutput,"%e\t\t%e\t\t%e\n",oTemp.GetX(),oTemp.GetY(),oTemp.GetZ());
		}
		
		fprintf(fpOutput,"VECTORS velocity double\n");
		for(i = 0 ; i < iNodesCount ; i++)
		{
			oTemp = ((FEMThermoMechanicalNode*)(m_poDataStructure->GetNode(i)))->GetVelocity();
			fprintf(fpOutput,"%e\t\t%e\t\t%e\n",oTemp.GetX(),oTemp.GetY(),oTemp.GetZ());
		}
		
		fprintf(fpOutput,"VECTORS acceleration double\n");
		for(i = 0 ; i < iNodesCount ; i++)
		{
			oTemp = ((FEMThermoMechanicalNode*)(m_poDataStructure->GetNode(i)))->GetAcceleration();
			fprintf(fpOutput,"%e\t\t%e\t\t%e\n",oTemp.GetX(),oTemp.GetY(),oTemp.GetZ());
		}
		
		fprintf(fpOutput,"VECTORS temperature double\n");
		for(i = 0 ; i < iNodesCount ; i++)
		{
			fprintf(fpOutput,"%e\t\t%e\t\t%e\n",((FEMThermoMechanicalNode*)(m_poDataStructure->GetNode(i)))->GetTemperature(),((FEMThermoMechanicalNode*)(m_poDataStructure->GetNode(i)))->GetHeatingRate(),((FEMThermoMechanicalNode*)(m_poDataStructure->GetNode(i)))->GetHeatingAcceleration());
		}
		
		fprintf(fpOutput,"TENSORS stress double\n");
		Matrix oStress;
		for(i = 0 ; i < iNodesCount ; i++)
		{
			oStress = ((FEMThermoMechanicalNode*)(m_poDataStructure->GetNode(i)))->GetStresses();
			fprintf(fpOutput,"%e\t\t%e\t\t%e\n",oStress.Get(1,1),oStress.Get(1,2),oStress.Get(1,3));
			fprintf(fpOutput,"%e\t\t%e\t\t%e\n",oStress.Get(2,1),oStress.Get(2,2),oStress.Get(2,3));
			fprintf(fpOutput,"%e\t\t%e\t\t%e\n",oStress.Get(3,1),oStress.Get(3,2),oStress.Get(3,3));
		}
		
		fprintf(fpOutput,"VECTORS flux double\n");
		for(i = 0 ; i < iNodesCount ; i++)
		{
			oTemp = ((FEMThermoMechanicalNode*)(m_poDataStructure->GetNode(i)))->GetFluxes();
			fprintf(fpOutput,"%e\t\t%e\t\t%e\n",oTemp.GetX(),oTemp.GetY(),oTemp.GetZ());
		}
		
		fclose(fpOutput);
	}
	void FEMImplicitDynamicsThermoMechanicalSolver::GenerateMatrices(const unsigned int& iDOFCount)
	{
		PrintOnScreen("Generating Matrices");
		vector<FEMElement*>* pvpoElements = m_poDataStructure->GetElements();
		SparseMatrix oStiffness(iDOFCount);
		SparseMatrix oDamping(iDOFCount);
		SparseMatrix oMass(iDOFCount);
		unsigned int iSize = pvpoElements->size();
		Matrix oTempStiffnessMatrix;
		Matrix oTempDampingMatrix;
		Matrix oTempMassMatrix;
		vector<unsigned int> viDOFIndices;
		unsigned int i = 0;
		unsigned int j = 0;
		unsigned int k = 0;
		unsigned int iElementDOFCount = 0;
		double dTemp = 0.0;
		for(i = 0; i < iSize ; i++)
		{
			viDOFIndices = pvpoElements->at(i)->GetDOFIndices();
			iElementDOFCount = (unsigned int)viDOFIndices.size();
			oTempStiffnessMatrix = ((FEMThermoMechanicalElement*)pvpoElements->at(i))->GetThermoMechanicalStiffnessMatrix();
			oTempDampingMatrix = ((FEMThermoMechanicalElement*)pvpoElements->at(i))->GetThermoMechanicalDampingMatrix();
			oTempMassMatrix = ((FEMThermoMechanicalElement*)pvpoElements->at(i))->GetThermoMechanicalMassMatrix();
			dTemp = 0.0;
			for(j = 0; j < iElementDOFCount ; j++)
			{
				for(k = 0; k < iElementDOFCount ; k++)
				{
					oStiffness.AddToEntry(viDOFIndices[j],viDOFIndices[k],oTempStiffnessMatrix.Get(j + 1,k + 1));
					oDamping.AddToEntry(viDOFIndices[j],viDOFIndices[k],oTempDampingMatrix.Get(j + 1,k + 1));
					oMass.AddToEntry(viDOFIndices[j],viDOFIndices[k],oTempMassMatrix.Get(j + 1,k + 1));
					dTemp = dTemp + oTempMassMatrix.Get(j + 1,k + 1);
				}
			}
		}
		PartitionMatrices(oStiffness,oDamping,oMass);
		PrintOnScreen("Done generating matrices");
	}
	void FEMImplicitDynamicsThermoMechanicalSolver::PartitionMatrices(const SparseMatrix& oStiffnessMatrix,const SparseMatrix& oDampingMatrix,const SparseMatrix& oMassMatrix)
	{
		unsigned int iSize1 = (unsigned int)m_viUnknownDisplacementsIndices.size();
		unsigned int iSize2 = (unsigned int)m_viKnownDisplacementsIndices.size();
		char cMessageToPrint[500];
		sprintf(cMessageToPrint,"Prtitioning global matrices, original size is %d and partition sizes are %d and %d",oStiffnessMatrix.GetRowsCount(),iSize1,iSize2);
		string sMessageToPrint = cMessageToPrint;
		PrintOnScreen(sMessageToPrint);
	
		m_oDynamicMatrix11.SetRowsCount(iSize1);
		m_oDynamicMatrix12.SetRowsCount(iSize1);
		m_oDynamicMatrix21.SetRowsCount(iSize2);
		m_oDynamicMatrix22.SetRowsCount(iSize2);
		
		m_oMassMatrix11.SetRowsCount(iSize1);
		m_oMassMatrix12.SetRowsCount(iSize1);
		m_oMassMatrix21.SetRowsCount(iSize2);
		m_oMassMatrix22.SetRowsCount(iSize2);
		
		m_oDampingMatrix11.SetRowsCount(iSize1);
		m_oDampingMatrix12.SetRowsCount(iSize1);
		m_oDampingMatrix21.SetRowsCount(iSize2);
		m_oDampingMatrix22.SetRowsCount(iSize2);
		
		double dStiffness = 0.0;
		double dDamping = 0.0;
		double dMass = 0.0;
		bool bStiffnessFound = false;
		bool bDampingFound = false;
		bool bMassFound = false;
		unsigned int i = 0;
		unsigned int j = 0;
		
		for(i = 0; i < iSize1 ; i++)
		{
			for(j = 0; j < iSize1 ; j++)
			{
				bStiffnessFound = false;
				bDampingFound = false;
				bMassFound = false;
				dMass = 0.0;
				dDamping = 0.0;
				dStiffness = 0.0;
				if(oMassMatrix.Get(m_viUnknownDisplacementsIndices[i],m_viUnknownDisplacementsIndices[j],dMass))
				{
					m_oMassMatrix11.Set(i + 1,j + 1,dMass);
					bMassFound = true;
				}
				if(oDampingMatrix.Get(m_viUnknownDisplacementsIndices[i],m_viUnknownDisplacementsIndices[j],dDamping))
				{
					m_oDampingMatrix11.Set(i + 1,j + 1,dDamping);
					bDampingFound = true;
				}
				if(oStiffnessMatrix.Get(m_viUnknownDisplacementsIndices[i],m_viUnknownDisplacementsIndices[j],dStiffness))
				{
					bStiffnessFound = true;
				}
				if(bMassFound || bStiffnessFound || bDampingFound)
				{
					m_oDynamicMatrix11.Set(i + 1,j + 1,dStiffness + m_dA0*dMass + m_dA1*dDamping);
				}
			}
		}
		for(i = 0; i < iSize1 ; i++)
		{
			for(j = 0; j < iSize2 ; j++)
			{
				bStiffnessFound = false;
				bDampingFound = false;
				bMassFound = false;
				dMass = 0.0;
				dDamping = 0.0;
				dStiffness = 0.0;
				if(oMassMatrix.Get(m_viUnknownDisplacementsIndices[i],m_viKnownDisplacementsIndices[j],dMass))
				{
					m_oMassMatrix12.Set(i + 1,j + 1,dMass);
					bMassFound = true;
				}
				if(oDampingMatrix.Get(m_viUnknownDisplacementsIndices[i],m_viKnownDisplacementsIndices[j],dDamping))
				{
					m_oDampingMatrix12.Set(i + 1,j + 1,dDamping);
					bDampingFound = true;
				}
				if(oStiffnessMatrix.Get(m_viUnknownDisplacementsIndices[i],m_viKnownDisplacementsIndices[j],dStiffness))
				{
					bStiffnessFound = true;
				}
				
				if(bMassFound || bStiffnessFound || bDampingFound)
				{
					m_oDynamicMatrix12.Set(i + 1,j + 1,dStiffness + m_dA0*dMass + m_dA1*dDamping);
				}
			}
		}
		for(i = 0; i < iSize2 ; i++)
		{
			for(j = 0; j < iSize1 ; j++)
			{
				bStiffnessFound = false;
				bDampingFound = false;
				bMassFound = false;
				dMass = 0.0;
				dDamping = 0.0;
				dStiffness = 0.0;
				if(oMassMatrix.Get(m_viKnownDisplacementsIndices[i],m_viUnknownDisplacementsIndices[j],dMass))
				{
					m_oMassMatrix21.Set(i + 1,j + 1,dMass);
					bMassFound = true;
				}
				if(oDampingMatrix.Get(m_viKnownDisplacementsIndices[i],m_viUnknownDisplacementsIndices[j],dDamping))
				{
					m_oDampingMatrix21.Set(i + 1,j + 1,dDamping);
					bDampingFound = true;
				}
				if(oStiffnessMatrix.Get(m_viKnownDisplacementsIndices[i],m_viUnknownDisplacementsIndices[j],dStiffness))
				{
					bStiffnessFound = true;
				}
				
				if(bMassFound || bStiffnessFound || bDampingFound)
				{
					m_oDynamicMatrix21.Set(i + 1,j + 1,dStiffness + m_dA0*dMass + m_dA1*dDamping);
				}
			}
		}
		for(i = 0; i < iSize2 ; i++)
		{
			for(j = 0; j < iSize2 ; j++)
			{
				bStiffnessFound = false;
				bDampingFound = false;
				bMassFound = false;
				dMass = 0.0;
				dDamping = 0.0;
				dStiffness = 0.0;
				if(oMassMatrix.Get(m_viKnownDisplacementsIndices[i],m_viKnownDisplacementsIndices[j],dMass))
				{
					m_oMassMatrix22.Set(i + 1,j + 1,dMass);
					bMassFound = true;
				}
				if(oDampingMatrix.Get(m_viKnownDisplacementsIndices[i],m_viKnownDisplacementsIndices[j],dDamping))
				{
					m_oDampingMatrix22.Set(i + 1,j + 1,dDamping);
					bDampingFound = true;
				}
				if(oStiffnessMatrix.Get(m_viKnownDisplacementsIndices[i],m_viKnownDisplacementsIndices[j],dStiffness))
				{
					bStiffnessFound = true;
				}
				
				if(bMassFound || bStiffnessFound || bDampingFound)
				{
					m_oDynamicMatrix22.Set(i + 1,j + 1,dStiffness + m_dA0*dMass + m_dA1*dDamping);
				}
			}
		}
		// now since the dynamic matrix is built, generate its preconditioner
		m_oDynamicMatrix11.InitializeGMRESMatrices();
		PrintOnScreen("Done partitioning matrices");
	}
	void FEMImplicitDynamicsThermoMechanicalSolver::GenerateNodalStateVectors()
	{
		// first, copy the current displacements to the previous displacements
		m_oPreviousDisplacements1 = m_oCurrentDisplacements1;
		m_oPreviousDisplacements2 = m_oCurrentDisplacements2;
		
		m_oPreviousVelocities1 = m_oCurrentVelocities1;
		m_oPreviousVelocities2 = m_oCurrentVelocities2;
		
		m_oPreviousAccelerations1 = m_oCurrentAccelerations1;
		m_oPreviousAccelerations2 = m_oCurrentAccelerations2;
		
		vector<FEMNode*>* pvpoNodes = m_poDataStructure->GetNodes();
		unsigned iSize = (unsigned int)pvpoNodes->size();
		unsigned int i = 0;
		FEMThermoMechanicalNode* poNode = NULL;
		FEMDegreeOfFreedom* poDOF = NULL;
		unsigned int iIndex = 0;
		Vector oDisplacement;
		Vector oForce;
		for(i = 0; i < iSize ; i++)
		{
			poNode = (FEMThermoMechanicalNode*)pvpoNodes->at(i);
			oDisplacement = poNode->GetDisplacement();
			oForce = poNode->GetForce();
			
			poDOF = poNode->GetXDOF();
			iIndex = poDOF->GetIndex() - 1;
			if(poDOF->GetCondition())
			{
				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
				m_oCurrentDisplacements2.Set(iIndex,1,oDisplacement.GetX());
			}
			else
			{
				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
				m_oForces1.Set(iIndex,1,oForce.GetX());
			}
	
			poDOF = poNode->GetYDOF();
			iIndex = poDOF->GetIndex() - 1;
			if(poDOF->GetCondition())
			{
				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
				m_oCurrentDisplacements2.Set(iIndex,1,oDisplacement.GetY());
			}
			else
			{
				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
				m_oForces1.Set(iIndex,1,oForce.GetY());
			}
	
			poDOF = poNode->GetZDOF();
			iIndex = poDOF->GetIndex() - 1;
			if(poDOF->GetCondition())
			{
				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
				m_oCurrentDisplacements2.Set(iIndex,1,oDisplacement.GetZ());
			}
			else
			{
				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
				m_oForces1.Set(iIndex,1,oForce.GetZ());
			}
			
			poDOF = poNode->GetTDOF();
			iIndex = poDOF->GetIndex() - 1;
			if(poDOF->GetCondition())
			{
				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
				m_oCurrentDisplacements2.Set(iIndex,1,poDOF->GetPrimaryValue());
			}
			else
			{
				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
				m_oForces1.Set(iIndex,1,poDOF->GetSecondaryValue());
			}
		}
	}
	void FEMImplicitDynamicsThermoMechanicalSolver::UpdateNodalStateVectors()
	{
		vector<FEMNode*>* pvpoNodes = m_poDataStructure->GetNodes();
		unsigned iSize = (unsigned int)pvpoNodes->size();
		unsigned int i = 0;
		FEMThermoMechanicalNode* poNode = NULL;
		FEMDegreeOfFreedom* poDOF = NULL;
		unsigned int iIndex = 0;
		double dVx = 0.0;
		double dVy = 0.0;
		double dVz = 0.0;
		double dVt = 0.0;
		double dAx = 0.0;
		double dAy = 0.0;
		double dAz = 0.0;
		double dAt = 0.0;
		for(i = 0; i < iSize ; i++)
		{
			poNode = (FEMThermoMechanicalNode*)pvpoNodes->at(i);
			
			poDOF = poNode->GetXDOF();
			iIndex = poDOF->GetIndex() - 1;
			if(poDOF->GetCondition())
			{
				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
				poDOF->SetSecondaryValue(m_oForces2.Get(iIndex,1));
				dVx = m_oCurrentVelocities2.Get(iIndex,1);
				dAx = m_oCurrentAccelerations2.Get(iIndex,1);
			}
			else
			{
				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
				poDOF->SetPrimaryValue(m_oCurrentDisplacements1.Get(iIndex,1));
				dVx = m_oCurrentVelocities1.Get(iIndex,1);
				dAx = m_oCurrentAccelerations1.Get(iIndex,1);
			}
	
			poDOF = poNode->GetYDOF();
			iIndex = poDOF->GetIndex() - 1;
			if(poDOF->GetCondition())
			{
				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
				poDOF->SetSecondaryValue(m_oForces2.Get(iIndex,1));
				dVy = m_oCurrentVelocities2.Get(iIndex,1);
				dAy = m_oCurrentAccelerations2.Get(iIndex,1);
			}
			else
			{
				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
				poDOF->SetPrimaryValue(m_oCurrentDisplacements1.Get(iIndex,1));
				dVy = m_oCurrentVelocities1.Get(iIndex,1);
				dAy = m_oCurrentAccelerations1.Get(iIndex,1);
			}
	
			poDOF = poNode->GetZDOF();
			iIndex = poDOF->GetIndex() - 1;
			if(poDOF->GetCondition())
			{
				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
				poDOF->SetSecondaryValue(m_oForces2.Get(iIndex,1));
				dVz = m_oCurrentVelocities2.Get(iIndex,1);
				dAz = m_oCurrentAccelerations2.Get(iIndex,1);
			}
			else
			{
				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
				poDOF->SetPrimaryValue(m_oCurrentDisplacements1.Get(iIndex,1));
				dVz = m_oCurrentVelocities1.Get(iIndex,1);
				dAz = m_oCurrentAccelerations1.Get(iIndex,1);
			}
			
			poDOF = poNode->GetTDOF();
			iIndex = poDOF->GetIndex() - 1;
			if(poDOF->GetCondition())
			{
				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
				poDOF->SetSecondaryValue(m_oForces2.Get(iIndex,1));
				dVt = m_oCurrentVelocities2.Get(iIndex,1);
				dAt = m_oCurrentAccelerations2.Get(iIndex,1);
			}
			else
			{
				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
				poDOF->SetPrimaryValue(m_oCurrentDisplacements1.Get(iIndex,1));
				dVt = m_oCurrentVelocities1.Get(iIndex,1);
				dAt = m_oCurrentAccelerations1.Get(iIndex,1);
			}
			
			poNode->SetVelocity(Vector(dVx,dVy,dVz));
			poNode->SetAcceleration(Vector(dAx,dAy,dAz));
			poNode->SetHeatingRate(dVt);
			poNode->SetHeatingAcceleration(dAt);
		}
	}
	void FEMImplicitDynamicsThermoMechanicalSolver::GenerateNodalStressesMatrix()
	{
		vector<FEMNode*>* pvpoNodes = m_poDataStructure->GetNodes();
		vector<FEMElement*>* pvpoElements = m_poDataStructure->GetElements();
		unsigned int iSize = pvpoNodes->size();
		unsigned int iStressesCount = 6;
		Matrix oNodalStresses(iSize,iStressesCount);
		Matrix oNodeRecurrence(iSize,1);
		unsigned int i = 0;
	
		iSize = pvpoElements->size();
		vector<FEMNode*>* pvpoElementNodes;
		Matrix oElementNodalStresses;
		unsigned int j = 0;
		unsigned int k = 0;
		unsigned int iNodesCount = 0;
		unsigned int iNodeID = 0;
	
		for(i = 0; i < iSize ; i++)
		{
			pvpoElementNodes = pvpoElements->at(i)->GetNodes();
			// update all Gauss points before getting the stresses
			pvpoElements->at(i)->UpdateGaussPoints();
			oElementNodalStresses = ((FEMSolidElement*)pvpoElements->at(i))->GetNodalStresses();
			iNodesCount = pvpoElementNodes->size();
			for(j = 0; j < iNodesCount ; j++)
			{
				iNodeID = pvpoElementNodes->at(j)->GetID();
				oNodeRecurrence.AddToEntry(iNodeID,1,1);
				for(k = 1; k <= iStressesCount ; k++)
				{
					oNodalStresses.AddToEntry(iNodeID,k,oElementNodalStresses.Get(j + 1,k));
				}
			}
		}
		
		iSize = pvpoNodes->size();
		char cMessageToPrint[500];
		string sMessageToPrint = "";
		unsigned int iRecurrence = 0;
		for(i = 1; i <= iSize ; i++)
		{
			iRecurrence = oNodeRecurrence.Get(i,1);
			if(iRecurrence == 0)
			{
				sprintf(cMessageToPrint,"node %d has never ocurred in stress distribution",i - 1);
				sMessageToPrint = cMessageToPrint;
				PrintOnScreen(sMessageToPrint);
				continue;
			}
			for(j = 1 ; j <= iStressesCount ; j++)
			{
				oNodalStresses.Set(i,j,oNodalStresses.Get(i,j)/(double)iRecurrence);
			}
		}
		
		Matrix oTempStress(3,3);
		for(i = 1; i <= iSize ; i++)
		{ 			
			oTempStress.Set(1,1,oNodalStresses.Get(i,1));
			oTempStress.Set(1,2,oNodalStresses.Get(i,4));
			oTempStress.Set(1,3,oNodalStresses.Get(i,6));
			oTempStress.Set(2,1,oNodalStresses.Get(i,4));
			oTempStress.Set(2,2,oNodalStresses.Get(i,2));
			oTempStress.Set(2,3,oNodalStresses.Get(i,5));
			oTempStress.Set(3,1,oNodalStresses.Get(i,6));
			oTempStress.Set(3,2,oNodalStresses.Get(i,5));
			oTempStress.Set(3,3,oNodalStresses.Get(i,3));
			((FEMThermoMechanicalNode*)pvpoNodes->at(i - 1))->SetStresses(oTempStress);
		}
	}
	void FEMImplicitDynamicsThermoMechanicalSolver::GenerateNodalFluxesMatrix()
	{
	 	vector<FEMNode*>* pvpoNodes = m_poDataStructure->GetNodes();
 		vector<FEMElement*>* pvpoElements = m_poDataStructure->GetElements();
 		unsigned int iSize = pvpoNodes->size();
 		unsigned int iFluxesCount = 3;
 		Matrix oNodalFluxes(iSize,iFluxesCount);
 		Matrix oNodeRecurrence(iSize,1);
 		unsigned int i = 0;
 
 		iSize = pvpoElements->size();
 		vector<FEMNode*>* pvpoElementNodes;
 		Matrix oElementNodalFluxes;
 		unsigned int j = 0;
 		unsigned int k = 0;
 		unsigned int iNodesCount = 0;
 		unsigned int iNodeID = 0;
 
 		for(i = 0; i < iSize ; i++)
 		{
 			pvpoElementNodes = pvpoElements->at(i)->GetNodes();
 			oElementNodalFluxes = ((FEMThermoMechanicalElement*)pvpoElements->at(i))->GetNodalFluxes();
 			iNodesCount = pvpoElementNodes->size();
 			for(j = 0; j < iNodesCount ; j++)
 			{
 				iNodeID = pvpoElementNodes->at(j)->GetID();
 				oNodeRecurrence.AddToEntry(iNodeID,1,1);
 				for(k = 1; k <= iFluxesCount ; k++)
 				{
 					oNodalFluxes.AddToEntry(iNodeID,k,oElementNodalFluxes.Get(j + 1,k));
 				}
 			}
 		}
 		
 		iSize = pvpoNodes->size();
 		char cMessageToPrint[500];
 		string sMessageToPrint = "";
		unsigned int iRecurrence = 0;
		for(i = 1; i <= iSize ; i++)
		{
			iRecurrence = oNodeRecurrence.Get(i,1);
			if(iRecurrence == 0)
			{
				sprintf(cMessageToPrint,"node %d has never ocurred in flux distribution",i - 1);
				sMessageToPrint = cMessageToPrint;
				PrintOnScreen(sMessageToPrint);
				continue;
			}
			for(j = 1 ; j <= iFluxesCount ; j++)
			{
				oNodalFluxes.Set(i,j,oNodalFluxes.Get(i,j)/(double)iRecurrence);
			}
		}
  		
 		Vector oTempFlux;
  		for(i = 1; i <= iSize ; i++)
  		{
  			oTempFlux.SetX(oNodalFluxes.Get(i,1));
  			oTempFlux.SetY(oNodalFluxes.Get(i,2));
  			oTempFlux.SetZ(oNodalFluxes.Get(i,3));
 			((FEMThermoMechanicalNode*)pvpoNodes->at(i - 1))->SetFluxes(oTempFlux);
  		}
	}
}



