#include "FEMExplicitDynamicsSolidSolver.h"
#include "FEMSolidElement.h"
#include "FEMSolidNode.h"
#include "Tools.h"
#include "cmath"

using namespace SupportSystem;

namespace FEMSystem
{
	FEMExplicitDynamicsSolidSolver::FEMExplicitDynamicsSolidSolver()
	{
		Initialize();
	}
	FEMExplicitDynamicsSolidSolver::FEMExplicitDynamicsSolidSolver(const FEMExplicitDynamicsSolidSolver& oSolver)
	{
		*this = oSolver;
	}
	FEMExplicitDynamicsSolidSolver::~FEMExplicitDynamicsSolidSolver()
	{
		Reset();
	}
	FEMExplicitDynamicsSolidSolver& FEMExplicitDynamicsSolidSolver::operator=(const FEMExplicitDynamicsSolidSolver& oSolver)
	{
		FEMSolver::operator=(oSolver);
		m_oForces1 = oSolver.m_oForces1;
		m_oForces2 = oSolver.m_oForces2;
		m_oNextDisplacements1 = oSolver.m_oNextDisplacements1;
		m_oNextDisplacements2 = oSolver.m_oNextDisplacements2;
		m_oCurrentDisplacements1 = oSolver.m_oCurrentDisplacements1;
		m_oCurrentDisplacements2 = oSolver.m_oCurrentDisplacements2;
		m_oPreviousDisplacements1 = oSolver.m_oPreviousDisplacements1;
		m_oPreviousDisplacements2 = oSolver.m_oPreviousDisplacements2;
		m_oVelocities1 = oSolver.m_oVelocities1;
		m_oVelocities2 = oSolver.m_oVelocities2;
		m_oAccelerations1 = oSolver.m_oAccelerations1;
		m_oAccelerations2 = oSolver.m_oAccelerations2;
		
		m_viUnknownDisplacementsIndices = oSolver.m_viUnknownDisplacementsIndices;
		m_viKnownDisplacementsIndices = oSolver.m_viKnownDisplacementsIndices;
		m_viUnknownDisplacementsIndicesReverseMap = oSolver.m_viUnknownDisplacementsIndicesReverseMap;
		m_viKnownDisplacementsIndicesReverseMap = oSolver.m_viKnownDisplacementsIndicesReverseMap;
		
		m_oDynamicMatrix1 = oSolver.m_oDynamicMatrix1;
		m_oDynamicMatrix2 = oSolver.m_oDynamicMatrix2;
		
		m_oT11Matrix = oSolver.m_oT11Matrix;
		m_oT12Matrix = oSolver.m_oT12Matrix;
		m_oT21Matrix = oSolver.m_oT21Matrix;
		m_oT22Matrix = oSolver.m_oT22Matrix;
		m_oS11Matrix = oSolver.m_oS11Matrix;
		m_oS12Matrix = oSolver.m_oS12Matrix;
		m_oS21Matrix = oSolver.m_oS21Matrix;
		m_oS22Matrix = oSolver.m_oS22Matrix;
		return *this;
	}
	void FEMExplicitDynamicsSolidSolver::Reset()
	{
		FEMSolver::Reset();
		Initialize();
	}
	void FEMExplicitDynamicsSolidSolver::Initialize()
	{ 		
		FEMSolver::Initialize();
		m_oForces1.Reset();
		m_oForces2.Reset();
		m_oNextDisplacements1.Reset();
		m_oNextDisplacements2.Reset();
		m_oCurrentDisplacements1.Reset();
		m_oCurrentDisplacements2.Reset();
		m_oPreviousDisplacements1.Reset();
		m_oPreviousDisplacements2.Reset();
		m_oVelocities1.Reset();
		m_oVelocities2.Reset();
		m_oAccelerations1.Reset();
		m_oAccelerations2.Reset();
		
		m_viUnknownDisplacementsIndices.clear();
		m_viKnownDisplacementsIndices.clear();
		m_viUnknownDisplacementsIndicesReverseMap.clear();
		m_viKnownDisplacementsIndicesReverseMap.clear();
		
		m_oDynamicMatrix1.Reset();
		m_oDynamicMatrix2.Reset();
		
		m_oT11Matrix.Reset();
		m_oT12Matrix.Reset();
		m_oT21Matrix.Reset();
		m_oT22Matrix.Reset();
		m_oS11Matrix.Reset();
		m_oS12Matrix.Reset();
		m_oS21Matrix.Reset();
		m_oS22Matrix.Reset();
	}
	Vector FEMExplicitDynamicsSolidSolver::GetDisplacement(Point* poPoint,unsigned int& iStatus)
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
			oDisplacement = ((FEMSolidElement*)poElement)->GetDisplacement(oNaturalCoordinates);
		}
		return oDisplacement;
	}
	Matrix FEMExplicitDynamicsSolidSolver::GetStress(Point* poPoint,unsigned int& iStatus)
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
			oStress = ((FEMSolidElement*)poElement)->GetStress(oNaturalCoordinates);
		}
		return oStress;
	}
	void FEMExplicitDynamicsSolidSolver::GenerateMatrices(const unsigned int& iDOFCount)
	{
		PrintOnScreen("Generating Matrices");
		vector<FEMElement*>* pvpoElements = m_poDataStructure->GetElements();
		SparseMatrix oTMatrix(iDOFCount);
		Matrix oMass(iDOFCount,1);
		unsigned int iSize = pvpoElements->size();
		Matrix oTempStiffnessMatrix;
		Matrix oTempMassMatrix;
		vector<unsigned int> viDOFIndices;
		unsigned int i = 0;
		unsigned int j = 0;
		unsigned int k = 0;
		unsigned int iElementDOFCount = 0;
		for(i = 0; i < iSize ; i++)
		{
			viDOFIndices = pvpoElements->at(i)->GetDOFIndices();
			iElementDOFCount = (unsigned int)viDOFIndices.size();
			oTempStiffnessMatrix = ((FEMSolidElement*)pvpoElements->at(i))->GetStiffnessMatrix();
			oTempMassMatrix = ((FEMSolidElement*)pvpoElements->at(i))->GetLumpedMassMatrix();
			for(j = 0; j < iElementDOFCount ; j++)
			{
				for(k = 0; k < iElementDOFCount ; k++)
				{
					oTMatrix.AddToEntry(viDOFIndices[j],viDOFIndices[k],oTempStiffnessMatrix.Get(j + 1,k + 1));
				}
				oMass.AddToEntry(viDOFIndices[j],1,oTempMassMatrix.Get(j + 1,1));
			}
		}
		// make sure that the time step is safe
		AdjustTimeStep(oTMatrix,oMass);
		// before partitioning, generate the T and S matrices
		// the S matrix is a vector, but it will be implemented as a matrix so that the same 
		// functions can be used in the simulations where damping is present. Notice that this does
		// not affect the performance since each map in the sparse matrix will end up having 1 elemnent
		// and the access time will be O(1)
		SparseMatrix oSMatrix(iDOFCount);
		double dMass = 0.0;
		for(i = 1; i <= iDOFCount ; i++)
		{
			dMass = oMass.Get(i,1);
			oTMatrix.AddToEntry(i,i,-m_dA2*dMass);
			oSMatrix.Set(i,i,m_dA0*dMass);
		}
		PartitionMatrices(oTMatrix,oSMatrix,oMass);
		PrintOnScreen("Done generating matrices");
	}
	void FEMExplicitDynamicsSolidSolver::PartitionMatrices(const SparseMatrix& oTMatrix,const SparseMatrix& oSMatrix,const Matrix& oMassMatrix)
	{
		unsigned int iSize1 = (unsigned int)m_viUnknownDisplacementsIndices.size();
		unsigned int iSize2 = (unsigned int)m_viKnownDisplacementsIndices.size();
		char cMessageToPrint[500];
		sprintf(cMessageToPrint,"Prtitioning global matrices, original size is %d and partition sizes are %d and %d",iSize1 + iSize2,iSize1,iSize2);
		string sMessageToPrint = cMessageToPrint;
		PrintOnScreen(sMessageToPrint);

		double dTemp = 0.0;
		unsigned int i = 0;
		unsigned int j = 0;
		
		m_oT11Matrix.SetRowsCount(iSize1);
		m_oT12Matrix.SetRowsCount(iSize1);
		m_oT21Matrix.SetRowsCount(iSize2);
		m_oT22Matrix.SetRowsCount(iSize2);
		m_oS11Matrix.SetRowsCount(iSize1);
		m_oS12Matrix.SetRowsCount(iSize1);
		m_oS21Matrix.SetRowsCount(iSize2);
		m_oS22Matrix.SetRowsCount(iSize2);
		m_oDynamicMatrix1.SetSize(iSize1,1);
		m_oDynamicMatrix2.SetSize(iSize2,1);
		
		for(i = 0; i < iSize1 ; i++)
		{
			for(j = 0; j < iSize1 ; j++)
			{
				if(oTMatrix.Get(m_viUnknownDisplacementsIndices[i],m_viUnknownDisplacementsIndices[j],dTemp))
				{
					m_oT11Matrix.Set(i + 1,j + 1,dTemp);
				}
				if(oSMatrix.Get(m_viUnknownDisplacementsIndices[i],m_viUnknownDisplacementsIndices[j],dTemp))
				{
					m_oS11Matrix.Set(i + 1,j + 1,dTemp);
				}
			}
			m_oDynamicMatrix1.Set(i + 1,1,m_dA0*oMassMatrix.Get(m_viUnknownDisplacementsIndices[i],1));
		}
		for(i = 0; i < iSize1 ; i++)
		{
			for(j = 0; j < iSize2 ; j++)
			{
				if(oTMatrix.Get(m_viUnknownDisplacementsIndices[i],m_viKnownDisplacementsIndices[j],dTemp))
				{
					m_oT12Matrix.Set(i + 1,j + 1,dTemp);
				}
				if(oSMatrix.Get(m_viUnknownDisplacementsIndices[i],m_viKnownDisplacementsIndices[j],dTemp))
				{
					m_oS12Matrix.Set(i + 1,j + 1,dTemp);
				}
			}
		}
		for(i = 0; i < iSize2 ; i++)
		{
			for(j = 0; j < iSize1 ; j++)
			{
				if(oTMatrix.Get(m_viKnownDisplacementsIndices[i],m_viUnknownDisplacementsIndices[j],dTemp))
				{
					m_oT21Matrix.Set(i + 1,j + 1,dTemp);
				}
				if(oSMatrix.Get(m_viKnownDisplacementsIndices[i],m_viUnknownDisplacementsIndices[j],dTemp))
				{
					m_oS21Matrix.Set(i + 1,j + 1,dTemp);
				}
			}
		}
		for(i = 0; i < iSize2 ; i++)
		{
			for(j = 0; j < iSize2 ; j++)
			{
				if(oTMatrix.Get(m_viKnownDisplacementsIndices[i],m_viKnownDisplacementsIndices[j],dTemp))
				{
					m_oT22Matrix.Set(i + 1,j + 1,dTemp);
				}
				if(oSMatrix.Get(m_viKnownDisplacementsIndices[i],m_viKnownDisplacementsIndices[j],dTemp))
				{
					m_oS22Matrix.Set(i + 1,j + 1,dTemp);
				}
			}
			m_oDynamicMatrix2.Set(i + 1,1,m_dA0*oMassMatrix.Get(m_viKnownDisplacementsIndices[i],1));
		}
		PrintOnScreen("Done partitioning matrices");
	}
	void FEMExplicitDynamicsSolidSolver::GenerateNodalStressesMatrix()
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
			((FEMSolidNode*)pvpoNodes->at(i - 1))->SetStresses(oTempStress);
		}
	}
	void FEMExplicitDynamicsSolidSolver::InitializeMatrices(const double& dTimeStep)
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
 		FEMSolidNode* poNode = NULL;
 		unsigned int iIndex = 0;
 		unsigned int iKnownIndicesCount = 0;
 		unsigned int iUnknownIndicesCount = 0;
 		for(i = 0; i < iSize ; i++)
 		{
 			poNode = (FEMSolidNode*)pvpoNodes->at(i);
 			
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
 		}

		// set the explicit method parameters
		SetMethodParameters(dTimeStep);
 		GenerateMatrices(iDOFCount);
 		// now resize the force, displacement, velocity and acceleration arrays
 		unsigned int iSize1 = (unsigned int)m_viUnknownDisplacementsIndices.size();
 		unsigned int iSize2 = (unsigned int)m_viKnownDisplacementsIndices.size();
 		m_oForces1.SetSize(iSize1,1);
 		m_oForces2.SetSize(iSize2,1);
 		m_oNextDisplacements1.SetSize(iSize1,1);
 		m_oNextDisplacements2.SetSize(iSize2,1);
		m_oCurrentDisplacements1.SetSize(iSize1,1);
		m_oCurrentDisplacements2.SetSize(iSize2,1);
		m_oPreviousDisplacements1.SetSize(iSize1,1);
		m_oPreviousDisplacements2.SetSize(iSize2,1);
		m_oVelocities1.SetSize(iSize1,1);
		m_oVelocities2.SetSize(iSize2,1);
		m_oAccelerations1.SetSize(iSize1,1);
		m_oAccelerations2.SetSize(iSize2,1);
		// no need to resize these matrices again, they should be directly used
 	 	PrintOnScreen("Done initializing");
	}
	void FEMExplicitDynamicsSolidSolver::Solve(const double& dTime)
	{
		printf("Solving the system .... ");
		GenerateNodalStateVectors();
		// get the unknown displacements
		Matrix oRHS1 = m_oT11Matrix*m_oCurrentDisplacements1 + m_oT12Matrix*m_oCurrentDisplacements2;
		Matrix oRHS2 = m_oS11Matrix*m_oPreviousDisplacements1 + m_oS12Matrix*m_oPreviousDisplacements2;
		Matrix oRHS = m_oForces1 - oRHS1 - oRHS2;
		// solve for the next displacements by looping over the force vector
		unsigned int i = 0;
		unsigned int iSize = m_oDynamicMatrix1.GetRowsCount();
		for(i = 1 ; i <= iSize ; i++)
		{
			m_oNextDisplacements1.Set(i,1,oRHS.Get(i,1)/m_oDynamicMatrix1.Get(i,1));
		}
		// get the unknown forces
		oRHS1 = m_oT21Matrix*m_oCurrentDisplacements1 + m_oT22Matrix*m_oCurrentDisplacements2;
		oRHS2 = m_oS21Matrix*m_oPreviousDisplacements1 + m_oS22Matrix*m_oPreviousDisplacements2;
		// get the product of the dynamic matrix and the displacements
		iSize = m_oDynamicMatrix2.GetRowsCount();
		oRHS.SetSize(iSize,1);
		for(i = 1 ; i <= iSize ; i++)
		{
			oRHS.Set(i,1,m_oDynamicMatrix2.Get(i,1)*m_oNextDisplacements2.Get(i,1));
		}
		m_oForces2 = oRHS + oRHS1 + oRHS2;
		// get the velocities
		m_oVelocities1 = (m_oNextDisplacements1 - m_oPreviousDisplacements1)*m_dA1;
		m_oVelocities2 = (m_oNextDisplacements2 - m_oPreviousDisplacements2)*m_dA1;
		// get the accelerations
		m_oAccelerations1 = (m_oNextDisplacements1 - (m_oCurrentDisplacements1*2.0) + m_oPreviousDisplacements1)*m_dA0;
		m_oAccelerations2 = (m_oNextDisplacements2 - (m_oCurrentDisplacements2*2.0) + m_oPreviousDisplacements2)*m_dA0;
		UpdateNodalStateVectors();
		GenerateNodalStressesMatrix();
		printf("Done\n");
	}
	void FEMExplicitDynamicsSolidSolver::GenerateNodalStateVectors()
	{
		vector<FEMNode*>* pvpoNodes = m_poDataStructure->GetNodes();
		unsigned iSize = (unsigned int)pvpoNodes->size();
		unsigned int i = 0;
		FEMSolidNode* poNode = NULL;
		FEMDegreeOfFreedom* poDOF = NULL;
		unsigned int iIndex = 0;
		Vector oDisplacement;
		Vector oVelocity;
		Vector oAcceleration;
		Vector oPreviousDisplacement;
		for(i = 0; i < iSize ; i++)
		{
			poNode = (FEMSolidNode*)pvpoNodes->at(i);
			oDisplacement = poNode->GetDisplacement();
			oVelocity = poNode->GetVelocity();
			oAcceleration = poNode->GetAcceleration();
			oPreviousDisplacement = oDisplacement - oVelocity*m_dTimeStep + oAcceleration*m_dA3;
			
			poDOF = poNode->GetXDOF();
			iIndex = poDOF->GetIndex() - 1;
			if(poDOF->GetCondition())
			{
				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
				m_oPreviousDisplacements2.Set(iIndex,1,oPreviousDisplacement.GetX());
				m_oCurrentDisplacements2.Set(iIndex,1,oDisplacement.GetX());
				m_oNextDisplacements2.Set(iIndex,1,poDOF->GetConstraintValue());
			}
			else
			{
				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
				m_oPreviousDisplacements1.Set(iIndex,1,oPreviousDisplacement.GetX());
				m_oCurrentDisplacements1.Set(iIndex,1,oDisplacement.GetX());				
				m_oForces1.Set(iIndex,1,poDOF->GetConstraintValue());
			}
	
			poDOF = poNode->GetYDOF();
			iIndex = poDOF->GetIndex() - 1;
			if(poDOF->GetCondition())
			{
				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
				m_oPreviousDisplacements2.Set(iIndex,1,oPreviousDisplacement.GetY());
				m_oCurrentDisplacements2.Set(iIndex,1,oDisplacement.GetY());
				m_oNextDisplacements2.Set(iIndex,1,poDOF->GetConstraintValue());
			}
			else
			{
				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
				m_oPreviousDisplacements1.Set(iIndex,1,oPreviousDisplacement.GetY());
				m_oCurrentDisplacements1.Set(iIndex,1,oDisplacement.GetY());	
				m_oForces1.Set(iIndex,1,poDOF->GetConstraintValue());
			}
	
			poDOF = poNode->GetZDOF();
			iIndex = poDOF->GetIndex() - 1;
			if(poDOF->GetCondition())
			{
				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
				m_oPreviousDisplacements2.Set(iIndex,1,oPreviousDisplacement.GetZ());
				m_oCurrentDisplacements2.Set(iIndex,1,oDisplacement.GetZ());
				m_oNextDisplacements2.Set(iIndex,1,poDOF->GetConstraintValue());
			}
			else
			{
				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
				m_oPreviousDisplacements1.Set(iIndex,1,oPreviousDisplacement.GetZ());
				m_oCurrentDisplacements1.Set(iIndex,1,oDisplacement.GetZ());
				m_oForces1.Set(iIndex,1,poDOF->GetConstraintValue());
			}
		}
	}
	void FEMExplicitDynamicsSolidSolver::UpdateNodalStateVectors()
	{
		vector<FEMNode*>* pvpoNodes = m_poDataStructure->GetNodes();
		unsigned iSize = (unsigned int)pvpoNodes->size();
		unsigned int i = 0;
		FEMSolidNode* poNode = NULL;
		FEMDegreeOfFreedom* poDOF = NULL;
		unsigned int iIndex = 0;
		double dVx = 0.0;
		double dVy = 0.0;
		double dVz = 0.0;
		double dAx = 0.0;
		double dAy = 0.0;
		double dAz = 0.0;
		for(i = 0; i < iSize ; i++)
		{
			poNode = (FEMSolidNode*)pvpoNodes->at(i);
			
			poDOF = poNode->GetXDOF();
			iIndex = poDOF->GetIndex() - 1;
			if(poDOF->GetCondition())
			{
				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
				poDOF->SetSecondaryValue(m_oForces2.Get(iIndex,1));
				dVx = m_oVelocities2.Get(iIndex,1);
				dAx = m_oAccelerations2.Get(iIndex,1);
			}
			else
			{
				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
				poDOF->SetPrimaryValue(m_oNextDisplacements1.Get(iIndex,1));
				dVx = m_oVelocities1.Get(iIndex,1);
				dAx = m_oAccelerations1.Get(iIndex,1);
			}
	
			poDOF = poNode->GetYDOF();
			iIndex = poDOF->GetIndex() - 1;
			if(poDOF->GetCondition())
			{
				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
				poDOF->SetSecondaryValue(m_oForces2.Get(iIndex,1));
				dVy = m_oVelocities2.Get(iIndex,1);
				dAy = m_oAccelerations2.Get(iIndex,1);
			}
			else
			{
				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
				poDOF->SetPrimaryValue(m_oNextDisplacements1.Get(iIndex,1));
				dVy = m_oVelocities1.Get(iIndex,1);
				dAy = m_oAccelerations1.Get(iIndex,1);
			}
	
			poDOF = poNode->GetZDOF();
			iIndex = poDOF->GetIndex() - 1;
			if(poDOF->GetCondition())
			{
				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
				poDOF->SetSecondaryValue(m_oForces2.Get(iIndex,1));
				dVz = m_oVelocities2.Get(iIndex,1);
				dAz = m_oAccelerations2.Get(iIndex,1);
			}
			else
			{
				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
				poDOF->SetPrimaryValue(m_oNextDisplacements1.Get(iIndex,1));
				dVz = m_oVelocities1.Get(iIndex,1);
				dAz = m_oAccelerations1.Get(iIndex,1);
			}
			
			poNode->SetVelocity(Vector(dVx,dVy,dVz));
			poNode->SetAcceleration(Vector(dAx,dAy,dAz));
		}
	}
	void FEMExplicitDynamicsSolidSolver::SetMethodParameters(const double& dTimeStep)
	{
		// set the parameters for the explicit dynamics solver
		m_dTimeStep = dTimeStep;
		m_dA0 = 1.0/dTimeStep/dTimeStep;
		m_dA1 = 0.5/dTimeStep;
		m_dA2 = 2.0*m_dA0;
		m_dA3 = 1.0/m_dA2;
	}
	void FEMExplicitDynamicsSolidSolver::AdjustTimeStep(const SparseMatrix& oStiffnessMatrix,const Matrix& oMassMatrix)
	{
		// this function makes sure that the time step input is consistent with the 
		// stability conditions for the explicit solver. if the time step is less than the
		// critical time step, no changes are made, otherwise, the time step is overriden and
		// all the parameters are modified accordingly
		unsigned int iSize = oStiffnessMatrix.GetRowsCount();
		SparseMatrix oDMatrix = oStiffnessMatrix;
		unsigned int i = 0;
		double dTemp = 0.0;
		for(i = 1 ; i <= iSize ; i++)
		{
			oDMatrix.MultiplyRow(i,1.0/oMassMatrix.Get(i,1));
		}
		double dCriticalTimeStep = 1.0/sqrt(oDMatrix.GetEigenValueUpperBound());
		if(dCriticalTimeStep < m_dTimeStep)
		{
			printf("warning: overriding the given time step %20.15f by %20.15f to ensure stability\n",m_dTimeStep,dCriticalTimeStep);
			SetMethodParameters(dCriticalTimeStep);
		}
	}
	void FEMExplicitDynamicsSolidSolver::WriteFEMSolutionToParaview(const unsigned int& iStep) const
	{
		double dMagnificationFactor = 100.0;
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
			oTemp = ((FEMSolidNode*)(m_poDataStructure->GetNode(i)))->GetDisplacement();
			fprintf(fpOutput,"%e\t\t%e\t\t%e\n",oTemp.GetX(),oTemp.GetY(),oTemp.GetZ());
		}

		fprintf(fpOutput,"VECTORS force double\n");
		for(i = 0 ; i < iNodesCount ; i++)
		{
			oTemp = ((FEMSolidNode*)(m_poDataStructure->GetNode(i)))->GetForce();
			fprintf(fpOutput,"%e\t\t%e\t\t%e\n",oTemp.GetX(),oTemp.GetY(),oTemp.GetZ());
		}
		
		fprintf(fpOutput,"VECTORS velocity double\n");
		for(i = 0 ; i < iNodesCount ; i++)
		{
			oTemp = ((FEMSolidNode*)(m_poDataStructure->GetNode(i)))->GetVelocity();
			fprintf(fpOutput,"%e\t\t%e\t\t%e\n",oTemp.GetX(),oTemp.GetY(),oTemp.GetZ());
		}
		
		fprintf(fpOutput,"VECTORS acceleration double\n");
		for(i = 0 ; i < iNodesCount ; i++)
		{
			oTemp = ((FEMSolidNode*)(m_poDataStructure->GetNode(i)))->GetAcceleration();
			fprintf(fpOutput,"%e\t\t%e\t\t%e\n",oTemp.GetX(),oTemp.GetY(),oTemp.GetZ());
		}
		
		fprintf(fpOutput,"TENSORS stress double\n");
		Matrix oStress;
		for(i = 0 ; i < iNodesCount ; i++)
		{
			oStress = ((FEMSolidNode*)(m_poDataStructure->GetNode(i)))->GetStresses();
			fprintf(fpOutput,"%e\t\t%e\t\t%e\n",oStress.Get(1,1),oStress.Get(1,2),oStress.Get(1,3));
			fprintf(fpOutput,"%e\t\t%e\t\t%e\n",oStress.Get(2,1),oStress.Get(2,2),oStress.Get(2,3));
			fprintf(fpOutput,"%e\t\t%e\t\t%e\n",oStress.Get(3,1),oStress.Get(3,2),oStress.Get(3,3));
		}
		fclose(fpOutput);		
	}
}






