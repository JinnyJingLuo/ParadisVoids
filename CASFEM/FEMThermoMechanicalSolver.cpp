// Ahmed M. Hussein
#include "FEMThermoMechanicalSolver.h"
#include "FEMThermoMechanicalElement.h"
#include "FEMThermoMechanicalNode.h"
#include "Tools.h"
#include "cmath"


namespace FEMSystem
{
	FEMThermoMechanicalSolver::FEMThermoMechanicalSolver()
	{
		Initialize();
	}
	FEMThermoMechanicalSolver::FEMThermoMechanicalSolver(const FEMThermoMechanicalSolver& oSolver)
	{
		*this = oSolver;
	}
	FEMThermoMechanicalSolver::~FEMThermoMechanicalSolver()
	{
		Reset();
	}
	FEMThermoMechanicalSolver& FEMThermoMechanicalSolver::operator=(const FEMThermoMechanicalSolver& oSolver)
	{
		FEMSolidSolver::operator=(oSolver);
		return *this;
	}
	void FEMThermoMechanicalSolver::Reset()
	{
		FEMSolidSolver::Reset();
	}
	void FEMThermoMechanicalSolver::Initialize()
	{
		FEMSolidSolver::Initialize();
	}
	void FEMThermoMechanicalSolver::InitializeMatrices(const double& dTimeStep)
	{
		PrintOnScreen("Initializing FEM problem");
		vector<FEMNode*>* pvpoNodes = m_poDataStructure->GetNodes();
 		unsigned int iSize = pvpoNodes->size();
 		unsigned int i = 0;
 		// get the exact node count
 		for(i = 0; i < iSize ; i++)
 		{
 			m_iDOFCount = m_iDOFCount + pvpoNodes->at(i)->GetDOFCount();
 		}
 		m_viUnknownDisplacementsIndices.reserve(m_iDOFCount);
 		m_viKnownDisplacementsIndices.reserve(m_iDOFCount);
 		// the maximum DOF index cannot be higher than the number of DOFs
 		m_viUnknownDisplacementsIndicesReverseMap.resize(m_iDOFCount);
 		m_viKnownDisplacementsIndicesReverseMap.resize(m_iDOFCount);
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
 	 	
		m_dTimeStep = dTimeStep;
		GenerateMatrices();
		// now resize the force and displacement arrays
		unsigned int iSize1 = (unsigned int)m_viUnknownDisplacementsIndices.size();
		unsigned int iSize2 = (unsigned int)m_viKnownDisplacementsIndices.size();
		m_oForces1.SetSize(iSize1,1);
		m_oForces2.SetSize(iSize2,1);
		m_oDisplacements1.SetSize(iSize1,1);
		m_oDisplacements2.SetSize(iSize2,1);
		// no need to resize these matrices again, they should be directly used
		PrintOnScreen("Done initializing");
	}
	void FEMThermoMechanicalSolver::Solve(const double& dTime)
	{
	 	PrintOnScreen("Generating nodal forces and displacements");
		GenerateNodalStateVectors();
 	 	PrintOnScreen("Done");
		PrintOnScreen("Solving the system");
		Matrix oRHS;
		Matrix oDeltaDisplacement;
		double dTolerance = 1.0E-2;
		double dError = 100.0*dTolerance;
		unsigned int i = 0;
		while(true)
		{
			GenerateMatrices();
			// get the unknown displacements
			oRHS = m_oForces1 - m_oStiffnessMatrix11*m_oDisplacements1 - m_oStiffnessMatrix12*m_oDisplacements2;
			oDeltaDisplacement = m_oStiffnessMatrix11.SolveGMRES(oRHS);
			i = i + 1;
			dError = oDeltaDisplacement.GetNorm()/m_oDisplacements1.GetNorm();
			m_oDisplacements1 = m_oDisplacements1 + oDeltaDisplacement;
			// get the unknown forces
			m_oForces2 = m_oStiffnessMatrix21*m_oDisplacements1 + m_oStiffnessMatrix22*m_oDisplacements2;
			UpdateNodalStateVectors();
			GenerateNodalStressesMatrix();
			GenerateNodalFluxesMatrix();
			printf("iteration : %d \t\t error : %e\n",i,dError);
			if(dError < dTolerance)
			{
				break;
			}
		}
 	 	PrintOnScreen("Done");
	}
	double FEMThermoMechanicalSolver::GetTemperature(Point* poPoint,unsigned int& iStatus)
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
	Vector FEMThermoMechanicalSolver::GetFlux(Point* poPoint,unsigned int& iStatus)
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
	void FEMThermoMechanicalSolver::WriteFEMSolutionToParaview(const unsigned int& iStep) const
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
	void FEMThermoMechanicalSolver::GenerateMatrices()
	{
		PrintOnScreen("Generating Matrices");
		vector<FEMElement*>* pvpoElements = m_poDataStructure->GetElements();
		SparseMatrix oStiffness(m_iDOFCount);
		unsigned int iSize = pvpoElements->size();
		Matrix oTempStiffnessMatrix;
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
			for(j = 0; j < iElementDOFCount ; j++)
			{
				for(k = 0; k < iElementDOFCount ; k++)
				{
				
					oStiffness.AddToEntry(viDOFIndices[j],viDOFIndices[k],oTempStiffnessMatrix.Get(j + 1,k + 1));
				}
			}
		}
		PartitionMatrices(oStiffness);
		PrintOnScreen("Done generating matrices");
	}
	void FEMThermoMechanicalSolver::PartitionMatrices(const SparseMatrix& oStiffnessMatrix)
	{
		FEMSolidSolver::PartitionMatrices(oStiffnessMatrix);
		m_oStiffnessMatrix11.InitializeGMRESMatrices();
	}
	void FEMThermoMechanicalSolver::GenerateNodalStateVectors()
	{
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
				m_oDisplacements2.Set(iIndex,1,oDisplacement.GetX());
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
				m_oDisplacements2.Set(iIndex,1,oDisplacement.GetY());
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
				m_oDisplacements2.Set(iIndex,1,oDisplacement.GetZ());
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
				m_oDisplacements2.Set(iIndex,1,poDOF->GetPrimaryValue());
			}
			else
			{
				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
				m_oForces1.Set(iIndex,1,poDOF->GetSecondaryValue());
			}
		}
	}
	void FEMThermoMechanicalSolver::UpdateNodalStateVectors()
	{
		vector<FEMNode*>* pvpoNodes = m_poDataStructure->GetNodes();
		unsigned iSize = (unsigned int)pvpoNodes->size();
		unsigned int i = 0;
		FEMThermoMechanicalNode* poNode = NULL;
		FEMDegreeOfFreedom* poDOF = NULL;
		unsigned int iIndex = 0;
		for(i = 0; i < iSize ; i++)
		{
			poNode = (FEMThermoMechanicalNode*)pvpoNodes->at(i);
			
			poDOF = poNode->GetXDOF();
			iIndex = poDOF->GetIndex() - 1;
			if(poDOF->GetCondition())
			{
				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
				poDOF->SetSecondaryValue(m_oForces2.Get(iIndex,1));
			}
			else
			{
				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
				poDOF->SetPrimaryValue(m_oDisplacements1.Get(iIndex,1));
			}
	
			poDOF = poNode->GetYDOF();
			iIndex = poDOF->GetIndex() - 1;
			if(poDOF->GetCondition())
			{
				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
				poDOF->SetSecondaryValue(m_oForces2.Get(iIndex,1));
			}
			else
			{
				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
				poDOF->SetPrimaryValue(m_oDisplacements1.Get(iIndex,1));
			}
	
			poDOF = poNode->GetZDOF();
			iIndex = poDOF->GetIndex() - 1;
			if(poDOF->GetCondition())
			{
				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
				poDOF->SetSecondaryValue(m_oForces2.Get(iIndex,1));
			}
			else
			{
				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
				poDOF->SetPrimaryValue(m_oDisplacements1.Get(iIndex,1));
			}
			
			poDOF = poNode->GetTDOF();
			iIndex = poDOF->GetIndex() - 1;
			if(poDOF->GetCondition())
			{
				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
				poDOF->SetSecondaryValue(m_oForces2.Get(iIndex,1));
			}
			else
			{
				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
				poDOF->SetPrimaryValue(m_oDisplacements1.Get(iIndex,1));
			}
		}
	}
	void FEMThermoMechanicalSolver::GenerateNodalFluxesMatrix()
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



