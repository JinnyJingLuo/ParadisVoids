 #include "FEMSolidSolver.h"
 #include "FEMSolidElement.h"
 #include "FEMSolidNode.h"
 #include "Tools.h"
 
 using namespace SupportSystem;
 
 namespace FEMSystem
 {
 	FEMSolidSolver::FEMSolidSolver()
 	{
 		Initialize();
 	}
 	FEMSolidSolver::FEMSolidSolver(const FEMSolidSolver& oSolver)
 	{
 		*this = oSolver;
 	}
 	FEMSolidSolver::~FEMSolidSolver()
 	{
 		Reset();
 	}
 	FEMSolidSolver& FEMSolidSolver::operator=(const FEMSolidSolver& oSolver)
 	{
 		FEMSolver::operator=(oSolver);
 		m_oForces1 = oSolver.m_oForces1;
 		m_oForces2 = oSolver.m_oForces2;
 		m_oDisplacements1 = oSolver.m_oDisplacements1;
 		m_oDisplacements2 = oSolver.m_oDisplacements2;

 		m_viUnknownDisplacementsIndices = oSolver.m_viUnknownDisplacementsIndices;
 		m_viKnownDisplacementsIndices = oSolver.m_viKnownDisplacementsIndices;
 		m_viUnknownDisplacementsIndicesReverseMap = oSolver.m_viUnknownDisplacementsIndicesReverseMap;
 		m_viKnownDisplacementsIndicesReverseMap = oSolver.m_viKnownDisplacementsIndicesReverseMap;
 		
 		m_oStiffnessMatrix11 = oSolver.m_oStiffnessMatrix11;
 		m_oStiffnessMatrix12 = oSolver.m_oStiffnessMatrix12;
 		m_oStiffnessMatrix21 = oSolver.m_oStiffnessMatrix21;
 		m_oStiffnessMatrix22 = oSolver.m_oStiffnessMatrix22;
 		return *this;
 	}
 	void FEMSolidSolver::Reset()
 	{
 		FEMSolver::Reset();
 		Initialize();
 	}
 	void FEMSolidSolver::Initialize()
 	{
 		FEMSolver::Initialize();
 		m_oForces1.Reset();
 		m_oForces2.Reset();
 		m_oDisplacements1.Reset();
 		m_oDisplacements2.Reset();

 		m_viUnknownDisplacementsIndices.clear();
 		m_viKnownDisplacementsIndices.clear();
 		m_viUnknownDisplacementsIndicesReverseMap.clear();
 		m_viKnownDisplacementsIndicesReverseMap.clear();
 		
 		m_oStiffnessMatrix11.Reset();
 		m_oStiffnessMatrix12.Reset();
 		m_oStiffnessMatrix21.Reset();
 		m_oStiffnessMatrix22.Reset();
 	}
 	Vector FEMSolidSolver::GetDisplacement(Point* poPoint,unsigned int& iStatus)
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
 	Matrix FEMSolidSolver::GetStress(Point* poPoint,unsigned int& iStatus)
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
 	void FEMSolidSolver::GenerateMatrices()
 	{
 		PrintOnScreen("Generating Matrices");
		vector<FEMNode*>* pvpoNodes = m_poDataStructure->GetNodes();
		vector<FEMElement*>* pvpoElements = m_poDataStructure->GetElements();
  		SparseMatrix oStiffness(m_iDOFCount);
  		unsigned int iSize = pvpoElements->size();
  		Matrix oTempStiffnessMatrix;
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
 	void FEMSolidSolver::PartitionMatrices(const SparseMatrix& oStiffnessMatrix)
 	{
 		unsigned int iSize1 = (unsigned int)m_viUnknownDisplacementsIndices.size();
 		unsigned int iSize2 = (unsigned int)m_viKnownDisplacementsIndices.size();
 
 		m_oStiffnessMatrix11.SetRowsCount(iSize1);
 		m_oStiffnessMatrix12.SetRowsCount(iSize1);
 		m_oStiffnessMatrix21.SetRowsCount(iSize2);
 		m_oStiffnessMatrix22.SetRowsCount(iSize2);
 		
 		double dStiffness = 0.0;
 		unsigned int i = 0;
 		unsigned int j = 0;
 		
 		for(i = 0; i < iSize1 ; i++)
 		{
 			for(j = 0; j < iSize1 ; j++)
 			{
				dStiffness = 0.0;
 				if(oStiffnessMatrix.Get(m_viUnknownDisplacementsIndices[i],m_viUnknownDisplacementsIndices[j],dStiffness))
 				{
 					m_oStiffnessMatrix11.Set(i + 1,j + 1,dStiffness);
 				}
 			}
 		}
 		for(i = 0; i < iSize1 ; i++)
 		{
 			for(j = 0; j < iSize2 ; j++)
 			{
				dStiffness = 0.0;
 				if(oStiffnessMatrix.Get(m_viUnknownDisplacementsIndices[i],m_viKnownDisplacementsIndices[j],dStiffness))
 				{
 					m_oStiffnessMatrix12.Set(i + 1,j + 1,dStiffness);
 				}
 			}
 		}
 		for(i = 0; i < iSize2 ; i++)
 		{
 			for(j = 0; j < iSize1 ; j++)
 			{
				dStiffness = 0.0;
 				if(oStiffnessMatrix.Get(m_viKnownDisplacementsIndices[i],m_viUnknownDisplacementsIndices[j],dStiffness))
 				{
 					m_oStiffnessMatrix21.Set(i + 1,j + 1,dStiffness);
 				}
 			}
 		}
 		for(i = 0; i < iSize2 ; i++)
 		{
 			for(j = 0; j < iSize2 ; j++)
 			{
				dStiffness = 0.0;
 				if(oStiffnessMatrix.Get(m_viKnownDisplacementsIndices[i],m_viKnownDisplacementsIndices[j],dStiffness))
 				{
 					m_oStiffnessMatrix22.Set(i + 1,j + 1,dStiffness);
 				}
 			}
 		}
 		m_oStiffnessMatrix11.BuildPreconditioner();
 	}
 	void FEMSolidSolver::GenerateNodalStressesMatrix()
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
 	void FEMSolidSolver::InitializeMatrices(const double& dTimeStep)
 	{
 	 	PrintOnScreen("Initializing FEM problem");			
 		vector<FEMNode*>* pvpoNodes = m_poDataStructure->GetNodes();
 		unsigned int iSize = pvpoNodes->size();
 		unsigned int i = 0;
 		// get the exact node count
 		m_iDOFCount = 0;
 		for(i = 0; i < iSize ; i++)
 		{
 			m_iDOFCount = m_iDOFCount + pvpoNodes->at(i)->GetDOFCount();
 		}
 		m_viUnknownDisplacementsIndices.reserve(m_iDOFCount);
 		m_viKnownDisplacementsIndices.reserve(m_iDOFCount);
 		// the maximum DOF index cannot be higher than the number of DOFs
 		m_viUnknownDisplacementsIndicesReverseMap.resize(m_iDOFCount);
 		m_viKnownDisplacementsIndicesReverseMap.resize(m_iDOFCount);
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
 	void FEMSolidSolver::Solve(const double& dTime)
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
			oDeltaDisplacement = m_oStiffnessMatrix11.SolveConjugateGradient(oRHS);
			i = i + 1;
			dError = oDeltaDisplacement.GetNorm()/m_oDisplacements1.GetNorm();
			m_oDisplacements1 = m_oDisplacements1 + oDeltaDisplacement;
			// get the unknown forces
			m_oForces2 = m_oStiffnessMatrix21*m_oDisplacements1 + m_oStiffnessMatrix22*m_oDisplacements2;
			UpdateNodalStateVectors();
			GenerateNodalStressesMatrix();
			printf("iteration : %d \t\t error : %e\n",i,dError);
			if(dError < dTolerance)
			{
				break;
			}
		}
 	 	PrintOnScreen("Done");
 	}
 	void FEMSolidSolver::GenerateNodalStateVectors()
 	{
		vector<FEMNode*>* pvpoNodes = m_poDataStructure->GetNodes();
 		unsigned iSize = (unsigned int)pvpoNodes->size();
 		unsigned int i = 0;
 		FEMSolidNode* poNode = NULL;
 		FEMDegreeOfFreedom* poDOF = NULL;
 		unsigned int iIndex = 0;
 		for(i = 0; i < iSize ; i++)
 		{
 			poNode = (FEMSolidNode*)pvpoNodes->at(i);
 			
 			poDOF = poNode->GetXDOF();
 			iIndex = poDOF->GetIndex() - 1;
 			if(poDOF->GetCondition())
 			{
 				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
 				m_oDisplacements2.Set(iIndex,1,poDOF->GetConstraintValue());
 			}
 			else
 			{
 				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
 				m_oForces1.Set(iIndex,1,poDOF->GetConstraintValue());
 			}
 
 			poDOF = poNode->GetYDOF();
 			iIndex = poDOF->GetIndex() - 1;
 			if(poDOF->GetCondition())
 			{
 				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
 				m_oDisplacements2.Set(iIndex,1,poDOF->GetConstraintValue());
 			}
 			else
 			{
 				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
 				m_oForces1.Set(iIndex,1,poDOF->GetConstraintValue());
 			}
 
 			poDOF = poNode->GetZDOF();
 			iIndex = poDOF->GetIndex() - 1;
 			if(poDOF->GetCondition())
 			{
 				iIndex = m_viKnownDisplacementsIndicesReverseMap[iIndex];
 				m_oDisplacements2.Set(iIndex,1,poDOF->GetConstraintValue());
 			}
 			else
 			{
 				iIndex = m_viUnknownDisplacementsIndicesReverseMap[iIndex];
 				m_oForces1.Set(iIndex,1,poDOF->GetConstraintValue());
 			}
 		}
 	}
 	void FEMSolidSolver::UpdateNodalStateVectors()
 	{
 		vector<FEMNode*>* pvpoNodes = m_poDataStructure->GetNodes();
 		unsigned iSize = (unsigned int)pvpoNodes->size();
 		unsigned int i = 0;
 		FEMSolidNode* poNode = NULL;
 		FEMDegreeOfFreedom* poDOF = NULL;
 		unsigned int iIndex = 0;
 		for(i = 0; i < iSize ; i++)
 		{
 			poNode = (FEMSolidNode*)pvpoNodes->at(i);
 			
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
 		}
 	}
 	void FEMSolidSolver::WriteFEMSolutionToParaview(const unsigned int& iStep) const
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




