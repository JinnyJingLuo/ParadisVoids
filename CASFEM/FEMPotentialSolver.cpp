#include "FEMPotentialSolver.h"



 #include "FEMPotentialSolver.h"
 #include "FEMPotentialElement.h"
 #include "FEMPotentialNode.h"
 #include "Tools.h"
 
 using namespace SupportSystem;
 
 namespace FEMSystem
 {
 	FEMPotentialSolver::FEMPotentialSolver()
 	{
 		Initialize();
 	}
 	FEMPotentialSolver::FEMPotentialSolver(const FEMPotentialSolver& oSolver)
 	{
 		*this = oSolver;
 	}
 	FEMPotentialSolver::~FEMPotentialSolver()
 	{
 		Reset();
 	}
	FEMPotentialSolver& FEMPotentialSolver::operator=(const FEMPotentialSolver& oSolver)
	{
	 	FEMSolver::operator=(oSolver);
 		m_oFluxes1 = oSolver.m_oFluxes1;
 		m_oFluxes2 = oSolver.m_oFluxes2;
 		m_oPotentials1 = oSolver.m_oPotentials1;
 		m_oPotentials2 = oSolver.m_oPotentials2;
 		
 		m_viUnknownPotentialsIndices = oSolver.m_viUnknownPotentialsIndices;
 		m_viKnownPotentialsIndices = oSolver.m_viKnownPotentialsIndices;
 		m_viUnknownPotentialsIndicesReverseMap = oSolver.m_viUnknownPotentialsIndicesReverseMap;
 		m_viKnownPotentialsIndicesReverseMap = oSolver.m_viKnownPotentialsIndicesReverseMap;
 		
 		m_oConductionMatrix11 = oSolver.m_oConductionMatrix11;
 		m_oConductionMatrix12 = oSolver.m_oConductionMatrix12;
 		m_oConductionMatrix21 = oSolver.m_oConductionMatrix21;
 		m_oConductionMatrix22 = oSolver.m_oConductionMatrix22;
 		return *this;
	}
	void FEMPotentialSolver::Reset()
	{
	 	FEMSolver::Reset();
 		Initialize();
	}
 	void FEMPotentialSolver::Initialize()
 	{ 		
 		FEMSolver::Initialize();
 		m_oFluxes1.Reset();
 		m_oFluxes2.Reset();
 		m_oPotentials1.Reset();
 		m_oPotentials2.Reset();
 		
 		m_viUnknownPotentialsIndices.clear();
 		m_viKnownPotentialsIndices.clear();
 		m_viUnknownPotentialsIndicesReverseMap.clear();
 		m_viKnownPotentialsIndicesReverseMap.clear();
 		
 		m_oConductionMatrix11.Reset();
 		m_oConductionMatrix12.Reset();
 		m_oConductionMatrix21.Reset();
 		m_oConductionMatrix22.Reset();
 	}
 	double FEMPotentialSolver::GetPotential(Point* poPoint,unsigned int& iStatus)
 	{
 		double dPotential = 0.0;
 		vector<double> vdCoordinates;
 		FEMElement* poElement = NULL;
 		iStatus = DeterminePointLocation(poPoint,poElement,vdCoordinates);
 		Matrix oNaturalCoordinates(3,1);
 		if(iStatus != 0)
 		{
 			oNaturalCoordinates.Set(1,1,vdCoordinates[0]);
 			oNaturalCoordinates.Set(2,1,vdCoordinates[1]);
 			oNaturalCoordinates.Set(3,1,vdCoordinates[2]);
 			dPotential = ((FEMPotentialElement*)poElement)->GetPotential(oNaturalCoordinates);
 		}
 		return dPotential;
 	}
 	Vector FEMPotentialSolver::GetFlux(Point* poPoint,unsigned int& iStatus)
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
 			oFlux = ((FEMPotentialElement*)poElement)->GetFlux(oNaturalCoordinates);
 		}
 		return oFlux;
 	}
 	void FEMPotentialSolver::GenerateMatrices(const unsigned int& iDOFCount)
 	{
 		PrintOnScreen("Generating Matrices");
		vector<FEMNode*>* pvpoNodes = m_poDataStructure->GetNodes();
		vector<FEMElement*>* pvpoElements = m_poDataStructure->GetElements();
  		SparseMatrix oConduction(iDOFCount);
  		unsigned int iSize = pvpoElements->size();
  		Matrix oTempConductionMatrix;
  		vector<unsigned int> viDOFIndices;
  		unsigned int i = 0;
  		unsigned int j = 0;
  		unsigned int k = 0;
		unsigned int iElementDOFCount = 0;
  		for(i = 0; i < iSize ; i++)
  		{
			viDOFIndices = pvpoElements->at(i)->GetDOFIndices();
			iElementDOFCount = (unsigned int)viDOFIndices.size();
  			oTempConductionMatrix = ((FEMPotentialElement*)pvpoElements->at(i))->GetConductionMatrix();
			for(j = 0; j < iElementDOFCount ; j++)
  			{
  				for(k = 0; k < iElementDOFCount ; k++)
  				{
  					oConduction.AddToEntry(viDOFIndices[j],viDOFIndices[k],oTempConductionMatrix.Get(j + 1,k + 1));
  				}
  			}
  		}
 		PartitionMatrices(oConduction);
		PrintOnScreen("Done generating matrices");
 	}
 	void FEMPotentialSolver::PartitionMatrices(const SparseMatrix& oConductionMatrix)
 	{
 		unsigned int iSize1 = (unsigned int)m_viUnknownPotentialsIndices.size();
 		unsigned int iSize2 = (unsigned int)m_viKnownPotentialsIndices.size();
 		char cMessageToPrint[500];
 		sprintf(cMessageToPrint,"Prtitioning global matrices, original size is %d and partition sizes are %d and %d",oConductionMatrix.GetRowsCount(),iSize1,iSize2);
 		string sMessageToPrint = cMessageToPrint;
 		PrintOnScreen(sMessageToPrint);
 
 		m_oConductionMatrix11.SetRowsCount(iSize1);
 		m_oConductionMatrix12.SetRowsCount(iSize1);
 		m_oConductionMatrix21.SetRowsCount(iSize2);
 		m_oConductionMatrix22.SetRowsCount(iSize2);
 		
 		double dConduction = 0.0;
 		unsigned int i = 0;
 		unsigned int j = 0;
 		
 		for(i = 0; i < iSize1 ; i++)
 		{
 			for(j = 0; j < iSize1 ; j++)
 			{
				dConduction = 0.0;
 				if(oConductionMatrix.Get(m_viUnknownPotentialsIndices[i],m_viUnknownPotentialsIndices[j],dConduction))
 				{
 					m_oConductionMatrix11.Set(i + 1,j + 1,dConduction);
 				}
 			}
 		}
 		for(i = 0; i < iSize1 ; i++)
 		{
 			for(j = 0; j < iSize2 ; j++)
 			{
				dConduction = 0.0;
 				if(oConductionMatrix.Get(m_viUnknownPotentialsIndices[i],m_viKnownPotentialsIndices[j],dConduction))
 				{
 					m_oConductionMatrix12.Set(i + 1,j + 1,dConduction);
 				}
 			}
 		}
 		for(i = 0; i < iSize2 ; i++)
 		{
 			for(j = 0; j < iSize1 ; j++)
 			{
				dConduction = 0.0;
 				if(oConductionMatrix.Get(m_viKnownPotentialsIndices[i],m_viUnknownPotentialsIndices[j],dConduction))
 				{
 					m_oConductionMatrix21.Set(i + 1,j + 1,dConduction);
 				}
 			}
 		}
 		for(i = 0; i < iSize2 ; i++)
 		{
 			for(j = 0; j < iSize2 ; j++)
 			{
				dConduction = 0.0;
 				if(oConductionMatrix.Get(m_viKnownPotentialsIndices[i],m_viKnownPotentialsIndices[j],dConduction))
 				{
 					m_oConductionMatrix22.Set(i + 1,j + 1,dConduction);
 				}
 			}
 		}
 		PrintOnScreen("Done partitioning matrices");
 	}
 	void FEMPotentialSolver::GenerateNodalFluxesMatrix()
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
 			oElementNodalFluxes = ((FEMPotentialElement*)pvpoElements->at(i))->GetNodalFluxes();
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
				sprintf(cMessageToPrint,"node %d has never ocurred in stress distribution",i - 1);
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
 			((FEMPotentialNode*)pvpoNodes->at(i - 1))->SetFluxes(oTempFlux);
  		}
 	}
 	void FEMPotentialSolver::InitializeMatrices(const double& dTimeStep)
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
 		m_viUnknownPotentialsIndices.reserve(iDOFCount);
 		m_viKnownPotentialsIndices.reserve(iDOFCount);
 		// the maximum DOF index cannot be higher than the number of DOFs
 		m_viUnknownPotentialsIndicesReverseMap.resize(iDOFCount);
 		m_viKnownPotentialsIndicesReverseMap.resize(iDOFCount);
 		FEMPotentialNode* poNode = NULL;
 		unsigned int iIndex = 0;
 		unsigned int iKnownIndicesCount = 0;
 		unsigned int iUnknownIndicesCount = 0;
 		for(i = 0; i < iSize ; i++)
 		{
 			poNode = (FEMPotentialNode*)pvpoNodes->at(i);
 			
 			iIndex = poNode->GetPotentialDOF()->GetIndex();
 			if(poNode->GetPotentialDOF()->GetCondition())
 			{
 				iKnownIndicesCount = iKnownIndicesCount + 1;
 				m_viKnownPotentialsIndices.push_back(iIndex);
 				m_viKnownPotentialsIndicesReverseMap[iIndex - 1] = iKnownIndicesCount;
 				m_viUnknownPotentialsIndicesReverseMap[iIndex - 1] = 0;
 			}
 			else
 			{
 				iUnknownIndicesCount = iUnknownIndicesCount + 1;
 				m_viUnknownPotentialsIndices.push_back(iIndex);
 				m_viUnknownPotentialsIndicesReverseMap[iIndex - 1] = iUnknownIndicesCount;
 				m_viKnownPotentialsIndicesReverseMap[iIndex - 1] = 0;
 			}
 		}

 		GenerateMatrices(iDOFCount);
 		// now resize the potential and flux arrays
 		unsigned int iSize1 = (unsigned int)m_viUnknownPotentialsIndices.size();
 		unsigned int iSize2 = (unsigned int)m_viKnownPotentialsIndices.size();
 		m_oFluxes1.SetSize(iSize1,1);
 		m_oFluxes2.SetSize(iSize2,1);
 		m_oPotentials1.SetSize(iSize1,1);
 		m_oPotentials2.SetSize(iSize2,1);
		// no need to resize these matrices again, they should be directly used
 	 	PrintOnScreen("Done initializing");
 	}
 	void FEMPotentialSolver::Solve(const double& dTime)
 	{
 	 	PrintOnScreen("Generating nodal forces and displacements");
		GenerateNodalStateVectors(dTime);
 	 	PrintOnScreen("Done");
		PrintOnScreen("Solving the system");	
		// get the unknown potentials
		Matrix oRHS = m_oFluxes1 - m_oConductionMatrix12*m_oPotentials2;
 		m_oPotentials1 = m_oConductionMatrix11.SolveConjugateGradient(oRHS);
 		// get the unknown fluxes
 		m_oFluxes2 = m_oConductionMatrix21*m_oPotentials1 + m_oConductionMatrix22*m_oPotentials2;
		UpdateNodalStateVectors();
 		PrintOnScreen("Done");	
 	 	PrintOnScreen("Calculating nodal stresses");
 		GenerateNodalFluxesMatrix();
 	 	PrintOnScreen("Done");
 	}
 	void FEMPotentialSolver::GenerateNodalStateVectors(const double& dTime)
 	{
		vector<FEMNode*>* pvpoNodes = m_poDataStructure->GetNodes();
 		unsigned iSize = (unsigned int)pvpoNodes->size();
 		unsigned int i = 0;
 		FEMPotentialNode* poNode = NULL;
 		FEMDegreeOfFreedom* poDOF = NULL;
 		unsigned int iIndex = 0;
 		for(i = 0; i < iSize ; i++)
 		{
 			poNode = (FEMPotentialNode*)pvpoNodes->at(i);
 			
 			poDOF = poNode->GetPotentialDOF();
 			iIndex = poDOF->GetIndex() - 1;
 			if(poDOF->GetCondition())
 			{
 				iIndex = m_viKnownPotentialsIndicesReverseMap[iIndex];
 				m_oPotentials2.Set(iIndex,1,poDOF->GetConstraintValue());
 			}
 			else
 			{
 				iIndex = m_viUnknownPotentialsIndicesReverseMap[iIndex];
 				m_oFluxes1.Set(iIndex,1,poDOF->GetConstraintValue());
 			}
 		}
 	}
 	void FEMPotentialSolver::UpdateNodalStateVectors()
 	{
 		vector<FEMNode*>* pvpoNodes = m_poDataStructure->GetNodes();
 		unsigned iSize = (unsigned int)pvpoNodes->size();
 		unsigned int i = 0;
 		FEMPotentialNode* poNode = NULL;
 		FEMDegreeOfFreedom* poDOF = NULL;
 		unsigned int iIndex = 0;
 		for(i = 0; i < iSize ; i++)
 		{
 			poNode = (FEMPotentialNode*)pvpoNodes->at(i);
 			
 			poDOF = poNode->GetPotentialDOF();
 			iIndex = poDOF->GetIndex() - 1;
 			if(poDOF->GetCondition())
 			{
 				iIndex = m_viKnownPotentialsIndicesReverseMap[iIndex];
 				poDOF->SetSecondaryValue(m_oFluxes2.Get(iIndex,1));
 			}
 			else
 			{
 				iIndex = m_viUnknownPotentialsIndicesReverseMap[iIndex];
 				poDOF->SetPrimaryValue(m_oPotentials1.Get(iIndex,1));
 			}
 		}
 	}
 	void FEMPotentialSolver::WriteFEMSolutionToParaview(const unsigned int& iStep) const
 	{
 	
 	}
 }




