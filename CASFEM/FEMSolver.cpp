// Ahmed M. Hussein

#include "FEMSolver.h"
#include "Tools.h"
#include "cmath"
#include "float.h"
#include "FEMPotentialSolver.h"
#include "FEMSolidSolver.h"
#include "FEMThermoMechanicalSolver.h"
#include "FEMExplicitDynamicsSolidSolver.h"
#include "FEMImplicitDynamicsSolidSolver.h"
#include "FEMImplicitDynamicsThermoMechanicalSolver.h"
#include "stack"
 
using namespace SupportSystem;
using namespace EZ;
using namespace std;

namespace FEMSystem
{
	FEMSolver::~FEMSolver()
	{
		Reset();
	}
	FEMSolver& FEMSolver::operator=(const FEMSolver& oSolver)
	{
		m_poDataStructure = oSolver.m_poDataStructure;
		m_iDOFCount = oSolver.m_iDOFCount;
		m_dTimeStep = oSolver.m_dTimeStep;
		return *this;
	}
	void FEMSolver::Reset()
	{
		Initialize();
	}
	void FEMSolver::Initialize()
	{
		m_poDataStructure = NULL;
		m_iDOFCount = 0;
		m_dTimeStep = 0.0;
	}
	void FEMSolver::SetDataStructure(MainDataStructure* poStructure)
	{
		m_poDataStructure = poStructure;
	}
	FEMSolver* FEMSolver::CreateSolverByPhysics(FEMSolverType eType,const ProblemType& eProblemType)
	{
		if(eType == NullFEMSolver)
		{
			return NULL;
		}
		else if(eType == PotentialFEMSolver)
		{
			return new FEMPotentialSolver;
		}
		else if(eType == SolidFEMSolver)
		{
			if(eProblemType == StaticProblem)
			{
				return new FEMSolidSolver;
			}
			else if(eProblemType == ExplicitDynamicsProblem)
			{
				return new FEMExplicitDynamicsSolidSolver;
			}
			else if(eProblemType == ImplicitDynamicsProblem)
			{
				return new FEMImplicitDynamicsSolidSolver;
			}
			return NULL;
		}
		else if(eType == ThermoMechanicalFEMSolver)
		{
			if(eProblemType == StaticProblem)
			{
				return new FEMThermoMechanicalSolver;
			}
			else
			{
				return new FEMImplicitDynamicsThermoMechanicalSolver;
			}
			return NULL;
		}
		return NULL;
	}
	FEMSolver* FEMSolver::CreateSolverByPhysicsIndex(const unsigned int& iIndex,const ProblemType& eProblemType)
	{
		if(iIndex == 1)
		{
			return CreateSolverByPhysics(PotentialFEMSolver,eProblemType);
		}
		if(iIndex == 2)
		{
			return CreateSolverByPhysics(SolidFEMSolver,eProblemType);
		}
		if(iIndex == 3)
		{
			return CreateSolverByPhysics(ThermoMechanicalFEMSolver,eProblemType);
		}
		return NULL;
	}
	double FEMSolver::GetTimeStep() const
	{
		return m_dTimeStep;
	}
	string FEMSolver::GetOutputFileName(const unsigned int& iStep) const
	{
		char cWrite[512];
		string sOutputBaseFileName = m_poDataStructure->GetOutputBaseFileName();
		sprintf(cWrite,"%s_%d",sOutputBaseFileName.c_str(),iStep);
		string sOutputFileName = cWrite;
		return sOutputFileName;
	}
	void FEMSolver::WriteFEMSolution(const unsigned int& iStep) const
	{
		unsigned int i = 0;
		unsigned int iSize = 0;
		string sOutputFileName = GetOutputFileName(iStep);
		FILE* fpFile = NULL;
		fpFile = fopen(sOutputFileName.c_str(),"w");
		// write the header, loads and materials
		m_poDataStructure->WriteHeader(fpFile);
		m_poDataStructure->WriteLoads(fpFile);
		m_poDataStructure->WriteMaterials(fpFile);
		// write the nodes
		iSize = m_poDataStructure->GetNodesCount();
		fprintf(fpFile,"* Nodes :\n");
		fprintf(fpFile,"%d\n",iSize);
		for(i = 0; i < iSize ; i++)
		{
			m_poDataStructure->GetNode(i)->WriteNode(fpFile);
		}
		
		// write the element connectivity
		iSize = m_poDataStructure->GetElementsCount();
		fprintf(fpFile,"* Elements :\n");
		fprintf(fpFile,"%d\n",iSize);
		for(i = 0 ; i < iSize ; i++)
		{
			m_poDataStructure->GetElement(i)->Write(fpFile);
		}
		fclose(fpFile);
	}
	int FEMSolver::DeterminePointLocationNaively(Point* poPoint,FEMElement*& poElement,vector<double>& vdNaturalCoordinates) const
	{
		unsigned int i = 0;
		vector<FEMElement*>* pvpoElements = m_poDataStructure->GetElements();
		unsigned int iSize = (unsigned int)pvpoElements->size();
		
		vector<double> vdXi;
		vector<double> vdEta;
		vector<double> vdZeta;
		vector<double> vdNaturalCoordinatesDistance;
		vector<FEMElement*> vpoCloseElements;
	
		vdXi.reserve(iSize);
		vdEta.reserve(iSize);
		vdZeta.reserve(iSize);
		vdNaturalCoordinatesDistance.reserve(iSize);
		vpoCloseElements.reserve(iSize);
		// initialize output
		vdNaturalCoordinates.clear();
		poElement = NULL;
		for(i = 0; i < iSize ; i++)
		{
			if(pvpoElements->at(i)->GetGeometry()->IsPointInside(poPoint,vdNaturalCoordinates))
			{
				vdXi.push_back(vdNaturalCoordinates[0]);
				vdEta.push_back(vdNaturalCoordinates[1]);
				vdZeta.push_back(vdNaturalCoordinates[2]);
				vdNaturalCoordinatesDistance.push_back(vdNaturalCoordinates[0]*vdNaturalCoordinates[0] + vdNaturalCoordinates[1]*vdNaturalCoordinates[1] + vdNaturalCoordinates[2]*vdNaturalCoordinates[2]);
				vpoCloseElements.push_back(pvpoElements->at(i));
			}
		}
	
		double dTemp = 0.0;
		double dMinimumDistance = DBL_MAX;
		unsigned int iElementIndex = 0;
		iSize = (unsigned int)vdNaturalCoordinatesDistance.size();
		int iStatus = 0;
		if(iSize != 0)
		{
			for(i = 0; i < iSize ; i++)
			{
				dTemp = vdNaturalCoordinatesDistance[i];
				if(dTemp < dMinimumDistance)
				{
					dMinimumDistance = dTemp;
					iElementIndex = i;
				}
			}
			// get the stress from the element with the least natural coordinates distance
			vdNaturalCoordinates.resize(3);
			vdNaturalCoordinates[0] = vdXi[iElementIndex];
			vdNaturalCoordinates[1] = vdEta[iElementIndex];
			vdNaturalCoordinates[2] = vdZeta[iElementIndex];
			iStatus = 1;
			poElement = vpoCloseElements[iElementIndex];
		}
	
		vdXi.clear();
		vdEta.clear();
		vdZeta.clear();
		vdNaturalCoordinatesDistance.clear();
		vpoCloseElements.clear();
		if(iStatus == 1)
		{
			return iStatus;
		}
		
		///////////////////////////////////////////////////
		// point couldn't be found in any of the elements
		// get an approximate stress
		///////////////////////////////////////////////////
		vdNaturalCoordinates.clear();
		poElement = NULL;
		iSize = (unsigned int)pvpoElements->size();
		vpoCloseElements.reserve(iSize);
		for(i = 0; i < iSize ; i++)
		{
			if(pvpoElements->at(i)->GetGeometry()->IsInAxisAlignedBoundingBox(poPoint))
			{
				vpoCloseElements.push_back(pvpoElements->at(i));
			}
		}
		iSize = (unsigned int)vpoCloseElements.size();
	
		if(iSize == 0)
		{
			iStatus = 0;
			return iStatus;
		}
	
		dMinimumDistance = DBL_MAX;
		for(i = 0; i < iSize ; i++)
		{
			dTemp = vpoCloseElements[i]->GetGeometry()->GetDistanceToElementCenterPoint(poPoint);
			if(dTemp < dMinimumDistance)
			{
				dMinimumDistance = dTemp;
				poElement = vpoCloseElements[i];
			}
		}
		vdNaturalCoordinates = poElement->GetGeometry()->GetNearestNodeNaturalCoordinates(poPoint);
		iStatus = 2;
		vpoCloseElements.clear();
		return iStatus;
	}
	int FEMSolver::DeterminePointLocation(Point* poPoint,FEMElement*& poElement,vector<double>& vdNaturalCoordinates) const
	{
		BSPTreeNode<FEMElement*>* poParentNode = m_poDataStructure->GetElementsBSPTree();
		stack< BSPTreeNode<FEMElement*>* > oNodesStack;
		oNodesStack.push(poParentNode);
		int iClassification = 0;
		list<FEMElement*>* plpoElements = NULL;
		list<FEMElement*>::iterator liElements;
		int iStatus = 0;
		while(!oNodesStack.empty())
		{
			poParentNode = oNodesStack.top();
			oNodesStack.pop();
			if(poParentNode->IsLeaf())
			{
				// if the node is a leaf, get the list of elements it has and look for the point
				// inside of them
				plpoElements = poParentNode->GetItems();
				for(liElements = plpoElements->begin() ; liElements != plpoElements->end() ; liElements++)
				{
					if((*liElements)->GetGeometry()->IsPointInside(poPoint,vdNaturalCoordinates))
					{
						poElement = (*liElements);
						iStatus = 1;
						return iStatus;
					}
				}
			}
			else
			{
				iClassification = poParentNode->Classify(poPoint);
				// push the nodes into the stack based on the point classification so that the 
				// most likely node is accessed first
				if(iClassification < 1)
				{
					oNodesStack.push(poParentNode->GetRightChild());
					oNodesStack.push(poParentNode->GetLeftChild());
				}
				else
				{
					oNodesStack.push(poParentNode->GetLeftChild());
					oNodesStack.push(poParentNode->GetRightChild());
				}
			}
		}
		return iStatus;
	}
}





