// Ahmed M. Hussein

#ifndef FEMNODE_H_
#define FEMNODE_H_

#include "GenericNode.h"
#include "FEMDegreeOfFreedom.h"
#include "iostream"
#include "vector"

using namespace std;
using namespace EZ;

namespace FEMSystem
{
	class FEMLoad;
	enum FEMNodeType
	{
		NULLFEMNode = 0,
		PotentialFEMNode = 1,
		SolidFEMNode = 2,
		ThermoMechanicalFEMNode = 3
	};
	
	class FEMNode : public GenericNode
	{
	public:
		FEMNode();
		virtual ~FEMNode();
		FEMNode(const double& dX,const double& dY,const double& dZ);
		FEMNode(const FEMNode& oNode);
		virtual FEMNode& operator=(const FEMNode& oNode);
		virtual void Reset();
		void SetOnSurface(const bool& bOnSurface = true);
		bool IsOnSurface() const;
		virtual FEMNode* Clone() const = 0;
		virtual unsigned int GetDOFCount() const = 0;
		virtual FEMNodeType GetType() const = 0;	
		virtual void ReadNode(FILE* fpFile,vector<FEMLoad*>* pvpoLoads) = 0;
		virtual void WriteNode(FILE* fpFile) const = 0;
		virtual void ResetLoads() = 0;
		virtual void ApplyLoads(const double& dTime) = 0;
		virtual bool IsConstrained() const = 0;
		static FEMNode* CreateNodeByType(FEMNodeType eType);
		static FEMNode* CreateNodeByTypeIndex(const unsigned int& iIndex);
		virtual unsigned int SetDOFIndices(const unsigned int& iCurrentDOFIndex) = 0;
		virtual vector<FEMDegreeOfFreedom*> GetDOFs() = 0;
		
	private:

	protected:
		virtual void Initialize();
		bool m_bOnSurface;
	};
}


#endif


