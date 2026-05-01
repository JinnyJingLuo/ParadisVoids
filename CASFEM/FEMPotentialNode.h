// Ahmed M. Hussein

#ifndef FEMPOTENTIALNODE_H_
#define FEMPOTENTIALNODE_H_

#include "FEMNode.h"
#include "vector"
#include "iostream"
#include "Matrix.h"
#include "FEMLoad.h"
#include "FEMDegreeOfFreedom.h"

using namespace std;
using namespace EZ;

namespace FEMSystem
{
	class FEMPotentialNode : public FEMNode
	{
	public:
		FEMPotentialNode();
		virtual ~FEMPotentialNode();
		FEMPotentialNode(const double& dX,const double& dY,const double& dZ);
		FEMPotentialNode(const FEMPotentialNode& oNode);
		virtual FEMPotentialNode& operator=(const FEMPotentialNode& oNode);
		virtual void Reset();
		
		unsigned int GetDOFCount() const;
		FEMNodeType GetType() const;
		
		FEMDegreeOfFreedom* GetPotentialDOF();
		
		void ResetLoads();
		FEMNode* Clone() const;
		void ApplyLoads(const double& dTime);
		void AddFlux(const double& dFlux);
		void SetLoad(FEMLoad* poLoad);
		FEMLoad* GetLoad() const;
		
		virtual void ReadNode(FILE* fpFile,vector<FEMLoad*>* pvpoLoads);
		virtual void WriteNode(FILE* fpFile) const;
		bool IsConstrained() const;
		double GetPotential() const;
		unsigned int SetDOFIndices(const unsigned int& iCurrentDOFIndex);
		virtual vector<FEMDegreeOfFreedom*> GetDOFs();
		void SetFluxes(const Vector& oFluxes);
		Vector GetFluxes() const;
		
	private:

	protected:
		virtual void Initialize();
		FEMDegreeOfFreedom m_oPotentialDOF;
		FEMLoad* m_poLoad;
		Vector m_oFluxes;
		static unsigned int PotentialDOFCount;
	};
}


#endif



