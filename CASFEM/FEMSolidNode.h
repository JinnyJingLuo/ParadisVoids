// Ahmed M. Hussein
#ifndef FEMSOLIDNODE_H_
#define FEMSOLIDNODE_H_

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
	class FEMSolidNode : public FEMNode
	{
	public:
		FEMSolidNode();
		virtual ~FEMSolidNode();
		FEMSolidNode(const double& dX,const double& dY,const double& dZ);
		FEMSolidNode(const FEMSolidNode& oNode);
		virtual FEMSolidNode& operator=(const FEMSolidNode& oNode);
		virtual void Reset();
		
		virtual unsigned int GetDOFCount() const;
		virtual FEMNodeType GetType() const;
		
		FEMDegreeOfFreedom* GetXDOF();
		FEMDegreeOfFreedom* GetYDOF();
		FEMDegreeOfFreedom* GetZDOF();
		
		virtual void ResetLoads();
		virtual FEMNode* Clone() const;
		virtual void ApplyLoads(const double& dTime);
		virtual void AddForce(Vector oForce);
		void SetStresses(const Matrix& oStresses);
		Matrix GetStresses() const;
		void SetXLoad(FEMLoad* poLoad);
		void SetYLoad(FEMLoad* poLoad);
		void SetZLoad(FEMLoad* poLoad);
		FEMLoad* GetXLoad() const;
		FEMLoad* GetYLoad() const;
		FEMLoad* GetZLoad() const;
		
		virtual void ReadNode(FILE* fpFile,vector<FEMLoad*>* pvpoLoads);
		virtual void WriteNode(FILE* fpFile) const;
		bool IsConstrained() const;
		
		Vector GetDisplacement() const;
		double GetXDisplacement() const;
		double GetYDisplacement() const;
		double GetZDisplacement() const;
		Vector GetForce() const;
		double GetXForce() const;
		double GetYForce() const;
		double GetZForce() const;
		virtual unsigned int SetDOFIndices(const unsigned int& iCurrentDOFIndex);
		virtual vector<FEMDegreeOfFreedom*> GetDOFs();
		
		void SetVelocity(const Vector& oVelocity);
		void SetAcceleration(const Vector& oAcceleration);
		Vector GetVelocity();
		Vector GetAcceleration();
		
	private:

	protected:
		virtual void Initialize();
		FEMDegreeOfFreedom m_oXDOF;
		FEMDegreeOfFreedom m_oYDOF;
		FEMDegreeOfFreedom m_oZDOF;
		Matrix m_oStresses;
		FEMLoad* m_poXLoad;
		FEMLoad* m_poYLoad;
		FEMLoad* m_poZLoad;
		Vector m_oVelocity;
		Vector m_oAcceleration;
		static unsigned int SolidDOFCount;
	};
}


#endif


