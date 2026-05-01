// Ahmed M. Hussein
#ifndef FEMTHERMOMECHANICALNODE_H_
#define FEMTHERMOMECHANICALNODE_H_

#include "FEMSolidNode.h"
#include "vector"
#include "iostream"
#include "Matrix.h"
#include "FEMLoad.h"
#include "FEMDegreeOfFreedom.h"

using namespace std;
using namespace EZ;

namespace FEMSystem
{
	class FEMThermoMechanicalNode : public FEMSolidNode
	{
	public:
		FEMThermoMechanicalNode();
		virtual ~FEMThermoMechanicalNode();
		FEMThermoMechanicalNode(const double& dX,const double& dY,const double& dZ);
		FEMThermoMechanicalNode(const FEMThermoMechanicalNode& oNode);
		virtual FEMThermoMechanicalNode& operator=(const FEMThermoMechanicalNode& oNode);
		virtual void Reset();
		
		virtual unsigned int GetDOFCount() const;
		virtual FEMNodeType GetType() const;
		
		FEMDegreeOfFreedom* GetTDOF();
		
		virtual void ResetLoads();
		virtual FEMNode* Clone() const;
		virtual void ApplyLoads(const double& dTime);
		void AddFlux(const double& dFlux);
		void SetTLoad(FEMLoad* poLoad);
		FEMLoad* GetTLoad() const;
		
		virtual void ReadNode(FILE* fpFile,vector<FEMLoad*>* pvpoLoads);	
		virtual void WriteNode(FILE* fpFile) const;
		
		double GetTemperature() const;
		virtual unsigned int SetDOFIndices(const unsigned int& iCurrentDOFIndex);
		virtual vector<FEMDegreeOfFreedom*> GetDOFs();
		void SetFluxes(const Vector& oFluxes);
		Vector GetFluxes() const;
		
		void SetHeatingRate(const double& dRate);
		void SetHeatingAcceleration(const double& dAcceleration);
		double GetHeatingRate();
		double GetHeatingAcceleration();
		
	private:

	protected:
		virtual void Initialize();
		FEMDegreeOfFreedom m_oTDOF;
		FEMLoad* m_poTLoad;
		Vector m_oFluxes;
		double m_dHeatingRate;
		double m_dHeatingAcceleration;
		static unsigned int ThermoMechanicalDOFCount;
	};
}


#endif


