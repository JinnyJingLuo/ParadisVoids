// Ahmed M. Hussein
#ifndef FEMTHERMOMECHANICALELEMENT_H_
#define FEMTHERMOMECHANICALELEMENT_H_

#include "FEMSolidElement.h"

using namespace EZ;
using namespace std;

namespace FEMSystem
{
	class FEMThermoMechanicalElement : public FEMSolidElement
	{
	public:		
		FEMThermoMechanicalElement();
		FEMThermoMechanicalElement(const FEMThermoMechanicalElement& oElement);
		virtual ~FEMThermoMechanicalElement();
		virtual FEMThermoMechanicalElement& operator=(const FEMThermoMechanicalElement& oElement);
		virtual void Reset();
		virtual void SetFacesLoadsSize();
		
		virtual Matrix GetConductionMatrix() const;	
		virtual Matrix GetStressTemperatureCouplingMatrix() const;
		virtual Matrix GetTemperatureStrainRateCouplingMatrix() const;
		virtual Matrix GetThermalDampingMatrix() const;
		virtual Matrix GetNodalFluxes() const;
		virtual double GetTemperature(const Matrix& oNaturalCoordinates) const;
		virtual Vector GetFlux(const Matrix& oNaturalCoordinates) const;
		void SetFacesLoads(vector<FEMLoad*>* pvpoXLoads,vector<FEMLoad*>* pvpoYLoads,vector<FEMLoad*>* pvpoZLoads,vector<FEMLoad*>* pvpoTLoads);
		void SetBodyLoads(FEMLoad* poXLoad,FEMLoad* poYLoad,FEMLoad* poZLoad,FEMLoad* poTLoad);

		virtual FEMElementType GetType() const;
		virtual unsigned int GetDegreesOfFreedomCount() const;
		virtual unsigned int GetThermalDegreesOfFreedomCount() const;
		virtual void ApplyLoads(const double& dTime);
		virtual void ReadLoads(FILE* fpFile,vector<FEMLoad*>* pvpoLoads);
		virtual void Write(FILE* fpFile) const;
		virtual vector<unsigned int> GetDOFIndices() const;
		
		virtual Matrix GetThermoMechanicalStiffnessMatrix();
		virtual Matrix GetThermoMechanicalDampingMatrix() const;
		virtual Matrix GetThermoMechanicalMassMatrix() const;
		
	private:
		
	protected:
		virtual void Initialize();
		Matrix GetThermalGradientTransformation(const double& dXi,const double& dEta,const double& dZeta,double& dJacobian) const;
		Matrix GetThermalShapeFunctions(const double& dXi,const double& dEta,const double& dZeta) const;
		virtual void ApplySolidLoads(const double& dTime);
		virtual void ApplyThermalLoads(const double& dTime);
		vector<FEMLoad*> m_vpoFacesTLoads;
		
		FEMLoad* m_poBodyTLoad;
		static unsigned int ThermoMechanicalDOFPerNode;
	};
}


#endif


