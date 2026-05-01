// Ahmed M. Hussein
#ifndef FEMPOTENTIALELEMENT_H_
#define FEMPOTENTIALELEMENT_H_


#include "FEMElement.h"
#include "vector"
#include "FEMLoad.h"

using namespace EZ;
using namespace std;

namespace FEMSystem
{
	class FEMPotentialElement : public FEMElement
	{
	public:
		FEMPotentialElement();
		FEMPotentialElement(const FEMPotentialElement& oElement);
		virtual ~FEMPotentialElement();
		virtual FEMPotentialElement& operator=(const FEMPotentialElement& oElement);
		virtual void Reset();
		virtual void SetFacesLoadsSize();
		
		virtual Matrix GetConductionMatrix() const;
		virtual Matrix GetNodalFluxes() const;
		virtual double GetPotential(const Matrix& oNaturalCoordinates) const;
		virtual Vector GetFlux(const Matrix& oNaturalCoordinates) const;
		virtual void SetFacesLoads(vector<FEMLoad*>* pvpoLoads);
		virtual void SetBodyLoad(FEMLoad* poLoad);

		virtual FEMElementType GetType() const;
		virtual unsigned int GetDegreesOfFreedomCount() const;
		virtual void ApplyLoads(const double& dTime);
		virtual void ReadLoads(FILE* fpFile,vector<FEMLoad*>* pvpoLoads);
		virtual void Write(FILE* fpFile) const;
		virtual vector<unsigned int> GetDOFIndices() const;

	private:
		
	protected:
		virtual void Initialize();
 		Matrix GetGradientTransformation(const double& dXi,const double& dEta,const double& dZeta,double& dJacobian) const;
		vector<FEMLoad*> m_vpoFacesLoads;
		FEMLoad* m_poBodyLoad;
		static unsigned int PotentialDOFPerNode;
	};
}


#endif




