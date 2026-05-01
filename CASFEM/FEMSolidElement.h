// Ahmed M. Hussein

#ifndef FEMSOLIDELEMENT_H_
#define FEMSOLIDELEMENT_H_

#include "FEMElement.h"
#include "vector"
#include "FEMLoad.h"

using namespace EZ;
using namespace std;

namespace FEMSystem
{
	class FEMSolidElement : public FEMElement
	{
	public:
		FEMSolidElement();
		FEMSolidElement(const FEMSolidElement& oElement);
		virtual ~FEMSolidElement();
		virtual FEMSolidElement& operator=(const FEMSolidElement& oElement);
		virtual void Reset();
		virtual void SetFacesLoadsSize();
		
		virtual Matrix GetStiffnessMatrix();
		virtual Matrix GetMassMatrix() const;
		virtual Matrix GetLumpedMassMatrix() const;
		virtual Matrix GetNodalStresses() const;
		virtual Vector GetDisplacement(const Matrix& oNaturalCoordinates) const;
		virtual Matrix GetStress(const Matrix& oNaturalCoordinates) const;
		void SetFacesLoads(vector<FEMLoad*>* pvpoXLoads,vector<FEMLoad*>* pvpoYLoads,vector<FEMLoad*>* pvpoZLoads);
		void SetBodyLoads(FEMLoad* poXLoad,FEMLoad* poYLoad,FEMLoad* poZLoad);

		virtual FEMElementType GetType() const;
		virtual unsigned int GetDegreesOfFreedomCount() const;
		virtual unsigned int GetSolidDegreesOfFreedomCount() const;
		virtual void ApplyLoads(const double& dTime);
		virtual void ReadLoads(FILE* fpFile,vector<FEMLoad*>* pvpoLoads);
		virtual void Write(FILE* fpFile) const;
		virtual vector<unsigned int> GetDOFIndices() const;
		Matrix GetStrainTransformation(const double& dXi,const double& dEta,const double& dZeta,double& dJacobian) const;

	private:
		
	protected:
		virtual void Initialize();
 		Matrix GetInertiaTransformation(const double& dXi,const double& dEta,const double& dZeta,double& dJacobian) const;
		vector<FEMLoad*> m_vpoFacesXLoads;
		vector<FEMLoad*> m_vpoFacesYLoads;
		vector<FEMLoad*> m_vpoFacesZLoads;
		
		FEMLoad* m_poBodyXLoad;
		FEMLoad* m_poBodyYLoad;
		FEMLoad* m_poBodyZLoad;
		static unsigned int SolidDOFPerNode;
	};
}


#endif




