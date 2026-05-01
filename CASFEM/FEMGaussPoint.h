#ifndef FEMGAUSSPOINT_H_
#define FEMGAUSSPOINT_H_

#include "Matrix.h"
#include "FEMPlasticityState.h"

using namespace EZ;

namespace FEMSystem
{
	class FEMElement;
	class FEMGaussPoint
	{
	public:
		FEMGaussPoint();
		FEMGaussPoint(const FEMGaussPoint& oPoint);
		~FEMGaussPoint();
		FEMGaussPoint& operator=(const FEMGaussPoint& oPoint);
		void Reset();
		void SetCoordinates(const double& dXi,const double& dEta,const double& dZeta);
		void Update();
		void SetElement(FEMElement* poElement);
		void SetWeight(const double& dWeight);
		Matrix GetMaterialStiffnessMatrix();
		Matrix GetStiffnessMatrix();
		Matrix GetStress() const;
		Matrix GetTotalStrain() const;
		void SetStress(const Matrix& oStress);
		FEMPlasticityState* GetPlasticityState();
		Matrix GetVectorizedStress() const;
		void SetVectorizedStress(const Matrix& oStress);

	private:

	protected:
		void Initialize();
		double m_dXi;
		double m_dEta;
		double m_dZeta;
		FEMElement* m_poElement;
		Matrix m_oStrainTransformation;
		double m_dWeight;
		double m_dJacobian;
		
		Matrix m_oStress;
		FEMPlasticityState m_oPlasticityState;
	};
}

#endif


