// Ahmed M. Hussein

#ifndef FEMJ2PLASTICISOTROPICMATERIAL_H_
#define FEMJ2PLASTICISOTROPICMATERIAL_H_

#include "FEMLinearElasticIsotropicMaterial.h"

namespace FEMSystem
{
	class FEMJ2PlasticIsotropicMaterial : public FEMLinearElasticIsotropicMaterial
	{
	public:
		FEMJ2PlasticIsotropicMaterial();
		~FEMJ2PlasticIsotropicMaterial();
		FEMJ2PlasticIsotropicMaterial(const FEMJ2PlasticIsotropicMaterial& oMaterial);
		FEMJ2PlasticIsotropicMaterial& operator=(const FEMJ2PlasticIsotropicMaterial& oMaterial);
		virtual void Reset();
		virtual void Read(FILE* fpFile);
		virtual void Write(FILE* fpFile) const;
		virtual Matrix GetStiffnessMatrix() const;
		virtual FEMMaterialType GetType() const;
		virtual bool IsPlastic() const;
		virtual Matrix UpdateGaussPoint(FEMGaussPoint* poPoint) const;

	private:

	protected:
		virtual void Initialize();
		virtual double GetGamma(const Matrix& oStress,const Matrix& oBackStress,const double& dAccumulatedStrain,Matrix& oN,double& dNormEta) const;
		static Matrix GetDeviatorOperator();
		double m_dInitialYieldStrength;
		double m_dIsotropicHardeningCoefficient;
		double m_dKinematicHardeningCoefficient;
	};
}

#endif


