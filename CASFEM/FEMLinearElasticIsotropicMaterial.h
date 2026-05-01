// Ahmed M. Hussein

#ifndef FEMLINEARELASTICISOTROPICMATERIAL_H_
#define FEMLINEARELASTICISOTROPICMATERIAL_H_

#include "FEMMaterial.h"

namespace FEMSystem
{
	class FEMLinearElasticIsotropicMaterial : public FEMMaterial
	{
	public:
		FEMLinearElasticIsotropicMaterial();
		virtual ~FEMLinearElasticIsotropicMaterial();
		FEMLinearElasticIsotropicMaterial(const FEMLinearElasticIsotropicMaterial& oMaterial);
		virtual FEMLinearElasticIsotropicMaterial& operator=(const FEMLinearElasticIsotropicMaterial& oMaterial);
		virtual void Reset();
		virtual void Read(FILE* fpFile);
		virtual void Write(FILE* fpFile) const;
		virtual Matrix GetStiffnessMatrix() const;
		virtual Matrix GetConductionMatrix() const;
		virtual double GetMassDensity() const;
		virtual Matrix GetStressTemperatureCoupling() const;
		virtual double GetSpecificHeatCapacity() const;
		virtual double GetReferenceTemperature() const;
		virtual FEMMaterialType GetType() const;
		virtual bool IsPlastic() const;
		virtual Matrix UpdateGaussPoint(FEMGaussPoint* poPoint) const;
		
	private:
	
	protected:
		virtual void Initialize();
		virtual void ConstructInitialMatrices();
		Matrix m_oStiffness;
		Matrix m_oConduction;
		Matrix m_oStressTemperatureCoupling;
		double m_dYoungsModulus;
		double m_dPoissonsRatio;
		double m_dThermalConduction;
		double m_dThermalExpansion;
		double m_dMassDensity;
		double m_dSpecificHeat;
		double m_dReferenceTemperature;
	};
}


#endif



