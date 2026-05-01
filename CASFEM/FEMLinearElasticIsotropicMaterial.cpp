// Ahmed M. Hussein

#include "FEMLinearElasticIsotropicMaterial.h"
#include "Tools.h"

using namespace SupportSystem;


namespace FEMSystem
{
	FEMLinearElasticIsotropicMaterial::FEMLinearElasticIsotropicMaterial()
	{
		Initialize();
	}
	FEMLinearElasticIsotropicMaterial::~FEMLinearElasticIsotropicMaterial()
	{
		Reset();
	}
	FEMLinearElasticIsotropicMaterial::FEMLinearElasticIsotropicMaterial(const FEMLinearElasticIsotropicMaterial& oMaterial)
	{
		*this = oMaterial;
	}
	FEMLinearElasticIsotropicMaterial& FEMLinearElasticIsotropicMaterial::operator=(const FEMLinearElasticIsotropicMaterial& oMaterial)
	{
		FEMMaterial::operator=(oMaterial);
		m_oStiffness = oMaterial.m_oStiffness;
		m_oConduction = oMaterial.m_oConduction;
		m_oStressTemperatureCoupling = oMaterial.m_oStressTemperatureCoupling;
		m_dYoungsModulus = oMaterial.m_dYoungsModulus;
		m_dPoissonsRatio = oMaterial.m_dPoissonsRatio;
		m_dThermalConduction = oMaterial.m_dThermalConduction;
		m_dThermalExpansion = oMaterial.m_dThermalExpansion;
		m_dMassDensity = oMaterial.m_dMassDensity;
		m_dSpecificHeat = oMaterial.m_dSpecificHeat;
		m_dReferenceTemperature = oMaterial.m_dReferenceTemperature;
		return *this;
	}
	void FEMLinearElasticIsotropicMaterial::Reset()
	{
		FEMMaterial::Reset();
		m_oStiffness.Reset();
		m_oConduction.Reset();
		m_oStressTemperatureCoupling.Reset();
	}
	void FEMLinearElasticIsotropicMaterial::Read(FILE* fpFile)
	{
		string sRead = GetRealString(500,fpFile);
		sscanf(sRead.c_str(),"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",&m_dYoungsModulus,&m_dPoissonsRatio,&m_dMassDensity,&m_dThermalConduction,&m_dThermalExpansion,&m_dSpecificHeat,&m_dReferenceTemperature);
		ConstructInitialMatrices();
	}
	void FEMLinearElasticIsotropicMaterial::Write(FILE* fpFile) const
	{
		fprintf(fpFile,"%d\n",GetType());
		fprintf(fpFile,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",m_dYoungsModulus,m_dPoissonsRatio,m_dMassDensity,m_dThermalConduction,m_dThermalExpansion,m_dSpecificHeat,m_dReferenceTemperature);
	}
	Matrix FEMLinearElasticIsotropicMaterial::GetStiffnessMatrix() const
	{
		return m_oStiffness;
	}
	Matrix FEMLinearElasticIsotropicMaterial::GetConductionMatrix() const
	{
		return m_oConduction;
	}
	double FEMLinearElasticIsotropicMaterial::GetMassDensity() const
	{
		return m_dMassDensity;
	}
	Matrix FEMLinearElasticIsotropicMaterial::GetStressTemperatureCoupling() const
	{
		return m_oStressTemperatureCoupling;
	}
	double FEMLinearElasticIsotropicMaterial::GetSpecificHeatCapacity() const
	{
		return m_dSpecificHeat;
	}
	double FEMLinearElasticIsotropicMaterial::GetReferenceTemperature() const
	{
		return m_dReferenceTemperature;
	}
	void FEMLinearElasticIsotropicMaterial::Initialize()
	{
		FEMMaterial::Initialize();
		m_oStiffness.Reset();
		m_oConduction.Reset();
		m_oStressTemperatureCoupling.Reset();
		m_dYoungsModulus = 0.0;
		m_dPoissonsRatio = 0.0;
		m_dThermalConduction = 0.0;
		m_dThermalExpansion = 0.0;
		m_dMassDensity = 0.0;
		m_dSpecificHeat = 0.0;
		m_dReferenceTemperature = 0.0;
	}
	FEMMaterialType FEMLinearElasticIsotropicMaterial::GetType() const
	{
		return LinearElasticIsotropicFEMMaterial;
	}
	void FEMLinearElasticIsotropicMaterial::ConstructInitialMatrices()
	{
		m_oStiffness.SetSize(6,6);
		double dFactor = m_dYoungsModulus/(1.0 + m_dPoissonsRatio)/(1.0 - 2.0*m_dPoissonsRatio);

		m_oStiffness.Set(1,1,dFactor*(1.0 - m_dPoissonsRatio));
		m_oStiffness.Set(1,2,dFactor*m_dPoissonsRatio);
		m_oStiffness.Set(1,3,dFactor*m_dPoissonsRatio);

		m_oStiffness.Set(2,1,dFactor*m_dPoissonsRatio);
		m_oStiffness.Set(2,2,dFactor*(1.0 - m_dPoissonsRatio));
		m_oStiffness.Set(2,3,dFactor*m_dPoissonsRatio);

		m_oStiffness.Set(3,1,dFactor*m_dPoissonsRatio);
		m_oStiffness.Set(3,2,dFactor*m_dPoissonsRatio);
		m_oStiffness.Set(3,3,dFactor*(1.0 - m_dPoissonsRatio));

		// tau = 2 mu epsilon_shear
		m_oStiffness.Set(4,4,dFactor*(1.0 - 2.0*m_dPoissonsRatio));
		m_oStiffness.Set(5,5,dFactor*(1.0 - 2.0*m_dPoissonsRatio));
		m_oStiffness.Set(6,6,dFactor*(1.0 - 2.0*m_dPoissonsRatio));

		m_oConduction.SetSize(3,3);
		m_oConduction.Set(1,1,m_dThermalConduction);
		m_oConduction.Set(1,2,0.0);
		m_oConduction.Set(1,3,0.0);

		m_oConduction.Set(2,1,0.0);
		m_oConduction.Set(2,2,m_dThermalConduction);
		m_oConduction.Set(2,3,0.0);

		m_oConduction.Set(3,1,0.0);
		m_oConduction.Set(3,2,0.0);
		m_oConduction.Set(3,3,m_dThermalConduction);

		m_oStressTemperatureCoupling.SetSize(6,1);
		dFactor = m_dThermalExpansion*m_dYoungsModulus/(1.0 - 2.0*m_dPoissonsRatio);
		m_oStressTemperatureCoupling.Set(1,1,dFactor);
		m_oStressTemperatureCoupling.Set(2,1,dFactor);
		m_oStressTemperatureCoupling.Set(3,1,dFactor);
		m_oStressTemperatureCoupling.Set(4,1,0.0);
		m_oStressTemperatureCoupling.Set(5,1,0.0);
		m_oStressTemperatureCoupling.Set(6,1,0.0);
	}
	bool FEMLinearElasticIsotropicMaterial::IsPlastic() const
	{
		return false;
	}
	Matrix FEMLinearElasticIsotropicMaterial::UpdateGaussPoint(FEMGaussPoint* poPoint) const
	{
		Matrix oWorkingStress = m_oStiffness*poPoint->GetTotalStrain();
		Matrix oStress(3,3);
		oStress.Set(1,1,oWorkingStress.Get(1,1));
		oStress.Set(1,2,oWorkingStress.Get(4,1));
		oStress.Set(1,3,oWorkingStress.Get(6,1));
		oStress.Set(2,1,oWorkingStress.Get(4,1));
		oStress.Set(2,2,oWorkingStress.Get(2,1));
		oStress.Set(2,3,oWorkingStress.Get(5,1));
		oStress.Set(3,1,oWorkingStress.Get(6,1));
		oStress.Set(3,2,oWorkingStress.Get(5,1));
		oStress.Set(3,3,oWorkingStress.Get(3,1));
		poPoint->SetStress(oStress);
		return m_oStiffness;
	}
}


