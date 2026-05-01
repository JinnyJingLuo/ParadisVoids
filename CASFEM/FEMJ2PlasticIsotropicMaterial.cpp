#include "FEMJ2PlasticIsotropicMaterial.h"
#include "Tools.h"
#include "cmath"

using namespace SupportSystem;

namespace FEMSystem
{
	FEMJ2PlasticIsotropicMaterial::FEMJ2PlasticIsotropicMaterial()
	{
		Initialize();
	}
	FEMJ2PlasticIsotropicMaterial::~FEMJ2PlasticIsotropicMaterial()
	{
		Reset();
	}
	FEMJ2PlasticIsotropicMaterial::FEMJ2PlasticIsotropicMaterial(const FEMJ2PlasticIsotropicMaterial& oMaterial)
	{
		*this = oMaterial;
	}
	FEMJ2PlasticIsotropicMaterial& FEMJ2PlasticIsotropicMaterial::operator=(const FEMJ2PlasticIsotropicMaterial& oMaterial)
	{
		FEMLinearElasticIsotropicMaterial::operator=(oMaterial);
		m_dInitialYieldStrength = oMaterial.m_dInitialYieldStrength;
		m_dIsotropicHardeningCoefficient = oMaterial.m_dIsotropicHardeningCoefficient;
		m_dKinematicHardeningCoefficient = oMaterial.m_dKinematicHardeningCoefficient;
		return *this;
	}
	void FEMJ2PlasticIsotropicMaterial::Reset()
	{
		FEMLinearElasticIsotropicMaterial::Reset();
	}
	void FEMJ2PlasticIsotropicMaterial::Read(FILE* fpFile)
	{
		string sRead = GetRealString(500,fpFile);
		sscanf(sRead.c_str(),"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",&m_dYoungsModulus,&m_dPoissonsRatio,&m_dMassDensity,&m_dThermalConduction,&m_dThermalExpansion,&m_dSpecificHeat,&m_dReferenceTemperature,&m_dInitialYieldStrength,&m_dIsotropicHardeningCoefficient,&m_dKinematicHardeningCoefficient);
		ConstructInitialMatrices();
	}
	void FEMJ2PlasticIsotropicMaterial::Write(FILE* fpFile) const
	{
		fprintf(fpFile,"%d\n",GetType());
		fprintf(fpFile,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",m_dYoungsModulus,m_dPoissonsRatio,m_dMassDensity,m_dThermalConduction,m_dThermalExpansion,m_dSpecificHeat,m_dReferenceTemperature,m_dInitialYieldStrength,m_dIsotropicHardeningCoefficient,m_dKinematicHardeningCoefficient);
	}
	Matrix FEMJ2PlasticIsotropicMaterial::GetStiffnessMatrix() const
	{
		return m_oStiffness;
	}
	FEMMaterialType FEMJ2PlasticIsotropicMaterial::GetType() const
	{
		return J2PlasticIsotropicFEMMaterial;
	}
	double FEMJ2PlasticIsotropicMaterial::GetGamma(const Matrix& oStress,const Matrix& oBackStress,const double& dAccumulatedStrain,Matrix& oN,double& dNormEta) const
	{
		// the flow rule expression has to be negative for elasticity
		// f = ||eta|| - sqrt(2/3) (sigma_y + k alpha)
		// eta = s - beta
		// s = deviator(sigma)
		// beta : back stress tensor
		// sigma_y : initial yield strength
		// k : isotropic hardening coefficient
		// all the input matrices are in vector form (6 x 1) 
		// tensor ordering : xx yy zz xy yz zx
		Matrix oEta = oStress;
		double dStressTrace = (oStress.Get(1,1) + oStress.Get(2,1) + oStress.Get(3,1))*(1.0/3.0);
		unsigned int i = 0;
		unsigned int j = 0;
		for(i = 1 ; i <= 3 ; i++)
		{
			oEta.Set(i,1,oEta.Get(i,1) - dStressTrace);
		}
		oEta = oEta - oBackStress;
		dNormEta = sqrt((oEta.GetTranspose()*oEta).Get(1,1));
		double dF = dNormEta - sqrt(2.0/3.0)*(m_dInitialYieldStrength - m_dIsotropicHardeningCoefficient*dAccumulatedStrain);
		// if the trial state is inside the current yield surface, gamma is zero
		if(dF < 0.0)
		{
			return 0.0;
		}
		// otherwise, get gamma, it will be used to update the plasticity state and the stress
		double dHPlusK = m_dIsotropicHardeningCoefficient + m_dKinematicHardeningCoefficient;
		double dMu = m_dYoungsModulus/2.0/(1.0 + m_dPoissonsRatio);
		double dDen = 2.0*dMu + 2.0/3.0*dHPlusK;
		double dNum = dNormEta - sqrt(2.0/3.0)*m_dInitialYieldStrength;
		double dGamma = dNum/dDen;
		oN = oEta*(1.0/dNormEta);
		return dGamma;
	}
	void FEMJ2PlasticIsotropicMaterial::Initialize()
	{
		FEMLinearElasticIsotropicMaterial::Initialize();
		m_dInitialYieldStrength = 0.0;
		m_dIsotropicHardeningCoefficient = 0.0;
		m_dKinematicHardeningCoefficient = 0.0;
	}
	bool FEMJ2PlasticIsotropicMaterial::IsPlastic() const
	{
		return true;
	}
	Matrix FEMJ2PlasticIsotropicMaterial::UpdateGaussPoint(FEMGaussPoint* poPoint) const
	{
		Matrix oTotalStrain = poPoint->GetTotalStrain();
		FEMPlasticityState* poState = poPoint->GetPlasticityState();
		Matrix oPreviousPlasticStrain = poState->GetVectorizedPlasticStrain();
		Matrix oPreviousBackStress = poState->GetVectorizedBackStress();
		double dPreviousAccumulatedStrain = poState->GetAccumulatedStrain();
		Matrix oTrialStress = m_oStiffness*(oTotalStrain - oPreviousPlasticStrain);
		// see if the state is elastic
		Matrix oN;
		double dNormEta = 0.0;
		double dGamma = GetGamma(oTrialStress,oPreviousBackStress,dPreviousAccumulatedStrain,oN,dNormEta);
		double dTolerance = 1.0E-6;
		if(dGamma < dTolerance)
		{
			// elastic state, update the stress and don't touch the plasticity state
			poPoint->SetVectorizedStress(oTrialStress);
			return m_oStiffness;
		}
		// plastic state, update everything given gamma
		Matrix oCurrentPlasticStrain = oPreviousPlasticStrain + oN*dGamma;
		Matrix oCurrentBackStress = oPreviousBackStress + oN*dGamma*(2.0/3.0*m_dKinematicHardeningCoefficient);
		double oCurrentAccumulatedStrain = dPreviousAccumulatedStrain + sqrt(2.0/3.0)*dGamma;
		
		poState->SetVectorizedPlasticStrain(oCurrentPlasticStrain);
		poState->SetVectorizedBackStress(oCurrentBackStress);
		poState->SetAccumulatedStrain(oCurrentAccumulatedStrain);
		
		Matrix oCurrentStress = m_oStiffness*(oTotalStrain - oCurrentPlasticStrain);
		poPoint->SetVectorizedStress(oCurrentStress);
		// finally, return the elastoplastic stiffness matrix, following Simo and Hughes 1998 (3.3), verified
		double dMu = m_dYoungsModulus/2.0/(1.0 + m_dPoissonsRatio);
		double dFactor1 = 2.0*dMu*dGamma/dNormEta;
		double dFactor2 = 1.0/(1.0 + (m_dIsotropicHardeningCoefficient + m_dKinematicHardeningCoefficient)/3.0/dMu) - dFactor1;
		Matrix oStiffness = m_oStiffness - GetDeviatorOperator()*(2.0*dMu*dFactor1) - (oN*oN.GetTranspose())*(2.0*dMu*dFactor2);
		return oStiffness;
	}
	Matrix FEMJ2PlasticIsotropicMaterial::GetDeviatorOperator()
	{
		double dOneOverThree = 1.0/3.0;
		Matrix oOperator(6,6);
		
		oOperator.Set(1,1,1.0 - dOneOverThree);
		oOperator.Set(1,2,-dOneOverThree);
		oOperator.Set(1,3,-dOneOverThree);
		
		oOperator.Set(2,1,-dOneOverThree);
		oOperator.Set(2,2,1.0 - dOneOverThree);
		oOperator.Set(2,3,-dOneOverThree);
		
		oOperator.Set(3,1,-dOneOverThree);
		oOperator.Set(3,2,-dOneOverThree);
		oOperator.Set(3,3,1.0 - dOneOverThree);

		oOperator.Set(4,4,1.0);		
		oOperator.Set(5,5,1.0);
		oOperator.Set(6,6,1.0);
		
		return oOperator;
	}
}


