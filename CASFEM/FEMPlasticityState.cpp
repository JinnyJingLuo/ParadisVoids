#include "FEMPlasticityState.h"

namespace FEMSystem
{
	FEMPlasticityState::FEMPlasticityState()
	{
		Initialize();
	}
	FEMPlasticityState::FEMPlasticityState(const FEMPlasticityState& oState)
	{
		*this = oState;
	}
	FEMPlasticityState::~FEMPlasticityState()
	{
		Reset();
	}
	FEMPlasticityState& FEMPlasticityState::operator=(const FEMPlasticityState& oState)
	{
		m_oPlasticStrain = oState.m_oPlasticStrain;
		m_oBackStress = oState.m_oBackStress;
		m_dAccumulatedStrain = oState.m_dAccumulatedStrain;
		return *this;
	}
	FEMPlasticityState FEMPlasticityState::operator*(const double& dFactor) const
	{
		FEMPlasticityState oResultState;
		oResultState.m_oPlasticStrain = m_oPlasticStrain*dFactor;
		oResultState.m_oBackStress = m_oBackStress*dFactor;
		oResultState.m_dAccumulatedStrain = m_dAccumulatedStrain*dFactor;
		return oResultState;
	}
	FEMPlasticityState FEMPlasticityState::operator+(const FEMPlasticityState& oState) const
	{
		FEMPlasticityState oResultState;
		oResultState.m_oPlasticStrain = m_oPlasticStrain + oState.m_oPlasticStrain;
		oResultState.m_oBackStress = m_oBackStress + oState.m_oBackStress;
		oResultState.m_dAccumulatedStrain = m_dAccumulatedStrain + oState.m_dAccumulatedStrain;
		return oResultState;
	}
	void FEMPlasticityState::Reset()
	{
		Initialize();
	}
	void FEMPlasticityState::SetPlasticStrain(const Matrix& oStrain)
	{
		m_oPlasticStrain = oStrain;
	}
	void FEMPlasticityState::SetBackStress(const Matrix& oBackStress)
	{
		m_oBackStress = oBackStress;
	}
	void FEMPlasticityState::SetAccumulatedStrain(const double& dStrain)
	{
		m_dAccumulatedStrain = dStrain;
	}
	const Matrix& FEMPlasticityState::GetPlasticStrain() const
	{
		return m_oPlasticStrain;
	}
	const Matrix& FEMPlasticityState::GetBackStress() const
	{
		return m_oBackStress;
	}
	double FEMPlasticityState::GetAccumulatedStrain() const
	{
		return m_dAccumulatedStrain;
	}
	Matrix FEMPlasticityState::GetVectorizedPlasticStrain() const
	{
		Matrix oResult(6,1);
		oResult.Set(1,1,m_oPlasticStrain.Get(1,1));
		oResult.Set(2,1,m_oPlasticStrain.Get(2,2));
		oResult.Set(3,1,m_oPlasticStrain.Get(3,3));
		oResult.Set(4,1,m_oPlasticStrain.Get(1,2));
		oResult.Set(5,1,m_oPlasticStrain.Get(2,3));
		oResult.Set(6,1,m_oPlasticStrain.Get(3,1));
		return oResult;
	}
	Matrix FEMPlasticityState::GetVectorizedBackStress() const
	{
		Matrix oResult(6,1);
		oResult.Set(1,1,m_oBackStress.Get(1,1));
		oResult.Set(2,1,m_oBackStress.Get(2,2));
		oResult.Set(3,1,m_oBackStress.Get(3,3));
		oResult.Set(4,1,m_oBackStress.Get(1,2));
		oResult.Set(5,1,m_oBackStress.Get(2,3));
		oResult.Set(6,1,m_oBackStress.Get(3,1));
		return oResult;
	}
	void FEMPlasticityState::SetVectorizedPlasticStrain(const Matrix& oStrain)
	{
		m_oPlasticStrain.Set(1,1,oStrain.Get(1,1));
		m_oPlasticStrain.Set(1,2,oStrain.Get(4,1));
		m_oPlasticStrain.Set(1,3,oStrain.Get(6,1));
		
		m_oPlasticStrain.Set(2,1,oStrain.Get(4,1));
		m_oPlasticStrain.Set(2,2,oStrain.Get(2,1));
		m_oPlasticStrain.Set(2,3,oStrain.Get(5,1));
		
		m_oPlasticStrain.Set(3,1,oStrain.Get(6,1));
		m_oPlasticStrain.Set(3,2,oStrain.Get(5,1));
		m_oPlasticStrain.Set(3,3,oStrain.Get(3,1));
	}
	void FEMPlasticityState::SetVectorizedBackStress(const Matrix& oBackStress)
	{
		m_oBackStress.Set(1,1,oBackStress.Get(1,1));
		m_oBackStress.Set(1,2,oBackStress.Get(4,1));
		m_oBackStress.Set(1,3,oBackStress.Get(6,1));
		
		m_oBackStress.Set(2,1,oBackStress.Get(4,1));
		m_oBackStress.Set(2,2,oBackStress.Get(2,1));
		m_oBackStress.Set(2,3,oBackStress.Get(5,1));
		
		m_oBackStress.Set(3,1,oBackStress.Get(6,1));
		m_oBackStress.Set(3,2,oBackStress.Get(5,1));
		m_oBackStress.Set(3,3,oBackStress.Get(3,1));
	}
	void FEMPlasticityState::Initialize()
	{
		m_oPlasticStrain.SetSize(3,3);
		m_oBackStress.SetSize(3,3);
		m_dAccumulatedStrain = 0.0;
	}
}


