// Ahmed M. Hussein

#include "FEMLoad.h"
#include "FEMConstantLoad.h"
#include "FEMTimeSinusoidalLoad.h"
#include "FEMXTorsionLoad.h"
#include "FEMYTorsionLoad.h"
#include "FEMTimeLinearLoad.h"
#include "FEMAxialTimeLinearTotalStrainLoad.h"
#include "FEMAxialTimeSinusoidalTotalStrainLoad.h"
#include "FEMStepLoad.h"
#include "FEMPulseTrainLoad.h"

namespace FEMSystem
{
	FEMLoad::~FEMLoad()
	{
		Reset();
	}
	FEMLoad& FEMLoad::operator=(const FEMLoad& oLoad)
	{
		m_iID = oLoad.m_iID;
		return *this;
	}
	void FEMLoad::Reset()
	{
		Initialize();
	}
	void FEMLoad::Initialize()
	{
		m_iID = 0;
	}
	FEMLoad* FEMLoad::GenerateLoadByType(const FEMLoadTypes& eType)
	{
		if(eType == NullFEMLoad)
		{
			return NULL;
		}
		else if(eType == ConstantFEMLoad)
		{
			return new FEMConstantLoad;
		}
		else if(eType == TimeSinusoidalFEMLoad)
		{
			return new FEMTimeSinusoidalLoad;
		}
		else if(eType == XTorsionFEMLoad)
		{
			return new FEMXTorsionLoad;
		}
		else if(eType == YTorsionFEMLoad)
		{
			return new FEMYTorsionLoad;
		}
		else if(eType == TimeLinearFEMLoad)
		{
			return new FEMTimeLinearLoad;
		}
		else if(eType == AxialTimeLinearTotalStrainFEMLoad)
		{
			return new FEMAxialTimeLinearTotalStrainLoad;
		}
		else if(eType == AxialTimeSinusoidalTotalStrainFEMLoad)
		{
			return new FEMAxialTimeSinusoidalTotalStrainLoad;
		}
		else if(eType == StepLoad)
		{
			return new FEMStepLoad;
		}
		else if(eType == PulseTrainLoad)
		{
			return new FEMPulseTrainLoad;
		}
		return NULL;
	}
	FEMLoad* FEMLoad::GenerateLoadByTypeIndex(const unsigned int& iIndex)
	{
		if(iIndex == 1)
		{
			return GenerateLoadByType(ConstantFEMLoad);
		}
		else if(iIndex == 2)
		{
			return GenerateLoadByType(TimeSinusoidalFEMLoad);
		}
		else if(iIndex == 3)
		{
			return GenerateLoadByType(XTorsionFEMLoad);
		}
		else if(iIndex == 4)
		{
			return GenerateLoadByType(YTorsionFEMLoad);
		}
		else if(iIndex == 5)
		{
			return GenerateLoadByType(TimeLinearFEMLoad);
		}
		else if(iIndex == 6)
		{
			return GenerateLoadByType(AxialTimeLinearTotalStrainFEMLoad);
		}
		else if(iIndex == 7)
		{
			return GenerateLoadByType(AxialTimeSinusoidalTotalStrainFEMLoad);
		}
		else if(iIndex == 8)
		{
			return GenerateLoadByType(StepLoad);
		}
		else if(iIndex == 9)
		{
			return GenerateLoadByType(PulseTrainLoad);
		}
		return NULL;
	}
	void FEMLoad::SetID(const unsigned int& iID)
	{
		m_iID = iID;
	}
	unsigned int FEMLoad::GetID()
	{
		return m_iID;
	}
}



