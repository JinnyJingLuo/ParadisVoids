// Ahmed M. Hussein

#include "FEMTimeSinusoidalLoad.h"
#include "Tools.h"
#include "cmath"

using namespace SupportSystem;

namespace FEMSystem
{
	FEMTimeSinusoidalLoad::FEMTimeSinusoidalLoad()
	{
		Initialize();
	}
	FEMTimeSinusoidalLoad::~FEMTimeSinusoidalLoad()
	{
		Reset();
	}
	FEMTimeSinusoidalLoad::FEMTimeSinusoidalLoad(const FEMTimeSinusoidalLoad& oLoad)
	{
		*this = oLoad;
	}
	FEMTimeSinusoidalLoad& FEMTimeSinusoidalLoad::operator=(const FEMTimeSinusoidalLoad& oLoad)
	{
		FEMLoad::operator=(oLoad);
		m_dAmplitude = oLoad.m_dAmplitude;
		m_dFrequency = oLoad.m_dFrequency;
		m_dPhaseShift = oLoad.m_dPhaseShift;
		return *this;
	}
	void FEMTimeSinusoidalLoad::Reset()
	{
		FEMLoad::Reset();
		Initialize();
	}
	void FEMTimeSinusoidalLoad::Initialize()
	{
		FEMLoad::Initialize();
		m_dAmplitude = 0.0;
		m_dFrequency = 0.0;
		m_dPhaseShift = 0.0;
	}
	void FEMTimeSinusoidalLoad::Set(const double& dAmplitude,const double& dFrequency,const double& dPhaseShift)
	{
		m_dAmplitude = dAmplitude;
		m_dFrequency = dFrequency;
		m_dPhaseShift = dPhaseShift;
	}
	double FEMTimeSinusoidalLoad::Get(const Point& oPoint,const double& dTime)
	{
		return (m_dAmplitude*sin(m_dFrequency*dTime + m_dPhaseShift));
	}
	FEMLoadTypes FEMTimeSinusoidalLoad::GetType() const
	{
		return TimeSinusoidalFEMLoad;
	}
	void FEMTimeSinusoidalLoad::Read(FILE* fpFile)
	{
		string sRead = GetRealString(500,fpFile);
		sscanf(sRead.c_str(),"%lf\t%lf\t%lf\n",&m_dAmplitude,&m_dFrequency,&m_dPhaseShift);
	}
	void FEMTimeSinusoidalLoad::Write(FILE* fpFile) const
	{
		fprintf(fpFile,"%d\n",GetType());
		fprintf(fpFile,"%e\t%e\t%e\n",m_dAmplitude,m_dFrequency,m_dPhaseShift);
	}
}




