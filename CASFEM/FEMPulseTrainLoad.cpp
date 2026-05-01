#include "FEMPulseTrainLoad.h"
#include "Tools.h"
#include "cmath"

using namespace SupportSystem;

namespace FEMSystem
{
	FEMPulseTrainLoad::FEMPulseTrainLoad()
	{
		Initialize();
	}
	FEMPulseTrainLoad::~FEMPulseTrainLoad()
	{
		Reset();
	}
	FEMPulseTrainLoad::FEMPulseTrainLoad(const FEMPulseTrainLoad& oLoad)
	{
		*this = oLoad;
	}
	FEMPulseTrainLoad& FEMPulseTrainLoad::operator=(const FEMPulseTrainLoad& oLoad)
	{
		FEMLoad::operator=(oLoad);
		m_dPeriod = oLoad.m_dPeriod;
		m_dDutyCycle = oLoad.m_dDutyCycle;
		m_dAmplitude = oLoad.m_dAmplitude;
		return *this;
	}
	void FEMPulseTrainLoad::Reset()
	{
		FEMLoad::Reset();
		Initialize();
	}
	void FEMPulseTrainLoad::Initialize()
	{
		FEMLoad::Initialize();
		m_dPeriod = 0.0;
		m_dDutyCycle = 0.0;
		m_dAmplitude = 0.0;
	}
	void FEMPulseTrainLoad::Set(const double& dPeriod,const double& dDutyCycle,const double& dAmplitude)
	{
		m_dPeriod = dPeriod;
		m_dDutyCycle = dDutyCycle;
		m_dAmplitude = dAmplitude;
	}
	double FEMPulseTrainLoad::Get(const Point& oPoint,const double& dTime)
	{
		double dTimeRatio = dTime/m_dPeriod;
		double dScaledTime = dTimeRatio - floor(dTimeRatio);
		if(dScaledTime <= m_dDutyCycle)
		{
			return m_dAmplitude;
		}
		return 0.0;
	}
	FEMLoadTypes FEMPulseTrainLoad::GetType() const
	{
		return PulseTrainLoad;
	}
	void FEMPulseTrainLoad::Read(FILE* fpFile)
	{
		string sRead = GetRealString(500,fpFile);
		sscanf(sRead.c_str(),"%lf\t%lf\t%lf\n",&m_dPeriod,&m_dDutyCycle,&m_dAmplitude);
	}
	void FEMPulseTrainLoad::Write(FILE* fpFile) const
	{
		fprintf(fpFile,"%d\n",GetType());
		fprintf(fpFile,"%e\t%e\t%e\n",m_dPeriod,m_dDutyCycle,m_dAmplitude);
	}
}


