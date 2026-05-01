#include "FEMStepLoad.h"
#include "Tools.h"

using namespace SupportSystem;

namespace FEMSystem
{
	FEMStepLoad::FEMStepLoad()
	{
		Initialize();
	}
	FEMStepLoad::~FEMStepLoad()
	{
		Reset();
	}
	FEMStepLoad::FEMStepLoad(const FEMStepLoad& oLoad)
	{
		*this = oLoad;
	}
	FEMStepLoad& FEMStepLoad::operator=(const FEMStepLoad& oLoad)
	{
		FEMLoad::operator=(oLoad);
		m_dStartTime = oLoad.m_dStartTime;
		m_dAmplitude = oLoad.m_dAmplitude;
		return *this;
	}
	void FEMStepLoad::Reset()
	{
		FEMLoad::Reset();
		Initialize();
	}
	void FEMStepLoad::Initialize()
	{
		FEMLoad::Initialize();
		m_dStartTime = 0.0;
		m_dAmplitude = 0.0;
	}
	void FEMStepLoad::Set(const double& dStartTime,const double& dAmplitude)
	{
		m_dStartTime = dStartTime;
		m_dAmplitude = dAmplitude;
	}
	double FEMStepLoad::Get(const Point& oPoint,const double& dTime)
	{
		if(dTime >= m_dStartTime)
		{
			return m_dAmplitude;
		}
		return 0.0;
	}
	FEMLoadTypes FEMStepLoad::GetType() const
	{
		return StepLoad;
	}
	void FEMStepLoad::Read(FILE* fpFile)
	{
		string sRead = GetRealString(500,fpFile);
		sscanf(sRead.c_str(),"%lf\t%lf\n",&m_dStartTime,&m_dAmplitude);
	}
	void FEMStepLoad::Write(FILE* fpFile) const
	{
		fprintf(fpFile,"%d\n",GetType());
		fprintf(fpFile,"%e\t%e\n",m_dStartTime,m_dAmplitude);
	}
}


