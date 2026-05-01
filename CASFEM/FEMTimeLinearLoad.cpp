// Ahmed M. Hussein

#include "FEMTimeLinearLoad.h"
#include "Tools.h"

using namespace SupportSystem;

namespace FEMSystem
{
	FEMTimeLinearLoad::FEMTimeLinearLoad()
	{
		Initialize();
	}
	FEMTimeLinearLoad::~FEMTimeLinearLoad()
	{
		Reset();
	}
	FEMTimeLinearLoad::FEMTimeLinearLoad(const FEMTimeLinearLoad& oLoad)
	{
		*this = oLoad;
	}
	FEMTimeLinearLoad& FEMTimeLinearLoad::operator=(const FEMTimeLinearLoad& oLoad)
	{
		FEMLoad::operator=(oLoad);
		m_dInitialLoad = oLoad.m_dInitialLoad;
		m_dRate = oLoad.m_dRate;
		return *this;
	}
	void FEMTimeLinearLoad::Reset()
	{
		FEMLoad::Reset();
		Initialize();
	}
	void FEMTimeLinearLoad::Initialize()
	{
		FEMLoad::Initialize();
		m_dInitialLoad = 0.0;
		m_dRate = 0.0;
	}
	void FEMTimeLinearLoad::Set(const double& dInitialLoad,const double& dRate)
	{
		m_dInitialLoad = dInitialLoad;
		m_dRate = dRate;
	}
	double FEMTimeLinearLoad::Get(const Point& oPoint,const double& dTime)
	{
		return (m_dInitialLoad + m_dRate*dTime);
	}
	FEMLoadTypes FEMTimeLinearLoad::GetType() const
	{
		return TimeLinearFEMLoad;
	}
	void FEMTimeLinearLoad::Read(FILE* fpFile)
	{
		string sRead = GetRealString(500,fpFile);
		sscanf(sRead.c_str(),"%lf\t%lf\n",&m_dInitialLoad,&m_dRate);
	}
	void FEMTimeLinearLoad::Write(FILE* fpFile) const
	{
		fprintf(fpFile,"%d\n",GetType());
		fprintf(fpFile,"%e\t%e\n",m_dInitialLoad,m_dRate);
	}
}


