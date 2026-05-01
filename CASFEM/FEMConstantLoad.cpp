// Ahmed M. Hussein

#include "FEMConstantLoad.h"
#include "Tools.h"

using namespace SupportSystem;

namespace FEMSystem
{
	FEMConstantLoad::FEMConstantLoad()
	{
		Initialize();
	}
	FEMConstantLoad::~FEMConstantLoad()
	{
		Reset();
	}
	FEMConstantLoad::FEMConstantLoad(const FEMConstantLoad& oLoad)
	{
		*this = oLoad;
	}
	FEMConstantLoad& FEMConstantLoad::operator=(const FEMConstantLoad& oLoad)
	{
		FEMLoad::operator=(oLoad);
		m_dLoad = oLoad.m_dLoad;
		return *this;
	}
	void FEMConstantLoad::Reset()
	{
		FEMLoad::Reset();
		m_dLoad = 0.0;
	}
	void FEMConstantLoad::Initialize()
	{
		FEMLoad::Initialize();
		m_dLoad = 0.0;
	}
	void FEMConstantLoad::Set(const double& dValue)
	{
		m_dLoad = dValue;
	}
	double FEMConstantLoad::Get(const Point& oPoint,const double& dTime)
	{
		return m_dLoad;
	}
	FEMLoadTypes FEMConstantLoad::GetType() const
	{
		return ConstantFEMLoad;
	}
	void FEMConstantLoad::Read(FILE* fpFile)
	{
		string sRead = GetRealString(500,fpFile);
		sscanf(sRead.c_str(),"%lf\n",&m_dLoad);
	}
	void FEMConstantLoad::Write(FILE* fpFile) const
	{
		fprintf(fpFile,"%d\n",GetType());
		fprintf(fpFile,"%e\n",m_dLoad);
	}
}

