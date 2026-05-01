// Ahmed M. Hussein

#include "FEMXTorsionLoad.h"
#include "Tools.h"
#include "cmath"

using namespace SupportSystem;

namespace FEMSystem
{
	FEMXTorsionLoad::FEMXTorsionLoad()
	{
		Initialize();	
	}
	FEMXTorsionLoad::~FEMXTorsionLoad()
	{
		Reset();
	}
	FEMXTorsionLoad::FEMXTorsionLoad(const FEMXTorsionLoad& oLoad)
	{
		*this = oLoad;
	}
	FEMXTorsionLoad& FEMXTorsionLoad::operator=(const FEMXTorsionLoad& oLoad)
	{
		FEMLoad::operator=(oLoad);
		m_dLoad = oLoad.m_dLoad;
		return *this;
	}
	void FEMXTorsionLoad::Reset()
	{
		FEMLoad::Reset();
		Initialize();
	}
	void FEMXTorsionLoad::Initialize()
	{
		FEMLoad::Initialize();
		m_dLoad = 0.0;
	}
	void FEMXTorsionLoad::Set(const double& dValue)
	{
		m_dLoad = dValue;
	}
	double FEMXTorsionLoad::Get(const Point& oPoint,const double& dTime)
	{
		double dX = oPoint.GetX();
		double dY = oPoint.GetY();
		double dR = sqrt(dX*dX + dY*dY);
		double dValue = -m_dLoad*dY/dR;
		return dValue;
	}
	FEMLoadTypes FEMXTorsionLoad::GetType() const
	{
		return XTorsionFEMLoad;
	}
	void FEMXTorsionLoad::Read(FILE* fpFile)
	{
		string sRead = GetRealString(500,fpFile);
		sscanf(sRead.c_str(),"%lf\n",&m_dLoad);
	}
	void FEMXTorsionLoad::Write(FILE* fpFile) const
	{
		fprintf(fpFile,"%d\n",GetType());
		fprintf(fpFile,"%e\n",m_dLoad);
	}
}



