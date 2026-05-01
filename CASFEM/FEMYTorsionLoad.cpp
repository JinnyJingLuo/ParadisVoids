// Ahmed M. Hussein

#include "FEMYTorsionLoad.h"
#include "Tools.h"
#include "cmath"

using namespace SupportSystem;

namespace FEMSystem
{
	FEMYTorsionLoad::FEMYTorsionLoad()
	{
		Initialize();
	}
	FEMYTorsionLoad::~FEMYTorsionLoad()
	{
		Reset();
	}
	FEMYTorsionLoad::FEMYTorsionLoad(const FEMYTorsionLoad& oLoad)
	{
		*this = oLoad;
	}
	FEMYTorsionLoad& FEMYTorsionLoad::operator=(const FEMYTorsionLoad& oLoad)
	{
		FEMLoad::operator=(oLoad);
		m_dLoad = oLoad.m_dLoad;
		return *this;
	}
	void FEMYTorsionLoad::Reset()
	{
		FEMLoad::Reset();
		Initialize();
	}
	void FEMYTorsionLoad::Initialize()
	{
		FEMLoad::Initialize();
		m_dLoad = 0.0;
	}
	void FEMYTorsionLoad::Set(const double& dValue)
	{
		m_dLoad = dValue;
	}
	double FEMYTorsionLoad::Get(const Point& oPoint,const double& dTime)
	{
		double dX = oPoint.GetX();
		double dY = oPoint.GetY();
		double dR = sqrt(dX*dX + dY*dY);
		double dValue = m_dLoad*dX/dR;
		return dValue;
	}
	FEMLoadTypes FEMYTorsionLoad::GetType() const
	{
		return YTorsionFEMLoad;
	}
	void FEMYTorsionLoad::Read(FILE* fpFile)
	{
		string sRead = GetRealString(500,fpFile);
		sscanf(sRead.c_str(),"%lf\n",&m_dLoad);
	}
	void FEMYTorsionLoad::Write(FILE* fpFile) const
	{
		fprintf(fpFile,"%d\n",GetType());
		fprintf(fpFile,"%e\n",m_dLoad);
	}
}



