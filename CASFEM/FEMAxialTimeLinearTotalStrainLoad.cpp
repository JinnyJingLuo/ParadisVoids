// Ahmed M. Hussein

#include "FEMAxialTimeLinearTotalStrainLoad.h"
#include "Tools.h"
#include "cmath"

using namespace SupportSystem;

namespace FEMSystem
{
	FEMAxialTimeLinearTotalStrainLoad::FEMAxialTimeLinearTotalStrainLoad()
	{
		Initialize();
	}
	FEMAxialTimeLinearTotalStrainLoad::~FEMAxialTimeLinearTotalStrainLoad()
	{
		Reset();
	}
	FEMAxialTimeLinearTotalStrainLoad::FEMAxialTimeLinearTotalStrainLoad(const FEMAxialTimeLinearTotalStrainLoad& oLoad)
	{
		*this = oLoad;
	}
	FEMAxialTimeLinearTotalStrainLoad& FEMAxialTimeLinearTotalStrainLoad::operator=(const FEMAxialTimeLinearTotalStrainLoad& oLoad)
	{
		FEMLoad::operator=(oLoad);
		m_dLength = oLoad.m_dLength;
		m_dTotalStrainRate = oLoad.m_dTotalStrainRate;
		m_dCurrentPlasticStrain = oLoad.m_dCurrentPlasticStrain;
		return *this;
	}
	void FEMAxialTimeLinearTotalStrainLoad::Reset()
	{
		FEMLoad::Reset();
		Initialize();
	}
	void FEMAxialTimeLinearTotalStrainLoad::Initialize()
	{
		FEMLoad::Initialize();
		m_dLength = 0.0;
		m_dTotalStrainRate = 0.0;
		m_dCurrentPlasticStrain = 0.0;
	}
	void FEMAxialTimeLinearTotalStrainLoad::Set(const double& dLength,const double& dTotalStrainRate)
	{
		m_dLength = dLength;
		m_dTotalStrainRate = dTotalStrainRate;
		m_dCurrentPlasticStrain = 0.0;
	}
	double FEMAxialTimeLinearTotalStrainLoad::Get(const Point& oPoint,const double& dTime)
	{
		// this returns the total axial displacement
		return m_dLength*(m_dTotalStrainRate*dTime - m_dCurrentPlasticStrain);
	}
	FEMLoadTypes FEMAxialTimeLinearTotalStrainLoad::GetType() const
	{
		return AxialTimeLinearTotalStrainFEMLoad;
	}
	void FEMAxialTimeLinearTotalStrainLoad::Read(FILE* fpFile)
	{
		string sRead = GetRealString(500,fpFile);
		sscanf(sRead.c_str(),"%lf\t%lf\n",&m_dLength,&m_dTotalStrainRate);
	}
	void FEMAxialTimeLinearTotalStrainLoad::Write(FILE* fpFile) const
	{
		fprintf(fpFile,"%d\n",GetType());
		fprintf(fpFile,"%e\t%e\n",m_dLength,m_dTotalStrainRate);
	}
	void FEMAxialTimeLinearTotalStrainLoad::SetCurrentPlasticStrain(const double& dPlasticStrain)
	{
		m_dCurrentPlasticStrain = dPlasticStrain;
	}
}




