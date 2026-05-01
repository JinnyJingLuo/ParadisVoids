#include "FEMAxialTimeSinusoidalTotalStrainLoad.h"
#include "Tools.h"
#include "cmath"

using namespace SupportSystem;

namespace FEMSystem
{
	FEMAxialTimeSinusoidalTotalStrainLoad::FEMAxialTimeSinusoidalTotalStrainLoad()
	{
		Initialize();
	}
	FEMAxialTimeSinusoidalTotalStrainLoad::~FEMAxialTimeSinusoidalTotalStrainLoad()
	{
		Reset();
	}
	FEMAxialTimeSinusoidalTotalStrainLoad::FEMAxialTimeSinusoidalTotalStrainLoad(const FEMAxialTimeSinusoidalTotalStrainLoad& oLoad)
	{
		*this = oLoad;
	}
	FEMAxialTimeSinusoidalTotalStrainLoad& FEMAxialTimeSinusoidalTotalStrainLoad::operator=(const FEMAxialTimeSinusoidalTotalStrainLoad& oLoad)
	{
		FEMLoad::operator=(oLoad);
		m_dLength = oLoad.m_dLength;
		m_dAmplitude = oLoad.m_dAmplitude;
		m_dFrequency = oLoad.m_dFrequency;
		m_dCurrentPlasticStrain = oLoad.m_dCurrentPlasticStrain;
		return *this;
	}
	void FEMAxialTimeSinusoidalTotalStrainLoad::Reset()
	{
		FEMLoad::Reset();
		Initialize();
	}
	void FEMAxialTimeSinusoidalTotalStrainLoad::Initialize()
	{
		FEMLoad::Initialize();
		m_dLength = 0.0;
		m_dAmplitude = 0.0;
		m_dFrequency = 0.0;
		m_dCurrentPlasticStrain = 0.0;
	}
	void FEMAxialTimeSinusoidalTotalStrainLoad::Set(const double& dLength,const double& dAmplitude,const double& dFrequency)
	{
		m_dLength = dLength;
		m_dAmplitude = dAmplitude;
		m_dFrequency = dFrequency;
		m_dCurrentPlasticStrain = 0.0;
	}
	double FEMAxialTimeSinusoidalTotalStrainLoad::Get(const Point& oPoint,const double& dTime)
	{
		// this returns the total axial displacement     
		return m_dLength*(m_dAmplitude*sin(m_dFrequency*dTime) - m_dCurrentPlasticStrain);
	}
	FEMLoadTypes FEMAxialTimeSinusoidalTotalStrainLoad::GetType() const
	{
		return AxialTimeSinusoidalTotalStrainFEMLoad;
	}
	void FEMAxialTimeSinusoidalTotalStrainLoad::Read(FILE* fpFile)
	{
		string sRead = GetRealString(500,fpFile);
		sscanf(sRead.c_str(),"%lf\t%lf\t%lf\n",&m_dLength,&m_dAmplitude,&m_dFrequency);
	}
	void FEMAxialTimeSinusoidalTotalStrainLoad::Write(FILE* fpFile) const
	{
		fprintf(fpFile,"%d\n",GetType());
		fprintf(fpFile,"%e\t%e\t%e\n",m_dLength,m_dAmplitude,m_dFrequency);
	}
	void FEMAxialTimeSinusoidalTotalStrainLoad::SetCurrentPlasticStrain(const double& dPlasticStrain)
	{
		m_dCurrentPlasticStrain = dPlasticStrain;
	}
}


