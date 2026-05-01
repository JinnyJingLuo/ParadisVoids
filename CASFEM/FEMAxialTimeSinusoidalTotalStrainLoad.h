#ifndef FEMAXIALTIMESINUSOIDALTOTALSTRAINLOAD_H_
#define FEMAXIALTIMESINUSOIDALTOTALSTRAINLOAD_H_

#include "FEMLoad.h"

namespace FEMSystem
{
	class FEMAxialTimeSinusoidalTotalStrainLoad : public FEMLoad
	{
	public:
		FEMAxialTimeSinusoidalTotalStrainLoad();
		~FEMAxialTimeSinusoidalTotalStrainLoad();
		FEMAxialTimeSinusoidalTotalStrainLoad(const FEMAxialTimeSinusoidalTotalStrainLoad& oLoad);
		FEMAxialTimeSinusoidalTotalStrainLoad& operator=(const FEMAxialTimeSinusoidalTotalStrainLoad& oLoad);
		void Reset();
		void Set(const double& dLength,const double& dAmplitude,const double& dFrequency);
		double Get(const Point& oPoint,const double& dTime);
		FEMLoadTypes GetType() const;
		virtual void Read(FILE* fpFile);
		virtual void Write(FILE* fpFile) const;
		void SetCurrentPlasticStrain(const double& dPlasticStrain);
		
	private:
	
	protected:
		void Initialize();
		double m_dLength;
		double m_dAmplitude;
		double m_dFrequency;
		double m_dCurrentPlasticStrain;
	};
}

#endif



