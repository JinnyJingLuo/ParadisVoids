// Ahmed M. Hussein

#ifndef FEMTIMESINUSOIDALLOAD_H_
#define FEMTIMESINUSOIDALLOAD_H_

#include "FEMLoad.h"

namespace FEMSystem
{
	class FEMTimeSinusoidalLoad : public FEMLoad
	{
	public:
		FEMTimeSinusoidalLoad();
		~FEMTimeSinusoidalLoad();
		FEMTimeSinusoidalLoad(const FEMTimeSinusoidalLoad& oLoad);
		FEMTimeSinusoidalLoad& operator=(const FEMTimeSinusoidalLoad& oLoad);
		void Reset();
		void Set(const double& dAmplitude,const double& dFrequency,const double& dPhaseShift);
		double Get(const Point& oPoint,const double& dTime);
		FEMLoadTypes GetType() const;
		virtual void Read(FILE* fpFile);
		virtual void Write(FILE* fpFile) const;
		
	private:
	
	protected:
		void Initialize();
		double m_dAmplitude;
		double m_dFrequency;
		double m_dPhaseShift;
	};
}

#endif



