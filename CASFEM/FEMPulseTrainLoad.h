#ifndef FEMPULSETRAINLOAD_H_
#define FEMPULSETRAINLOAD_H_

#include "FEMLoad.h"

namespace FEMSystem
{
	class FEMPulseTrainLoad : public FEMLoad
	{
	public:
		FEMPulseTrainLoad();
		~FEMPulseTrainLoad();
		FEMPulseTrainLoad(const FEMPulseTrainLoad& oLoad);
		FEMPulseTrainLoad& operator=(const FEMPulseTrainLoad& oLoad);
		void Reset();
 		void Set(const double& dPeriod,const double& dDutyCycle,const double& dAmplitude);
		double Get(const Point& oPoint,const double& dTime);
		FEMLoadTypes GetType() const;
		virtual void Read(FILE* fpFile);
		virtual void Write(FILE* fpFile) const;
		
	private:
	
	protected:
		void Initialize();
		double m_dPeriod;
		double m_dDutyCycle;
		double m_dAmplitude;
	};
}

#endif


