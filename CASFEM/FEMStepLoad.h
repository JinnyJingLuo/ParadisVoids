#ifndef FEMSTEPLOAD_H_
#define FEMSTEPLOAD_H_

#include "FEMLoad.h"

namespace FEMSystem
{
	class FEMStepLoad : public FEMLoad
	{
	public:
		FEMStepLoad();
		~FEMStepLoad();
		FEMStepLoad(const FEMStepLoad& oLoad);
		FEMStepLoad& operator=(const FEMStepLoad& oLoad);
		void Reset();
		void Set(const double& dStartTime,const double& dAmplitude);
		double Get(const Point& oPoint,const double& dTime);
		FEMLoadTypes GetType() const;
		virtual void Read(FILE* fpFile);
		virtual void Write(FILE* fpFile) const;
		
	private:
	
	protected:
		void Initialize();
		double m_dStartTime;
		double m_dAmplitude;
	};
}


#endif


