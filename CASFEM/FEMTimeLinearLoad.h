// Ahmed M. Hussein

#ifndef FEMTIMELINEARLOAD_H_
#define FEMTIMELINEARLOAD_H_

#include "FEMLoad.h"

namespace FEMSystem
{
	class FEMTimeLinearLoad : public FEMLoad
	{
	public:
		FEMTimeLinearLoad();
		~FEMTimeLinearLoad();
		FEMTimeLinearLoad(const FEMTimeLinearLoad& oLoad);
		FEMTimeLinearLoad& operator=(const FEMTimeLinearLoad& oLoad);
		void Reset();
		void Set(const double& dInitialLoad,const double& dRate);
		double Get(const Point& oPoint,const double& dTime);
		FEMLoadTypes GetType() const;
		virtual void Read(FILE* fpFile);
		virtual void Write(FILE* fpFile) const;
		
	private:
	
	protected:
		void Initialize();
		double m_dInitialLoad;
		double m_dRate;
	};
}

#endif



