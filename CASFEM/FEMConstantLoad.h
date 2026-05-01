// Ahmed M. Hussein

#ifndef FEMCONSTANTLOAD_H_
#define FEMCONSTANTLOAD_H_

#include "FEMLoad.h"


namespace FEMSystem
{
	class FEMConstantLoad : public FEMLoad
	{
	public:
		FEMConstantLoad();
		~FEMConstantLoad();
		FEMConstantLoad(const FEMConstantLoad& oLoad);
		FEMConstantLoad& operator=(const FEMConstantLoad& oLoad);
		void Reset();
		void Set(const double& dValue);
		FEMLoadTypes GetType() const;
		double Get(const Point& oPoint,const double& dTime);
		virtual void Read(FILE* fpFile);
		virtual void Write(FILE* fpFile) const;
	private:
	
	protected:
		void Initialize();
		double m_dLoad;
	};
}

#endif


