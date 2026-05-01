// Ahmed M. Hussein

#ifndef FEMXTORSIONLOAD_H_
#define FEMXTORSIONLOAD_H_

#include "FEMLoad.h"


namespace FEMSystem
{
	class FEMXTorsionLoad : public FEMLoad
	{
	public:
		FEMXTorsionLoad();
		~FEMXTorsionLoad();
		FEMXTorsionLoad(const FEMXTorsionLoad& oLoad);
		FEMXTorsionLoad& operator=(const FEMXTorsionLoad& oLoad);
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


