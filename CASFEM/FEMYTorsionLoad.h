// Ahmed M. Hussein

#ifndef FEMYTORSIONLOAD_H_
#define FEMYTORSIONLOAD_H_

#include "FEMLoad.h"


namespace FEMSystem
{
	class FEMYTorsionLoad : public FEMLoad
	{
	public:
		FEMYTorsionLoad();
		~FEMYTorsionLoad();
		FEMYTorsionLoad(const FEMYTorsionLoad& oLoad);
		FEMYTorsionLoad& operator=(const FEMYTorsionLoad& oLoad);
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


