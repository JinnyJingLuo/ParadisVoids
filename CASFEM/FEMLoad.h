// Ahmed M. Hussein

#ifndef FEMLOAD_H_
#define FEMLOAD_H_

#include "iostream"
#include "Point.h"

using namespace std;
using namespace EZ;

namespace FEMSystem
{
	enum FEMLoadTypes
	{
		NullFEMLoad = 0,
		ConstantFEMLoad = 1,
		TimeSinusoidalFEMLoad = 2,
		XTorsionFEMLoad = 3,
		YTorsionFEMLoad = 4,
		TimeLinearFEMLoad = 5,
		AxialTimeLinearTotalStrainFEMLoad = 6,
		AxialTimeSinusoidalTotalStrainFEMLoad = 7,
		StepLoad = 8,
		PulseTrainLoad = 9
	};
	
	class FEMLoad
	{
	public:
		virtual ~FEMLoad();
		virtual FEMLoad& operator=(const FEMLoad& oLoad);
		virtual void Reset();
		virtual double Get(const Point& oPoint,const double& dTime) = 0;
		virtual void Read(FILE* fpFile) = 0;
		virtual void Write(FILE* fpFile) const = 0;
		virtual FEMLoadTypes GetType() const = 0;
		static FEMLoad* GenerateLoadByType(const FEMLoadTypes& eType);
		static FEMLoad* GenerateLoadByTypeIndex(const unsigned int& iIndex);
		void SetID(const unsigned int& iID);
		unsigned int GetID();

	private:
	
	protected:
		virtual void Initialize();
		unsigned int m_iID;
	};
}


#endif


