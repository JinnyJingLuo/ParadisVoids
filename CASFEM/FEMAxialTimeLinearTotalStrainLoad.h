#ifndef FEMAXIALTIMELINEARTOTALSTRAINLOAD_H_
#define FEMAXIALTIMELINEARTOTALSTRAINLOAD_H_

#include "FEMLoad.h"

namespace FEMSystem
{
	class FEMAxialTimeLinearTotalStrainLoad : public FEMLoad
	{
	public:
		FEMAxialTimeLinearTotalStrainLoad();
		~FEMAxialTimeLinearTotalStrainLoad();
		FEMAxialTimeLinearTotalStrainLoad(const FEMAxialTimeLinearTotalStrainLoad& oLoad);
		FEMAxialTimeLinearTotalStrainLoad& operator=(const FEMAxialTimeLinearTotalStrainLoad& oLoad);
		void Reset();
		void Set(const double& dLength,const double& dTotalStrainRate);
		double Get(const Point& oPoint,const double& dTime);
		FEMLoadTypes GetType() const;
		virtual void Read(FILE* fpFile);
		virtual void Write(FILE* fpFile) const;
		void SetCurrentPlasticStrain(const double& dPlasticStrain);
		
	private:
	
	protected:
		void Initialize();
		double m_dLength;
		double m_dTotalStrainRate;
		double m_dCurrentPlasticStrain;
	};
}


#endif


