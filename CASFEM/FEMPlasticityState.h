// Ahmed M. Hussein

#ifndef FEMPLASTICITYSTATE_H_
#define FEMPLASTICITYSTATE_H_

#include "Matrix.h"
#include "Vector.h"

using namespace EZ;

namespace FEMSystem
{
	class FEMPlasticityState
	{
	public:
		FEMPlasticityState();
		FEMPlasticityState(const FEMPlasticityState& oState);
		~FEMPlasticityState();
		FEMPlasticityState& operator=(const FEMPlasticityState& oState);
		FEMPlasticityState operator*(const double& dFactor) const;
		FEMPlasticityState operator+(const FEMPlasticityState& oState) const;
		void Reset();
		void SetPlasticStrain(const Matrix& oStrain);
		void SetBackStress(const Matrix& oBackStress);
		void SetAccumulatedStrain(const double& dStrain);
		const Matrix& GetPlasticStrain() const;
		const Matrix& GetBackStress() const;
		double GetAccumulatedStrain() const;
		Matrix GetVectorizedPlasticStrain() const;
		Matrix GetVectorizedBackStress() const;
		void SetVectorizedPlasticStrain(const Matrix& oStrain);
		void SetVectorizedBackStress(const Matrix& oBackStress);

	private:

	protected:
		void Initialize();
		Matrix m_oPlasticStrain;
		Matrix m_oBackStress;
		double m_dAccumulatedStrain;
	};
}

#endif


