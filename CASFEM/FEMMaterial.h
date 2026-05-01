// Ahmed M. Hussein

#ifndef FEMMATERIAL_H_
#define FEMMATERIAL_H_

#include "Matrix.h"
#include "FEMPlasticityState.h"
#include "FEMGaussPoint.h"

using namespace EZ;

namespace FEMSystem
{
	enum FEMMaterialType
	{
		NullFEMMaterial = 0,
		LinearElasticIsotropicFEMMaterial = 1,
		J2PlasticIsotropicFEMMaterial = 2
	};
	
	class FEMMaterial
	{
	public:
		virtual ~FEMMaterial();
		virtual FEMMaterial& operator=(const FEMMaterial& oMaterial);
		virtual void Reset();
		virtual Matrix GetStiffnessMatrix() const = 0;
		virtual Matrix GetConductionMatrix() const = 0;
		virtual double GetMassDensity() const = 0;
		virtual Matrix GetStressTemperatureCoupling() const = 0;
		virtual double GetSpecificHeatCapacity() const = 0;
		virtual double GetReferenceTemperature() const = 0;
		virtual Matrix UpdateGaussPoint(FEMGaussPoint* poPoint) const = 0;
		static FEMMaterial* GenerateMaterialByType(const FEMMaterialType& eType);
		static FEMMaterial* GenerateMaterialByTypeIndex(const unsigned int& iIndex);
		virtual void Read(FILE* fpFile) = 0;
		virtual void Write(FILE* fpFile) const = 0;
		virtual FEMMaterialType GetType() const = 0;
		void SetID(const unsigned int& iID);
		unsigned int GetID();
		virtual bool IsPlastic() const = 0;
		
	private:
	
	protected:
		virtual void Initialize();
		unsigned int m_iID;
	};
}

#endif


