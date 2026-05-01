// Ahmed M. Hussein

#include "FEMMaterial.h"
#include "FEMLinearElasticIsotropicMaterial.h"
#include "FEMJ2PlasticIsotropicMaterial.h"

namespace FEMSystem
{
	FEMMaterial::~FEMMaterial()
	{
		Reset();
	}
	FEMMaterial& FEMMaterial::operator=(const FEMMaterial& oMaterial)
	{
		m_iID = oMaterial.m_iID;
		return *this;
	}
	void FEMMaterial::Reset()
	{
		Initialize();
	}
	void FEMMaterial::Initialize()
	{
		m_iID = 0;
	}
	FEMMaterial* FEMMaterial::GenerateMaterialByType(const FEMMaterialType& eType)
	{
		if(eType == NullFEMMaterial)
		{
			return NULL;
		}
		else if(eType == LinearElasticIsotropicFEMMaterial)
		{
			return new FEMLinearElasticIsotropicMaterial;
		}
		else if(eType == J2PlasticIsotropicFEMMaterial)
		{
			return new FEMJ2PlasticIsotropicMaterial;
		}
		return NULL;
	}
	FEMMaterial* FEMMaterial::GenerateMaterialByTypeIndex(const unsigned int& iIndex)
	{
		if(iIndex == 1)
		{
			return GenerateMaterialByType(LinearElasticIsotropicFEMMaterial);
		}
		if(iIndex == 2)
		{
			return GenerateMaterialByType(J2PlasticIsotropicFEMMaterial);
		}
		return NULL;
	}
	void FEMMaterial::SetID(const unsigned int& iID)
	{
		m_iID = iID;
	}
	unsigned int FEMMaterial::GetID()
	{
		return m_iID;
	}
}


