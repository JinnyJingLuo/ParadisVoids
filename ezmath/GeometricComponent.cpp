#include "GeometricComponent.h"


namespace GeometrySystem
{
	GeometricComponent::GeometricComponent()
	{
		Initialize();
	}
	GeometricComponent::GeometricComponent(const GeometricComponent& oComponent)
	{
		*this = oComponent;
	}
	GeometricComponent::~GeometricComponent()
	{
		Reset();
	}
	GeometricComponent& GeometricComponent::operator=(const GeometricComponent& oComponent)
	{
		m_iID = oComponent.m_iID;
		m_iCategory = oComponent.m_iCategory;
		return *this;
	}
	void GeometricComponent::Reset()
	{
		Initialize();
	}
	void GeometricComponent::SetID(const unsigned int& iID)
	{
		m_iID = iID;
	}
	unsigned int GeometricComponent::GetID() const
	{
		return m_iID;
	}
	void GeometricComponent::SetCategory(const int& iCategory)
	{
		m_iCategory = iCategory;
	}
	int GeometricComponent::GetCategory() const
	{
		return m_iCategory;
	}
	void GeometricComponent::Initialize()
	{
		m_iID = 0;
		m_iCategory = 0;
	}
}



