#include "DislocationNodeTag.h"



namespace DislocationSystem
{
	DislocationNodeTag::DislocationNodeTag()
	{
		Initialize();
	}
	DislocationNodeTag::DislocationNodeTag(const DislocationNodeTag& oTag)
	{
		*this = oTag;
	}
	DislocationNodeTag::~DislocationNodeTag()
	{
		
	}
	DislocationNodeTag& DislocationNodeTag::operator=(const DislocationNodeTag& oTag)
	{
		m_iDomainID = oTag.m_iDomainID;
		m_iNodeID = oTag.m_iNodeID;
		return *this;
	}
	void DislocationNodeTag::Reset()
	{
		m_iDomainID = 0;
		m_iNodeID = 0;
	}
	void DislocationNodeTag::SetDomainID(const unsigned int& iID)
	{
		m_iDomainID = iID;
	}
	void DislocationNodeTag::SetNodeID(const unsigned int& iID)
	{
		m_iNodeID = iID;
	}
	unsigned int DislocationNodeTag::GetDomainID() const
	{
		return m_iDomainID;
	}
	unsigned int DislocationNodeTag::GetNodeID() const
	{
		return m_iNodeID;
	}
	void DislocationNodeTag::Initialize()
	{
		m_iDomainID = 0;
		m_iNodeID = 0;
	}
}



