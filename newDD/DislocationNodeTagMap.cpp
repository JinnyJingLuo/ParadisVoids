#include "DislocationNodeTagMap.h"



namespace DislocationSystem
{
	DislocationNodeTagMap::DislocationNodeTagMap()
	{
		Initialize();
	}
	DislocationNodeTagMap::DislocationNodeTagMap(const DislocationNodeTagMap& oTag)
	{
		*this = oTag;
	}
	DislocationNodeTagMap::~DislocationNodeTagMap()
	{
		Reset();
	}
	DislocationNodeTagMap& DislocationNodeTagMap::operator=(const DislocationNodeTagMap& oTag)
	{
		m_iSourceDomainID = oTag.m_iSourceDomainID;
		m_iSourceNodeID = oTag.m_iSourceNodeID;
		m_iTargetDomainID = oTag.m_iTargetDomainID;
		m_iTargetNodeID = oTag.m_iTargetNodeID;
		return *this;
	}
	void DislocationNodeTagMap::Reset()
	{
	
	}
	void DislocationNodeTagMap::SetSourceDomainID(const unsigned int& iID)
	{
		m_iSourceDomainID = iID;
	}
	void DislocationNodeTagMap::SetSourceNodeID(const unsigned int& iID)
	{
		m_iSourceNodeID = iID;
	}
	void DislocationNodeTagMap::SetTargetDomainID(const unsigned int& iID)
	{
		m_iTargetDomainID = iID;
	}
	void DislocationNodeTagMap::SetTargetNodeID(const unsigned int& iID)
	{
		m_iTargetNodeID = iID;
	}
	unsigned int DislocationNodeTagMap::GetSourceDomainID() const
	{
		return m_iSourceDomainID;
	}
	unsigned int DislocationNodeTagMap::GetSourceNodeID() const
	{
		return m_iSourceNodeID;
	}
	unsigned int DislocationNodeTagMap::GetTargetDomainID() const
	{
		return m_iTargetDomainID;
	}
	unsigned int DislocationNodeTagMap::GetTargetNodeID() const
	{
		return m_iTargetNodeID;
	}
	void DislocationNodeTagMap::Initialize()
	{
		m_iSourceDomainID = 0;
		m_iSourceNodeID = 0;
		m_iTargetDomainID = 0;
		m_iTargetNodeID = 0;
	}
	void DislocationNodeTagMap::Pack(list<unsigned int>* pliData) const
	{
		pliData->push_back(m_iSourceDomainID);
		pliData->push_back(m_iSourceNodeID);
		pliData->push_back(m_iTargetDomainID);
		pliData->push_back(m_iTargetNodeID);
	}
	void DislocationNodeTagMap::Unpack(list<unsigned int>* pliData)
	{
		m_iSourceDomainID = pliData->front();
		pliData->pop_front();
		m_iSourceNodeID = pliData->front();
		pliData->pop_front();
		m_iTargetDomainID = pliData->front();
		pliData->pop_front();
		m_iTargetNodeID = pliData->front();
		pliData->pop_front();
	}
}


