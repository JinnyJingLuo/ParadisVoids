#include "ChangeConnectionOperation.h"

namespace TopologySystem
{
	ChangeConnectionOperation::ChangeConnectionOperation()
	{
		Initialize();
	}
	ChangeConnectionOperation::ChangeConnectionOperation(const ChangeConnectionOperation& oOperation)
	{
		*this = oOperation;
	}
	ChangeConnectionOperation::~ChangeConnectionOperation()
	{
		Reset();
	}
	ChangeConnectionOperation& ChangeConnectionOperation::operator=(const ChangeConnectionOperation& oOperation)
	{
		m_iResponsibleDomain = oOperation.m_iResponsibleDomain;
		m_iSourceDomainID = oOperation.m_iSourceDomainID;
		m_iSourceNodeID = oOperation.m_iSourceNodeID;
		m_iOldTargetDomainID = oOperation.m_iOldTargetDomainID;
		m_iOldTargetNodeID = oOperation.m_iOldTargetNodeID;
		m_iNewTargetDomainID = oOperation.m_iNewTargetDomainID;
		m_iNewTargetNodeID = oOperation.m_iNewTargetNodeID;
		return *this;
	}
	void ChangeConnectionOperation::Reset()
	{
	
	}
	TopologicalOperationType ChangeConnectionOperation::GetType() const
	{
		return ChangeConnection;
	}
	void ChangeConnectionOperation::Pack(list<double>* pldData) const
	{
		pldData->push_back((double)ChangeConnection);
		pldData->push_back((double)m_iSourceDomainID);
		pldData->push_back((double)m_iSourceNodeID);
		pldData->push_back((double)m_iOldTargetDomainID);
		pldData->push_back((double)m_iOldTargetNodeID);
		pldData->push_back((double)m_iNewTargetDomainID);
		pldData->push_back((double)m_iNewTargetNodeID);
	}
	void ChangeConnectionOperation::Unpack(list<double>* pldData)
	{
		m_iSourceDomainID = (unsigned int)pldData->front();
		pldData->pop_front();
		m_iSourceNodeID = (unsigned int)pldData->front();
		pldData->pop_front();
		m_iOldTargetDomainID = (unsigned int)pldData->front();
		pldData->pop_front();
		m_iOldTargetNodeID = (unsigned int)pldData->front();
		pldData->pop_front();
		m_iNewTargetDomainID = (unsigned int)pldData->front();
		pldData->pop_front();
		m_iNewTargetNodeID = (unsigned int)pldData->front();
		pldData->pop_front();
	}
	void ChangeConnectionOperation::SetSourceData(const unsigned int& iDomainID,const unsigned int& iNodeID)
	{
		m_iSourceDomainID = iDomainID;
		m_iSourceNodeID = iNodeID;
	}
	void ChangeConnectionOperation::SetOldTargetData(const unsigned int& iDomainID,const unsigned int& iNodeID)
	{
		m_iOldTargetDomainID = iDomainID;
		m_iOldTargetNodeID = iNodeID;
	}
	void ChangeConnectionOperation::SetNewTargetData(const unsigned int& iDomainID,const unsigned int& iNodeID)
	{
		m_iNewTargetDomainID = iDomainID;
		m_iNewTargetNodeID = iNodeID;
	}
	unsigned int ChangeConnectionOperation::GetSourceDomainID() const
	{
		return m_iSourceDomainID;
	}
	unsigned int ChangeConnectionOperation::GetSourceNodeID() const
	{
		return m_iSourceNodeID;
	}
	unsigned int ChangeConnectionOperation::GetOldTargetDomainID() const
	{
		return m_iOldTargetDomainID;
	}
	unsigned int ChangeConnectionOperation::GetOldTargetNodeID() const
	{
		return m_iOldTargetNodeID;
	}
	unsigned int ChangeConnectionOperation::GetNewTargetDomainID() const
	{
		return m_iNewTargetDomainID;
	}
	unsigned int ChangeConnectionOperation::GetNewTargetNodeID() const
	{
		return m_iNewTargetNodeID;
	}
	void ChangeConnectionOperation::Initialize()
	{
		m_iSourceDomainID = 0;
		m_iSourceNodeID = 0;
		m_iOldTargetDomainID = 0;
		m_iOldTargetNodeID = 0;
		m_iNewTargetDomainID = 0;
		m_iNewTargetNodeID = 0;
	}
}


