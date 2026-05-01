#include "TopologicalOperation.h"
#include "ChangeConnectionOperation.h"

namespace TopologySystem
{
	TopologicalOperation::~TopologicalOperation()
	{
		Reset();
	}
	void TopologicalOperation::Reset()
	{
	
	}
	void TopologicalOperation::SetResponsibleDomainID(const unsigned int& iID)
	{
		m_iResponsibleDomain = iID;
	}
	unsigned int TopologicalOperation::GetResponsibleDomainID() const
	{
		return m_iResponsibleDomain;
	}
	TopologicalOperation* TopologicalOperation::CreateAndUnpackOperation(list<double>* pldData)
	{
		int iOperationType = (TopologicalOperationType)(int)pldData->front();
		pldData->pop_front();
		TopologicalOperation* poOperation = NULL;
		if(iOperationType == ChangeConnection)
		{
			poOperation = new ChangeConnectionOperation;
			poOperation->Unpack(pldData);
		}
		return poOperation;
	}
	void TopologicalOperation::Initialize()
	{
		m_iResponsibleDomain = 0;
	}
}


