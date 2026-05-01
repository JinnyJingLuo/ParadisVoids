#ifndef CHANGECONNECTIONOPERATION_H_
#define CHANGECONNECTIONOPERATION_H_

#include "TopologicalOperation.h"

namespace TopologySystem
{
	class ChangeConnectionOperation : public TopologicalOperation
	{
	public:
		ChangeConnectionOperation();
		ChangeConnectionOperation(const ChangeConnectionOperation& oOperation);
		~ChangeConnectionOperation();
		ChangeConnectionOperation& operator=(const ChangeConnectionOperation& oOperation);
		void Reset();
		TopologicalOperationType GetType() const;
		void Pack(list<double>* pldData) const;
		void Unpack(list<double>* pldData);
		void SetSourceData(const unsigned int& iDomainID,const unsigned int& iNodeID);
		void SetOldTargetData(const unsigned int& iDomainID,const unsigned int& iNodeID);
		void SetNewTargetData(const unsigned int& iDomainID,const unsigned int& iNodeID);
		unsigned int GetSourceDomainID() const;
		unsigned int GetSourceNodeID() const;
		unsigned int GetOldTargetDomainID() const;
		unsigned int GetOldTargetNodeID() const;
		unsigned int GetNewTargetDomainID() const;
		unsigned int GetNewTargetNodeID() const;
		
	private:
	
	protected:
		void Initialize();
		unsigned int m_iSourceDomainID;
		unsigned int m_iSourceNodeID;
		unsigned int m_iOldTargetDomainID;
		unsigned int m_iOldTargetNodeID;
		unsigned int m_iNewTargetDomainID;
		unsigned int m_iNewTargetNodeID;
	};
}


#endif


