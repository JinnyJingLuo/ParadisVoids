#ifndef TOPOLOGICALOPERATION_H_
#define TOPOLOGICALOPERATION_H_

#include "list"

using namespace std;

namespace TopologySystem
{
	enum TopologicalOperationType
	{
		NullOperation = 0,
		ChangeConnection = 1
	};
	
	class TopologicalOperation
	{
	public:
		virtual ~TopologicalOperation();
		virtual void Reset();
		virtual TopologicalOperationType GetType() const = 0;
		virtual void Pack(list<double>* pldData) const = 0;
		virtual void Unpack(list<double>* pldData) = 0;
		void SetResponsibleDomainID(const unsigned int& iID);
		unsigned int GetResponsibleDomainID() const;
		static TopologicalOperation* CreateAndUnpackOperation(list<double>* pldData);
		
	private:
	
	protected:
		virtual void Initialize();
		unsigned int m_iResponsibleDomain;
	};
}

#endif


