#ifndef DISLOCATIONNODETAGMAP_H_
#define DISLOCATIONNODETAGMAP_H_

#include "list"

using namespace std;

namespace DislocationSystem
{
	class DislocationNodeTagMap
	{
	public:
		DislocationNodeTagMap();
		DislocationNodeTagMap(const DislocationNodeTagMap& oTag);
		~DislocationNodeTagMap();
		DislocationNodeTagMap& operator=(const DislocationNodeTagMap& oTag);
		void Reset();
		void SetSourceDomainID(const unsigned int& iID);
		void SetSourceNodeID(const unsigned int& iID);
		void SetTargetDomainID(const unsigned int& iID);
		void SetTargetNodeID(const unsigned int& iID);
		unsigned int GetSourceDomainID() const;
		unsigned int GetSourceNodeID() const;
		unsigned int GetTargetDomainID() const;
		unsigned int GetTargetNodeID() const;
		void Pack(list<unsigned int>* pliData) const;
		void Unpack(list<unsigned int>* pliData);
		
	private:
	
	protected:
		void Initialize();
		unsigned int m_iSourceDomainID;
		unsigned int m_iSourceNodeID;
		unsigned int m_iTargetDomainID;
		unsigned int m_iTargetNodeID;
	};
}

#endif


