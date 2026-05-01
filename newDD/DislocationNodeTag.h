#ifndef DISLOCATIONNOETAG_H_
#define DISLOCATIONNODETAG_H_

namespace DislocationSystem
{
	class DislocationNodeTag
	{
	public:
		DislocationNodeTag();
		DislocationNodeTag(const DislocationNodeTag& oTag);
		~DislocationNodeTag();
		DislocationNodeTag& operator=(const DislocationNodeTag& oTag);
		void Reset();
		void SetDomainID(const unsigned int& iID);
		void SetNodeID(const unsigned int& iID);
		unsigned int GetDomainID() const;
		unsigned int GetNodeID() const;
		
	private:
	
	protected:
		void Initialize();
		unsigned int m_iDomainID;
		unsigned int m_iNodeID;
	};
}

#endif

