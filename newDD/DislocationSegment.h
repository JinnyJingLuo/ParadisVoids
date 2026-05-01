// Ahmed M. Hussein

#ifndef DISLOCATIONSEGMENT_H_
#define DISLOCATIONSEGMENT_H_

#include "Vector.h"
#include "DislocationNode.h"

using namespace EZ;

namespace DislocationSystem
{	
	class DislocationSegment
	{
	public:
		DislocationSegment();
		DislocationSegment(const DislocationSegment& oSegment);
		~DislocationSegment();
		DislocationSegment& operator=(const DislocationSegment& oSegment);
		void Reset();
		void SetEndDomainID(const unsigned int iID);
		void SetEndNodeID(const unsigned int iID);
		unsigned int GetEndDomainID() const;
		unsigned int GetEndNodeID() const;
		DislocationNodeTag* GetEndTag();
		void SetSlipPlaneNormal(const Vector& oNormal);
		void SetBurgersVector(const Vector& oBurgersVector);
		Vector GetSlipPlaneNormal() const;
		Vector GetBurgersVector() const;
		void Write(FILE* fpFile) const;
		void Pack(list<double>* pldData) const;
		void Unpack(list<double>* pldData);
		
	private:
	
	protected:
		void Initialize();
		DislocationNodeTag m_oEndTag;
		Vector m_oSlipPlaneNormal;
		Vector m_oBurgersVector;
	};
}

#endif


