// Ahmed M. Hussein

#include "DislocationSegment.h"

namespace DislocationSystem
{
	DislocationSegment::DislocationSegment()
	{
		Initialize();
	}
	DislocationSegment::DislocationSegment(const DislocationSegment& oSegment)
	{
		*this = oSegment;
	}
	DislocationSegment::~DislocationSegment()
	{
		
	}
	DislocationSegment& DislocationSegment::operator=(const DislocationSegment& oSegment)
	{
		Reset();
		m_oEndTag = oSegment.m_oEndTag;
		m_oSlipPlaneNormal = oSegment.m_oSlipPlaneNormal;
		m_oBurgersVector = oSegment.m_oBurgersVector;
		return *this;
	}
	void DislocationSegment::Reset()
	{
		m_oEndTag.Reset();
		m_oSlipPlaneNormal.Set(0.0,0.0,0.0);
		m_oBurgersVector.Set(0.0,0.0,0.0);
	}
	void DislocationSegment::SetEndDomainID(const unsigned int iID)
	{
		m_oEndTag.SetDomainID(iID);
	}
	void DislocationSegment::SetEndNodeID(const unsigned int iID)
	{
		m_oEndTag.SetNodeID(iID);
	}
	unsigned int DislocationSegment::GetEndDomainID() const
	{
		return m_oEndTag.GetDomainID();
	}
	unsigned int DislocationSegment::GetEndNodeID() const
	{
		return m_oEndTag.GetNodeID();
	}
	DislocationNodeTag* DislocationSegment::GetEndTag()
	{
		return &m_oEndTag;
	}
	void DislocationSegment::SetSlipPlaneNormal(const Vector& oNormal)
	{
		m_oSlipPlaneNormal = oNormal;
		m_oSlipPlaneNormal.Normalize();
	}
	void DislocationSegment::SetBurgersVector(const Vector& oBurgersVector)
	{
		m_oBurgersVector = oBurgersVector;
	}
	Vector DislocationSegment::GetSlipPlaneNormal() const
	{
		return m_oSlipPlaneNormal;
	}
	Vector DislocationSegment::GetBurgersVector() const
	{
		return m_oBurgersVector;
	}
	void DislocationSegment::Write(FILE* fpFile) const
	{
		fprintf(fpFile,"\t%d,%d\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",m_oEndTag.GetDomainID(),m_oEndTag.GetNodeID(),m_oBurgersVector.GetX(),m_oBurgersVector.GetY(),m_oBurgersVector.GetZ(),m_oSlipPlaneNormal.GetX(),m_oSlipPlaneNormal.GetY(),m_oSlipPlaneNormal.GetZ());
	}
	void DislocationSegment::Pack(list<double>* pldData) const
	{
		pldData->push_back((double)m_oEndTag.GetDomainID());
		pldData->push_back((double)m_oEndTag.GetNodeID());
		pldData->push_back(m_oSlipPlaneNormal.GetX());
		pldData->push_back(m_oSlipPlaneNormal.GetY());
		pldData->push_back(m_oSlipPlaneNormal.GetZ());
		pldData->push_back(m_oBurgersVector.GetX());
		pldData->push_back(m_oBurgersVector.GetY());
		pldData->push_back(m_oBurgersVector.GetZ());
	}
	void DislocationSegment::Unpack(list<double>* pldData)
	{
		m_oEndTag.SetDomainID((unsigned int)pldData->front());
		pldData->pop_front();
		m_oEndTag.SetNodeID((unsigned int)pldData->front());
		pldData->pop_front();
		m_oSlipPlaneNormal.SetX(pldData->front());
		pldData->pop_front();
		m_oSlipPlaneNormal.SetY(pldData->front());
		pldData->pop_front();
		m_oSlipPlaneNormal.SetZ(pldData->front());
		pldData->pop_front();
		m_oBurgersVector.SetX(pldData->front());
		pldData->pop_front();
		m_oBurgersVector.SetY(pldData->front());
		pldData->pop_front();
		m_oBurgersVector.SetZ(pldData->front());
		pldData->pop_front();
	}
	void DislocationSegment::Initialize()
	{
		m_oEndTag.SetDomainID(0);
		m_oEndTag.SetNodeID(0);
		m_oSlipPlaneNormal.Set(0.0,0.0,0.0);
		m_oBurgersVector.Set(0.0,0.0,0.0);
	}
}



