#include "Edge.h"
#include "Face.h"

namespace GeometrySystem
{
	Edge::Edge()
	{
		Initialize();
	}
	Edge::Edge(const Edge& oEdge)
	{
		*this = oEdge;
	}
	Edge::~Edge()
	{
		Reset();
	}
	Edge& Edge::operator=(const Edge& oEdge)
	{
		GeometricComponent::operator=(oEdge);
		m_poStartPoint = oEdge.m_poStartPoint;
		m_poEndPoint = oEdge.m_poEndPoint;
		m_poRightFace = oEdge.m_poRightFace;
		m_poLeftFace = oEdge.m_poLeftFace;
		return *this;
	}
	void Edge::Reset()
	{
		Initialize();
	}
	void Edge::SetEndPoints(GenericNode* poStartPoint,GenericNode* poEndPoint)
	{
		m_poStartPoint = poStartPoint;
		m_poEndPoint = poEndPoint;
	}
	void Edge::SetRightFace(Face* poFace)
	{
		m_poRightFace = poFace;
	}
	void Edge::SetLeftFace(Face* poFace)
	{
		m_poLeftFace = poFace;
	}
	GenericNode* Edge::GetStartPoint() const
	{
		return m_poStartPoint;
	}
	GenericNode* Edge::GetEndPoint() const
	{
		return m_poEndPoint;
	}
	Face* Edge::GetRightFace() const
	{
		return m_poRightFace;
	}
	Face* Edge::GetLeftFace() const
	{
		return m_poLeftFace;
	}
	Face* Edge::GetOtherFace(Face* poFace) const
	{
		if(poFace == m_poRightFace)
		{
			return m_poLeftFace;
		}
		if(poFace == m_poLeftFace)
		{
			return m_poRightFace;
		}
		return NULL;
	}
	string Edge::ToString() const
	{
		char cString[512];
		sprintf(cString,"%d : %d -> %d\t\tL : %d\t\tR : %d",m_iID,m_poStartPoint->GetID(),m_poEndPoint->GetID(),m_poLeftFace->GetID(),m_poRightFace->GetID());
		return string(cString);
	}
	void Edge::ReplaceFace(Face* poOldFace,Face* poNewFace)
	{
		if(poOldFace == m_poRightFace)
		{
			m_poRightFace = poNewFace;
		}
		if(poOldFace == m_poLeftFace)
		{
			m_poLeftFace = poNewFace;
		}
	}
	void Edge::RemoveFace(Face* poFace)
	{
		ReplaceFace(poFace,NULL);
	}
	Vector Edge::GetVector() const
	{
		return Vector(*m_poStartPoint,*m_poEndPoint);
	}
	double Edge::GetLength() const
	{
		return m_poStartPoint->Distance(*m_poEndPoint);
	}
	void Edge::Initialize()
	{
		GeometricComponent::Initialize();
		m_poStartPoint = NULL;
		m_poEndPoint = NULL;
		m_poRightFace = NULL;
		m_poLeftFace = NULL;
	}
}



