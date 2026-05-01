#include "PolyhedronFace.h"
#include "math.h"
#include "float.h"
#include "Segment.h"

namespace GeometrySystem
{
	PolyhedronFace::PolyhedronFace()
	{
		Initialize();
	}
	PolyhedronFace::PolyhedronFace(const PolyhedronFace& oFace)
	{
		*this = oFace;
	}
	PolyhedronFace::~PolyhedronFace()
	{
		Reset();
	}
	PolyhedronFace& PolyhedronFace::operator=(const PolyhedronFace& oFace)
	{
		Face::operator=(oFace);
		TriPatch::operator=(oFace);
		m_poEdge1 = oFace.m_poEdge1;
		m_poEdge2 = oFace.m_poEdge2;
		m_poEdge3 = oFace.m_poEdge3;
		return *this;
	}
	void PolyhedronFace::Reset()
	{
		Face::Reset();
		TriPatch::Reset();
		Initialize();
	}
	void PolyhedronFace::SetPoints(GenericNode* poPoint1,GenericNode* poPoint2,GenericNode* poPoint3)
	{
		Set(poPoint1,poPoint2,poPoint3);
	}
	void PolyhedronFace::SetEdges(Edge* poEdge1,Edge* poEdge2,Edge* poEdge3)
	{
		m_poEdge1 = poEdge1;
		m_poEdge2 = poEdge2;
		m_poEdge3 = poEdge3;
	}
	bool PolyhedronFace::IsPointOnPlane(const Point& oPoint)
	{
		double dTolerance = 1.0E-6;
		if(fabs(Vector(*m_poPoint1,oPoint)*m_oNormal) < dTolerance)
		{
			return true;
		}
		return false;
	}
	Edge* PolyhedronFace::GetEdge1() const
	{
		return m_poEdge1;
	}
	Edge* PolyhedronFace::GetEdge2() const
	{
		return m_poEdge2;
	}
	Edge* PolyhedronFace::GetEdge3() const
	{
		return m_poEdge3;
	}
	bool PolyhedronFace::IsVisibleFromPoint(const Point& oPoint) const
	{
		Vector oTestVector = Vector(*m_poPoint1,oPoint);
		oTestVector.Normalize();
		if(oTestVector*m_oNormal > 0.0)
		{
			return true;
		}
		return false;
	}
	string PolyhedronFace::ToString() const
	{
		char cString[512];
		sprintf(cString,"%d : %d -> %d -> %d : %d -> %d -> %d",m_iID,m_poPoint1->GetID(),m_poPoint2->GetID(),m_poPoint3->GetID(),m_poEdge1->GetID(),m_poEdge2->GetID(),m_poEdge3->GetID());
		return string(cString);
	}
	bool PolyhedronFace::IsEdgeAligned(Edge* poEdge) const
	{
		if((poEdge->GetStartPoint() == m_poPoint1) && (poEdge->GetEndPoint() == m_poPoint2))
		{
			return true;
		}
		if((poEdge->GetStartPoint() == m_poPoint2) && (poEdge->GetEndPoint() == m_poPoint3))
		{
			return true;
		}
		if((poEdge->GetStartPoint() == m_poPoint3) && (poEdge->GetEndPoint() == m_poPoint1))
		{
			return true;
		}
		return false;
	}
	void PolyhedronFace::Initialize()
	{
		Face::Initialize();
		TriPatch::Initialize();
		m_poEdge1 = NULL;
		m_poEdge2 = NULL;
		m_poEdge3 = NULL;
	}
}


