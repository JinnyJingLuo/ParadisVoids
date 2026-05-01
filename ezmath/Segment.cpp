// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#include "Segment.h"
#include "Line.h"

namespace GeometrySystem
{
	double Segment::m_dToleranceFactor = 1E-10;
	Segment::Segment()
	{
		Initialize();
	}
	Segment::Segment(const Segment& oSegment)
	{
		*this = oSegment;
	}
	Segment::Segment(GenericNode* poPoint1,GenericNode* poPoint2,bool bIsOwner)
	{
		Initialize();
		Set(poPoint1,poPoint2,bIsOwner);
	}
	Segment::~Segment()
	{
		Reset();
	}
	Segment& Segment::operator=(const Segment& oSegment)
	{
		Set(oSegment.m_poPoint1,oSegment.m_poPoint2,oSegment.m_bIsOwner);
		return *this;
	}
	void Segment::Reset()
	{
		if(m_bIsOwner)
		{
			if(m_poPoint1 != NULL)		delete m_poPoint1;
			if(m_poPoint2 != NULL)		delete m_poPoint2;
		}
		Initialize();
	}
	void Segment::Set(GenericNode* poPoint1,GenericNode* poPoint2,bool bIsOwner)
	{
		Reset();
		m_poPoint1 = poPoint1;
		m_poPoint2 = poPoint2;
		m_oDirection.SetByPoints((const Point&)(*m_poPoint1),(const Point&)(*m_poPoint2));
		m_dLength = m_oDirection.Length();
		m_oDirection.Normalize();
		m_bIsOwner = bIsOwner;
	}
	GenericNode* Segment::GetPoint1() const
	{
		return m_poPoint1;
	}
	GenericNode* Segment::GetPoint2() const
	{
		return m_poPoint2;
	}
	Vector Segment::GetDirection() const
	{
		return m_oDirection;
	}
	bool Segment::IsPointOnSegment(const Point& oPoint) const
	{
		Line oL(m_oDirection,*m_poPoint1);
		if(!oL.IsPointOnLine(oPoint))
		{
			return false;
		}
		Vector oV(*m_poPoint1,oPoint);
		double dTemp = oV*m_oDirection;
		if((dTemp >= 0.0) && (dTemp <= m_dLength))
		{
			return true;
		}
		return false;
	}
	Point Segment::GetClosestPoint(const Point& oPoint,double& dDistance) const
	{
		Line oL(m_oDirection,*m_poPoint1);
		Point oProjection = oL.GetPointProjection(oPoint);
		Point ClosesetPoint;
		if(IsPointOnSegment(oProjection))
		{
			ClosesetPoint = oProjection;
		}
		else
		{
			dDistance = m_poPoint1->Distance(oPoint);
			ClosesetPoint = *m_poPoint1;
			if(m_poPoint2->Distance(oPoint) < dDistance)
			{
				ClosesetPoint = *m_poPoint2;
			}
		}
		dDistance = oPoint.Distance(ClosesetPoint);
		return ClosesetPoint;
	}
	bool Segment::GetIntersectionPoint(const Line& oLine,Point& oPoint) const
	{
	 	Line oSegmentLine(m_oDirection,*m_poPoint1);
	 	if(oSegmentLine.GetIntersectionPoint(oLine,oPoint))
	 	{
	 		if(IsPointOnSegment(oPoint))
	 		{
	 			return true;
	 		}
	 		return false;
	 	}
		return false;
	}
	bool Segment::GetIntersectionPoint(const Segment& oSegment,Point& oPoint) const
	{
	 	Line oSegmentLine(m_oDirection,*m_poPoint1);
	 	if(oSegment.GetIntersectionPoint(oSegmentLine,oPoint))
	 	{
	 		if(IsPointOnSegment(oPoint))
	 		{
	 			return true;
	 		}
	 		return false;
	 	}
	 	return false;
	}
	void Segment::Flip()
	{
		GenericNode* poTemp = m_poPoint1;
		m_poPoint1 = m_poPoint2;
		m_poPoint2 = poTemp;
		m_oDirection.Reverse();
	}
	bool Segment::GetIntersectionPoints(const AxisAlignedBoundingBox& oBox,Point& oStartPoint,Point& oEndPoint) const
	{
		double dTolerance = 1.0E-10;
		
		double dStart = m_poPoint1->GetX();
		double dXi1 = 0.0;
		double dXi2 = m_dLength;
		double dGradient = m_oDirection.GetX();
		if(dGradient > dTolerance)
		{
			dXi1 = (oBox.GetXMin() - dStart)/dGradient;
			dXi2 = (oBox.GetXMax() - dStart)/dGradient;
		}
		else if(dGradient < -dTolerance)
		{
			dXi2 = (oBox.GetXMin() - dStart)/dGradient;
			dXi1 = (oBox.GetXMax() - dStart)/dGradient;
		}
		if((dXi1 < 0.0) && (dXi2 < 0.0))
		{
			return false;
		}
		if((dXi1 > m_dLength) && (dXi2 > m_dLength))
		{
			return false;
		}
		if(dXi1 < 0.0)					dXi1 = 0.0;
		if(dXi2 < 0.0)					dXi2 = 0.0;
		if(dXi1 > m_dLength)			dXi1 = m_dLength;
		if(dXi2 > m_dLength)			dXi2 = m_dLength;
		
		dStart = m_poPoint1->GetY();
		double dEta1 = 0.0;
		double dEta2 = m_dLength;
		dGradient = m_oDirection.GetY();
		if(dGradient > dTolerance)
		{
			dEta1 = (oBox.GetYMin() - dStart)/dGradient;
			dEta2 = (oBox.GetYMax() - dStart)/dGradient;
		}
		else if(dGradient < -dTolerance)
		{
			dEta2 = (oBox.GetZMin() - dStart)/dGradient;
			dEta1 = (oBox.GetZMax() - dStart)/dGradient;
		}
		if((dEta1 < 0.0) && (dEta2 < 0.0))
		{
			return false;
		}
		if((dEta1 > m_dLength) && (dEta2 > m_dLength))
		{
			return false;
		}
		if(dEta1 < 0.0)					dEta1 = 0.0;
		if(dEta2 < 0.0)					dEta2 = 0.0;
		if(dEta1 > m_dLength)			dEta1 = m_dLength;
		if(dEta2 > m_dLength)			dEta2 = m_dLength;
		
		dStart = m_poPoint1->GetZ();
		double dZeta1 = 0.0;
		double dZeta2 = m_dLength;
		dGradient = m_oDirection.GetZ();
		if(dGradient > dTolerance)
		{
			dZeta1 = (oBox.GetZMin() - dStart)/dGradient;
			dZeta2 = (oBox.GetZMax() - dStart)/dGradient;
		}
		else if(dGradient < -dTolerance)
		{
			dZeta2 = (oBox.GetZMin() - dStart)/dGradient;
			dZeta1 = (oBox.GetZMax() - dStart)/dGradient;
		}
		if((dZeta1 < 0.0) && (dZeta2 < 0.0))
		{
			return false;
		}
		if((dZeta1 > m_dLength) && (dZeta2 > m_dLength))
		{
			return false;
		}
		if(dZeta1 < 0.0)				dZeta1 = 0.0;
		if(dZeta2 < 0.0)				dZeta2 = 0.0;
		if(dZeta1 > m_dLength)			dZeta1 = m_dLength;
		if(dZeta2 > m_dLength)			dZeta2 = m_dLength;
		
		// get max min and min max
		double dTIn = max(max(dXi1,dEta1),dZeta1);
		double dTOut = min(min(dXi2,dEta2),dZeta2);
		
		oStartPoint = *m_poPoint1 + m_oDirection*dTIn;
		oEndPoint = *m_poPoint1 + m_oDirection*dTOut;
		
		return true;
	}
	void Segment::Initialize()
	{
		m_poPoint1 = NULL;
		m_poPoint2 = NULL;
		m_oDirection.Set(0.0,0.0,0.0);
		m_dLength = 0.0;
		m_bIsOwner = false;
	}
}


