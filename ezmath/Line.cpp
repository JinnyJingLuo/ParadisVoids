// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#include "Line.h"
#include "math.h"
#include "Matrix.h"

namespace GeometrySystem
{
	Line::Line()
	{
		m_oDirection.Set(0.0,0.0,1.0);
		m_oPoint.Set(0.0,0.0,0.0);
	}
	Line::Line(const Line& oLine)
	{
		*this = oLine;
	}
	Line::Line(const Vector& oDirection,const Point& oPoint)
	{
		Set(oDirection,oPoint);
	}
	Line::Line(const Point& oPoint1,const Point& oPoint2)
	{
		Set(oPoint1,oPoint2);
	}
	Line::~Line()
	{
	
	}
	Line& Line::operator=(const Line& oLine)
	{
		Set(oLine.m_oDirection,oLine.m_oPoint);
		return *this;
	}
	void Line::Set(const Vector& oDirection,const Point& oPoint)
	{
		SetDirection(oDirection);
		SetPoint(oPoint);
	}
	void Line::Set(const Point& oPoint1,const Point& oPoint2)
	{
		Vector oDirection(oPoint1,oPoint2);
		SetDirection(oDirection);
		SetPoint(oPoint1);
	}
	void Line::SetDirection(const Vector& oDirection)
	{
		m_oDirection = oDirection;
		m_oDirection.Normalize();
	}
	void Line::SetPoint(const Point& oPoint)
	{
		m_oPoint = oPoint;
	}
	Vector Line::GetDirection() const
	{
		return m_oDirection;
	}
	Point Line::GetPoint() const
	{
		return m_oPoint;
	}
	bool Line::IsPointOnLine(const Point& oPoint) const
	{
		Vector oV(m_oPoint,oPoint);
		oV.Normalize();
		double dTolerance = 1E-6;
		if(fabs(fabs(oV*m_oDirection) - 1.0) < dTolerance)
		{
			return true;
		}
		return false;
	}
	void Line::ReverseDirection()
	{
		m_oDirection.Reverse();
	}
	Point Line::GetPointProjection(const Point& oPoint) const
	{
		Vector oV(m_oPoint,oPoint);
		double dT = oV*m_oDirection;
		Point oProjection = m_oPoint + m_oDirection*dT;
		return oProjection;
	}
	bool Line::GetNearestPoints(const Line& oLine,Point& oPoint1,Point& oPoint2) const
	{
		double dOffDiagonalEntry = m_oDirection*oLine.m_oDirection;
		double dTolerance = 1E-6;
		if(fabs(fabs(dOffDiagonalEntry) - 1.0) < dTolerance)
		{
			return false;
		}
		Matrix oSystem(2,2);
		Matrix oRHS(2,1);
		oSystem.Set(1,1,1.0);
		oSystem.Set(1,2,-dOffDiagonalEntry);
		oSystem.Set(2,1,-dOffDiagonalEntry);
		oSystem.Set(2,2,1.0);
		Point oPointDifference = m_oPoint - oLine.m_oPoint;
		oRHS.Set(1,1,-(m_oDirection*oPointDifference));
		oRHS.Set(2,1,(oLine.m_oDirection*oPointDifference));
		Matrix oSolution = Matrix::Solve2x2System(oSystem,oRHS);
		oPoint1 = m_oPoint + m_oDirection*oSolution.Get(1,1);
		oPoint2 = oLine.m_oPoint + oLine.m_oDirection*oSolution.Get(2,1);
		return true;
	}
	bool Line::GetIntersectionPoint(const Line& oLine,Point& oPoint,const double& dTolerance) const
	{
		Point oP1;
		Point oP2;
		if(GetNearestPoints(oLine,oP1,oP2))
		{
			double dDistance = oP1.Distance(oP2);
			if(dDistance <= dTolerance)
			{
				oPoint = oP1;
				return true;
			}
			return false;
		}
		return false;
	}
}


