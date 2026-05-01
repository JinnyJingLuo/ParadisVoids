// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef LINE_H_
#define LINE_H_

#include "Vector.h"

using namespace EZ;

namespace GeometrySystem
{
	class Line
	{
	public:
		Line();
		Line(const Line& oLine);
		Line(const Vector& oDirection,const Point& oPoint);
		Line(const Point& oPoint1,const Point& oPoint2);
		~Line();
		Line& operator=(const Line& oLine);
		void Set(const Vector& oDirection,const Point& oPoint);
		void Set(const Point& oPoint1,const Point& oPoint2);
		void SetDirection(const Vector& oDirection);
		void SetPoint(const Point& oPoint);
		Vector GetDirection() const;
		Point GetPoint() const;
		bool IsPointOnLine(const Point& oPoint) const;
		void ReverseDirection();
		Point GetPointProjection(const Point& oPoint) const;
		bool GetNearestPoints(const Line& oLine,Point& oPoint1,Point& oPoint2) const;
		bool GetIntersectionPoint(const Line& oLine,Point& oPoint,const double& dTolerance = 1E-6) const;
	private:

	protected:
		Vector m_oDirection;
		Point m_oPoint;
	};
}


#endif


