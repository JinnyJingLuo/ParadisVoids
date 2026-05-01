// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef PLANE_H_
#define PLANE_H_

#include "Vector.h"
#include "Line.h"
#include "AxisAlignedBoundingBox.h"

using namespace EZ;

namespace GeometrySystem
{
	class Plane
	{
	public:
		Plane();
		Plane(const Plane& oPlane);
		Plane(const Vector& oNormal,const Point& oPoint);
		~Plane();
		Plane& operator=(const Plane& oPlane);
		void Reset();
		void Set(const Vector& oNormal,const Point& oPoint);
		void SetNormal(const Vector& oNormal);
		void SetPoint(const Point& oPoint);
		Vector GetNormal() const;
		Point GetPoint() const;
		Vector* GetNormal();
		Point* GetPoint();
		int ClassifyPoint(const Point& oPoint) const;
		int ClassifyBox(AxisAlignedBoundingBox* poBox) const;
		bool IsPointOnOrAbove(const Point& oPoint) const;
		void ReverseNormal();
		bool IsParallel(const Plane& oPlane) const;
		Point GetPointProjection(const Point& oPoint) const;
		bool GetIntersectionLine(const Plane& oPlane,Line& oIntersectionLine) const;
		bool GetLineIntersection(const Line& oLine,Point& oIntersectionPoint) const;
		bool GetSegmentIntersection(const Line& oLine,const double& dLength,Point& oIntersectionPoint) const;
		bool GetSegmentIntersection(const Point& oStartPoint,const Point& oEndPoint,Point& oIntersectionPoint) const;
		double GetOriginSignedDistance() const;
		double GetPointDistance(const Point& oPoint) const;
		
	private:

	protected:
		void Initialize();
		Vector m_oNormal;
		Point m_oPoint;
	};
}

#endif

