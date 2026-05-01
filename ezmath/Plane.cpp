// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#include "Plane.h"
#include "Matrix.h"
#include "math.h"

namespace GeometrySystem
{
	Plane::Plane()
	{
		Initialize();
	}
	Plane::Plane(const Plane& oPlane)
	{
		*this = oPlane;
	}
	Plane::Plane(const Vector& oNormal,const Point& oPoint)
	{
		Set(oNormal,oPoint);
	}
	Plane::~Plane()
	{

	}
	Plane& Plane::operator=(const Plane& oPlane)
	{
		Set(oPlane.m_oNormal,oPlane.m_oPoint);
		return *this;
	}
	void Plane::Reset()
	{
		Initialize();
	}
	void Plane::Initialize()
	{
		m_oNormal.Set(0.0,0.0,1.0);
		m_oPoint.Set(0.0,0.0,0.0);
	}
	void Plane::Set(const Vector& oNormal,const Point& oPoint)
	{
		SetNormal(oNormal);
		SetPoint(oPoint);
	}
	void Plane::SetNormal(const Vector& oNormal)
	{
		if(oNormal.Length() < 1E-14)
		{
			return;
		}
		m_oNormal = oNormal;
		m_oNormal.Normalize();
	}
	void Plane::SetPoint(const Point& oPoint)
	{
		m_oPoint = oPoint;
	}
	Vector Plane::GetNormal() const
	{
		return m_oNormal;
	}
	Point Plane::GetPoint() const
	{
		return m_oPoint;
	}
	Vector* Plane::GetNormal()
	{
		return &m_oNormal;
	}
	Point* Plane::GetPoint()
	{
		return &m_oPoint;
	}
	int Plane::ClassifyPoint(const Point& oPoint) const
	{
		Vector oV(m_oPoint,oPoint);
		double dProduct = oV*m_oNormal;
		if(dProduct > 0.0)
		{
			return 1;
		}
		else if(dProduct < 0.0)
		{
			return -1;
		}
		return 0;
	}
	int Plane::ClassifyBox(AxisAlignedBoundingBox* poBox) const
	{
		double dXMin = poBox->GetXMin();
		double dXMax = poBox->GetXMax();
		double dYMin = poBox->GetYMin();
		double dYMax = poBox->GetYMax();
		double dZMin = poBox->GetZMin();
		double dZMax = poBox->GetZMax();
		int iSum = 0;
		iSum = iSum + ClassifyPoint(Point(dXMin,dYMin,dZMin));
		iSum = iSum + ClassifyPoint(Point(dXMax,dYMin,dZMin));
		iSum = iSum + ClassifyPoint(Point(dXMax,dYMax,dZMin));
		iSum = iSum + ClassifyPoint(Point(dXMin,dYMax,dZMin));
		iSum = iSum + ClassifyPoint(Point(dXMin,dYMin,dZMax));
		iSum = iSum + ClassifyPoint(Point(dXMax,dYMin,dZMax));
		iSum = iSum + ClassifyPoint(Point(dXMax,dYMax,dZMax));
		iSum = iSum + ClassifyPoint(Point(dXMin,dYMax,dZMax));
		
		if(iSum == 8)
		{
			return 1;
		}
		else if(iSum == -8)
		{
			return -1;
		}
		return 0;
	}
	bool Plane::IsPointOnOrAbove(const Point& oPoint) const
	{
		if(ClassifyPoint(oPoint) >= 0)
		{
			return true;
		}
		return false;
	}
	void Plane::ReverseNormal()
	{
		m_oNormal.Reverse();
	}
	bool Plane::IsParallel(const Plane& oPlane) const
	{
		double dTolerance = 1E-6;
		if(fabs(fabs(m_oNormal*oPlane.m_oNormal) - 1.0) < dTolerance)
		{
			return true;
		}
		return false;
	}
	Point Plane::GetPointProjection(const Point& oPoint) const
	{
		Vector oV(m_oPoint,oPoint);
		double dT = oV*m_oNormal;
		Point oProjection = oPoint - m_oNormal*dT;
		return oProjection;
	}
	bool Plane::GetIntersectionLine(const Plane& oPlane,Line& oIntersectionLine) const
	{
		if(IsParallel(oPlane))
		{
			return false;
		}
		// get the plane equations for both planes, we already
		// know the coefficients (from the plane normals) of the 
		// x,y and z terms, we are missing the constant d
		// the form of the equation is ax + by + cz = d
		double dD1 = m_oNormal*m_oPoint;
		double dD2 = oPlane.m_oNormal*oPlane.m_oPoint;
		Vector oLineDirection = m_oNormal^oPlane.m_oNormal;
		oLineDirection.Normalize();
		double dTolernace = 1E-6;
		// the system matrix and rhs
		Matrix oSystem(2,2);
		Matrix oRHS(2,1);
		
		oRHS.Set(1,1,dD1);
		oRHS.Set(2,1,dD2);
		unsigned int iZeroComponent = 0;
		if(fabs(oLineDirection.GetZ()) > dTolernace)
		{
			oSystem.Set(1,1,m_oNormal.GetX());
			oSystem.Set(1,2,m_oNormal.GetY());
			oSystem.Set(2,1,oPlane.m_oNormal.GetX());
			oSystem.Set(2,2,oPlane.m_oNormal.GetY());
			iZeroComponent = 3;
		}
		else if(fabs(oLineDirection.GetY()) > dTolernace)
		{
			oSystem.Set(1,1,m_oNormal.GetX());
			oSystem.Set(1,2,m_oNormal.GetZ());
			oSystem.Set(2,1,oPlane.m_oNormal.GetX());
			oSystem.Set(2,2,oPlane.m_oNormal.GetZ());
			iZeroComponent = 2;
		}
		else
		{
			oSystem.Set(1,1,m_oNormal.GetY());
			oSystem.Set(1,2,m_oNormal.GetZ());
			oSystem.Set(2,1,oPlane.m_oNormal.GetY());
			oSystem.Set(2,2,oPlane.m_oNormal.GetZ());
			iZeroComponent = 1;
		}
		Matrix oSolution = Matrix::Solve2x2System(oSystem,oRHS);
		Point oLinePoint;
		if(iZeroComponent == 1)
		{
			oLinePoint.SetX(0.0);
			oLinePoint.SetY(oSolution.Get(1,1));
			oLinePoint.SetZ(oSolution.Get(2,1));
		}
		else if(iZeroComponent == 2)
		{
			oLinePoint.SetX(oSolution.Get(1,1));
			oLinePoint.SetY(0.0);
			oLinePoint.SetZ(oSolution.Get(2,1));
		}
		else
		{
			oLinePoint.SetX(oSolution.Get(1,1));
			oLinePoint.SetY(oSolution.Get(2,1));
			oLinePoint.SetZ(0.0);
		}
		oIntersectionLine = Line(oLineDirection,oLinePoint);
		return true;
	}
	bool Plane::GetLineIntersection(const Line& oLine,Point& oIntersectionPoint) const
	{
		double dTolerance = 1E-6;
		Vector oLineDirection = oLine.GetDirection();
		double dDen = m_oNormal*oLineDirection;
		if(fabs(dDen) < dTolerance)
		{
			return false;
		}
		Point oLinePoint = oLine.GetPoint();
		double dNum = m_oNormal*oLinePoint;
		dNum = dNum - m_oNormal*m_oPoint;
		double dT = -dNum/dDen;
		oIntersectionPoint = oLinePoint + oLineDirection*dT;
		return true;
	}
	bool Plane::GetSegmentIntersection(const Line& oLine,const double& dLength,Point& oIntersectionPoint) const
	{
		double dTolerance = 1E-6;
		Vector oLineDirection = oLine.GetDirection();
		double dDen = m_oNormal*oLineDirection;
		if(fabs(dDen) < dTolerance)
		{
			return false;
		}
		Point oLinePoint = oLine.GetPoint();
		double dNum = m_oNormal*oLinePoint;
		dNum = dNum - m_oNormal*m_oPoint;
		double dT = -dNum/dDen;
		if(dT < 0.0)
		{
			return false;
		}
		if(dT > dLength)
		{
			return false;
		}
		oIntersectionPoint = oLinePoint + oLineDirection*dT;
		return true;
	}
	bool Plane::GetSegmentIntersection(const Point& oStartPoint,const Point& oEndPoint,Point& oIntersectionPoint) const
	{
		Line oLine(oStartPoint,oEndPoint);
		double dLength = oStartPoint.Distance(oEndPoint);
		return GetSegmentIntersection(oLine,dLength,oIntersectionPoint);
	}
	double Plane::GetOriginSignedDistance() const
	{
		Vector oTemp = m_oPoint*(-1.0);
		return oTemp*m_oNormal;
	}
	double Plane::GetPointDistance(const Point& oPoint) const
	{
		Vector oV(m_oPoint,oPoint);
		return fabs(oV*m_oNormal);
	}
}


