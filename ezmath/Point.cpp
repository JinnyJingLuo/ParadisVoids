// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#include "Point.h"
#include "math.h"
#include "Tools.h"
#include "string"

using namespace std;
using namespace SupportSystem;

namespace EZ
{
	Point::Point()
	{
		Reset();
	}
	Point::Point(const Point& oPoint)
	{
		*this = oPoint;
	}
	Point::Point(const double& dX,const double& dY,const double& dZ)
	{
		Set(dX,dY,dZ);
	}
	Point::~Point()
	{

	}
	Point& Point::operator=(const Point& oPoint)
	{
		Set(oPoint.m_dX,oPoint.m_dY,oPoint.m_dZ);
		return *this;
	}
	double Point::GetX() const
	{
		return m_dX;
	}
	double Point::GetY() const
	{
		return m_dY;
	}
	double Point::GetZ() const
	{
		return m_dZ;
	}
	void Point::SetX(const double& dValue)
	{
		m_dX = dValue;
	}
	void Point::SetY(const double& dValue)
	{
		m_dY = dValue;
	}
	void Point::SetZ(const double& dValue)
	{
		m_dZ = dValue;
	}
	void Point::Set(const double& dX,const double& dY,const double& dZ)
	{
		m_dX = dX;
		m_dY = dY;
		m_dZ = dZ;
	}
	void Point::Reset()
	{
		m_dX = 0.0;
		m_dY = 0.0;
		m_dZ = 0.0;
	}
	Point Point::operator+(const Point& oPoint) const
	{
		return Point(m_dX + oPoint.m_dX,m_dY + oPoint.m_dY,m_dZ + oPoint.m_dZ);
	}
	Point Point::operator-(const Point& oPoint) const
	{
		return Point(m_dX - oPoint.m_dX,m_dY - oPoint.m_dY,m_dZ - oPoint.m_dZ);
	}
	Point Point::operator*(const double& dFactor) const
	{
		return Point(dFactor*m_dX,dFactor*m_dY,dFactor*m_dZ);
	}
	double Point::Distance() const
	{
		return sqrt(m_dX*m_dX + m_dY*m_dY + m_dZ*m_dZ);
	}
	double Point::Distance(const Point& oPoint) const
	{
		Point oRelative = *this - oPoint;
		return oRelative.Distance();
	}
	double Point::GetDistanceSquared() const
	{
		return (m_dX*m_dX + m_dY*m_dY + m_dZ*m_dZ);
	}
	double Point::GetDistanceSquared(const Point& oPoint) const
	{
		double dX = m_dX - oPoint.m_dX;
		double dY = m_dY - oPoint.m_dY;
		double dZ = m_dZ - oPoint.m_dZ;
		return (dX*dX + dY*dY + dZ*dZ);
	}
	double Point::GetRX(const Point& oPoint) const
	{
		Point oRelative = oPoint - (*this);
		return (oRelative.m_dX/oRelative.Distance());
	}
	double Point::GetRY(const Point& oPoint) const
	{
		Point oRelative = oPoint - (*this);
		return (oRelative.m_dY/oRelative.Distance());
	}
	double Point::GetRZ(const Point& oPoint) const
	{
		Point oRelative = oPoint - (*this);
		return (oRelative.m_dZ/oRelative.Distance());
	}
	double Point::GetRXX(const Point& oPoint) const
	{
		Point oRelative = oPoint - (*this);
		double dDistance = oRelative.Distance();
		double dDistanceSquared = dDistance*dDistance;
		double dDistanceCubed = dDistanceSquared*dDistance;
		return ((dDistanceSquared - oRelative.m_dX*oRelative.m_dX)/dDistanceCubed);
	}
	double Point::GetRXY(const Point& oPoint) const
	{
		Point oRelative = oPoint - (*this);
		double dDistance = oRelative.Distance();
		double dDistanceCubed = dDistance*dDistance*dDistance;
		return ((- oRelative.m_dX*oRelative.m_dY)/dDistanceCubed);
	}
	double Point::GetRXZ(const Point& oPoint) const
	{
		Point oRelative = oPoint - (*this);
		double dDistance = oRelative.Distance();
		double dDistanceCubed = dDistance*dDistance*dDistance;
		return ((- oRelative.m_dX*oRelative.m_dZ)/dDistanceCubed);
	}
	double Point::GetRYX(const Point& oPoint) const
	{
		return GetRXY(oPoint);
	}
	double Point::GetRYY(const Point& oPoint) const
	{
		Point oRelative = oPoint - (*this);
		double dDistance = oRelative.Distance();
		double dDistanceSquared = dDistance*dDistance;
		double dDistanceCubed = dDistanceSquared*dDistance;
		return ((dDistanceSquared - oRelative.m_dY*oRelative.m_dY)/dDistanceCubed);
	}
	double Point::GetRYZ(const Point& oPoint) const
	{
		Point oRelative = oPoint - (*this);
		double dDistance = oRelative.Distance();
		double dDistanceCubed = dDistance*dDistance*dDistance;
		return ((- oRelative.m_dY*oRelative.m_dZ)/dDistanceCubed);
	}
	double Point::GetRZX(const Point& oPoint) const
	{
		return GetRXZ(oPoint);
	}
	double Point::GetRZY(const Point& oPoint) const
	{
		return GetRYZ(oPoint);
	}
	double Point::GetRZZ(const Point& oPoint) const
	{
		Point oRelative = oPoint - (*this);
		double dDistance = oRelative.Distance();
		double dDistanceSquared = dDistance*dDistance;
		double dDistanceCubed = dDistanceSquared*dDistance;
		return ((dDistanceSquared - oRelative.m_dZ*oRelative.m_dZ)/dDistanceCubed);
	}
	void Point::Shift(const double& dX,const double& dY,const double& dZ)
	{
		m_dX = m_dX + dX;
		m_dY = m_dY + dY;
		m_dZ = m_dZ + dZ;
	}
	void Point::Read(FILE* fpFile)
	{
		string sRead = GetRealString(500,fpFile);
		sscanf(sRead.c_str(),"%lf\t\t%lf\t\t%lf",&m_dX,&m_dY,&m_dZ);
	}
	void Point::Write(FILE* fpFile) const
	{
		char cWrite[500];
		sprintf(cWrite,"%E\t\t%E\t\t%E",m_dX,m_dY,m_dZ);
		fputs(cWrite,fpFile);
	}
	double Point::GetComponent(const unsigned int& iIndex) const
	{
		if(iIndex == 1)
		{
			return m_dX;
		}
		if(iIndex == 2)
		{
			return m_dY;
		}
		if(iIndex == 3)
		{
			return m_dZ;
		}
		return 0.0;
	}
	void Point::SetComponent(const unsigned int& iIndex,const double& dValue)
	{
		if(iIndex == 1)
		{
			m_dX = dValue;
		}
		if(iIndex == 2)
		{
			m_dY = dValue;
		}
		if(iIndex == 3)
		{
			m_dZ = dValue;
		}
	}
	double Point::GetComponentZeroBased(const unsigned int& iIndex) const
	{
		return GetComponent(iIndex + 1);
	}
	void Point::SetComponentZeroBased(const unsigned int& iIndex,const double& dValue)
	{
		SetComponent(iIndex + 1,dValue);
	}
	bool Point::IsSame(const Point& oPoint,const double& dTolerance) const
	{
		if(fabs(m_dX - oPoint.m_dX) < dTolerance)
		{
			if(fabs(m_dY - oPoint.m_dY) < dTolerance)
			{
				if(fabs(m_dZ - oPoint.m_dZ) < dTolerance)
				{
					return true;
				}
			}
		}
		return false;
	}
	bool Point::IsOpposite(const Point& oPoint,const double& dTolerance) const
	{
		if(fabs(m_dX + oPoint.m_dX) < dTolerance)
		{
			if(fabs(m_dY + oPoint.m_dY) < dTolerance)
			{
				if(fabs(m_dZ + oPoint.m_dZ) < dTolerance)
				{
					return true;
				}
			}
		}
		return false;
	}
	bool Point::IsSimilar(const Point& oPoint,const double& dTolerance) const
	{
		return (IsSame(oPoint,dTolerance) || IsOpposite(oPoint,dTolerance));
	}
	void Point::Set(const Point& oPoint)
	{
		m_dX = oPoint.m_dX;
		m_dY = oPoint.m_dY;
		m_dZ = oPoint.m_dZ;
	}
	string Point::ToString() const
	{
		char cString[256];
		sprintf(cString,"(%25.20f,%25.20f,%25.20f)",m_dX,m_dY,m_dZ);
		string sString = cString;
		return sString;
	}
}

