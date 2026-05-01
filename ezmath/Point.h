// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef POINT_H_
#define POINT_H_

#include "stdio.h"
#include "string"

const double PI = 3.14159265358979323846264338327950288419716939937510582097;

using namespace std;

namespace EZ
{
	class Point
	{
	public:

		Point();
		Point(const Point& oPoint);
		Point(const double& dX,const double& dY,const double& dZ);
		virtual ~Point();
		Point& operator=(const Point& oPoint);
		double GetX() const;
		double GetY() const;
		double GetZ() const;
		void SetX(const double& dValue);
		void SetY(const double& dValue);
		void SetZ(const double& dValue);
		void Set(const double& dX,const double& dY,const double& dZ);
		void Set(const Point& oPoint);
		virtual void Reset();
		Point operator+(const Point& oPoint) const;
		Point operator-(const Point& oPoint) const;
		Point operator*(const double& dFactor) const;
		double Distance() const;
		double Distance(const Point& oPoint) const;
		double GetDistanceSquared() const;
		double GetDistanceSquared(const Point& oPoint) const;
		double GetRX(const Point& oPoint) const;
		double GetRY(const Point& oPoint) const;
		double GetRZ(const Point& oPoint) const;
		double GetRXX(const Point& oPoint) const;
		double GetRXY(const Point& oPoint) const;
		double GetRXZ(const Point& oPoint) const;
		double GetRYX(const Point& oPoint) const;
		double GetRYY(const Point& oPoint) const;
		double GetRYZ(const Point& oPoint) const;
		double GetRZX(const Point& oPoint) const;
		double GetRZY(const Point& oPoint) const;
		double GetRZZ(const Point& oPoint) const;
		double GetComponent(const unsigned int& iIndex) const;
		void SetComponent(const unsigned int& iIndex,const double& dValue);
		double GetComponentZeroBased(const unsigned int& iIndex) const;
		void SetComponentZeroBased(const unsigned int& iIndex,const double& dValue);
		void Shift(const double& dX,const double& dY,const double& dZ);
		bool IsSame(const Point& oPoint,const double& dTolerance = 1.0E-6) const;
		bool IsOpposite(const Point& oPoint,const double& dTolerance = 1.0E-6) const;
		bool IsSimilar(const Point& oPoint,const double& dTolerance = 1.0E-6) const;
		virtual void Read(FILE* fpFile);
		virtual void Write(FILE* fpFile) const;
		string ToString() const;
		

	private:

	protected:
		double m_dX;
		double m_dY;
		double m_dZ;
	};
}

#endif

