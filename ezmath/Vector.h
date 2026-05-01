// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef VECTOR_H_
#define VECTOR_H_

// this is a general purpose 3D vector class. it supports most of the common vector
// operations, inherited from the point class

// by Ahmed M. Hussein
// August 2009

#include "Point.h"
#include "string"

using namespace std;

namespace EZ
{
	class IntegerVector
	{
	public:
		IntegerVector();
		IntegerVector(const IntegerVector& oVector);
		IntegerVector(int iX,int iY,int iZ);
		~IntegerVector();
		IntegerVector& operator=(const IntegerVector& oVector);
		IntegerVector operator+(const IntegerVector& oVector);
		IntegerVector operator-(const IntegerVector& oVector);
		int operator*(const IntegerVector& oVector);
		IntegerVector operator*(const int& iFactor);
		IntegerVector operator^(const IntegerVector& oVector);
		void Reset();
		void Set(int iX,int iY,int iZ);
		void SetX(const int& iValue);
		void SetY(const int& iValue);
		void SetZ(const int& iValue);
		int GetX() const;
		int GetY() const;
		int GetZ() const;
		bool IsZero() const;
		void Reverse();
		int GetComponent(const unsigned int& iIndex) const;
		string ToString() const;
		bool IsSimilar(const IntegerVector& oVector) const;
		void FlipX();
		void FlipY();
		void FlipZ();
		double GetNorm() const;
		bool Is111() const;

	private:

	protected:
		void Initialize();
		int m_iX;
		int m_iY;
		int m_iZ;
	};

	class Vector : public Point
	{
	public:
		Vector();
		Vector(const Vector& oVector);
		Vector(const IntegerVector& oVector);
		Vector(const double& dX,const double& dY,const double& dZ);
		Vector(const Point& oPoint);
		Vector(const Point& oPoint1,const Point& oPoint2);
		~Vector();
		Vector& operator=(const Vector& oVector);
		void SetByPoints(const Point& oPoint1,const Point& oPoint2);
		Vector operator+(const Vector& oVector) const;
		Vector operator-(const Vector& oVector) const;
		double operator*(const Vector& oVector) const;
		Vector operator*(const double& dFactor) const;
		Vector operator^(const Vector& oVector) const;
		double Length() const;
		double LengthSquared() const;
		double Magnitude() const;
		void Normalize();
		double GetAngle(const Vector& oVector) const;
		Vector GetDirection() const;
		double GetXAngle() const;
		double GetYAngle() const;
		double GetZAngle() const;
		void Reverse();
		bool IsSameDirection(const Vector& oVector) const;
		bool IsCollinear(const Vector& oVector) const;
		static Vector GenerateRandomNormalizedVector();
		static Vector GenerateRandomNormalizedVectorNormalTo(const Vector oNormal);
		IntegerVector GetEquivalentIntegerVector() const;

	private:

	protected:

	};
}

#endif

