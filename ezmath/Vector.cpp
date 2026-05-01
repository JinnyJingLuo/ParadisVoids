// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#include "Vector.h"
#include "math.h"
#include "Randomizer.h"

namespace EZ
{
	IntegerVector::IntegerVector()
	{
		Initialize();
	}
	IntegerVector::IntegerVector(const IntegerVector& oVector)
	{
		*this = oVector;
	}
	IntegerVector::IntegerVector(int iX,int iY,int iZ)
	{
		Set(iX,iY,iZ);
	}
	IntegerVector::~IntegerVector()
	{
		Reset();
	}
	IntegerVector& IntegerVector::operator=(const IntegerVector& oVector)
	{
		m_iX = oVector.m_iX;
		m_iY = oVector.m_iY;
		m_iZ = oVector.m_iZ;
		return *this;
	}
	IntegerVector IntegerVector::operator+(const IntegerVector& oVector)
	{
		IntegerVector oResult;
		oResult.Set(m_iX + oVector.m_iX,m_iY + oVector.m_iY,m_iZ + oVector.m_iZ);
		return oResult;
	}
	IntegerVector IntegerVector::operator-(const IntegerVector& oVector)
	{
		IntegerVector oResult;
		oResult.Set(m_iX - oVector.m_iX,m_iY - oVector.m_iY,m_iZ - oVector.m_iZ);
		return oResult;
	}
	int IntegerVector::operator*(const IntegerVector& oVector)
	{
		return (m_iX*oVector.m_iX + m_iY*oVector.m_iY + m_iZ*oVector.m_iZ);
	}
	IntegerVector IntegerVector::operator*(const int& iFactor)
	{
		IntegerVector oResult;
		oResult.Set(m_iX*iFactor,m_iY*iFactor,m_iZ*iFactor);
		return oResult;
	}
	IntegerVector IntegerVector::operator^(const IntegerVector& oVector)
	{
		IntegerVector oResult;
		oResult.m_iX = m_iY*oVector.m_iZ - m_iZ*oVector.m_iY;
		oResult.m_iY = m_iZ*oVector.m_iX - m_iX*oVector.m_iZ;
		oResult.m_iZ = m_iX*oVector.m_iY - m_iY*oVector.m_iX;
		return oResult;
	}
	void IntegerVector::Reset()
	{
		m_iX = 0;
		m_iY = 0;
		m_iZ = 0;
	}
	void IntegerVector::Set(int iX,int iY,int iZ)
	{
		m_iX = iX;
		m_iY = iY;
		m_iZ = iZ;
	}
	void IntegerVector::SetX(const int& iValue)
	{
		m_iX = iValue;
	}
	void IntegerVector::SetY(const int& iValue)
	{
		m_iY = iValue;
	}
	void IntegerVector::SetZ(const int& iValue)
	{
		m_iZ = iValue;
	}
	int IntegerVector::GetX() const
	{
		return m_iX;
	}
	int IntegerVector::GetY() const
	{
		return m_iY;
	}
	int IntegerVector::GetZ() const
	{
		return m_iZ;
	}
	bool IntegerVector::IsZero() const
	{
		if(m_iX == 0)
		{
			if(m_iY == 0)
			{
				if(m_iZ == 0)
				{
					return true;
				}
			}
		}
		return false;
	}
	void IntegerVector::Reverse()
	{
		m_iX = -m_iX;
		m_iY = -m_iY;
		m_iZ = -m_iZ;
	}
	int IntegerVector::GetComponent(const unsigned int& iIndex) const
	{
		if(iIndex == 1)
		{
			return m_iX;
		}
		if(iIndex == 2)
		{
			return m_iY;
		}
		if(iIndex == 3)
		{
			return m_iZ;
		}
		return 0;
	}
	string IntegerVector::ToString() const
	{
		char cWrite[512];
		sprintf(cWrite,"%d,%d,%d",m_iX,m_iY,m_iZ);
		string sString = cWrite;
		return sString;
	}
	void IntegerVector::Initialize()
	{
		m_iX = 0;
		m_iY = 0;
		m_iZ = 0;
	}
	bool IntegerVector::IsSimilar(const IntegerVector& oVector) const
	{
		if((m_iX == oVector.m_iX) && (m_iY == oVector.m_iY) && (m_iZ == oVector.m_iZ))
		{
			return true;
		}
		if((m_iX == -oVector.m_iX) && (m_iY == -oVector.m_iY) && (m_iZ == -oVector.m_iZ))
		{
			return true;
		}
		return false;
	}
	void IntegerVector::FlipX()
	{
		m_iX = -m_iX;
	}
	void IntegerVector::FlipY()
	{
		m_iY = -m_iY;
	}
	void IntegerVector::FlipZ()
	{
		m_iZ = -m_iZ;
	}
	double IntegerVector::GetNorm() const
	{
		double dNorm = m_iX*m_iX + m_iY*m_iY + m_iZ*m_iZ;
		dNorm = sqrt(dNorm);
		return dNorm;
	}
	bool IntegerVector::Is111() const
	{
		if((m_iX == m_iY) || (m_iX == -m_iY))
		{
			if((m_iX == m_iZ) || (m_iX == -m_iZ))
			{
				return true;
			}
		}
		return false;
	}


	Vector::Vector():Point()
	{

	}
	Vector::Vector(const Vector& oVector)
	{
		*this = oVector;
	}
	Vector::Vector(const IntegerVector& oVector)
	{
		Set(oVector.GetX(),oVector.GetY(),oVector.GetZ());
	}
	Vector::Vector(const double& dX,const double& dY,const double& dZ)
	{
		Set(dX,dY,dZ);
	}
	Vector::Vector(const Point& oPoint)
	{
		Set(oPoint.GetX(),oPoint.GetY(),oPoint.GetZ());
	}
	Vector::Vector(const Point& oPoint1,const Point& oPoint2)
	{
		SetByPoints(oPoint1,oPoint2);
	}
	Vector::~Vector()
	{

	}
	Vector& Vector::operator=(const Vector& oVector)
	{
		Set(oVector.m_dX,oVector.m_dY,oVector.m_dZ);
		return *this;
	}
	void Vector::SetByPoints(const Point& oPoint1,const Point& oPoint2)
	{
		Point oRelative = oPoint2 - oPoint1;
		Set(oRelative.GetX(),oRelative.GetY(),oRelative.GetZ());
	}
	Vector Vector::operator+(const Vector& oVector) const
	{
		return Vector(m_dX + oVector.m_dX,m_dY + oVector.m_dY,m_dZ + oVector.m_dZ);
	}
	Vector Vector::operator-(const Vector& oVector) const
	{
		return Vector(m_dX - oVector.m_dX,m_dY - oVector.m_dY,m_dZ - oVector.m_dZ);
	}
	double Vector::operator*(const Vector& oVector) const
	{
		return (m_dX*oVector.m_dX + m_dY*oVector.m_dY + m_dZ*oVector.m_dZ);
	}
	Vector Vector::operator*(const double& dFactor) const
	{
		return Vector(dFactor*m_dX,dFactor*m_dY,dFactor*m_dZ);
	}
	Vector Vector::operator^(const Vector& oVector) const
	{
		Vector oResult;
		oResult.SetX(m_dY*oVector.m_dZ - m_dZ*oVector.m_dY);
		oResult.SetY(m_dZ*oVector.m_dX - m_dX*oVector.m_dZ);
		oResult.SetZ(m_dX*oVector.m_dY - m_dY*oVector.m_dX);
		return oResult;
	}
	double Vector::Length() const
	{
		return Distance();
	}
	double Vector::LengthSquared() const
	{
		return (m_dX*m_dX + m_dY*m_dY + m_dZ*m_dZ);
	}
	void Vector::Normalize()
	{
		double dLength = Length();
		if(dLength < 1E-6)
		{
			m_dX = 0.0;
			m_dY = 0.0;
			m_dZ = 0.0;
		}
		else
		{
			m_dX = m_dX/dLength;
			m_dY = m_dY/dLength;
			m_dZ = m_dZ/dLength;
		}
	}
	double Vector::GetAngle(const Vector& oVector) const
	{
		Vector oV1 = *this;
		Vector oV2 = oVector;
		oV1.Normalize();
		oV2.Normalize();
		return acos(oV1*oV2);
	}
	Vector Vector::GetDirection() const
	{
		Vector oDirection = *this;
		oDirection.Normalize();
		return oDirection;
	}
	double Vector::GetXAngle() const
	{
		Vector oV = *this;
		oV.Normalize();
		return acos(oV.GetX());
	}
	double Vector::GetYAngle() const
	{
		Vector oV = *this;
		oV.Normalize();
		return acos(oV.GetY());
	}
	double Vector::GetZAngle() const
	{
		Vector oV = *this;
		oV.Normalize();
		return acos(oV.GetZ());
	}
	void Vector::Reverse()
	{
		m_dX = -m_dX;
		m_dY = -m_dY;
		m_dZ = -m_dZ;
	}
	double Vector::Magnitude() const
	{
		return Length();
	}
	bool Vector::IsSameDirection(const Vector& oVector) const
	{
		Vector oCrossProduct = *this^oVector;
		double dTolerance = 1.0E-8;
		if(oCrossProduct.LengthSquared() > dTolerance)
		{
			return false;
		}
		if(*this*oVector > 0.0)
		{
			return true;
		}
		return false;
	}
	Vector Vector::GenerateRandomNormalizedVector()
	{
		double dX = 0.0;
		double dY = 0.0;
		double dZ = 0.0;
		double dLengthSquared = 0.0;
		double dToleranceSquared = 1.0E-12;
		double dLength = 0.0;
		while(true)
		{
			dX = Randomizer::Random(-1.0,1.0);
			dY = Randomizer::Random(-1.0,1.0);
			dZ = Randomizer::Random(-1.0,1.0);
			dLengthSquared = dX*dX + dY*dY + dZ*dZ;
			if(dLengthSquared > dToleranceSquared)
			{
				dLength = sqrt(dLengthSquared);
				dX = dX/dLength;
				dY = dY/dLength;
				dZ = dZ/dLength;
				break;
			}
		}
		return Vector(dX,dY,dZ);
	}
	Vector Vector::GenerateRandomNormalizedVectorNormalTo(const Vector oNormal)
	{
		Vector oVector;
		Vector oWorkingNormal = oNormal;
		oWorkingNormal.Normalize();
		double dProjection = 0.0;
		double dTolerance = 1.0E-6;
		while(true)
		{
			oVector = GenerateRandomNormalizedVector();
			dProjection = oVector*oWorkingNormal;
			oVector = oVector - oWorkingNormal*dProjection;
			if(oVector.Length() > dTolerance)
			{
				oVector.Normalize();
				break;
			}
		}
		return oVector;
	}
	IntegerVector Vector::GetEquivalentIntegerVector() const
	{
		double dTolerance = 1.0e-9;
		double dX = fabs(m_dX);
		double dY = fabs(m_dY);
		double dZ = fabs(m_dZ);
		IntegerVector oVector(0,0,0);
		if(dX > dTolerance)
		{
			if((dY > dX) && (dZ > dX))
			{
				// this component is the smallest nonzero
				if(m_dX > 0.0)
				{
					oVector.SetX(1);
				}
				else
				{
					oVector.SetX(-1);
				}
				oVector.SetY((int)floor(m_dY/m_dX + 0.5));
				oVector.SetZ((int)floor(m_dZ/m_dX + 0.5));
			}
		}
		else if(dY > dTolerance)
		{
			if((dX > dY) && (dZ > dY))
			{
				// this component is the smallest nonzero
				if(m_dY > 0.0)
				{
					oVector.SetY(1);
				}
				else
				{
					oVector.SetY(-1);
				}
				oVector.SetX((int)floor(m_dX/m_dY + 0.5));
				oVector.SetZ((int)floor(m_dZ/m_dY + 0.5));
			}
		}
		else if(dZ > dTolerance)
		{
			if((dX > dZ) && (dY > dZ))
			{
				// this component is the smallest nonzero
				if(m_dZ > 0.0)
				{
					oVector.SetZ(1);
				}
				else
				{
					oVector.SetZ(-1);
				}
				oVector.SetX((int)floor(m_dX/m_dZ + 0.5));
				oVector.SetY((int)floor(m_dY/m_dZ + 0.5));
			}
		}
		return oVector;
	}
}


