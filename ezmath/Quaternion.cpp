// Paradis Processor Project
// Ahmed M. Hussein (mailto : am.hussin@gmail.com)
// June 2012

#include "Quaternion.h"
#include "math.h"

namespace EZ
{
	Quaternion::Quaternion()
	{
		m_dScalar = 0.0;
		m_oVector.Set(0.0,0.0,0.0);
	}
	Quaternion::Quaternion(const Quaternion& oQuaternion)
	{
		*this = oQuaternion;
	}
	Quaternion::Quaternion(const double& dA,const double& dB,const double& dC,const double& dD)
	{
		Set(dA,dB,dC,dD);
	}
	Quaternion::Quaternion(const double& dScalar,const Vector& oVector)
	{
		Set(dScalar,oVector);
	}
	Quaternion::~Quaternion()
	{

	}
	void Quaternion::Set(const double& dA,const double& dB,const double& dC,const double& dD)
	{
		m_dScalar = dA;
		m_oVector.Set(dB,dC,dD);
	}
	void Quaternion::SetByAxisAndAngle(const double& dAngle,const Vector& oAxis)
	{
		m_dScalar = cos(0.5*dAngle);
		m_oVector = oAxis.GetDirection()*sin(0.5*dAngle);
	}
	void Quaternion::Set(const double& dScalar,const Vector& oVector)
	{
		m_dScalar = dScalar;
		m_oVector = oVector;
	}
	Quaternion& Quaternion::operator=(const Quaternion& oQuaternion)
	{
		m_dScalar = oQuaternion.m_dScalar;
		m_oVector = oQuaternion.m_oVector;
		return *this;
	}
	Quaternion Quaternion::operator+(const Quaternion& oQuaternion) const
	{
		Quaternion oResult;
		oResult.m_dScalar = m_dScalar + oQuaternion.m_dScalar;
		oResult.m_oVector = m_oVector + oQuaternion.m_oVector;
		return oResult;
	}
	Quaternion Quaternion::operator-(const Quaternion& oQuaternion) const
	{
		Quaternion oResult;
		oResult.m_dScalar = m_dScalar - oQuaternion.m_dScalar;
		oResult.m_oVector = m_oVector - oQuaternion.m_oVector;
		return oResult;
	}
	Quaternion Quaternion::operator*(const Quaternion& oQuaternion) const
	{
		Quaternion oResult;
		double dA1 = m_dScalar;
		double dB1 = m_oVector.GetX();
		double dC1 = m_oVector.GetY();
		double dD1 = m_oVector.GetZ();
		double dA2 = oQuaternion.m_dScalar;
		double dB2 = oQuaternion.m_oVector.GetX();
		double dC2 = oQuaternion.m_oVector.GetY();
		double dD2 = oQuaternion.m_oVector.GetZ();
		oResult.m_dScalar = dA1*dA2 - dB1*dB2 - dC1*dC2 - dD1*dD2;
		oResult.m_oVector.SetX(dA1*dB2 + dA2*dB1 + dC1*dD2 - dC2*dD1);
		oResult.m_oVector.SetY(dA1*dC2 + dA2*dC1 + dD1*dB2 - dD2*dB1);
		oResult.m_oVector.SetZ(dA1*dD2 + dA2*dD1 + dB1*dC2 - dB2*dC1);
		return oResult;
	}
	Quaternion Quaternion::operator*(const double& dFactor) const
	{
		Quaternion oResult;
		oResult.m_dScalar = m_dScalar*dFactor;
		oResult.m_oVector = m_oVector*dFactor;
		return oResult;
	}
	Quaternion Quaternion::Conjugate() const
	{
		Quaternion oResult;
		oResult.m_dScalar = m_dScalar;
		oResult.m_oVector = m_oVector*(-1.0);
		return oResult;
	}
	double Quaternion::GetScalar() const
	{
		return m_dScalar;
	}
	Vector Quaternion::GetVector() const
	{
		return m_oVector;
	}
	double Quaternion::GetAngle() const
	{
		return 2.0*acos(m_dScalar);
	}
	Vector Quaternion::GetAxis() const
	{
		Vector oAxis = m_oVector;
		oAxis.Normalize();
		return oAxis;
	}
	void Quaternion::Normalize()
	{
		m_oVector.Normalize();
	}
}

