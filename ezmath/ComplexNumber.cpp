#include "ComplexNumber.h"
#include "cmath"
#include "stdio.h"

namespace EZ
{
	ComplexNumber::ComplexNumber()
	{
		Initialize();
	}
	ComplexNumber::ComplexNumber(const ComplexNumber& oNumber)
	{
		*this = oNumber;
	}
	ComplexNumber::ComplexNumber(const double& dReal,const double& dImaginary)
	{
		Set(dReal,dImaginary);
	}
	ComplexNumber::~ComplexNumber()
	{
	
	}
	ComplexNumber& ComplexNumber::operator=(const ComplexNumber& oNumber)
	{
		m_dReal = oNumber.m_dReal;
		m_dImaginary = oNumber.m_dImaginary;
		return *this;
	}
	void ComplexNumber::Reset()
	{
		Initialize();
	}
	void ComplexNumber::Set(const double& dReal,const double& dImaginary)
	{
		m_dReal = dReal;
		m_dImaginary = dImaginary;
	}
	void ComplexNumber::SetAmplitudeAndPhase(const double& dAmplitude,const double& dPhase)
	{
		m_dReal = dAmplitude*cos(dPhase);
		m_dImaginary = dAmplitude*sin(dPhase);
	}
	void ComplexNumber::SetReal(const double& dValue)
	{
		m_dReal = dValue;
	}
	void ComplexNumber::SetImaginary(const double& dValue)
	{
		m_dImaginary = dValue;
	}
	double ComplexNumber::GetReal() const
	{
		return m_dReal;
	}
	double ComplexNumber::GetImaginary() const
	{
		return m_dImaginary;
	}
	ComplexNumber ComplexNumber::operator+(const ComplexNumber& oNumber) const
	{
		return ComplexNumber(m_dReal + oNumber.m_dReal,m_dImaginary + oNumber.m_dImaginary);
	}
	ComplexNumber ComplexNumber::operator-(const ComplexNumber& oNumber) const
	{
		return ComplexNumber(m_dReal - oNumber.m_dReal,m_dImaginary - oNumber.m_dImaginary);
	}
	ComplexNumber ComplexNumber::operator*(const ComplexNumber& oNumber) const
	{
		ComplexNumber oResult;
		oResult.SetReal(m_dReal*oNumber.m_dReal - m_dImaginary*oNumber.m_dImaginary);
		oResult.SetImaginary(m_dReal*oNumber.m_dImaginary + m_dImaginary*oNumber.m_dReal);
		return oResult;
	}
	ComplexNumber ComplexNumber::operator*(const double& dValue) const
	{
		return ComplexNumber(m_dReal*dValue,m_dImaginary*dValue);
	}
	ComplexNumber ComplexNumber::operator/(const ComplexNumber& oNumber) const
	{
		double dDenominator = oNumber.m_dReal*oNumber.m_dReal + oNumber.m_dImaginary*oNumber.m_dImaginary;
		ComplexNumber oResult = (*this*oNumber.Conjugate())/dDenominator;
		return oResult;
	}
	ComplexNumber ComplexNumber::operator/(const double& dValue) const
	{
		return ComplexNumber(m_dReal/dValue,m_dImaginary/dValue);
	}
	ComplexNumber ComplexNumber::Conjugate() const
	{
		return ComplexNumber(m_dReal,-m_dImaginary);
	}
	string ComplexNumber::ToString() const
	{
		char cString[128];
		sprintf(cString,"%lf + %lf i",m_dReal,m_dImaginary);
		return string(cString);
	}
	double ComplexNumber::GetAmplitude() const
	{
		return sqrt(GetSquareAmplitude());
	}
	double ComplexNumber::GetSquareAmplitude() const
	{
		return (m_dReal*m_dReal + m_dImaginary*m_dImaginary);
	}
	void ComplexNumber::Initialize()
	{
		m_dReal = 0.0;
		m_dImaginary = 0.0;
	}
}



