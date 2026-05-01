#ifndef COMPLEXNUMBER_H_
#define COMPLEXNUMBER_H_

#include "string"

using namespace std;

namespace EZ
{
	class ComplexNumber
	{
	public:
		ComplexNumber();
		ComplexNumber(const ComplexNumber& oNumber);
		ComplexNumber(const double& dReal,const double& dImaginary);
		~ComplexNumber();
		ComplexNumber& operator=(const ComplexNumber& oNumber);
		void Reset();
		void Set(const double& dReal,const double& dImaginary);
		void SetAmplitudeAndPhase(const double& dAmplitude,const double& dPhase);
		void SetReal(const double& dValue);
		void SetImaginary(const double& dValue);
		double GetReal() const;
		double GetImaginary() const;
		ComplexNumber operator+(const ComplexNumber& oNumber) const;
		ComplexNumber operator-(const ComplexNumber& oNumber) const;
		ComplexNumber operator*(const ComplexNumber& oNumber) const;
		ComplexNumber operator*(const double& dValue) const;
		ComplexNumber operator/(const ComplexNumber& oNumber) const;
		ComplexNumber operator/(const double& dValue) const;
		ComplexNumber Conjugate() const;
		string ToString() const;
		double GetAmplitude() const;
		double GetSquareAmplitude() const;
		
	private:
	
	protected:
		void Initialize();
		double m_dReal;
		double m_dImaginary;
	};
}

#endif

