// Paradis Processor Project
// Ahmed M. Hussein (mailto : am.hussin@gmail.com)
// June 2012

#ifndef QUATERNION_H_
#define QUATERNION_H_

#include "Vector.h"

namespace EZ
{
	class Quaternion
	{
	public:
		Quaternion();
		Quaternion(const Quaternion& oQuaternion);
		Quaternion(const double& dA,const double& dB,const double& dC,const double& dD);
		Quaternion(const double& dScalar,const Vector& oVector);
		~Quaternion();
		void Set(const double& dA,const double& dB,const double& dC,const double& dD);
		void Set(const double& dScalar,const Vector& oVector);
		void SetByAxisAndAngle(const double& dAngle,const Vector& oAxis);
		Quaternion& operator=(const Quaternion& oQuaternion);
		Quaternion operator+(const Quaternion& oQuaternion) const;
		Quaternion operator-(const Quaternion& oQuaternion) const;
		Quaternion operator*(const Quaternion& oQuaternion) const;
		Quaternion operator*(const double& dFactor) const;
		Quaternion Conjugate() const;
		double GetScalar() const;
		Vector GetVector() const;
		double GetAngle() const;
		Vector GetAxis() const;
		void Normalize();
	protected:

	private:
		double m_dScalar;
		Vector m_oVector;
	};
}

#endif

