// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef CARTESIANORTHOGONALCOORDINATESYSTEM_H_
#define CARTESIANORTHOGONALCOORDINATESYSTEM_H_

#include "Vector.h"
#include "Matrix.h"

using namespace EZ;

namespace GeometrySystem
{
	class CartesianOrthogonalCoordinateSystem
	{
	public:
		CartesianOrthogonalCoordinateSystem();
		CartesianOrthogonalCoordinateSystem(const CartesianOrthogonalCoordinateSystem& oSystem);
		~CartesianOrthogonalCoordinateSystem();
		CartesianOrthogonalCoordinateSystem& operator=(const CartesianOrthogonalCoordinateSystem& oSystem);
		void SetOrigin(const Point& oPoint);
		void SetXY(const Vector& oX,const Vector& oY);
		void SetYZ(const Vector& oY,const Vector& oZ);
		void SetXZ(const Vector& oX,const Vector& oZ);
		Point GetOrigin() const;
		Vector GetX() const;
		Vector GetY() const;
		Vector GetZ() const;
		void Reset();
		Point GetInGlobalCoordinates(const Point& oPoint) const;
		Point GetInLocalCoordinates(const Point& oPoint) const;
		Vector GetInGlobalCoordinates(const Vector& oVector) const;
		Vector GetInLocalCoordinates(const Vector& oVector) const;
		Matrix GetInGlobalCoordinates(const Matrix& oMatrix) const;
		Matrix GetInLocalCoordinates(const Matrix& oMatrix) const;
		void Move(const Vector& oStep);
		void Move(const double& dX,const double& dY,const double& dZ);
	private:

	protected:
		Matrix GetRotationMatrix() const;
		Matrix GetInverseRotationMatrix() const;
		Point m_oOrigin;
		Vector m_oX;
		Vector m_oY;
		Vector m_oZ;
		Matrix m_oRotationMatrix;
	};
}

#endif

