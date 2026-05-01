// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#include "CartesianOrthogonalCoordinateSystem.h"
#include "math.h"

using namespace EZ;

namespace GeometrySystem
{
	CartesianOrthogonalCoordinateSystem::CartesianOrthogonalCoordinateSystem()
	{
		Reset();
	}
	CartesianOrthogonalCoordinateSystem::CartesianOrthogonalCoordinateSystem(const CartesianOrthogonalCoordinateSystem& oSystem)
	{
		*this = oSystem;
	}
	CartesianOrthogonalCoordinateSystem::~CartesianOrthogonalCoordinateSystem()
	{

	}
	CartesianOrthogonalCoordinateSystem& CartesianOrthogonalCoordinateSystem::operator=(const CartesianOrthogonalCoordinateSystem& oSystem)
	{
		SetOrigin(oSystem.m_oOrigin);
		SetXY(oSystem.GetX(),oSystem.GetY());
		return *this;
	}
	void CartesianOrthogonalCoordinateSystem::SetOrigin(const Point& oPoint)
	{
		m_oOrigin = oPoint;
	}
	void CartesianOrthogonalCoordinateSystem::SetXY(const Vector& oX,const Vector& oY)
	{
		m_oX = oX;
		m_oY = oY;
		m_oX.Normalize();
		m_oY.Normalize();
		m_oZ = oX^oY;
		m_oZ.Normalize();
		m_oRotationMatrix = GetRotationMatrix();
	}
	void CartesianOrthogonalCoordinateSystem::SetYZ(const Vector& oY,const Vector& oZ)
	{
		m_oY = oY;
		m_oZ = oZ;
		m_oY.Normalize();
		m_oZ.Normalize();
		m_oX = oY^oZ;
		m_oX.Normalize();
		m_oRotationMatrix = GetRotationMatrix();
	}
	void CartesianOrthogonalCoordinateSystem::SetXZ(const Vector& oX,const Vector& oZ)
	{
		m_oX = oX;
		m_oZ = oZ;
		m_oX.Normalize();
		m_oZ.Normalize();
		m_oY = oZ^oX;
		m_oY.Normalize();
		m_oRotationMatrix = GetRotationMatrix();
	}
	Point CartesianOrthogonalCoordinateSystem::GetOrigin() const
	{
		return m_oOrigin;
	}
	Vector CartesianOrthogonalCoordinateSystem::GetX() const
	{
		return m_oX;
	}
	Vector CartesianOrthogonalCoordinateSystem::GetY() const
	{
		return m_oY;
	}
	Vector CartesianOrthogonalCoordinateSystem::GetZ() const
	{
		return m_oZ;
	}
	void CartesianOrthogonalCoordinateSystem::Reset()
	{
		m_oOrigin.Set(0.0,0.0,0.0);
		m_oX.Set(1.0,0.0,0.0);
		m_oY.Set(0.0,1.0,0.0);
		m_oZ.Set(0.0,0.0,1.0);
		m_oRotationMatrix.SetSize(3,3);
		m_oRotationMatrix.Set(1,1,1.0);
		m_oRotationMatrix.Set(2,2,1.0);
		m_oRotationMatrix.Set(3,3,1.0);
	}
	Point CartesianOrthogonalCoordinateSystem::GetInGlobalCoordinates(const Point& oPoint) const
	{
		//Matrix oT = GetRotationMatrix();
		Matrix oP(3,1);
		oP.Set(1,1,oPoint.GetX());
		oP.Set(2,1,oPoint.GetY());
		oP.Set(3,1,oPoint.GetZ());
		//oP = oT*oP;
		oP = m_oRotationMatrix*oP;
		oP.Set(1,1,oP.Get(1,1) + m_oOrigin.GetX());
		oP.Set(2,1,oP.Get(2,1) + m_oOrigin.GetY());
		oP.Set(3,1,oP.Get(3,1) + m_oOrigin.GetZ());
		double dX = oP.Get(1,1);
		double dY = oP.Get(2,1);
		double dZ = oP.Get(3,1);
		return Point(dX,dY,dZ);
	}
	Point CartesianOrthogonalCoordinateSystem::GetInLocalCoordinates(const Point& oPoint) const
	{
		//Matrix oT = GetInverseRotationMatrix();
		Matrix oP(3,1);
		oP.Set(1,1,oPoint.GetX() - m_oOrigin.GetX());
		oP.Set(2,1,oPoint.GetY() - m_oOrigin.GetY());
		oP.Set(3,1,oPoint.GetZ() - m_oOrigin.GetZ());
		//oP = oT*oP;
		oP = m_oRotationMatrix.GetTranspose()*oP;
		double dX = oP.Get(1,1);
		double dY = oP.Get(2,1);
		double dZ = oP.Get(3,1);
		return Point(dX,dY,dZ);
	}
	Vector CartesianOrthogonalCoordinateSystem::GetInGlobalCoordinates(const Vector& oVector) const
	{
		//Matrix oT = GetRotationMatrix();
		Matrix oP(3,1);
		oP.Set(1,1,oVector.GetX());
		oP.Set(2,1,oVector.GetY());
		oP.Set(3,1,oVector.GetZ());
		//oP = oT*oP;
		oP = m_oRotationMatrix*oP;
		double dX = oP.Get(1,1);
		double dY = oP.Get(2,1);
		double dZ = oP.Get(3,1);
		return Vector(dX,dY,dZ);
	}
	Vector CartesianOrthogonalCoordinateSystem::GetInLocalCoordinates(const Vector& oVector) const
	{
		//Matrix oT = GetInverseRotationMatrix();
		Matrix oP(3,1);
		oP.Set(1,1,oVector.GetX());
		oP.Set(2,1,oVector.GetY());
		oP.Set(3,1,oVector.GetZ());
		//oP = oT*oP;
		oP = m_oRotationMatrix.GetTranspose()*oP;
		double dX = oP.Get(1,1);
		double dY = oP.Get(2,1);
		double dZ = oP.Get(3,1);
		return Vector(dX,dY,dZ);
	}
	Matrix CartesianOrthogonalCoordinateSystem::GetInGlobalCoordinates(const Matrix& oMatrix) const
	{
		Matrix oP = m_oRotationMatrix*oMatrix*m_oRotationMatrix.GetTranspose();
		return oP;
	}
	Matrix CartesianOrthogonalCoordinateSystem::GetInLocalCoordinates(const Matrix& oMatrix) const
	{
		Matrix oP = m_oRotationMatrix.GetTranspose()*oMatrix*m_oRotationMatrix;
		return oP;
	}
	Matrix CartesianOrthogonalCoordinateSystem::GetRotationMatrix() const
	{
		Matrix oT(3,3);
		oT.Set(1,1,cos(m_oX.GetXAngle()));
		oT.Set(2,1,cos(m_oX.GetYAngle()));
		oT.Set(3,1,cos(m_oX.GetZAngle()));

		oT.Set(1,2,cos(m_oY.GetXAngle()));
		oT.Set(2,2,cos(m_oY.GetYAngle()));
		oT.Set(3,2,cos(m_oY.GetZAngle()));

		oT.Set(1,3,cos(m_oZ.GetXAngle()));
		oT.Set(2,3,cos(m_oZ.GetYAngle()));
		oT.Set(3,3,cos(m_oZ.GetZAngle()));

		return oT;
	}
	Matrix CartesianOrthogonalCoordinateSystem::GetInverseRotationMatrix() const
	{
		return GetRotationMatrix().GetTranspose();
	}
	void CartesianOrthogonalCoordinateSystem::Move(const Vector& oStep)
	{
		m_oOrigin.Shift(oStep.GetX(),oStep.GetY(),oStep.GetZ());
	}
	void CartesianOrthogonalCoordinateSystem::Move(const double& dX,const double& dY,const double& dZ)
	{
		m_oOrigin.Shift(dX,dY,dZ);
	}
}

