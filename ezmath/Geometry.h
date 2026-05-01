// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include "Point.h"
#include "CartesianOrthogonalCoordinateSystem.h"

using namespace EZ;

namespace GeometrySystem
{
	class CartesianOrthogonalCoordinateSystem;
	class Geometry
	{
	public:
		virtual ~Geometry();
		virtual Geometry* Clone() = 0;
		virtual double GetVolume() const = 0;
		virtual bool IsPointInside(const Point& oPoint,const double& dToleranceFactor = 1.0E-6) const = 0;
		CartesianOrthogonalCoordinateSystem* GetSystem();
		void SetSystem(const CartesianOrthogonalCoordinateSystem& oSystem);

	private:

	protected:
		CartesianOrthogonalCoordinateSystem m_oSystem;
	};
}

#endif


