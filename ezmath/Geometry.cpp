// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#include "Geometry.h"
#include "CartesianOrthogonalCoordinateSystem.h"

namespace GeometrySystem
{
	Geometry::~Geometry()
	{

	}
	CartesianOrthogonalCoordinateSystem* Geometry::GetSystem()
	{
		return &m_oSystem;
	}
	void Geometry::SetSystem(const CartesianOrthogonalCoordinateSystem& oSystem)
	{
		m_oSystem = oSystem;
	}
}



