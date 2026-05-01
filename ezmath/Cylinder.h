// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef _CYLINDER_H_
#define _CYLINDER_H_

#include "Geometry.h"

using namespace EZ;

namespace GeometrySystem
{
	class Cylinder : public Geometry
	{
	public:
		Cylinder();
		~Cylinder();
		Cylinder(const Cylinder& oCylinder);
		Cylinder& operator=(const Cylinder& oCylinder);
		void Reset();
		bool IsOnSurface(const Point& oPoint) const;
		double GetRadius() const;
		double GetLength() const;
		void SetRadius(const double& dRadius);
		void SetLength(const double&  dLength);
		void SetResolution(const unsigned int& iRadialResolution,const unsigned int& iCircumferentialResolution,const unsigned int& iAxialResolution);
		unsigned int GetRadialResolution() const;
		unsigned int GetCircumferentialResolution() const;
		unsigned int GetAxialResolution() const;
		bool IsOnLowerFace(const Point& oPoint) const;
		bool IsOnUpperFace(const Point& oPoint) const;
		bool IsOnLateralFace(const Point& oPoint) const;
		bool IsCenterOfLowerFace(const Point& oPoint) const;
		bool IsCenterOfUpperFace(const Point& oPoint) const;
		virtual Geometry* Clone();
		double GetVolume() const;
		bool IsPointInside(const Point& oPoint,const double& dTolerance = 1.0E-6) const;
		
	private:

	protected:
		double m_dRadius;
		double m_dLength;
		unsigned int m_iRadialResolution;
		unsigned int m_iCircumferentialResolution;
		unsigned int m_iAxialResolution;
	};
}

#endif



