// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef AXISALIGNEDBOUNDINGBOX_H_
#define AXISALIGNEDBOUNDINGBOX_H_

#include "Point.h"
#include "list"
#include "vector"
#include "Geometry.h"

using namespace EZ;
using namespace std;

namespace GeometrySystem
{
	class AxisAlignedBoundingBox : public Geometry
	{
	public:
		AxisAlignedBoundingBox();
		AxisAlignedBoundingBox(const AxisAlignedBoundingBox& oBox);
		~AxisAlignedBoundingBox();
		AxisAlignedBoundingBox& operator=(const AxisAlignedBoundingBox& oBox);
		void Reset();
		void SetXMin(const double& dValue);
		void SetXMax(const double& dValue);
		void SetYMin(const double& dValue);
		void SetYMax(const double& dValue);
		void SetZMin(const double& dValue);
		void SetZMax(const double& dValue);
		double GetXMin() const;
		double GetXMax() const;
		double GetYMin() const;
		double GetYMax() const;
		double GetZMin() const;
		double GetZMax() const;
		double GetMinimumDimension() const;
		double GetMaximumDimension() const;
		list<AxisAlignedBoundingBox*> UniformPartition(const unsigned int& iXSlicesCount,const unsigned int& iYSlicesCount,const unsigned int& iZSlicesCount) const;
		list<AxisAlignedBoundingBox*> NonUniformPartition(const unsigned int& iXSlicesCount,const unsigned int& iYSlicesCount,const unsigned int& iZSlicesCount) const;
		list<AxisAlignedBoundingBox*> PointBasedPartition(list<Point>* ploPoints,const unsigned int& iXSlicesCount,const unsigned int& iYSlicesCount,const unsigned int& iZSlicesCount) const;
		void UniformXSlice(const unsigned int& iSlicesCount,list<AxisAlignedBoundingBox*>* plpoBoxes) const;
		void UniformYSlice(const unsigned int& iSlicesCount,list<AxisAlignedBoundingBox*>* plpoBoxes) const;
		void UniformZSlice(const unsigned int& iSlicesCount,list<AxisAlignedBoundingBox*>* plpoBoxes) const;
		void NonUniformXSlice(const unsigned int& iSlicesCount,list<AxisAlignedBoundingBox*>* plpoBoxes) const;
		void NonUniformYSlice(const unsigned int& iSlicesCount,list<AxisAlignedBoundingBox*>* plpoBoxes) const;
		void NonUniformZSlice(const unsigned int& iSlicesCount,list<AxisAlignedBoundingBox*>* plpoBoxes) const;
		void PointBasedXSlice(list<Point>* ploPoints,const unsigned int& iSlicesCount,list<AxisAlignedBoundingBox*>* plpoBoxes) const;
		void PointBasedYSlice(list<Point>* ploPoints,const unsigned int& iSlicesCount,list<AxisAlignedBoundingBox*>* plpoBoxes) const;
		void PointBasedZSlice(list<Point>* ploPoints,const unsigned int& iSlicesCount,list<AxisAlignedBoundingBox*>* plpoBoxes) const;
		double GetVolume() const;
		bool IsPointInside(const Point& oPoint,const double& dToleranceFactor = 1.0E-6) const;
		Point GenerateRandomPoint(const double& dMargin = 0.0) const;
		Point GenerateRandomExternalPoint(const double& dSpacingFactor = 2.0) const;
		void CenterAt(const Point& oPoint);
		void ExpandToContain(const Point& oPoint);
		int GetMaximumDimensionIndex() const;
		Point GetCenter() const;
		Geometry* Clone();
		string ToString() const;
		Point GetRelativePosition(const Point& oPoint) const;
		AxisAlignedBoundingBox GetSubVolume(const double& dVolumeFraction) const;

		
	private:
	
	protected:
		virtual void Initialize();
		double m_dXMin;
		double m_dXMax;
		double m_dYMin;
		double m_dYMax;
		double m_dZMin;
		double m_dZMax;
	};
}

#endif


