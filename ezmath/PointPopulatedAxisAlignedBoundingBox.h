// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef POINTPOPULATEDAXISALIGNEDBOUNDINGBOX_H_
#define POINTPOPULATEDAXISALIGNEDBOUNDINGBOX_H_

#include "AxisAlignedBoundingBox.h"
#include "list"

using namespace EZ;
using namespace std;

namespace GeometrySystem
{
	class PointPopulatedAxisAlignedBoundingBox : public AxisAlignedBoundingBox
	{
	public:
		PointPopulatedAxisAlignedBoundingBox();
		PointPopulatedAxisAlignedBoundingBox(const AxisAlignedBoundingBox& oBox);
		PointPopulatedAxisAlignedBoundingBox(const PointPopulatedAxisAlignedBoundingBox& oBox);
		~PointPopulatedAxisAlignedBoundingBox();
		PointPopulatedAxisAlignedBoundingBox& operator=(const PointPopulatedAxisAlignedBoundingBox& oBox);
		list<Point*>* GetPoints();
		void AddPoint(Point* poPoint);
		void SetPoints(list<Point*>* plpoPoints);
		AxisAlignedBoundingBox GetBox() const;
		list<AxisAlignedBoundingBox*> PointBasedPartition(const unsigned int& iXSlicesCount,const unsigned int& iYSlicesCount,const unsigned int& iZSlicesCount);
		void PointBasedXSlice(const unsigned int& iSlicesCount,list<PointPopulatedAxisAlignedBoundingBox*>* plpoBoxes);
		void PointBasedYSlice(const unsigned int& iSlicesCount,list<PointPopulatedAxisAlignedBoundingBox*>* plpoBoxes);
		void PointBasedZSlice(const unsigned int& iSlicesCount,list<PointPopulatedAxisAlignedBoundingBox*>* plpoBoxes);
	private:
	
	protected:
		virtual void Initialize();
		list<Point*> m_lpoPoints;
	};
}

#endif


