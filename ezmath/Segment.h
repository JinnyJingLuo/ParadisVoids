// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef SEGMENT_H_
#define SEGMENT_H_

#include "GenericNode.h"
#include "Vector.h"
#include "Line.h"
#include "AxisAlignedBoundingBox.h"

using namespace EZ;
using namespace std;

namespace GeometrySystem
{
	class Segment
	{
	public:
		Segment();
		Segment(const Segment& oSegment);
		Segment(GenericNode* poPoint1,GenericNode* poPoint2,bool bIsOwner = false);
		~Segment();
		Segment& operator=(const Segment& oSegment);
		void Reset();
		void Set(GenericNode* poPoint1,GenericNode* poPoint2,bool bIsOwner = false);
		GenericNode* GetPoint1() const;
		GenericNode* GetPoint2() const;
		Vector GetDirection() const;
		bool IsPointOnSegment(const Point& oPoint) const;
        Point GetClosestPoint(const Point& oPoint,double& dDistance) const;
        bool GetIntersectionPoint(const Line& oLine,Point& oPoint) const;
        bool GetIntersectionPoint(const Segment& oSegment,Point& oPoint) const;
        void Flip();
        bool GetIntersectionPoints(const AxisAlignedBoundingBox& oBox,Point& oStartPoint,Point& oEndPoint) const;
        
	private:

	protected:
		void Initialize();
		GenericNode* m_poPoint1;
		GenericNode* m_poPoint2;
		double m_dLength;
		Vector m_oDirection;
		static double m_dToleranceFactor;
		bool m_bIsOwner;
	};
}


#endif


