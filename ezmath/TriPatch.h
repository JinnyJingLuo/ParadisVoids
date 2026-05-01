// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef TRIPATCH_H_
#define TRIPATCH_H_

#include "GenericNode.h"
#include "Plane.h"
#include "list"
#include "Patch.h"
#include "Segment.h"

using namespace EZ;
using namespace std;

namespace GeometrySystem
{
	class TriPatch : public Patch
	{
	public:
		TriPatch();
		TriPatch(const TriPatch& oPatch);
		TriPatch(GenericNode* poPoint1,GenericNode* poPoint2,GenericNode* poPoint3);
		TriPatch(GenericNode* poPoint1,GenericNode* poPoint2,GenericNode* poPoint3,const bool& bIsConstrained);
		virtual ~TriPatch();
		virtual TriPatch& operator=(const TriPatch& oPatch);
		virtual void Reset();
		void Set(GenericNode* poPoint1,GenericNode* poPoint2,GenericNode* poPoint3);
		GenericNode* GetPoint1() const;
		GenericNode* GetPoint2() const;
		GenericNode* GetPoint3() const;
		Vector GetNormal() const;
		bool GetLineIntersection(const Point& oLineStart,const Point& oLineEnd,Point& oIntersectionPoint,double& dIntersectionParameter) const;
		bool IsIntersecting(const Plane& oPlane) const;
		bool IsPointInTriangle(const Point& oPoint) const;
		double GetMinSideLength() const;
		void Refine(list<GenericNode*>& lpoPoints,list<TriPatch*>& lpoTriangles,const double& dTolerance) const;
        vector<Patch*> GenerateTriangulation() const;
        Point GetNearestPoint(const Point& oPoint,double& dDistance) const;
        bool GetNearestPointOnPlane(const Point& oPoint,const Plane& oPlane,Point& oNearestPoint,double& dDistance) const;
        bool GetNearestPointOnLine(const Point& oPoint,const Line& oLine,Point& oNearestPoint,double& dDistance) const;
        void GeneratePoints(const unsigned int& iResolution,list<Point*>* plpoPoints) const;
        void SetConstraint(const bool& bIsConstrained = true);
        bool IsConstrained() const;
        double GetNormalDistance(const Point& oPoint) const;
        double GetArea() const;
        void UpdateProperties();
        bool DoesSegmentIntersect(const Point& oStartPoint,const Point& oEndPoint);
        bool DoesPlaneCut(const Plane& oPlane) const;
        bool GetPlaneIntersection(const Plane& oPlane,Point& oPoint1,Point& oPoint2) const;
		Point GetCentroid() const;
		Point GetCircumcenter() const;
		
	private:

	protected:
		virtual void Initialize();
		GenericNode* m_poPoint1;
		GenericNode* m_poPoint2;
		GenericNode* m_poPoint3;
		bool m_bIsConstrained;
		double m_dTolerance;
		double m_dArea;
		Vector m_oNormal;
	};
}

#endif

