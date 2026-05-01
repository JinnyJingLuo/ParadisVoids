// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#include "TriPatch.h"
#include "Vector.h"
#include "iostream"
#include "float.h"
#include "math.h"
#include "vector"
#include "Segment.h"
#include "Matrix.h"

using namespace EZ;
using namespace std;

namespace GeometrySystem
{
	TriPatch::TriPatch()
	{
		Initialize();
	}
	TriPatch::TriPatch(const TriPatch& oPatch)
	{
		*this = oPatch;
	}
	TriPatch::TriPatch(GenericNode* poPoint1,GenericNode* poPoint2,GenericNode* poPoint3)
	{
		Set(poPoint1,poPoint2,poPoint3);
	}
   	TriPatch::TriPatch(GenericNode* poPoint1,GenericNode* poPoint2,GenericNode* poPoint3,const bool& bIsConstrained)
   	{
		Set(poPoint1,poPoint2,poPoint3);
		m_bIsConstrained = bIsConstrained;
   	}
	TriPatch::~TriPatch()
	{
		Reset();
	}
	TriPatch& TriPatch::operator=(const TriPatch& oPatch)
	{
		Set(oPatch.m_poPoint1,oPatch.m_poPoint2,oPatch.m_poPoint3);
		m_bIsConstrained = oPatch.m_bIsConstrained;
		return *this;
	}
	void TriPatch::Reset()
	{
		Initialize();
	}
	void TriPatch::Set(GenericNode* poPoint1,GenericNode* poPoint2,GenericNode* poPoint3)
	{
		m_poPoint1 = poPoint1;
		m_poPoint2 = poPoint2;
		m_poPoint3 = poPoint3;
		UpdateProperties();
	}
	GenericNode* TriPatch::GetPoint1() const
	{
		return m_poPoint1;
	}
	GenericNode* TriPatch::GetPoint2() const
	{
		return m_poPoint2;
	}
	GenericNode* TriPatch::GetPoint3() const
	{
		return m_poPoint3;
	}
	Vector TriPatch::GetNormal() const
	{
		return m_oNormal;
	}
	bool TriPatch::GetLineIntersection(const Point& oLineStart,const Point& oLineEnd,Point& oIntersectionPoint,double& dIntersectionParameter) const
	{
		dIntersectionParameter = 0.0;
		Vector oV(oLineStart,oLineEnd);
		Vector oNormal = m_oNormal*-1.0;
		
		double dDen = oV*oNormal;
		if(fabs(dDen) < m_dTolerance)
		{
			return false;
		}
		Vector oTest(oLineStart,*m_poPoint1);
		dIntersectionParameter = oNormal*oTest/dDen;
		if(dIntersectionParameter < -m_dTolerance || dIntersectionParameter > 1.0 + m_dTolerance)
		{
			return false;
		}
		oIntersectionPoint = Point(oLineStart + oV*dIntersectionParameter);
		return IsPointInTriangle(oIntersectionPoint);
	}
	bool TriPatch::IsIntersecting(const Plane& oPlane) const
	{
		int	iPoint1Classification = oPlane.ClassifyPoint(*m_poPoint1);
		int	iPoint2Classification = oPlane.ClassifyPoint(*m_poPoint2);
		int	iPoint3Classification = oPlane.ClassifyPoint(*m_poPoint3);
		if(iPoint1Classification == 1 && iPoint2Classification == 1 && iPoint3Classification == 1)
		{
			return false;
		}
		if(iPoint1Classification == -1 && iPoint2Classification == -1 && iPoint3Classification == -1)
		{
			return false;
		}
		return true;
	}
	double TriPatch::GetMinSideLength() const
	{
		double dLength1 = m_poPoint1->Distance(*m_poPoint2);
		double dLength2 = m_poPoint2->Distance(*m_poPoint3);
		double dLength3 = m_poPoint3->Distance(*m_poPoint1);
		return min(min(dLength1,dLength2),dLength3);
	}
	void TriPatch::Refine(list<GenericNode*>& lpoPoints,list<TriPatch*>& lpoTriangles,const double& dTolerance) const
	{
		vector<GenericNode> voMidPoints;
		voMidPoints.resize(3);
		voMidPoints[0] = (*m_poPoint1 + *m_poPoint2)*0.5;
		voMidPoints[1] = (*m_poPoint2 + *m_poPoint3)*0.5;
		voMidPoints[2] = (*m_poPoint3 + *m_poPoint1)*0.5;
		
		vector<GenericNode*> vpoTempPoints;
		vpoTempPoints.resize(3);
		unsigned int i = 0;
		for(i = 0; i < 3 ; i++)
		{
			vpoTempPoints[i] = NULL;
		}
		list<GenericNode*>::iterator liPoints;
		for(liPoints = lpoPoints.begin() ; liPoints != lpoPoints.end() ; liPoints++)
		{
			for(i = 0; i < 3 ; i++)
			{
				if(vpoTempPoints[i] == NULL)
				{
					if(voMidPoints[i].Distance(*(*liPoints)) < dTolerance)
					{
						vpoTempPoints[i] = *liPoints;
						break;
					}
				}
			}

			if(vpoTempPoints[0] != NULL && vpoTempPoints[1] != NULL && vpoTempPoints[2] != NULL)
			{
				break;
			}
		}
		for(i = 0; i < 3 ; i++)
		{
			if(vpoTempPoints[i] == NULL)
			{
				lpoPoints.push_back(new GenericNode(voMidPoints[i]));
				vpoTempPoints[i] = lpoPoints.back();
			}
		}
		lpoTriangles.push_back(new TriPatch(m_poPoint1,vpoTempPoints[0],vpoTempPoints[2]));
		lpoTriangles.push_back(new TriPatch(m_poPoint2,vpoTempPoints[1],vpoTempPoints[0]));
		lpoTriangles.push_back(new TriPatch(m_poPoint3,vpoTempPoints[2],vpoTempPoints[1]));
		lpoTriangles.push_back(new TriPatch(vpoTempPoints[0],vpoTempPoints[1],vpoTempPoints[2]));
	}
	bool TriPatch::IsPointInTriangle(const Point& oPoint) const
	{
		Vector oArm1(oPoint,*m_poPoint1);
		Vector oArm2(oPoint,*m_poPoint2);
		Vector oArm3(oPoint,*m_poPoint3);
		double dTotalArea = 0.5*((oArm1^oArm2).Length() + (oArm2^oArm3).Length() + (oArm3^oArm1).Length());
		if(fabs(dTotalArea - m_dArea) > m_dTolerance)
		{
			return false;
		}
		return true;
	}
	vector<Patch*> TriPatch::GenerateTriangulation() const
	{
		vector<Patch*> vpoTriangles;
		vpoTriangles.resize(1);
		vpoTriangles[0] = new TriPatch(m_poPoint1,m_poPoint2,m_poPoint3);
		return vpoTriangles;
	}
	Point TriPatch::GetNearestPoint(const Point& oPoint,double& dDistance) const
	{
		// get point projection on the triangle plane
    	Plane oTrianglePlane = Plane(m_oNormal,*m_poPoint1);
		Point oProjection = oTrianglePlane.GetPointProjection(oPoint);
		// see if the projection lies inside the triangle
		bool bIsInside = IsPointInTriangle(oProjection);
		Point oClosestPoint;
		double dT = 0.0;
		if(bIsInside)
		{
			oClosestPoint = oProjection;
		}
		else
		{
			Vector oV13 = Vector(*m_poPoint1,*m_poPoint3);
			Vector oV32 = Vector(*m_poPoint3,*m_poPoint2);
			Vector oV21 = Vector(*m_poPoint2,*m_poPoint1);
			Vector oN13 = m_oNormal^oV13;
			Vector oN32 = m_oNormal^oV32;
			Vector oN21 = m_oNormal^oV21;
			oN13.Normalize();
			oN32.Normalize();
			oN21.Normalize();
			Vector oV1p = Vector(*m_poPoint1,oProjection);
			Vector oV2p = Vector(*m_poPoint2,oProjection);
			Vector oV3p = Vector(*m_poPoint3,oProjection);
			bool bIsEdgeRegion = false;
			// check if it is in the edge 1 region
			if(oV13*oV1p >= 0)
			{
				if(oV13*oV3p <= 0)
				{
					dT = oN13*oV1p;
					if(dT >= 0)
					{
						bIsEdgeRegion = true;
						oClosestPoint = oProjection - oN13*dT;
					}
				}
			}
			// check if it is in the edge 2 region
			if(!bIsEdgeRegion)
			{
				if(oV32*oV3p >= 0)
				{
					if(oV32*oV2p <= 0)
					{
						dT = oN32*oV3p;
						if(dT >= 0)
						{
							bIsEdgeRegion = true;
							oClosestPoint = oProjection - oN32*dT;
						}
					}
				}
			}
			// check if it is in the edge 3 region
			if(!bIsEdgeRegion)
			{
				if(oV21*oV2p >= 0)
				{
					if(oV21*oV1p <= 0)
					{
						dT = oN21*oV2p;
						if(dT >= 0)
						{
							bIsEdgeRegion = true;
							oClosestPoint = oProjection - oN21*dT;
						}
					}
				}
			}
			// check to see if it is in a vertex region
			if(!bIsEdgeRegion)
			{
				dDistance = m_poPoint1->Distance(oProjection);
				oClosestPoint = *m_poPoint1;
				dT = m_poPoint2->Distance(oProjection);
				if(dT < dDistance)
				{
					dDistance = dT;
					oClosestPoint = *m_poPoint2;
				}
				dT = m_poPoint3->Distance(oProjection);
				if(dT < dDistance)
				{
					dDistance = dT;
					oClosestPoint = *m_poPoint3;
				}
			}
		}
		dDistance = oPoint.Distance(oClosestPoint);
		return oClosestPoint;
	}
	bool TriPatch::GetNearestPointOnPlane(const Point& oPoint,const Plane& oPlane,Point& oNearestPoint,double& dDistance) const
	{
		// initialize return values
		oNearestPoint.Set(0.0,0.0,0.0);
		dDistance = DBL_MAX;
		// if the plane is not intersecting the traingle, there is no way
		// to satisfy all the constraints, return false
		if(!IsIntersecting(oPlane))
		{
			return false;
		}
		// otherwise, see if the plane is coinciding with the triangle plane
		double dTemp = fabs(m_oNormal*oPlane.GetNormal());
		double dTolerance = 1E-6;
		// this is the product of normalized vectors, no need for relative tolerance
		if(fabs(dTemp - 1.0) < dTolerance)
		{
			// the planes are parallel, since the plane intersects the
			// triangle, then they are guaranteed to be coinciding
			// in that case, all of the points on and inside the triange boundary
			// satisfy the plane constraint, the problem reduces to finding the
			// closest point to the triangle without considering the plane
			// constraint
			oNearestPoint = GetNearestPoint(oPoint,dDistance);
			return true;
		}
		// the plane and the triangle do NOT coincide, get the intersection
		// line between the triangle plane and the constraint plane
		Plane oTrianglePlane(m_oNormal,*m_poPoint1);
		Line oIntersectionLine;
		if(!oPlane.GetIntersectionLine(oTrianglePlane,oIntersectionLine))
		{
			return false;
		}
		GenericNode oIntersectionPoint1(0.0,0.0,0.0);
		GenericNode oIntersectionPoint2(0.0,0.0,0.0);
		GenericNode oIntersectionPoint3(0.0,0.0,0.0);
		int iAvoidEdge = 0;
		if(!Segment(m_poPoint1,m_poPoint2).GetIntersectionPoint(oIntersectionLine,oIntersectionPoint1))
		{
			iAvoidEdge = 1;
		}
		if(!Segment(m_poPoint2,m_poPoint3).GetIntersectionPoint(oIntersectionLine,oIntersectionPoint2))
		{
			iAvoidEdge = 2;
		}
		if(!Segment(m_poPoint3,m_poPoint1).GetIntersectionPoint(oIntersectionLine,oIntersectionPoint3))
		{
			iAvoidEdge = 3;
		}
		// now, use two of the edges to get the closest point but make sure that we are avoiding an edge
		if(iAvoidEdge == 0)
		{
			// the line intersects all of the 3 edges, get the maximum length segment and use it
			if(oIntersectionPoint1.Distance(oIntersectionPoint2) > oIntersectionPoint1.Distance(oIntersectionPoint3))
			{
				iAvoidEdge = 3;
			}
			else
			{
				iAvoidEdge = 2;
			}
		}
		
		if(iAvoidEdge == 1)
		{
			oNearestPoint = Segment(&oIntersectionPoint2,&oIntersectionPoint3).GetClosestPoint(oPoint,dDistance);
			return true;
		}
		else if(iAvoidEdge == 2)
		{
			oNearestPoint = Segment(&oIntersectionPoint3,&oIntersectionPoint1).GetClosestPoint(oPoint,dDistance);
			return true;
		}
		else if(iAvoidEdge == 3)
		{
			oNearestPoint = Segment(&oIntersectionPoint1,&oIntersectionPoint2).GetClosestPoint(oPoint,dDistance);
			return true;
		}
		return false;
	}
	bool TriPatch::GetNearestPointOnLine(const Point& oPoint,const Line& oLine,Point& oNearestPoint,double& dDistance) const
	{
		Plane oTrianglePlane(m_oNormal,*m_poPoint1);
		if(!oTrianglePlane.GetLineIntersection(oLine,oNearestPoint))
		{
			return false;
		}
		if(!IsPointInTriangle(oNearestPoint))
		{
			return false;
		}
		dDistance = oPoint.Distance(oNearestPoint);
		return true;
	}
	void TriPatch::GeneratePoints(const unsigned int& iResolution,list<Point*>* plpoPoints) const
	{
		unsigned int iWorkingResolution = iResolution;
		if(iWorkingResolution < 2)
		{
			iWorkingResolution = 2;
		}
		double dEpsilon = 1.0/2.0/((double)(iWorkingResolution));
		double dFactor = 2*dEpsilon;
		unsigned int i = 0;
		unsigned int j = 0;
		double dXi = 0.0;
		double dEta = 0.0;
		double dX = 0.0;
		double dY = 0.0;
		double dZ = 0.0;
		for(i = 0 ; i < iWorkingResolution ; i++)
		{
			dXi = dEpsilon + dFactor*i;
			for(j = 0 ; j < iWorkingResolution ; j++)
			{
				dEta = dEpsilon + dFactor*j;
				if((dXi + dEta) > 1.0)
				{
					break;
				}
				dX = (1.0 - dXi - dEta)*m_poPoint1->GetX() + dXi*m_poPoint2->GetX() + dEta*m_poPoint3->GetX();
				dY = (1.0 - dXi - dEta)*m_poPoint1->GetY() + dXi*m_poPoint2->GetY() + dEta*m_poPoint3->GetY();
				dZ = (1.0 - dXi - dEta)*m_poPoint1->GetZ() + dXi*m_poPoint2->GetZ() + dEta*m_poPoint3->GetZ();
				plpoPoints->push_back(new Point(dX,dY,dZ));
			}
		}
	}
	double TriPatch::GetNormalDistance(const Point& oPoint) const
	{
		Vector oV(*m_poPoint1,oPoint);
		return fabs(oV*m_oNormal);
	}
	double TriPatch::GetArea() const
	{
		return m_dArea;
	}
	void TriPatch::SetConstraint(const bool& bIsConstrained)
    {
       m_bIsConstrained = bIsConstrained;
    }
    bool TriPatch::IsConstrained() const
    {
       return m_bIsConstrained;
    }
    void TriPatch::UpdateProperties()
    {
    	Vector oE1(*m_poPoint1,*m_poPoint2);
		Vector oE2(*m_poPoint1,*m_poPoint3);
		Vector oE3(*m_poPoint2,*m_poPoint3);
		m_dTolerance = 1E-6*min(min(oE1.Length(),oE2.Length()),oE3.Length());
		m_oNormal = oE1^oE2;
		m_dArea = 0.5*(m_oNormal.Length());
		m_oNormal.Normalize();
    }
    bool TriPatch::DoesSegmentIntersect(const Point& oStartPoint,const Point& oEndPoint)
    {
    	Point oIntersectionPoint;
    	double dIntersectionParameter = 0.0;
    	return GetLineIntersection(oStartPoint,oEndPoint,oIntersectionPoint,dIntersectionParameter);
    }
    bool TriPatch::DoesPlaneCut(const Plane& oPlane) const
    {
    	int iP1Class = oPlane.ClassifyPoint(*m_poPoint1);
    	int iP2Class = oPlane.ClassifyPoint(*m_poPoint2);
    	int iP3Class = oPlane.ClassifyPoint(*m_poPoint3);
    	if((iP1Class == iP2Class) && (iP1Class == iP3Class))		return false;
    	return true;
    }
    bool TriPatch::GetPlaneIntersection(const Plane& oPlane,Point& oPoint1,Point& oPoint2) const
    {
    	oPoint1.Reset();
    	oPoint2.Reset();
    	if(!DoesPlaneCut(oPlane))			return false;
    	Point oL1Intersection;
    	Point oL2Intersection;
    	Point oL3Intersection;
    	bool bL1Intersecting = oPlane.GetSegmentIntersection(*m_poPoint1,*m_poPoint2,oL1Intersection);
    	bool bL2Intersecting = oPlane.GetSegmentIntersection(*m_poPoint2,*m_poPoint3,oL2Intersection);
    	bool bL3Intersecting = oPlane.GetSegmentIntersection(*m_poPoint3,*m_poPoint1,oL3Intersection);
    	// make sure that not all 3 segments are intersecting or not intersecting
    	if(bL1Intersecting && bL2Intersecting && bL3Intersecting)				return false;
    	if((!bL1Intersecting) && (!bL2Intersecting) && (!bL3Intersecting))		return false;
    	if(bL1Intersecting && bL2Intersecting)
    	{
    		oPoint1 = oL1Intersection;
    		oPoint2 = oL2Intersection;
    	}
    	if(bL2Intersecting && bL3Intersecting)
    	{
    	    oPoint1 = oL2Intersection;
    		oPoint2 = oL3Intersection;
    	}
    	if(bL3Intersecting && bL1Intersecting)
    	{
    	    oPoint1 = oL3Intersection;
    		oPoint2 = oL1Intersection;
    	}
    	return true;
    }
    void TriPatch::Initialize()
    {
    	m_poPoint1 = NULL;
		m_poPoint2 = NULL;
		m_poPoint3 = NULL;
		m_dTolerance = 0.0;
		m_dArea = 0.0;
		m_oNormal.Set(0.0,0.0,0.0);
		m_bIsConstrained = false;
    }
    Point TriPatch::GetCentroid() const
    {
    	Point Centroid = *m_poPoint1 + *m_poPoint2 + *m_poPoint3;
    	Centroid = Centroid*(1.0/3.0);
    	return Centroid;
    }
    Point TriPatch::GetCircumcenter() const
    {
    	Vector oM = Vector(*m_poPoint1,*m_poPoint2)^m_oNormal;
    	Vector oN = Vector(*m_poPoint2,*m_poPoint3)^m_oNormal;
    	oM.Normalize();
    	oN.Normalize();
    	Point oP = (*m_poPoint1 + *m_poPoint2)*0.5;
    	Point oQ = (*m_poPoint2 + *m_poPoint3)*0.5;
    	// equate the x and y coordinates of the two lines, the z coordinate will be equal 
    	// automatically because the 2 lines lie in the plane of the triangle so that for a given
    	// x and y, there is only 1 possible z on that plane
    	// as long as this is a triangle, there is always a solution to the linear system of equations
    	Matrix oMat(2,2);
    	oMat.Set(1,1,oM.GetX());
    	oMat.Set(1,2,-oN.GetX());
    	oMat.Set(2,1,oM.GetY());
    	oMat.Set(2,2,-oN.GetY());
    	Matrix oRHS(2,1);
    	oRHS.Set(1,1,oQ.GetX() - oP.GetX());
    	oRHS.Set(2,1,oQ.GetY() - oP.GetY());
    	
    	Matrix oX = oMat.Solve(oRHS);
    	// get the point
    	Point oResult = oP + oM*oX.Get(1,1);
    	double dD1 = oResult.Distance(*m_poPoint1);
    	double dD2 = oResult.Distance(*m_poPoint2);
    	double dD3 = oResult.Distance(*m_poPoint3);
    	double dError = fabs(dD1 - dD2);
    	if(fabs(dD1 - dD3) > dError)		dError = fabs(dD1 - dD3);
    	if(fabs(dD2 - dD3) > dError)		dError = fabs(dD2 - dD3);
    	return oResult;
    }
}



