#ifndef POLYHEDRON_H_
#define POLYHEDRON_H_

#include "PolyhedronFace.h"
#include "CircularLinkedList.h"
#include "AxisAlignedBoundingBox.h"
#include "Curve.h"

#define FRESH_NODE				1
#define PROCESSED_NODE			2
#define INTERNAL_EDGE			3
#define BOUNDARY_EDGE			4
#define INVISIBLE_EDGE			5
#define VISIBLE_FACE			6
#define INVISIBLE_FACE			7


namespace GeometrySystem
{
	class Polyhedron
	{
	public:
		Polyhedron();
		Polyhedron(const Polyhedron& oPolyhedron);
		~Polyhedron();
		Polyhedron& operator=(const Polyhedron& oPolyhedron);
		void Reset();
		void AddNode(GenericNode* poNode);
		void AddEdge(Edge* poEdge);
		void AddFace(PolyhedronFace* poFace);
		CircularLinkedList< GenericNode* >* GetNodes();
		list< Edge* >* GetEdges();
		list< PolyhedronFace* >* GetFaces();
		string ToString();
		bool CreateAsHullFromPoints(list< Point >* ploPoints);
		void WriteNodes(FILE* fpFile);
		void WriteParaview(const string& sFileName);
		bool IsPointInside(const Point& oPoint) const;
		void SetID(const unsigned int& iID);
		unsigned int GetID() const;
		double GetDistance(const Point& oPoint) const;
		double GetSignedDistance(const Point& oPoint) const;
		Point GetNearestPoint(const Point& oPoint,double& dDistance) const;
		bool GetNearestPointOnPlane(const Point& oPoint,const Plane& oPlane,Point& oNearestPoint,double& dDistance) const;
		bool GetNearestPointOnLine(const Point& oPoint,const Line& oLine,Point& oNearestPoint,double& dDistance) const;
		bool DoesSegmentPierce(const Point& oStartPoint,const Point& oEndPoint) const;
		const AxisAlignedBoundingBox* GetBox() const;
		Curve* GetIntersectionCurve(const Plane& oPlane) const;
		bool IsIntersecting(Polyhedron& oPolyhedron);
		bool IsIntersecting(const Plane& oPlane) const;
		bool Contains(Polyhedron& oPolyhedron);
		double GetVolume() const;
		Polyhedron* PlaneCut(const Plane& oPlane,const bool& bPerturb = false,const double& dPerturbationDistance = 0.0);
		Point GenerateInternalPoint() const;
		int ClassifyPlane(const Plane& oPlane) const;
		void UpdateIDs();
		
	private:

	protected:
		void Initialize();
		bool CreateDihedron();
		bool ConstructHull();
		bool ProcessHullPoint(GenericNode* poPoint);
		void CleanUp();
		void CleanUpNodes();
		
		CircularLinkedList< GenericNode* > m_lpoNodes;
		list< Edge* > m_lpoEdges;
		list< PolyhedronFace* > m_lpoFaces;
		unsigned int m_iLastNodeID;
		unsigned int m_iLastEdgeID;
		unsigned int m_iLastFaceID;
		bool m_bIsConvex;
		bool m_bHullConstructionError;
		unsigned int m_iID;
		AxisAlignedBoundingBox m_oBox;
	};
}

#endif

