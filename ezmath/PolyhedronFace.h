#ifndef POLYHEDRONFACE_H_
#define POLYHEDRONFACE_H_

// disable the override intended warning
#pragma warning disable 1125

#include "Edge.h"
#include "Vector.h"
#include "TriPatch.h"
#include "Face.h"

using namespace EZ;

namespace GeometrySystem
{
	class PolyhedronFace : public Face, public TriPatch
	{
	public:
		PolyhedronFace();
		PolyhedronFace(const PolyhedronFace& oFace);
		virtual ~PolyhedronFace();
		virtual PolyhedronFace& operator=(const PolyhedronFace& oFace);
		virtual void Reset();
		void SetPoints(GenericNode* poPoint1,GenericNode* poPoint2,GenericNode* poPoint3);
		void SetEdges(Edge* poEdge1,Edge* poEdge2,Edge* poEdge3);
		bool IsPointOnPlane(const Point& oPoint);
		Edge* GetEdge1() const;
		Edge* GetEdge2() const;
		Edge* GetEdge3() const;
		bool IsVisibleFromPoint(const Point& oPoint) const;
		string ToString() const;
		bool IsEdgeAligned(Edge* poEdge) const;

	private:

	protected:
		virtual void Initialize();
		Edge* m_poEdge1;
		Edge* m_poEdge2;
		Edge* m_poEdge3;
	};
}

#endif



