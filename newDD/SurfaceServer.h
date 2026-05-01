// Ahmed M. Hussein

#ifndef SURFACESERVER_H_
#define SURFACESERVER_H_

#include "MainDataStructure.h"
#include "Vector.h"
#include "TriPatch.h"
#include "list"

using namespace std;
using namespace EZ;
using namespace GeometrySystem;

class SurfaceServer
{
public:
	static SurfaceServer* GetInstance();
	~SurfaceServer();
	void SetDataStructure(MainDataStructure* poDataStructure);
	void CheckNodes();
	void RemoveSurfaceArms();
	
private:

protected:
	static SurfaceServer* m_poInstance;
	SurfaceServer();
	void Initialize();
	void Reset();
	TriPatch* GetNearestTriangle(const Point& oPoint,Point& oNearestPoint) const;
	TriPatch* GetNearestTriangleOnPlane(const Plane& oPlane,const Point& oPoint,Point& oNearestPoint) const;
	TriPatch* GetNearestTriangleOnLine(const Line& oLine,const Point& oPoint,Point& oNearestPoint) const;
	bool IsPointInBox(const Point& oPoint);
	bool IsPointInside(const Point& oPoint);
	bool IsPointOnSurface(const Point& oPoint) const;
	double GetLeastSurfaceDistance(const Point& oPoint) const;
	void GenerateTriangulations();
	MainDataStructure* m_poDataStructure;
	list<GenericNode*> m_lpoPoints;
	list<TriPatch*> m_lpoTriangles;
	double m_dXMin;
	double m_dXMax;
	double m_dYMin;
	double m_dYMax;
	double m_dZMin;
	double m_dZMax;
};


#endif


