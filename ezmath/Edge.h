#ifndef EDGE_H_
#define EDGE_H_

#include "GenericNode.h"
#include "GeometricComponent.h"
#include "Vector.h"
#include "Face.h"

using namespace EZ;

namespace GeometrySystem
{
class Edge : public GeometricComponent
{
public:
	Edge();
	Edge(const Edge& oEdge);
	~Edge();
	Edge& operator=(const Edge& oEdge);
	void Reset();
	void SetEndPoints(GenericNode* poStartPoint,GenericNode* poEndPoint);
	void SetRightFace(Face* poFace);
	void SetLeftFace(Face* poFace);
	GenericNode* GetStartPoint() const;
	GenericNode* GetEndPoint() const;
	Face* GetRightFace() const;
	Face* GetLeftFace() const;
	Face* GetOtherFace(Face* poFace) const;
	void ReplaceFace(Face* poOldFace,Face* poNewFace);
	string ToString() const;
	void RemoveFace(Face* poFace);
	Vector GetVector() const;
	double GetLength() const;

private:

protected:
	void Initialize();
	GenericNode* m_poStartPoint;
	GenericNode* m_poEndPoint;
	Face* m_poRightFace;
	Face* m_poLeftFace;
};
}

#endif


