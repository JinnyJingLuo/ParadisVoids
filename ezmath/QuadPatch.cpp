// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#include "QuadPatch.h"


namespace GeometrySystem
{
	QuadPatch::QuadPatch()
	{
		m_vpoPoints.clear();
	}
	QuadPatch::QuadPatch(const QuadPatch& oElement)
	{
		*this = oElement;
	}
	QuadPatch::QuadPatch(vector<GenericNode*> vpoPoints)
	{
		Set(vpoPoints);
	}
	QuadPatch::~QuadPatch()
	{
		m_vpoPoints.clear();
	}
	QuadPatch& QuadPatch::operator=(const QuadPatch& oElement)
	{
		Set(oElement.m_vpoPoints);
		return *this;
	}
	void QuadPatch::Set(vector<GenericNode*> vpoPoints)
	{
		if(vpoPoints.size() != PointsPerQuadPatch)
		{
			return;
		}
		unsigned int i = 0;
		m_vpoPoints.resize(PointsPerQuadPatch);
		for(i = 0; i < PointsPerQuadPatch ; i++)
		{
			m_vpoPoints[i] = vpoPoints[i];
		}
	}
	vector<Patch*> QuadPatch::GenerateTriangulation() const
	{
		vector<Patch*> vpoTriangles;
		vpoTriangles.resize(8);
		vpoTriangles[0] = new TriPatch(m_vpoPoints[0],m_vpoPoints[1],m_vpoPoints[7]);
		vpoTriangles[1] = new TriPatch(m_vpoPoints[1],m_vpoPoints[8],m_vpoPoints[7]);
		vpoTriangles[2] = new TriPatch(m_vpoPoints[1],m_vpoPoints[2],m_vpoPoints[3]);
		vpoTriangles[3] = new TriPatch(m_vpoPoints[1],m_vpoPoints[3],m_vpoPoints[8]);
		vpoTriangles[4] = new TriPatch(m_vpoPoints[7],m_vpoPoints[8],m_vpoPoints[5]);
		vpoTriangles[5] = new TriPatch(m_vpoPoints[7],m_vpoPoints[5],m_vpoPoints[6]);
		vpoTriangles[6] = new TriPatch(m_vpoPoints[8],m_vpoPoints[3],m_vpoPoints[5]);
		vpoTriangles[7] = new TriPatch(m_vpoPoints[3],m_vpoPoints[4],m_vpoPoints[5]);
		return vpoTriangles;
	}
}

