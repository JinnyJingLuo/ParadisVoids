// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef QUADPATCH_H_
#define QUADPATCH_H_

#include "vector"
#include "GenericNode.h"
#include "TriPatch.h"

using namespace std;
using namespace EZ;

namespace GeometrySystem
{
	const unsigned int PointsPerQuadPatch = 9;
	class QuadPatch : public Patch
	{
	public:
		QuadPatch();
		QuadPatch(const QuadPatch& oPatch);
		QuadPatch(vector<GenericNode*> vpoPoints);
		~QuadPatch();
		QuadPatch& operator=(const QuadPatch& oElement);
		void Set(vector<GenericNode*> vpoPoints);
		vector<Patch*> GenerateTriangulation() const;
	private:

	protected:
		vector<GenericNode*> m_vpoPoints;	// the patch mid point is the last point
	};
}

#endif

