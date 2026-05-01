#ifndef FEMBOUNDARYELEMENTFACE_H_
#define FEMBOUNDARYELEMENTFACE_H_

#include "vector"
#include "FEMElementGeometry.h"

using namespace std;

namespace FEMSystem
{
	class FEMBoundaryElementFace
	{
	public:
		FEMBoundaryElementFace();
		FEMBoundaryElementFace(const FEMBoundaryElementFace& oFace);
		~FEMBoundaryElementFace();
		FEMBoundaryElementFace& operator=(const FEMBoundaryElementFace& oFace);
		void Reset();
		void Set(FEMElementGeometry* poGeometry,const unsigned int& iFaceIndex,const unsigned int& iPointsCount);
		void Print() const;
		vector<Point>* GetGaussPoints();
		void SetGaussPointStress(const unsigned int& iIndex,const Matrix& oStress);
		unsigned int GetGaussPointsCount() const;
		void ApplyNodalForcesFromStresses();		// this function assumes a solid mechanics or a thermomechanical problem only
		
	private:
	
	protected:
		void Initialize();
		// an object of this class represents a face in an element that is on the domain's boundary
		// it is mainly used to add extra forces, typically image forces from DD, to the nodes
		// directly without modeling the process using the more restrictive FEMLoad classes.
		// for this calculation to run efficiently, we need to store the following
		// 1. the owning element's nodes
		// 2. the face Gauss points
		// 3. the face normals at the Gauss points
		// 4. the face Gauss points weights (for area integration)
		// 5. the shape function of each node at each Gauss point (for fast load distribution over nodes)
		vector<FEMNode*> m_vpoNodes;
		vector<Point> m_voGaussPoints;
		vector<Vector> m_voGaussPointsNormals;
		vector<double> m_vdGaussPointsAreas;
		Matrix m_oShapeFunctions;
		vector<Matrix> m_voGaussPointsStresses;
	};
}


#endif


