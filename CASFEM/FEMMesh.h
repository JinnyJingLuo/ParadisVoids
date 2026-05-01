// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef FEMMESH_H_
#define FEMMESH_H_

#include "FEMElement.h"
#include "string"
#include "Block.h"
#include "Cylinder.h"

using namespace std;
using namespace FEMSystem;
using namespace GeometrySystem;

namespace FEMSystem
{
	class FEMMesh
	{
	public:
		// imp mesh
		FEMMesh();
		FEMMesh(const FEMMesh& oMesh);
		virtual ~FEMMesh();
		virtual FEMMesh& operator=(const FEMMesh& oMesh);
		virtual void Reset();
		static void GenerateMeshFromFile(FILE* fpFile,vector<FEMNode*>* pvpoNodes,vector<FEMElement*>* pvpoElements,vector<FEMLoad*>* pvpoLoads,vector<FEMMaterial*>* pvpoMaterials);
		static void GenerateMeshFromCylinder(Cylinder* poCylinder,vector<FEMNode*>& vpoNodes,vector<FEMElement*>& vpoElements,FEMNodeType eNodeType,FEMElementType eElementType);
		static void GenerateMeshFromBlock(Block* poBlock,vector<FEMNode*>& vpoNodes,vector<FEMElement*>& vpoElements,FEMNodeType eNodeType,FEMElementType eElementType);

	private:

	protected:
		virtual void Initialize();
	};
}
 
#endif
 


