#ifndef MAINDATASTRUCTURE_H_
#define MAINDATASTRUCTURE_H_

#include "SparseMatrix.h"
#include "TriPatch.h"
#include "AxisAlignedBoundingBox.h"
#include "FEMLoad.h"
#include "FEMNode.h"
#include "FEMElement.h"
#include "GenericNode.h"
#include "BSPTreeNode.h"
#include "FEMBoundaryElementFace.h"

using namespace FEMSystem;
using namespace EZ;
using namespace GeometrySystem;

enum ProblemType
{
	NullProblemType = 0,
	StaticProblem = 1,
	ExplicitDynamicsProblem = 2,
	ImplicitDynamicsProblem = 3
};
	
	
// this class is a singleton
class MainDataStructure
{
public:
	static MainDataStructure* CreateInstance();
 	~MainDataStructure();
 	void Reset();
	
	// file names handling
	void SetInputFileName(const string& sFileName);
 	void ReadInput();
    string GetInputFileName() const;
    string GetOutputBaseFileName() const;

	// parameters handling
    ProblemType GetProblemType() const;
    double GetCurrentTime() const;
    double GetTargetTime() const;
    double GetTimeStep() const;
	unsigned int GetProblemPhysics() const;
	unsigned int GetCurrentOutputCount() const;
	void IncrementCurrentOutputCount();
	void UpdateTime(const double& dTime);
	void WriteHeader(FILE* fpFile) const;
	void WriteLoads(FILE* fpFile) const;
	void WriteMaterials(FILE* fpFile) const;

	// data communication
 	FEMNode* GetNode(const unsigned int& iID) const;
 	FEMElement* GetElement(const unsigned int& iID) const;
 	unsigned int GetNodesCount() const;
 	unsigned int GetElementsCount() const;
 	vector<FEMNode*>* GetNodes();
 	vector<FEMElement*>* GetElements();
 	vector<FEMDegreeOfFreedom*>* GetDOFs();
 	void SetNodes(const vector<FEMNode*>& vpoNodes);
 	void SetElements(const vector<FEMElement*>& vpoElements);
 	const vector<FEMMaterial*>* GetMaterials() const;

 	// load handling
 	void ApplyLoads(const double& dTime);
 	void ResetLoads();
 	vector<FEMLoad*>* GetLoads();
 	void ApplyImageForces();		// this function assumes a solid mechanics or a thermomechanical problem only

 	double GetModelVolume() const;
     
    // surface triangulations
    list<TriPatch*>* GetSurfaceTriangulation();
  	list<GenericNode*>* GetSurfaceTriangulationPoints();
  	void RefineSurfaceTriangulation(const unsigned int& iRefinementsCount = 1);
  	void ClearSurfaceTriangulation();
  	void ReadSurfaceTriangulationFromFile(const string& sFileName);
  	void WriteSurfaceTriangulation(const string& sFileName);
  	void WriteSurfaceTriangulationToParaview(const string& sFileName);
  	
  	// element tree
  	BSPTreeNode<FEMElement*>* GetElementsBSPTree();
  	void GenerateBoundaryFaces();
  	list<FEMBoundaryElementFace*>* GetBoundaryFaces();
 
 private:
 
 protected:
	static MainDataStructure* MainDataStructureInstance;
	MainDataStructure();
 	void Initialize();
 	void ReadModel(FILE* fpFile);
 	void GenerateSurfaceTriangulation();
 	void RefineTriangulation();
 	void RegisterDOFs();
 	void BuildElementsBSPTree();
 	void BuildElementsBSPTreeNode(list<FEMElement*>* plpoElements,BSPTreeNode<FEMElement*>* poTreeNode);
	
	// file names
 	string m_sInputFileName;
 	string m_sOutputBaseFileName;

	// solution parameters
 	ProblemType m_eProblemType;
	unsigned int m_iProblemPhysics;		// 1 - potential, 2 - solid mechanics
 	unsigned int m_iCurrentOutputCount;
 	double m_dCurrentTime;
 	double m_dTimeStep;
 	double m_dTargetTime;

	// model data
 	vector<FEMElement*> m_vpoElements;
 	vector<FEMNode*> m_vpoNodes;
 	vector<FEMLoad*> m_vpoLoads;
 	vector<FEMDegreeOfFreedom*> m_vpoDOFs;
 	vector<FEMMaterial*> m_vpoMaterials;

	// data frequently used by external agents
 	list<TriPatch*> m_lpoSurfaceTriangles;
 	list<GenericNode*> m_lpoTriangulationPoints;
	BSPTreeNode<FEMElement*> m_oElementsBSPTree;
	list<FEMBoundaryElementFace*> m_lpoBoundaryFaces;
};

#endif



