 #include "MainDataStructure.h"
 #include "QuadPatch.h"
 #include "float.h"
 #include "Tools.h"
 #include "FEMMesh.h"
 
 using namespace GeometrySystem;
 using namespace SupportSystem;
 
MainDataStructure* MainDataStructure::MainDataStructureInstance = NULL;
MainDataStructure* MainDataStructure::CreateInstance()
{
	if(MainDataStructureInstance == NULL)
	{
		MainDataStructureInstance = new MainDataStructure;
	}
	return MainDataStructureInstance;
}
MainDataStructure::MainDataStructure()
{
 	Initialize();
}
MainDataStructure::~MainDataStructure()
{
 	Reset();
}
void MainDataStructure::ReadInput()
{
 	if(m_sInputFileName == "")
 	{
 		printf("Input file not specified\n");
 		return;
 	}
 	FILE* fpFile = fopen(m_sInputFileName.c_str(),"r");
 	string sRead;
 	m_sOutputBaseFileName = GetRealString(500,fpFile);
 	// get problem type, total time, number of steps, and physics
 	sRead = GetRealString(500,fpFile);
 	sscanf(sRead.c_str(),"%d,%d,%lf,%lf,%lf,%d\n",&m_eProblemType,&m_iProblemPhysics,&m_dCurrentTime,&m_dTimeStep,&m_dTargetTime,&m_iCurrentOutputCount);
 	// read the mesh, initial and boundary conditions
 	ReadModel(fpFile);
 	fclose(fpFile);
}
void MainDataStructure::ReadModel(FILE* fpFile)
{
	// read loads
	string sRead = GetRealString(500,fpFile);
	unsigned int iLoadsCount = 0;
 	sscanf(sRead.c_str(),"%d\n",&iLoadsCount);
 	unsigned int i = 0;
 	FEMLoad* poLoad = NULL;
 	m_vpoLoads.resize(iLoadsCount);
 	unsigned int iTemp = 0;
 	for(i = 0 ; i < iLoadsCount ; i++)
 	{
 		sRead = GetRealString(500,fpFile);
 		sscanf(sRead.c_str(),"%d\n",&iTemp);
 		poLoad = FEMLoad::GenerateLoadByTypeIndex(iTemp);
 		if(poLoad == NULL)
 		{
 			printf("error: unknown load type %d\n",iTemp);
 		}
 		poLoad->Read(fpFile);
 		poLoad->SetID(i + 1);
 		m_vpoLoads[i] = poLoad;
 	}
 	// read materials
	sRead = GetRealString(500,fpFile);
	unsigned int iMaterialsCount = 0;
 	sscanf(sRead.c_str(),"%d\n",&iMaterialsCount);
 	FEMMaterial* poMaterial = NULL;
 	m_vpoMaterials.resize(iMaterialsCount);
 	for(i = 0 ; i < iMaterialsCount ; i++)
 	{
 		sRead = GetRealString(500,fpFile);
 		sscanf(sRead.c_str(),"%d\n",&iTemp);
 		poMaterial = FEMMaterial::GenerateMaterialByTypeIndex(iTemp);
 		if(poMaterial == NULL)
 		{
 			printf("error: unknown material type %d\n",iTemp);
 		}
 		poMaterial->Read(fpFile);
 		poMaterial->SetID(i + 1);
 		m_vpoMaterials[i] = poMaterial;
 	}
 	// read mesh
	char cWrite[500];
 	vector<FEMNode*> vpoNodes;
 	vector<FEMElement*> vpoElements;
 	vpoNodes.clear();
 	vpoElements.clear();
  	FEMMesh::GenerateMeshFromFile(fpFile,&vpoNodes,&vpoElements,&m_vpoLoads,&m_vpoMaterials);
	// set the nodes and the elements
  	SetNodes(vpoNodes);
	SetElements(vpoElements);
	RegisterDOFs();
 }
 void MainDataStructure::RegisterDOFs()
 {
 	// this function is called after everything is read and it simply
 	// registers the DOFs
 	// first, set their indices
 	unsigned int iNodesCount = (unsigned int)m_vpoNodes.size();
 	unsigned int i = 0;
 	unsigned int iCurrentIndex = 1;
 	for(i = 0 ; i < iNodesCount ; i++)
 	{
		iCurrentIndex = m_vpoNodes[i]->SetDOFIndices(iCurrentIndex);
 	}
 	// second, get them and keep them in record
 	m_vpoDOFs.resize(iCurrentIndex);
 	vector<FEMDegreeOfFreedom*> vpoDOFs;
 	unsigned int iDOFCount = 0;
 	unsigned int j = 0;
 	for(i = 0 ; i < iNodesCount ; i++)
 	{
 		vpoDOFs = m_vpoNodes[i]->GetDOFs();
 		iDOFCount = (unsigned int)vpoDOFs.size();
 		for(j = 0 ; j < iDOFCount ; j++)
 		{
 			m_vpoDOFs[vpoDOFs[j]->GetIndex() - 1] = vpoDOFs[j];
 		}
 	}
 }
 void MainDataStructure::Reset()
 {
 	unsigned int i = 0;
 
 	unsigned int iSize = m_vpoElements.size();
 	for(i = 0; i < iSize ; i++)
 	{
 		if(m_vpoElements[i] != NULL)
 		{
 			delete m_vpoElements[i];
 		}
 	}
 	m_vpoElements.clear();
 	// clear the DOF pointers, DOF objects are owned by the nodes
 	m_vpoDOFs.clear();
 	
 	iSize = m_vpoNodes.size();
 	for(i = 0; i < iSize ; i++)
 	{
 		if(m_vpoNodes[i] != NULL)
 		{
 			delete m_vpoNodes[i];
 		}
 	}
 	m_vpoNodes.clear();
 	
 	iSize = m_vpoLoads.size();
 	for(i = 0; i < iSize ; i++)
 	{
 		if(m_vpoLoads[i] != NULL)
 		{
 			delete m_vpoLoads[i];
 		}
 	}
 	m_vpoLoads.clear();
 	
 	iSize = m_vpoMaterials.size();
 	for(i = 0; i < iSize ; i++)
 	{
 		if(m_vpoMaterials[i] != NULL)
 		{
 			delete m_vpoMaterials[i];
 		}
 	}
 	m_vpoMaterials.clear();
 	
 	list<TriPatch*>::iterator liSurfaceTriangles;
 	for(liSurfaceTriangles = m_lpoSurfaceTriangles.begin() ; liSurfaceTriangles != m_lpoSurfaceTriangles.end() ; liSurfaceTriangles++)
 	{
 		if(*liSurfaceTriangles != NULL)
 		{
 			delete *liSurfaceTriangles;
 		}
 	}
 	m_lpoSurfaceTriangles.clear();
 
 	list<GenericNode*>::iterator liPoints;
 	for(liPoints = m_lpoTriangulationPoints.begin() ; liPoints != m_lpoTriangulationPoints.end() ; liPoints++)
 	{
 		if(*liPoints != NULL)
 		{
 			delete *liPoints;
 		}
 	}
 	
 	list<FEMBoundaryElementFace*>::iterator liBoundaryFaces;
 	for(liBoundaryFaces = m_lpoBoundaryFaces.begin() ; liBoundaryFaces != m_lpoBoundaryFaces.end() ; liBoundaryFaces++)
 	{
 		if(*liBoundaryFaces != NULL)
 		{
 			delete *liBoundaryFaces;
 		}
 	}
 	m_lpoBoundaryFaces.clear();
 	
 	m_lpoTriangulationPoints.clear();
 	m_sInputFileName = "";
 	m_sOutputBaseFileName = "";
 	m_eProblemType = NullProblemType;
 	m_iCurrentOutputCount = 0;
 	m_dCurrentTime = 0.0;
 	m_dTimeStep = 0.0;
 	m_dTargetTime = 0.0;
	m_iProblemPhysics = 0;
  	m_oElementsBSPTree.Collapse();
 }
 void MainDataStructure::SetInputFileName(const string& sFileName)
 {
 	m_sInputFileName = sFileName;
 }
 string MainDataStructure::GetInputFileName() const
 {
 	return m_sInputFileName;
 }
 string MainDataStructure::GetOutputBaseFileName() const
 {
 	return m_sOutputBaseFileName;
 }
 ProblemType MainDataStructure::GetProblemType() const
 {
 	return m_eProblemType;
 }
 unsigned int MainDataStructure::GetProblemPhysics() const
 {
 	return m_iProblemPhysics;
 }
unsigned int MainDataStructure::GetCurrentOutputCount() const
{
	return m_iCurrentOutputCount;
}
void MainDataStructure::IncrementCurrentOutputCount()
{
	m_iCurrentOutputCount = m_iCurrentOutputCount + 1;
}
void MainDataStructure::UpdateTime(const double& dTime)
{
	m_dCurrentTime = dTime;
}
void MainDataStructure::WriteHeader(FILE* fpFile) const
{
	fprintf(fpFile,"* output base file name\n");
	fprintf(fpFile,"%s\n",m_sOutputBaseFileName.c_str());
	fprintf(fpFile,"* problem type, physics, current time, time step, target time, this output count\n");
	fprintf(fpFile,"%d,%d,%e,%e,%e,%d\n",m_eProblemType,m_iProblemPhysics,m_dCurrentTime,m_dTimeStep,m_dTargetTime,m_iCurrentOutputCount);
}
void MainDataStructure::WriteLoads(FILE* fpFile) const
{
	fprintf(fpFile,"* Loads :\n");
	unsigned int iSize = (unsigned int)m_vpoLoads.size();
	unsigned int i = 0;
	fprintf(fpFile,"%d\n",iSize);
	for(i = 0 ; i < iSize ; i++)
	{
		m_vpoLoads[i]->Write(fpFile);
	}
}
void MainDataStructure::WriteMaterials(FILE* fpFile) const
{
	fprintf(fpFile,"* Materials :\n");
	unsigned int iSize = (unsigned int)m_vpoMaterials.size();
	unsigned int i = 0;
	fprintf(fpFile,"%d\n",iSize);
	for(i = 0 ; i < iSize ; i++)
	{
		m_vpoMaterials[i]->Write(fpFile);
	}
}
double MainDataStructure::GetCurrentTime() const
{
	return m_dCurrentTime;
}
double MainDataStructure::GetTargetTime() const
{
	return m_dTargetTime;
}
double MainDataStructure::GetTimeStep() const
{
	return m_dTimeStep;
}
 FEMNode* MainDataStructure::GetNode(const unsigned int& iID) const
 {
 	if(iID >= m_vpoNodes.size())
 	{
 		return NULL;
 	}
 	return m_vpoNodes[iID];
 }
 FEMElement* MainDataStructure::GetElement(const unsigned int& iID) const
 {
 	if(iID >= m_vpoElements.size())
 	{
 		return NULL;
 	}
 	return m_vpoElements[iID];
 }
 unsigned int MainDataStructure::GetNodesCount() const
 {
 	return (unsigned int)m_vpoNodes.size();
 }
 unsigned int MainDataStructure::GetElementsCount() const
 {
 	return (unsigned int)m_vpoElements.size();
 }
 vector<FEMNode*>* MainDataStructure::GetNodes()
 {
 	return &m_vpoNodes;
 }
 vector<FEMElement*>* MainDataStructure::GetElements()
 {
 	return &m_vpoElements;
 }
 vector<FEMDegreeOfFreedom*>* MainDataStructure::GetDOFs()
 {
 	return &m_vpoDOFs;
 }
 void MainDataStructure::SetNodes(const vector<FEMNode*>& vpoNodes)
 {
 	unsigned int i = 0;
 	unsigned iSize = m_vpoNodes.size();
 	for(i = 0; i < iSize ; i++)
 	{
 		if(m_vpoNodes[i] != NULL)
 		{
 			delete m_vpoNodes[i];
 		}
 	}
 	m_vpoNodes.clear();
 	iSize = vpoNodes.size();
 	m_vpoNodes.resize(iSize);
 	for(i = 0; i < iSize ; i++)
 	{
 		m_vpoNodes[i] = vpoNodes[i];
 	}
 }
 void MainDataStructure::SetElements(const vector<FEMElement*>& vpoElements)
 {
 	unsigned int i = 0;
 	unsigned iSize = m_vpoElements.size();
 	for(i = 0; i < iSize ; i++)
 	{
 		if(m_vpoElements[i] != NULL)
 		{
 			delete m_vpoElements[i];
 		}
 	}
 	m_vpoElements.clear();
 	iSize = vpoElements.size();
 	m_vpoElements.resize(iSize);
 	for(i = 0; i < iSize ; i++)
 	{
 		m_vpoElements[i] = vpoElements[i];
 	}
 	GenerateSurfaceTriangulation();
 	BuildElementsBSPTree();
 }
 void MainDataStructure::ClearSurfaceTriangulation()
 {
	list<TriPatch*>::iterator liSurfaceTriangles;
 	for(liSurfaceTriangles = m_lpoSurfaceTriangles.begin() ; liSurfaceTriangles != m_lpoSurfaceTriangles.end() ; liSurfaceTriangles++)
 	{
 		if(*liSurfaceTriangles != NULL)
 		{
 			delete *liSurfaceTriangles;
 		}
 	}
 	m_lpoSurfaceTriangles.clear();
 
 	list<GenericNode*>::iterator liPoints;
 	for(liPoints = m_lpoTriangulationPoints.begin() ; liPoints != m_lpoTriangulationPoints.end() ; liPoints++)
 	{
 		if(*liPoints != NULL)
 		{
 			delete *liPoints;
 		}
 	}
 	m_lpoTriangulationPoints.clear();
 }
 void MainDataStructure::ReadSurfaceTriangulationFromFile(const string& sFileName)
 {
	 ClearSurfaceTriangulation();
	 FILE* fpFile = fopen(sFileName.c_str(),"r");
	 unsigned int iPointsCount = 0;
	 
	 string sRead;
	 sRead = GetRealString(500,fpFile);
	 sscanf(sRead.c_str(),"%d\n",&iPointsCount);
	 unsigned int i = 0;
	 double dX = 0.0;
	 double dY = 0.0;
	 double dZ = 0.0;
	 vector<GenericNode*> vpoPoints;
	 vpoPoints.resize(iPointsCount);
	 for(i = 0 ; i < iPointsCount ; i++)
	 {
		 sRead = GetRealString(500,fpFile);
		 sscanf(sRead.c_str(),"%lf\t\t%lf\t\t%lf\n",&dX,&dY,&dZ);
		 vpoPoints[i] = new GenericNode;
		 vpoPoints[i]->Set(dX,dY,dZ);
		 vpoPoints[i]->SetID(i + 1);
		 m_lpoTriangulationPoints.push_back(vpoPoints[i]);
	 }
	 
	 unsigned int iTrianglesCount = 0;
	 sRead = GetRealString(500,fpFile);
	 sscanf(sRead.c_str(),"%d\n",&iTrianglesCount);
	 unsigned int iPoint1ID = 0;
	 unsigned int iPoint2ID = 0;
	 unsigned int iPoint3ID = 0;
	 TriPatch* poTriangle = NULL;
	 for(i = 0 ; i < iTrianglesCount ; i++)
	 {
		 sRead = GetRealString(500,fpFile);
		 sscanf(sRead.c_str(),"%d\t\t%d\t\t%d\n",&iPoint1ID,&iPoint2ID,&iPoint3ID);
		 poTriangle = new TriPatch;
		 poTriangle->Set(vpoPoints[iPoint1ID - 1],vpoPoints[iPoint2ID - 1],vpoPoints[iPoint3ID - 1]);
		 m_lpoSurfaceTriangles.push_back(poTriangle);
	 }
	 vpoPoints.clear();
	 fclose(fpFile);
 }
 void MainDataStructure::ApplyLoads(const double& dTime)
 {
 	unsigned int i = 0;
 	unsigned int iSize = m_vpoNodes.size();
 	for(i = 0; i < iSize ; i++)
 	{
 		m_vpoNodes[i]->ApplyLoads(dTime);
 	}
 	iSize = m_vpoElements.size();
 	for(i = 0; i < iSize ; i++)
 	{	
 		m_vpoElements[i]->ApplyLoads(dTime);
 	}
 }
 void MainDataStructure::ResetLoads()
 {
 	unsigned int i = 0;
 	unsigned int iSize = m_vpoNodes.size();
 	for(i = 0; i < iSize ; i++)
 	{
 		m_vpoNodes[i]->ResetLoads();
 	}
 }
 vector<FEMLoad*>* MainDataStructure::GetLoads()
 {
 	return &m_vpoLoads;
 }
 void MainDataStructure::Initialize()
 {
 	m_vpoElements.clear();
 	m_vpoDOFs.clear();
 	m_vpoNodes.clear();
 	m_vpoLoads.clear();
 	m_vpoMaterials.clear();
 	m_lpoSurfaceTriangles.clear();
 	m_lpoTriangulationPoints.clear();
 	m_lpoBoundaryFaces.clear();
 	m_sInputFileName = "";
 	m_sOutputBaseFileName = "";
 	m_eProblemType = NullProblemType;
 	m_iCurrentOutputCount = 0;
 	m_dCurrentTime = 0.0;
 	m_dTimeStep = 0.0;
 	m_dTargetTime = 0.0;
	m_iProblemPhysics = 0;
 	m_oElementsBSPTree.Initialize();
 }
 double MainDataStructure::GetModelVolume() const
 {
 	unsigned int i = 0;
 	unsigned int iSize = (unsigned int)m_vpoElements.size();
 	double dVolume = 0.0;
 	for(i = 0; i < iSize ; i++)
 	{
 		dVolume = dVolume + m_vpoElements[i]->GetVolume();
 	}
 	return dVolume;
 }
 list<TriPatch*>* MainDataStructure::GetSurfaceTriangulation()
 {
 	return &m_lpoSurfaceTriangles;
 }
 list<GenericNode*>* MainDataStructure::GetSurfaceTriangulationPoints()
 {
 	return &m_lpoTriangulationPoints;
 }
void MainDataStructure::GenerateSurfaceTriangulation()
{
 	unsigned int i = 0;
 	unsigned int iNodesCount = (unsigned int)m_vpoNodes.size();
 	unsigned int iElementsCount = (unsigned int)m_vpoElements.size();
 	vector<GenericNode*> vpoPoints;
 	vector<bool> vbIsUsed;
 	vpoPoints.reserve(iNodesCount + 2*iElementsCount);      // the 2 factor is just an estimate
 	vbIsUsed.reserve(iNodesCount + 2*iElementsCount);
 	GenericNode* poNode = NULL;
 	for(i = 0; i < iNodesCount ; i++)
 	{
		poNode = new GenericNode;
		poNode->Set(m_vpoNodes[i]->GetX(),m_vpoNodes[i]->GetY(),m_vpoNodes[i]->GetZ());
		poNode->SetID(m_vpoNodes[i]->GetID());
 		vpoPoints.push_back(poNode);
 		vbIsUsed.push_back(false);
 	}
 	
 	vector< vector<unsigned int> > vviNodesIndices;
 	vector<GenericNode*> vpoMidPoints;
 	unsigned int j = 0;
 	unsigned int k = 0;
 	unsigned int iPatchesCount = 0;
 	vector<GenericNode*> vpoTempPoints;
 	vector<QuadPatch> voPatches;
 	voPatches.reserve(2*iElementsCount);
 	
 	// loop over the actual elements and get the quad patches
 	for(i = 0; i < iElementsCount ; i++)
 	{
 		m_vpoElements[i]->GetGeometry()->GenerateSurfacePatches(vviNodesIndices,vpoMidPoints);
 		iPatchesCount = vviNodesIndices.size();
 		for(j = 0; j < iPatchesCount ; j++)
 		{
 			if(vviNodesIndices[j].size() != PointsPerQuadPatch - 1 || vpoMidPoints[j] == NULL)
 			{
 				continue;
 			}
 			vpoPoints.push_back(vpoMidPoints[j]);
 			vbIsUsed.push_back(false);
 			vpoTempPoints.resize(PointsPerQuadPatch);
 			for(k = 0; k < PointsPerQuadPatch - 1 ; k++)
 			{
 				vpoTempPoints[k] = vpoPoints[vviNodesIndices[j].at(k) - 1];
 				vbIsUsed[vviNodesIndices[j].at(k) - 1] = true;
 			}
 			vpoTempPoints[PointsPerQuadPatch - 1] = vpoPoints.back();
 			vbIsUsed.back() = true;
 			voPatches.push_back(QuadPatch(vpoTempPoints));
 		}
 	}
 	
 	// now all the used points and the quad patches are ready
 	unsigned int iSize = vpoPoints.size();
 	for(i = 0; i < iSize ; i++)
 	{
 		if(vbIsUsed[i])
 		{
 			m_lpoTriangulationPoints.push_back(vpoPoints[i]);
 		}
 		else
 		{
 			delete vpoPoints[i];
 		}
 	}
 	vpoPoints.clear();
 
 	list<TriPatch*>::iterator liSurfaceTriangles;
 	for(liSurfaceTriangles = m_lpoSurfaceTriangles.begin() ; liSurfaceTriangles != m_lpoSurfaceTriangles.end() ; liSurfaceTriangles++)
 	{
 		if(*liSurfaceTriangles != NULL)
 		{
 			delete *liSurfaceTriangles;
 		}
 	}
 	m_lpoSurfaceTriangles.clear();
 
 	iSize = voPatches.size();
 	vector<Patch*> vpoTempTriangles;
 	unsigned int iTrianglesCount = 0;
 	for(i = 0; i < iSize ; i++)
 	{
 		vpoTempTriangles = voPatches[i].GenerateTriangulation();
 		iTrianglesCount = vpoTempTriangles.size();
 		for(j = 0 ; j < iTrianglesCount ; j++)
 		{
 			m_lpoSurfaceTriangles.push_back((TriPatch*)vpoTempTriangles[j]);
 		}
 	}
 	unsigned int iIndex = 0;
 	list<GenericNode*>::iterator liPoints;
 	for(liPoints = m_lpoTriangulationPoints.begin() ; liPoints != m_lpoTriangulationPoints.end() ; liPoints++)
 	{
 		(*liPoints)->SetID(iIndex);
 		iIndex = iIndex + 1;
 	}
}
 void MainDataStructure::RefineTriangulation()
 {
 	list<TriPatch*>::iterator liTriangles;
 	double dMinSideLength = DBL_MAX;
 	double dTemp = 0.0;
 	for(liTriangles = m_lpoSurfaceTriangles.begin() ; liTriangles != m_lpoSurfaceTriangles.end() ; liTriangles++)
 	{
 		dTemp = (*liTriangles)->GetMinSideLength();
 		if(dTemp < dMinSideLength)
 		{
 			dMinSideLength = dTemp;
 		}
 	}
 	double dTolerance = 1E-6*dMinSideLength;
 	list<GenericNode*> lpoNewPoints;
 	lpoNewPoints.clear();
 	list<TriPatch*> lpoNewTriangles;
 	lpoNewTriangles.clear();
 	for(liTriangles = m_lpoSurfaceTriangles.begin() ; liTriangles != m_lpoSurfaceTriangles.end() ; liTriangles++)
 	{
 		(*liTriangles)->Refine(lpoNewPoints,lpoNewTriangles,dTolerance);
 		delete *liTriangles;
 	}
 	m_lpoSurfaceTriangles.clear();
 	list<GenericNode*>::iterator liPoints;
 	// add the new points to the old ones
 	for(liPoints = lpoNewPoints.begin() ; liPoints != lpoNewPoints.end() ; liPoints++)
 	{
 		m_lpoTriangulationPoints.push_back(*liPoints);
 	}
 	lpoNewPoints.clear();
 	// add the new triangles to the now empty list
 	for(liTriangles = lpoNewTriangles.begin() ; liTriangles != lpoNewTriangles.end() ; liTriangles++)
 	{
 		m_lpoSurfaceTriangles.push_back(*liTriangles);
 	}
         
 	unsigned int iIndex = 0;
 	for(liPoints = m_lpoTriangulationPoints.begin() ; liPoints != m_lpoTriangulationPoints.end() ; liPoints++)
 	{
 		(*liPoints)->SetID(iIndex);
 		iIndex = iIndex + 1;
 	}
 }
void MainDataStructure::RefineSurfaceTriangulation(const unsigned int& iRefinementsCount)
{
 	unsigned int i = 0;
 	for(i = 0; i < iRefinementsCount ; i++)
 	{
 		RefineTriangulation();
 	}
}
BSPTreeNode<FEMElement*>* MainDataStructure::GetElementsBSPTree()
{
	return &m_oElementsBSPTree;
}
void MainDataStructure::WriteSurfaceTriangulation(const string& sFileName)
{
	FILE* fpFile = fopen(sFileName.c_str(),"w");
	
	list<GenericNode*>::iterator liPoints;
	fprintf(fpFile,"%d\n",(unsigned int)m_lpoTriangulationPoints.size());
	for(liPoints = m_lpoTriangulationPoints.begin() ; liPoints != m_lpoTriangulationPoints.end() ; liPoints++)
	{
		fprintf(fpFile,"%lf\t\t%lf\t\t%lf\n",(*liPoints)->GetX(),(*liPoints)->GetY(),(*liPoints)->GetZ());
	}
	
	list<TriPatch*>::iterator liTriangles;
	fprintf(fpFile,"%d\n",(unsigned int)m_lpoSurfaceTriangles.size());
	for(liTriangles = m_lpoSurfaceTriangles.begin() ; liTriangles != m_lpoSurfaceTriangles.end() ; liTriangles++)
	{
		fprintf(fpFile,"%d\t\t%d\t\t%d\n",(*liTriangles)->GetPoint1()->GetID(),(*liTriangles)->GetPoint2()->GetID(),(*liTriangles)->GetPoint3()->GetID());
	}
	
	fclose(fpFile);
}
void MainDataStructure::WriteSurfaceTriangulationToParaview(const string& sFileName)
{
	FILE* fpFile = fopen(sFileName.c_str(),"w");
	
	// write file header
	fprintf(fpFile,"# vtk DataFile Version 2.0\n");
	fprintf(fpFile,"surface triagulation from FEM mesh\n");
	fprintf(fpFile,"ASCII\n");
	fprintf(fpFile,"DATASET UNSTRUCTURED_GRID\n");
	
	// write the nodes
	list<GenericNode*>::iterator liPoints;
	fprintf(fpFile,"POINTS %d float\n",(unsigned int)m_lpoTriangulationPoints.size());
	for(liPoints = m_lpoTriangulationPoints.begin() ; liPoints != m_lpoTriangulationPoints.end() ; liPoints++)
	{
		fprintf(fpFile,"%e\t\t%e\t\t%e\n",(*liPoints)->GetX(),(*liPoints)->GetY(),(*liPoints)->GetZ());
	}
	
	list<TriPatch*>::iterator liTriangles;
	unsigned int iTrianglesCount = (unsigned int)m_lpoSurfaceTriangles.size();
	fprintf(fpFile,"CELLS %d %d\n",iTrianglesCount,4*iTrianglesCount);
	for(liTriangles = m_lpoSurfaceTriangles.begin() ; liTriangles != m_lpoSurfaceTriangles.end() ; liTriangles++)
	{
		fprintf(fpFile,"3\t\t%d\t\t%d\t\t%d\n",(*liTriangles)->GetPoint1()->GetID(),(*liTriangles)->GetPoint2()->GetID(),(*liTriangles)->GetPoint3()->GetID());
	}
	
	fprintf(fpFile,"CELL_TYPES %d\n",iTrianglesCount);
	for(liTriangles = m_lpoSurfaceTriangles.begin() ; liTriangles != m_lpoSurfaceTriangles.end() ; liTriangles++)
	{
		fprintf(fpFile,"5\n");
	}
	
	fprintf(fpFile,"POINT_DATA %d\n",(unsigned int)m_lpoTriangulationPoints.size());
	fprintf(fpFile,"SCALARS node_id float 1\n");
	fprintf(fpFile,"LOOKUP_TABLE default\n");
	unsigned int i = 0;
	for(liPoints = m_lpoTriangulationPoints.begin() ; liPoints != m_lpoTriangulationPoints.end() ; liPoints++)
	{
		fprintf(fpFile,"%f\n",(float)i);
		i = i + 1;
	}
	
	fclose(fpFile);
}
void MainDataStructure::BuildElementsBSPTree()
{
	list<FEMElement*> lpoElements;
	lpoElements.clear();
	unsigned int iSize = (unsigned int)m_vpoElements.size();
	unsigned int i = 0;
	for(i = 0; i < iSize ; i++)
	{
		lpoElements.push_back(m_vpoElements[i]);
	}
	BuildElementsBSPTreeNode(&lpoElements,&m_oElementsBSPTree);
	lpoElements.clear();
}
void MainDataStructure::BuildElementsBSPTreeNode(list<FEMElement*>* plpoElements,BSPTreeNode<FEMElement*>* poTreeNode)
{
	// if the number of elements is too small, do not split the node space any further, 
	// just add the elements to the node
	unsigned int iMaxNodeElementsCount = 2;
	list<FEMElement*>::iterator liElements;
	unsigned int iSize = (unsigned int)plpoElements->size();
	if(iSize <= iMaxNodeElementsCount)
	{
		for(liElements = plpoElements->begin() ; liElements != plpoElements->end() ; liElements++)
		{
			poTreeNode->AddItem((*liElements));
		}
		return;
	}
	// if there are many elements, split the node space into two spaces along a random dimension
	list<FEMElement*> lpoLeftList;
	list<FEMElement*> lpoRightList;
	lpoLeftList.clear();
	lpoRightList.clear();
	vector<Point> voCenters;
	voCenters.reserve(iSize);
	// get the element centroids and the bounding box of the centroids
	Point oCenter;
	AxisAlignedBoundingBox oBox;
	oBox.CenterAt(plpoElements->front()->GetGeometry()->GetCenterPoint());
	for(liElements = plpoElements->begin() ; liElements != plpoElements->end() ; liElements++)
	{
		oCenter = (*liElements)->GetGeometry()->GetCenterPoint();
		voCenters.push_back(oCenter);
		oBox.ExpandToContain(oCenter);
	}
	// split along the maximum dimension
	int iDimension = oBox.GetMaximumDimensionIndex();
	// set the splitting plane
	poTreeNode->SetSeparator(&voCenters,iDimension);
	int iClassification = 0;
	// classify the elements, balance the tree as much as possible
	unsigned int iLeftTreeElementsCount = 0;
	unsigned int iRightTreeElementsCount = 0;
	for(liElements = plpoElements->begin() ; liElements != plpoElements->end() ; liElements++)
	{
		iClassification = poTreeNode->ClassifyBox((*liElements)->GetGeometry()->GetBox());
		//iClassification = poTreeNode->Classify((*liElements)->GetGeometry()->GetCenterPoint());
		if(iClassification == -1)
		{
			lpoLeftList.push_back((*liElements));
			iLeftTreeElementsCount = iLeftTreeElementsCount + 1;
		}
		else if(iClassification == 1)
		{
			lpoRightList.push_back((*liElements));
			iRightTreeElementsCount = iRightTreeElementsCount + 1;
		}
		else
		{
			if(iLeftTreeElementsCount < iRightTreeElementsCount)
			{
				lpoLeftList.push_back((*liElements));
				iLeftTreeElementsCount = iLeftTreeElementsCount + 1;
			}
			else
			{
				lpoRightList.push_back((*liElements));
				iRightTreeElementsCount = iRightTreeElementsCount + 1;
			}
		}
	}
	// now generate the two child nodes based on the elements list
	BSPTreeNode<FEMElement*>* poLeftNode = new BSPTreeNode<FEMElement*>;
	BSPTreeNode<FEMElement*>* poRightNode = new BSPTreeNode<FEMElement*>;
	BuildElementsBSPTreeNode(&lpoLeftList,poLeftNode);
	BuildElementsBSPTreeNode(&lpoRightList,poRightNode);
	poTreeNode->SetLeftChild(poLeftNode);
	poTreeNode->SetRightChild(poRightNode);
	lpoLeftList.clear();
	lpoRightList.clear();
}
void MainDataStructure::GenerateBoundaryFaces()
{
	unsigned int i = 0;
	unsigned int j = 0;
	unsigned int iSize = (unsigned int)m_vpoElements.size();
	unsigned int iFacesCount = 0;
	vector<unsigned int> viSurfaceFacesIndices;
	FEMElementGeometry* poGeometry = NULL;
	FEMBoundaryElementFace* poFace = NULL;
	unsigned int iFacePointResolution = 4;
	for(i = 0 ; i < iSize ; i++)
	{
		poGeometry = m_vpoElements[i]->GetGeometry();
		viSurfaceFacesIndices = poGeometry->GetSurfaceFacesIndices();
		iFacesCount = (unsigned int)viSurfaceFacesIndices.size();
		for(j = 0 ; j < iFacesCount ; j++)
		{
			poFace = new FEMBoundaryElementFace;
			poFace->Set(poGeometry,viSurfaceFacesIndices[j],iFacePointResolution);
			m_lpoBoundaryFaces.push_back(poFace);
		}
	}
}
list<FEMBoundaryElementFace*>* MainDataStructure::GetBoundaryFaces()
{
	return &m_lpoBoundaryFaces;
}
void MainDataStructure::ApplyImageForces()
{
	list<FEMBoundaryElementFace*>::iterator liFaces;
	for(liFaces = m_lpoBoundaryFaces.begin() ; liFaces != m_lpoBoundaryFaces.end() ; liFaces++)
	{
		(*liFaces)->ApplyNodalForcesFromStresses();
	}
}
const vector<FEMMaterial*>* MainDataStructure::GetMaterials() const
{
	return &m_vpoMaterials;
}





