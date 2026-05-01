//// Ahmed M. Hussein
//
//#include "LoopFilter.h"
//#include "string"
//#include "Tools.h"
//#include "map"
//
// using namespace std;
// using namespace SupportSystem;
//
// LoopFilter* LoopFilter::m_poLoopFilterInstance = NULL;
// LoopFilter* LoopFilter::CreateInstance()
//{
//	if(m_poLoopFilterInstance == NULL)
//	{
//		m_poLoopFilterInstance = new LoopFilter;
//	}
//	return m_poLoopFilterInstance;
//}
// LoopFilter::LoopFilter()
//{
//	Initialize();
//}
// LoopFilter::~LoopFilter()
//{
//	Reset();
//}
// void LoopFilter::Reset()
//{
//	list< DislocationChain* >::iterator liChains;
//	for(liChains = m_lpoSurvivingChains.begin() ; liChains !=
// m_lpoSurvivingChains.end() ; liChains++)
//	{
//		if((*liChains) != NULL)
//		{
//			(*liChains)->Collapse();
//			delete (*liChains);
//		}
//	}
//	m_lpoSurvivingChains.clear();
//	m_oNetwork.Collapse();
//	m_oDualGraph.Collapse();
//	m_oBox.Reset();
//}
// void LoopFilter::Initialize()
//{
//	m_oNetwork.Reset();
//	m_oDualGraph.Reset();
//	m_lpoSurvivingChains.clear();
//	m_oBox.Reset();
//}
// void LoopFilter::Set(const string& sInputFileName,const unsigned int&
// iProcessesCount)
//{
//	Reset();
//	SetBoundingBox(sInputFileName);
//	FILE* fpFile = fopen(sInputFileName.c_str(),"r");
//	SeekRestartFileToData(fpFile);
//	// now loop and read the node data, generate a vector of maps for this
//	string sRead = "";
//	vector< map<unsigned int,DislocationNetworkNode*> > vpmoNodes;
//	vpmoNodes.resize(iProcessesCount);
//	unsigned int iDomainIndex = 0;
//	unsigned int iNodeIndex = 0;
//	double dX = 0.0;
//	double dY = 0.0;
//	double dZ = 0.0;
//	unsigned int iNeighboursCount = 0;
//	unsigned int iConstraint = 0;
//	DislocationNetworkNode* poGraphNode = NULL;
//	DislocationNode oNode;
//	unsigned int iGlobalIndex = 0;
//	unsigned int i = 0;
//	list<DislocationNetworkNode*> lpoNodes;
//	while(!feof(fpFile))
//	{
//		sRead = GetRealString(500,fpFile);
//		if(sRead.empty())
//		{
//			continue;
//		}
//		sscanf(sRead.c_str(),"%d,%d\t%lf\t%lf\t%lf\t%d\t%d\n",&iDomainIndex,&iNodeIndex,&dX,&dY,&dZ,&iNeighboursCount,&iConstraint);
//		iGlobalIndex = iGlobalIndex + 1;
//		poGraphNode = new DislocationNetworkNode;
//		oNode.Set(dX,dY,dZ);
//		poGraphNode->Set(&oNode);
//		poGraphNode->GetDataPointer()->SetCategory(iConstraint);
//		// read the surface normal
//		sRead = GetRealString(500,fpFile);
//		sscanf(sRead.c_str(),"%lf\t%lf\t%lf\n",&dX,&dY,&dZ);
//		poGraphNode->GetDataPointer()->SetSurfaceNormal(Vector(dX,dY,dZ));
//		poGraphNode->SetID(iGlobalIndex);
//		vpmoNodes[iDomainIndex].operator[](iNodeIndex) = poGraphNode;
//		// skip the arms data for now
//		for(i = 0 ; i < iNeighboursCount ; i++)
//		{
//			sRead = GetRealString(500,fpFile);
//			sRead = GetRealString(500,fpFile);
//		}
//		lpoNodes.push_back(poGraphNode);
//	}
//
//	unsigned int iMaxMapIndex = 0;
//	unsigned int iTemp = 0;
//	for(i = 0 ; i < iProcessesCount ; i++)
//	{
//		if(vpmoNodes[i].empty())
//		{
//			continue;
//		}
//		iTemp = (*vpmoNodes[i].rbegin()).first;
//		if(iTemp > iMaxMapIndex)
//		{
//			iMaxMapIndex = iTemp;
//		}
//	}
//	// set the global indices for the nodes
//	unsigned int iIndexBinSize = (unsigned
// int)pow(10,(floor(log((double)iMaxMapIndex)/log(10.0)) + 1.0));
// map<unsigned int,DislocationNetworkNode*>::iterator miNodes; 	for(i = 0 ;
// i <
// iProcessesCount ; i++)
//	{
//		for(miNodes = vpmoNodes[i].begin() ; miNodes !=
// vpmoNodes[i].end() ; miNodes++)
//		{
//			iGlobalIndex = i*iIndexBinSize + (*miNodes).first;
//			(*miNodes).second->GetDataPointer()->SetID(iGlobalIndex);
//		}
//	}
//	// now read and generate the arms
//	SeekRestartFileToData(fpFile);
//	DislocationNetworkNode* poGraphNeighbour = NULL;
//	DislocationNetworkArm* poEdge = NULL;
//	unsigned int iNeighbourDomain = 0;
//	unsigned int iNeighbourIndex = 0;
//	DislocationSegment oSegment;
//	list<DislocationNetworkArm*> lpoArms;
//	while(!feof(fpFile))
//	{
//		sRead = GetRealString(500,fpFile);
//		if(sRead.empty())
//		{
//			continue;
//		}
//		sscanf(sRead.c_str(),"%d,%d\t%*f\t%*f\t%*f\t%d\t%*d\n",&iDomainIndex,&iNodeIndex,&iNeighboursCount);
//		poGraphNode = vpmoNodes[iDomainIndex].operator[](iNodeIndex);
//		// skip the surface normal data
//		sRead = GetRealString(500,fpFile);
//		// read the arms data
//		for(i = 0 ; i < iNeighboursCount ; i++)
//		{
//			sRead = GetRealString(500,fpFile);
//			sscanf(sRead.c_str(),"%d,%d\t%lf\t%lf\t%lf\n",&iNeighbourDomain,&iNeighbourIndex,&dX,&dY,&dZ);
//			poGraphNeighbour =
// vpmoNodes[iNeighbourDomain].operator[](iNeighbourIndex);
//			// prepare segment data
//			oSegment.SetBurgersVector(Vector(dX,dY,dZ));
//			sRead = GetRealString(500,fpFile);
//			sscanf(sRead.c_str(),"%lf\t%lf\t%lf\n",&dX,&dY,&dZ);
//			oSegment.SetSlipPlaneNormal(Vector(dX,dY,dZ));
//			if(poGraphNode->ConnectTo(poGraphNeighbour,poEdge))
//			{
//				poEdge->SetData(oSegment);
//				lpoArms.push_back(poEdge);
//			}
//		}
//	}
//	fclose(fpFile);
//	m_oNetwork.Set(lpoNodes,lpoArms);
//}
// void LoopFilter::ExtractLoops(const string& sOutputFileName)
//{
//	// generate chains based on the network topology
//	GenerateDislocationChains();
//	// extract loops
//	list< DislocationChain* >::iterator liChains;
//	list< DislocationChain* > lpoLoops;
//	liChains = m_lpoSurvivingChains.begin();
//	while(liChains != m_lpoSurvivingChains.end())
//	{
//		if(IsSurfaceLoop((*liChains)))
//		{
//			FixLoopBurgersVector((*liChains));
//			lpoLoops.push_back((*liChains));
//			liChains = m_lpoSurvivingChains.erase(liChains);
//		}
//		else
//		{
//			liChains++;
//		}
//	}
//
//	FILE* fpFile = fopen(sOutputFileName.c_str(),"a");
//	fprintf(fpFile,"%d\n",(unsigned int)lpoLoops.size());
//	for(liChains = lpoLoops.begin() ; liChains != lpoLoops.end() ;
// liChains++)
//	{
//		WriteLoop((*liChains),fpFile);
//	}
//	fclose(fpFile);
//
//	for(liChains = lpoLoops.begin() ; liChains != lpoLoops.end() ;
// liChains++)
//	{
//		if((*liChains) != NULL)
//		{
//			delete (*liChains);
//		}
//	}
//}
// void LoopFilter::GenerateModifiedRestartFile(const string& sOutputFileName)
//{
//	list<DislocationNetworkNode*>* plpoNodes = m_oNetwork.GetNodes();
//	list<DislocationNetworkNode*>::iterator liNodes;
//	FILE* fpFile = fopen(sOutputFileName.c_str(),"w");
//	DislocationNetworkNode* poNetworkNode = NULL;
//	DislocationNode* poNode = NULL;
//	// write the header section
//	fprintf(fpFile,"dataFileVersion =   5\n");
//	fprintf(fpFile,"minCoordinates = [\n");
//	fprintf(fpFile,"%lf\n",m_oBox.GetXMin());
//	fprintf(fpFile,"%lf\n",m_oBox.GetYMin());
//	fprintf(fpFile,"%lf\n",m_oBox.GetZMin());
//	fprintf(fpFile,"]\n");
//	fprintf(fpFile,"maxCoordinates = [\n");
//	fprintf(fpFile,"%lf\n",m_oBox.GetXMax());
//	fprintf(fpFile,"%lf\n",m_oBox.GetYMax());
//	fprintf(fpFile,"%lf\n",m_oBox.GetZMax());
//	fprintf(fpFile,"]\n");
//
//	// write the node data
//	list< DislocationNetworkArm* >* plpoArms = NULL;
//	list< DislocationNetworkArm* >::iterator liArms;
//	Vector oTempVector;
//	// count the surviving nodes
//	list<DislocationChain*>::iterator liChains;
//	unsigned int i = 0;
//	unsigned int iSize = (unsigned int)plpoNodes->size();
//	vector<bool> vbHasNodeSurvived;
//	vbHasNodeSurvived.resize(iSize);
//	for(i = 0 ; i < iSize ; i++)
//	{
//		vbHasNodeSurvived[i] = false;
//	}
//
//	unsigned int iCount = 0;
//	unsigned int iIndex = 0;
//	for(liChains = m_lpoSurvivingChains.begin() ; liChains !=
// m_lpoSurvivingChains.end() ; liChains++)
//	{
//		(*liChains)->ResetIterator();
//		iSize = (*liChains)->GetSize();
//		for(i = 0 ; i < iSize ; i++)
//		{
//			iIndex = (*liChains)->GetCurrentNode()->GetID() - 1;
//			if(!vbHasNodeSurvived[iIndex])
//			{
//				vbHasNodeSurvived[iIndex] = true;
//				iCount = iCount + 1;
//			}
//			(*liChains)->IncrementIterator();
//		}
//	}
//	fprintf(fpFile,"nodeCount =   %d\n",iCount);
//	fprintf(fpFile,"nodalData =\n");
//	fprintf(fpFile,"#  Primary lines: node_tag, x, y, z, num_arms,
// constraint\n"); 	fprintf(fpFile,"#  Secondary lines: arm_tag, burgx,
// burgy, burgz, nx, ny, nz\n"); 	unsigned int iNeighbourIndex = 0;
//	UpdateNodeIndices(vbHasNodeSurvived);
//	DislocationNetworkNode* poNetworkNeighbourNode = NULL;
//	unsigned int iNeighboursCount = 0;
//	for(liNodes = plpoNodes->begin() ; liNodes != plpoNodes->end() ;
// liNodes++)
//	{
//		poNetworkNode = (*liNodes);
//		if(!vbHasNodeSurvived[poNetworkNode->GetID() - 1])
//		{
//			continue;
//		}
//		poNode = poNetworkNode->GetDataPointer();
//		plpoArms = poNetworkNode->GetEdges();
//		// count the surviving arms
//		iNeighboursCount = 0;
//		for(liArms = plpoArms->begin() ; liArms != plpoArms->end() ;
// liArms++)
//		{
//			poNetworkNeighbourNode =
//(*liArms)->GetOther(poNetworkNode);
//			if(vbHasNodeSurvived[poNetworkNeighbourNode->GetID() -
// 1])
//			{
//				iNeighboursCount = iNeighboursCount + 1;
//			}
//		}
//		fprintf(fpFile,"0,%d\t\t%25.20f\t\t%25.20f\t\t%25.20f\t\t%d\t\t%d\n",poNetworkNode->GetDataPointer()->GetID(),poNode->GetX(),poNode->GetY(),poNode->GetZ(),iNeighboursCount,poNode->GetCategory());
//		oTempVector = poNode->GetSurfaceNormal();
//		fprintf(fpFile,"\t\t\t%25.20f\t\t%25.20f\t\t%25.20f\n",oTempVector.GetX(),oTempVector.GetY(),oTempVector.GetZ());
//		for(liArms = plpoArms->begin() ; liArms != plpoArms->end() ;
// liArms++)
//		{
//			poNetworkNeighbourNode =
//(*liArms)->GetOther(poNetworkNode); 			iNeighbourIndex =
// poNetworkNeighbourNode->GetDataPointer()->GetID();
//			if(!vbHasNodeSurvived[poNetworkNeighbourNode->GetID() -
// 1])
//			{
//				continue;
//			}
//			oTempVector =
//(*liArms)->GetDataPointer()->GetBurgersVector();
// if((*liArms)->GetFirst() != poNetworkNode)
//			{
//				oTempVector = oTempVector*(-1.0);
//			}
//			fprintf(fpFile,"\t\t0,%d\t%25.20f\t\t%25.20f\t\t%25.20f\n",iNeighbourIndex,oTempVector.GetX(),oTempVector.GetY(),oTempVector.GetZ());
//			oTempVector =
//(*liArms)->GetDataPointer()->GetSlipPlaneNormal();
//			fprintf(fpFile,"\t\t\t%25.20f\t\t%25.20f\t\t%25.20f\n",oTempVector.GetX(),oTempVector.GetY(),oTempVector.GetZ());
//		}
//	}
//	fclose(fpFile);
//}
// void LoopFilter::WriteLoop(DislocationChain* poChain,FILE* fpFile)
//{
//	unsigned int iSize = poChain->GetSize();
//	Vector oBurgersVector = poChain->GetBurgersVector();
//	fprintf(fpFile,"%d\t\t%e,%e,%e\n",iSize,oBurgersVector.GetX(),oBurgersVector.GetY(),oBurgersVector.GetZ());
//	poChain->ResetIterator();
//	unsigned int i = 0;
//	DislocationNode* poCurrentNode = NULL;
//	DislocationNetworkNode* poCurrentNetworkNode = NULL;
//	DislocationNetworkNode* poNextNetworkNode = NULL;
//	DislocationNetworkArm* poNetworkEdge = NULL;
//	Vector oPlaneNormal;
//	Vector oSurfaceNormal;
//	for(i = 0 ; i < iSize - 1 ; i++)
//	{
//		poCurrentNetworkNode = poChain->GetCurrentNode();
//		poNextNetworkNode = poChain->GetNextNode();
//		poCurrentNode = poCurrentNetworkNode->GetDataPointer();
//		if(poCurrentNetworkNode->IsConnected(poNextNetworkNode,poNetworkEdge))
//		{
//			oPlaneNormal =
// poNetworkEdge->GetDataPointer()->GetSlipPlaneNormal();
//		}
//		else
//		{
//			printf("error: inconsistent loop\n");
//			return;
//		}
//		oSurfaceNormal = poCurrentNode->GetSurfaceNormal();
//		fprintf(fpFile,"%e,%e,%e,%e,%e,%e,%e,%e,%e\n",poCurrentNode->GetX(),poCurrentNode->GetY(),poCurrentNode->GetZ(),oPlaneNormal.GetX(),oPlaneNormal.GetY(),oPlaneNormal.GetZ(),oSurfaceNormal.GetX(),oSurfaceNormal.GetY(),oSurfaceNormal.GetZ());
//		poChain->IncrementIterator();
//	}
//}
// bool LoopFilter::IsSurfaceLoop(DislocationChain* poChain)
//{
//	if(!poChain->IsLoop())
//	{
//		return false;
//	}
//
//	unsigned int iSize = poChain->GetSize();
//	poChain->ResetIterator();
//	unsigned int i = 0;
//	for(i = 0 ; i < iSize - 1 ; i++)
//	{
//		if(poChain->GetCurrentNode()->GetDataPointer()->GetCategory() !=
// SurfaceNode)
//		{
//			return false;
//		}
//		poChain->IncrementIterator();
//	}
//	return true;
//}
// void LoopFilter::FixLoopBurgersVector(DislocationChain* poChain)
//{
//	poChain->ResetIterator();
//	DislocationNetworkNode* poNode1 = poChain->GetCurrentNode();
//	DislocationNetworkNode* poNode2 = poChain->GetNextNode();
//	DislocationNetworkArm* poArm = NULL;
//	Vector oBurgersVector;
//	if(poNode1->IsConnected(poNode2,poArm))
//	{
//		oBurgersVector = poArm->GetDataPointer()->GetBurgersVector();
//		if(poNode1 != poArm->GetFirst())
//		{
//			oBurgersVector = oBurgersVector*(-1.0);
//		}
//	}
//	else
//	{
//		printf("error: first loop arm missing\n");
//	}
//	poChain->OverrideBurgersVector(oBurgersVector);
//}
// void LoopFilter::SetBoundingBox(const string& sInputFileName)
//{
//	FILE* fpFile = fopen(sInputFileName.c_str(),"r");
//	fseek(fpFile,0,SEEK_SET);
//	string sRead = "";
//	double dTemp = 0.0;
//	// seek till the first node data
//	while(!feof(fpFile))
//	{
//		sRead = GetRealString(500,fpFile);
//		if(sRead.compare("minCoordinates = [") == 0)
//		{
//			sRead = GetRealString(500,fpFile);
//			sscanf(sRead.c_str(),"%lf",&dTemp);
//			m_oBox.SetXMin(dTemp);
//
//			sRead = GetRealString(500,fpFile);
//			sscanf(sRead.c_str(),"%lf",&dTemp);
//			m_oBox.SetYMin(dTemp);
//
//			sRead = GetRealString(500,fpFile);
//			sscanf(sRead.c_str(),"%lf",&dTemp);
//			m_oBox.SetZMin(dTemp);
//		}
//		if(sRead.compare("maxCoordinates = [") == 0)
//		{
//			sRead = GetRealString(500,fpFile);
//			sscanf(sRead.c_str(),"%lf",&dTemp);
//			m_oBox.SetXMax(dTemp);
//
//			sRead = GetRealString(500,fpFile);
//			sscanf(sRead.c_str(),"%lf",&dTemp);
//			m_oBox.SetYMax(dTemp);
//
//			sRead = GetRealString(500,fpFile);
//			sscanf(sRead.c_str(),"%lf",&dTemp);
//			m_oBox.SetZMax(dTemp);
//		}
//	}
//	fclose(fpFile);
//}
// void LoopFilter::SeekRestartFileToData(FILE* fpFile)
//{
//	fseek(fpFile,0,SEEK_SET);
//	string sRead = "";
//	// seek till the first node data
//	while(!feof(fpFile))
//	{
//		sRead = GetRealString(500,fpFile);
//		if((sRead.compare("nodalData =") == 0) ||
//(sRead.compare("nodalData = ") == 0))
//		{
//			break;
//		}
//	}
//}
// void LoopFilter::GenerateDislocationChains()
//{
//	// generate chains based on the network topology
//	list< GraphChain<DislocationNode,DislocationSegment>* > lpoGraphChains =
// m_oNetwork.GenerateGraphChains(); 	list<
// GraphChain<DislocationNode,DislocationSegment>* >::iterator liGraphChains;
//	for(liGraphChains = lpoGraphChains.begin() ; liGraphChains !=
// lpoGraphChains.end() ; liGraphChains++)
//	{
//		m_lpoSurvivingChains.push_back(new
// DislocationChain((*liGraphChains)));
//	}
//
//	// split chains based on the node type
//	DislocationSystem::NodeCategory eInitialCategory =
// DislocationSystem::UnconstrainedNode; 	DislocationSystem::NodeCategory
// eCategory = DislocationSystem::UnconstrainedNode; 	DislocationChain*
// poNewChain = NULL; 	unsigned int i = 0; 	unsigned int iChainLength = 0;
//	list<DislocationChain*>::iterator liChains;
//	DislocationChain* poChain = NULL;
//	// split chains based on changes in node types
//	for(liChains = m_lpoSurvivingChains.begin() ; liChains !=
// m_lpoSurvivingChains.end() ; liChains++)
//	{
//		poChain = (*liChains);
//		eInitialCategory =
//(DislocationSystem::NodeCategory)(poChain->GetCurrentNode()->GetDataPointer()->GetCategory());
//		iChainLength = poChain->GetSize();
//		poChain->ResetIterator();
//		for(i = 0 ; i < iChainLength ; i++)
//		{
//			eCategory =
//(DislocationSystem::NodeCategory)(poChain->GetCurrentNode()->GetDataPointer()->GetCategory());
//			if(eCategory != eInitialCategory)
//			{
//				if((eInitialCategory ==
// DislocationSystem::SurfaceNode) || (eCategory ==
// DislocationSystem::SurfaceNode))
//				{
//					poNewChain = poChain->Split();
//					if(poNewChain == NULL)
//					{
//						poChain->IncrementIterator();
//						continue;
//					}
//					else
//					{
//						m_lpoSurvivingChains.push_back(poNewChain);
//						break;
//					}
//				}
//			}
//			poChain->IncrementIterator();
//		}
//	}
//	Vector oInitialVector;
//	Vector oCurrentVector;
//	DislocationNetworkNode* poCurrentNode = NULL;
//	DislocationNetworkNode* poNextNode = NULL;
//	DislocationNetworkArm* poEdge = NULL;
//	double dTolerance = 1E-6;
//	/*// split chains based on the node slip plane
//	for(liChains = m_lpoSurvivingChains.begin() ; liChains !=
// m_lpoSurvivingChains.end() ; liChains++)
//	{
//		poChain = (*liChains);
//		oInitialVector = poChain->GetSlipPlaneNormal();
//		iChainLength = poChain->GetSize();
//		poChain->ResetIterator();
//		for(i = 0 ; i < iChainLength ; i++)
//		{
//			poCurrentNode = poChain->GetCurrentNode();
//			poNextNode = poChain->GetNextNode();
//			if(poNextNode == NULL)
//			{
//				break;
//			}
//			if(poCurrentNode->IsConnected(poNextNode,poEdge))
//			{
//				oCurrentVector =
// poEdge->GetDataPointer()->GetSlipPlaneNormal();
//			}
//			else
//			{
//				break;		// this should be unreachable
//			}
//			// check to see if the normals are the same
//			if((oCurrentVector^oInitialVector).Length() >
// dTolerance)
//			{
//				poNewChain = poChain->Split();
//				if(poNewChain == NULL)
//				{
//					poChain->IncrementIterator();
//					continue;
//				}
//				else
//				{
//					m_lpoSurvivingChains.push_back(poNewChain);
//					break;
//				}
//			}
//			poChain->IncrementIterator();
//		}
//	}*/
//	// split chains based on the node Burgers vector
//	for(liChains = m_lpoSurvivingChains.begin() ; liChains !=
// m_lpoSurvivingChains.end() ; liChains++)
//	{
//		poChain = (*liChains);
//		oInitialVector = poChain->GetBurgersVector();
//		iChainLength = poChain->GetSize();
//		poChain->ResetIterator();
//		for(i = 0 ; i < iChainLength ; i++)
//		{
//			poCurrentNode = poChain->GetCurrentNode();
//			poNextNode = poChain->GetNextNode();
//			if(poNextNode == NULL)
//			{
//				break;
//			}
//			if(poCurrentNode->IsConnected(poNextNode,poEdge))
//			{
//				oCurrentVector =
// poEdge->GetDataPointer()->GetBurgersVector();
//			}
//			else
//			{
//				break;		// this should be unreachable
//			}
//			// check to see if the burgers vectors are the same
//			if((oCurrentVector^oInitialVector).Length() >
// dTolerance)
//			{
//				poNewChain = poChain->Split();
//				if(poNewChain == NULL)
//				{
//					poChain->IncrementIterator();
//					continue;
//				}
//				else
//				{
//					m_lpoSurvivingChains.push_back(poNewChain);
//					break;
//				}
//			}
//			poChain->IncrementIterator();
//		}
//	}
//
//	Fire();
//}
// void LoopFilter::SeparateInternalChains()
//{
//	// split chains based on the node type
//	unsigned int i = 0;
//	unsigned int iChainLength = 0;
//	list<DislocationChain*>::iterator liChains;
//	DislocationChain* poChain = NULL;
//	bool bRemoveChain = false;
//	liChains = m_lpoSurvivingChains.begin();
//	while(liChains != m_lpoSurvivingChains.end())
//	{
//		poChain = (*liChains);
//		iChainLength = poChain->GetSize();
//		poChain->ResetIterator();
//		bRemoveChain = true;
//		for(i = 0 ; i < iChainLength ; i++)
//		{
//			if(poChain->GetCurrentNode()->GetDataPointer()->GetCategory()
//== DislocationSystem::SurfaceNode)
//			{
//				bRemoveChain = false;
//				break;
//			}
//			poChain->IncrementIterator();
//		}
//		if(bRemoveChain)
//		{
//			(*liChains)->Collapse();
//			delete (*liChains);
//			liChains = m_lpoSurvivingChains.erase(liChains);
//		}
//		else
//		{
//			liChains++;
//		}
//	}
//	Fire();
//}
// void LoopFilter::BuildDualGraph()
//{
//	list<DislocationChain*>::iterator liChains;
//	DualGraphNode* poGraphNode = NULL;
//	unsigned int i = 0;
//	for(liChains = m_lpoSurvivingChains.begin() ; liChains !=
// m_lpoSurvivingChains.end(); liChains++)
//	{
//		poGraphNode = new DualGraphNode;
//		poGraphNode->Set((*liChains));
//		i = i + 1;
//		poGraphNode->SetID(i);
//		m_oDualGraph.AddNode(poGraphNode);
//	}
//
//	list<DualGraphNode*>* plpoGraphNodes = m_oDualGraph.GetNodes();
//	list<DualGraphNode*>::iterator liOuterNodes;
//	list<DualGraphNode*>::iterator liInnerNodes;
//	DislocationNetworkNode* poLocalFirstNode = NULL;
//	DislocationNetworkNode* poLocalLastNode = NULL;
//	DislocationNetworkNode* poRemoteFirstNode = NULL;
//	DislocationNetworkNode* poRemoteLastNode = NULL;
//	DislocationNetworkNode* poConnectionNode = NULL;
//	DualGraphNode* poLocalChain = NULL;
//	DualGraphNode* poRemoteChain = NULL;
//	DualGraphEdge* poEdge = NULL;
//	for(liOuterNodes = plpoGraphNodes->begin() ; liOuterNodes !=
// plpoGraphNodes->end() ; liOuterNodes++)
//	{
//		// get the end points of this chain
//		poLocalChain = (*liOuterNodes);
//		poLocalFirstNode = poLocalChain->GetData()->GetFirst();
//		poLocalLastNode = poLocalChain->GetData()->GetLast();
//		// loop over all of the other chains
//		for(liInnerNodes = plpoGraphNodes->begin() ; liInnerNodes !=
// plpoGraphNodes->end() ; liInnerNodes++)
//		{
//			poRemoteChain = (*liInnerNodes);
//			// skip the chain in question
//			if(poLocalChain == poRemoteChain)
//			{
//				continue;
//			}
//			poRemoteFirstNode =
// poRemoteChain->GetData()->GetFirst(); 			poRemoteLastNode
// =
// poRemoteChain->GetData()->GetLast(); 			poConnectionNode
// =
// NULL; 			if((poLocalFirstNode == poRemoteFirstNode) ||
// (poLocalFirstNode ==
// poRemoteLastNode))
//			{
//				poConnectionNode = poLocalFirstNode;
//			}
//			else if((poLocalLastNode == poRemoteFirstNode) ||
//(poLocalLastNode == poRemoteLastNode))
//			{
//				poConnectionNode = poLocalLastNode;
//			}
//
//			if(poConnectionNode != NULL)
//			{
//				if(!poLocalChain->IsConnected(poRemoteChain,poEdge))
//				{
//					// if the chains are not connected,
//connect
// them with an undirected edge
// poLocalChain->ConnectTo(poRemoteChain,poEdge);
//					poEdge->SetData(poConnectionNode);
//					m_oDualGraph.AddEdge(poEdge);
//				}
//			}
//		}
//	}
//}
// void LoopFilter::ExtractChains()
//{
//	list<DualGraphNode*>* plpoGraphNodes;
//	list<DualGraphNode*>::iterator liInnerNodes;
//	list<DualGraphNode*>::iterator liOuterNodes;
//	DislocationNetworkNode* poConnectionNode = NULL;
//	DualGraphNode* poLocalChain = NULL;
//	DualGraphNode* poRemoteChain = NULL;
//	Vector oTargetBurgersVector;
//	Vector oTestBurgersVector;
//	list<DualGraphNode*> lpoPath;
//	list<DualGraphNode*> lpoTempPath;
//	while(m_oDualGraph.GetNodesCount() > 0)
//	{
//		plpoGraphNodes = m_oDualGraph.GetNodes();
//		lpoPath.clear();
//		for(liOuterNodes = plpoGraphNodes->begin() ; liOuterNodes !=
// plpoGraphNodes->end() ; liOuterNodes++)
//		{
//			poLocalChain = (*liOuterNodes);
//			oTargetBurgersVector =
// poLocalChain->GetData()->GetBurgersVector();
// for(liInnerNodes
// = plpoGraphNodes->begin() ; liInnerNodes != plpoGraphNodes->end() ;
// liInnerNodes++)
//			{
//				poRemoteChain = (*liInnerNodes);
//				if(poLocalChain == poRemoteChain)
//				{
//					continue;
//				}
//				oTestBurgersVector =
// poRemoteChain->GetData()->GetBurgersVector();
//				if(oTargetBurgersVector.IsSimilar(oTestBurgersVector))
//				{
//					// try to get a path between the two
//					lpoTempPath =
// m_oDualGraph.GetShortestPath(poLocalChain,poRemoteChain);
//					if(lpoTempPath.size() == 0)
//					{
//						continue;
//					}
//					if(lpoTempPath.size() > lpoPath.size())
//					{
//						lpoPath = lpoTempPath;
//					}
//				}
//			}
//		}
//		// now if the longest path is longer than 1 step, which means
// that the chain is actually connected
//		// to other chains, take it, if not, take the chain to be the
// path (in this case, the path
//		// will just consist of 1 chain
//		if(lpoPath.size() > 1)
//		{
//			printf("%p longest path is %d steps
// long\n",poLocalChain,lpoPath.size());
// ExtractPath(&lpoPath);
// printf("%d nodes remaining\n",m_oDualGraph.GetNodesCount());
// break;
//		}
//		else
//		{
//			printf("%p is an isolated
// chain\n",poLocalChain,lpoPath.size()); 			break;
//		}
//	}
//	Fire();
//}
// void LoopFilter::ExtractPath(list<DualGraphNode*>* plpoPath)
//{
//	list<DualGraphNode*>::iterator liChains = plpoPath->begin();
//	DualGraphNode* poCurrentChain = (*liChains);
//	DualGraphNode* poNextChain = NULL;
//	// get the path's Burgers vector
//	Vector oTargetBurgersVector =
// poCurrentChain->GetData()->GetBurgersVector(); 	Vector
// oTempBurgersVector;
//
//	int iCode = 0;
//	liChains++;
//	DislocationChain* poNewChain = NULL;
//	DualGraphNode* poNewGraphNode = NULL;
//	DualGraphEdge* poGraphEdge = NULL;
//	list<DualGraphNode*> lpoNeighbours;
//	list<DualGraphNode*>::iterator liNeighbours;
//	while(liChains != plpoPath->end())
//	{
//		poNextChain = (*liChains);
//		iCode =
// EncodeConnection(poCurrentChain->GetData(),poNextChain->GetData());
//		printf("code is %d\n",iCode);
//		if(iCode == 8)
//		{
//			poNextChain->GetData()->Reverse();
//		}
//
//		oTempBurgersVector = poNextChain->GetData()->GetBurgersVector();
//		// if the burgers vectors are not the same, dissociate this
// chain into two chains
// if(!oTargetBurgersVector.IsSame(oTempBurgersVector))
//		{
//			poNewChain =
// poNextChain->GetData()->Dissociate(oTargetBurgersVector);
//			// now we need to do 3 things
//			// 1. add the new chains to the surviving chains list,
// no need to add it to the
//			// network because the chains don't own their nodes
//			m_lpoSurvivingChains.push_back(poNewChain);
//			// 2. add this chain to the dual graph, now the number
// of nodes in the dual graph will
//			// increase by one for each added chain, we need this so
// that we can include this
//			// chain in the future chain extractions
//			poNewGraphNode = new DualGraphNode;
//			poNewGraphNode->Set(poNewChain);
//			// here, we need to get an ID, since we can't do that
// now, we will update the IDs at the
//			// end of the loop for all the newly added chains, also
// we can't do this now because
//			// the nodes of this path will be removed from the dual
// graph anyway 			poNewGraphNode->SetID(0);
// m_oDualGraph.AddNode(poNewGraphNode);
//			// 3. connect this chain properly in the dual graph, for
// that, we need to get all
//			// of the neighbours of the original chain and connect
// them to this chain as well 			lpoNeighbours =
// poNextChain->GetNeighbours(); 			for(liNeighbours =
// lpoNeighbours.begin() ;
// liNeighbours != lpoNeighbours.end() ; liNeighbours++)
//			{
//				// get the connected edge
//				if(!poNextChain->IsConnected((*liNeighbours),poGraphEdge))
//				{
//					// this condition is unreachable
//					printf("error: chain extraction
// connectivity problem - edge not found\n");
//				}
//				// now switch the edge's connection from the old
// node to the new one, this step is
//				// in fact much better the disconnecting and
// reconnecting because
//				// A. we don't have to remove the edge
//connecting
// the old node, this step involves a list search O(N)
//				// B. we don't have to create a new edge
//				// C. we don't have to worry about copying the
// data from the old edge to the
//				// new edge
//				// 4. the node gets automatically isolated at
//the
// end of the process so it is
//				// safe to remove it
//				if(!poGraphEdge->SwitchConnection(poNextChain,poNewGraphNode))
//				{
//					// this should be unreachable
//					printf("error: chain extraction
// connectivity problem - node not found\n");
//				}
//			}
//			// just a quick check, make sure that the old node has
// been isolated by now 			if(!poNextChain->IsIsolated())
//			{
//				printf("error: chain extraction node isolation
// problem\n");
//			}
//			printf("dissociated\n");
//			fflush(NULL);
//			// at this point, we have a new chain, and its
// corresponding dual graph node, in
//			// place, we can proceed to other chains in the path
//		}
//		// move on to the next chain
//		poCurrentChain = poNextChain;
//		liChains++;
//	}
//	// now we need to do 3 things
//	// 1. write the extracted path
//	WritePath(plpoPath);
//	// 2. remove all the path nodes (dislocation chains) from the dual graph
//	liChains = plpoPath->begin();
//	while(liChains != plpoPath->end())
//	{
//		m_oDualGraph.RemoveNode((*liChains),true);
//		liChains++;
//	}
//	// 3. reassign dual graph node IDs
//	m_oDualGraph.GenerateSequentialNodeIDs();
//}
// int LoopFilter::EncodeConnection(DislocationChain*
// poCurrentChain,DislocationChain* poNextChain) const
//{
//	// 0 - no connection
//	// 1 - CF->NF		2 - CF->NL		4 - CL->NF 8
//- CL->NL
//	// 3 - 1&2		5 - 1&3		9 - 1&4		6 - 2&3
// 10 - 2&4		12 - 3&4
//	// 7 - 1&2&3		11 - 1&2&4		13 - 1&3&4
//14
//- 2&3&4
//	// 15 - 1&2&3&4
//	int iCode = 0;
//	DislocationNetworkNode* poCurrentFirstNode = poCurrentChain->GetFirst();
//	DislocationNetworkNode* poCurrentLastNode = poCurrentChain->GetLast();
//	DislocationNetworkNode* poNextFirstNode = poNextChain->GetFirst();
//	DislocationNetworkNode* poNextLastNode = poNextChain->GetLast();
//	if(poCurrentFirstNode == poNextFirstNode)
//	{
//		iCode |= 0x01;
//	}
//	if(poCurrentFirstNode == poNextLastNode)
//	{
//		iCode |= 0x02;
//	}
//	if(poCurrentLastNode == poNextFirstNode)
//	{
//		iCode |= 0x04;
//	}
//	if(poCurrentLastNode == poNextLastNode)
//	{
//		iCode |= 0x08;
//	}
//	return iCode;
//}
// void LoopFilter::WritePath(list<DualGraphNode*>* plpoPath) const
//{
//	FILE* fpFile = fopen("paths.txt","a");
//	list<DualGraphNode*>::iterator liChains;
//	DislocationChain* poChain = NULL;
//	unsigned int i = 0;
//	unsigned int iChainSize = 0;
//	DislocationNode* poNode;
//	Vector oSurfaceNormal;
//	for(liChains = plpoPath->begin() ; liChains != plpoPath->end() ;
// liChains++)
//	{
//		poChain = (*liChains)->GetData();
//		iChainSize = poChain->GetSize();
//		fprintf(fpFile,"size %d\n",iChainSize);
//		for(i = 0 ; i < iChainSize ; i++)
//		{
//			poNode = poChain->GetCurrentNode()->GetDataPointer();
//			oSurfaceNormal = poNode->GetSurfaceNormal();
//			fprintf(fpFile,"%25.20f\t%25.20f\t%25.20f\t%e\t%e\t%e\t%d\n",poNode->GetX(),poNode->GetY(),poNode->GetZ(),oSurfaceNormal.GetX(),oSurfaceNormal.GetY(),oSurfaceNormal.GetZ(),poNode->GetCategory());
//			poChain->IncrementIterator();
//		}
//	}
//	fclose(fpFile);
//}
// void LoopFilter::UpdateNodeIndices(const vector<bool>& vbHasNodeSurvived)
//{
//	list<DislocationNetworkNode*>* plpoNodes = m_oNetwork.GetNodes();
//	list<DislocationNetworkNode*>::iterator liNodes;
//	unsigned int iCurrentIndex = 1;
//	for(liNodes = plpoNodes->begin() ; liNodes != plpoNodes->end() ;
// liNodes++)
//	{
//		if(vbHasNodeSurvived[(*liNodes)->GetID() - 1])
//		{
//			(*liNodes)->GetDataPointer()->SetID(iCurrentIndex);
//			iCurrentIndex = iCurrentIndex + 1;
//		}
//	}
//}
//
//
// void LoopFilter::Fire()
//{
//	// coverage check
//	list<DislocationNetworkNode*>* plpoNodes = m_oNetwork.GetNodes();
//	unsigned int iNodesCount = (unsigned int)plpoNodes->size();
//	list<DislocationNetworkNode*>::iterator liNodes;
//	vector<bool> vbIsUsed;
//	vbIsUsed.resize(iNodesCount);
//	printf("nodes count %d\n",iNodesCount);
//	printf("chains count %d\n",(unsigned int)m_lpoSurvivingChains.size());
//	for(liNodes = plpoNodes->begin() ; liNodes != plpoNodes->end() ;
// liNodes++)
//	{
//		vbIsUsed[(*liNodes)->GetID() - 1] = false;
//	}
//
//	unsigned int iChainSize = 0;
//	DislocationChain* poChain = NULL;
//	list< DislocationChain* >::iterator liChains;
//	unsigned int i = 0;
//	for(liChains = m_lpoSurvivingChains.begin() ; liChains !=
// m_lpoSurvivingChains.end() ; liChains++)
//	{
//		poChain = (*liChains);
//		iChainSize = poChain->GetSize();
//		poChain->ResetIterator();
//		for(i = 0 ; i < iChainSize ; i++)
//		{
//			vbIsUsed[poChain->GetCurrentNode()->GetID() - 1] = true;
//			poChain->IncrementIterator();
//		}
//	}
//
//	for(i = 0 ; i < iNodesCount ; i++)
//	{
//		if(!vbIsUsed[i])
//		{
//			//printf("node %d was not used !!!\n",i);
//		}
//	}
//}
//
//
//// int main(int argc,char** argv)
//// {
//// 	if(argc < 3)
//// 	{
//// 		printf("incorrent arguments count\n");
//// 		fflush(NULL);
//// 		return 1;
//// 	}
//// 	string sInputFileName = argv[1];
//// 	unsigned int iProcessesCount = 0;
//// 	sscanf(argv[2],"%d",&iProcessesCount);
//// 	LoopFilter* poLoopFilter = LoopFilter::CreateInstance();
//// 	poLoopFilter->Set(sInputFileName,iProcessesCount);
//// 	poLoopFilter->GenerateDislocationChains();
//// 	poLoopFilter->SeparateInternalChains();
//// 	poLoopFilter->BuildDualGraph();
//// // 	poLoopFilter->ExtractLoops("out_loops.txt");
////  	poLoopFilter->GenerateModifiedRestartFile("new_restart.txt");
////  	poLoopFilter->ExtractChains();
////  	//poLoopFilter->GenerateModifiedRestartFile("new_restart2.txt");
//// 	delete poLoopFilter;
//// 	return 0;
//// }
//
//
