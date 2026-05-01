#include "MainDataStructure.h"
#include "Tools.h"
#include "vector"
#include "map"
#include "mpi.h"
#include "stdio.h"
#include "cmath"

using namespace std;
using namespace SupportSystem;

MainDataStructure* MainDataStructure::m_poInstance = NULL;
MainDataStructure* MainDataStructure::GetInstance()
{
	if(m_poInstance == NULL)
	{
		m_poInstance = new MainDataStructure;
	}
	return m_poInstance;
}
MainDataStructure::~MainDataStructure()
{
	Reset();
}
MainDataStructure::MainDataStructure()
{
	Initialize();
}
void MainDataStructure::Reset()
{
	list<DislocationNode*>::iterator liNodes;
	for(liNodes = m_lpoLocalNodes.begin() ; liNodes != m_lpoLocalNodes.end() ; liNodes++)
	{
		if((*liNodes) != NULL)
		{
			delete (*liNodes);
		}
	}
	
	for(liNodes = m_lpoRemoteNodes.begin() ; liNodes != m_lpoRemoteNodes.end() ; liNodes++)
	{
		if((*liNodes) != NULL)
		{
			delete (*liNodes);
		}
	}
	
	list<AxisAlignedBoundingBox*>::iterator liDomains;
	unsigned int iSize = (unsigned int)m_vpoDomains.size();
	unsigned int i = 0;
	for(i = 0 ; i < iSize ; i++)
	{
		if(m_vpoDomains[i] != NULL)
		{
			delete m_vpoDomains[i];
		}
	}
	
	iSize = (unsigned int)m_vvpoNodes.size();
	for(i = 0 ; i < iSize ; i++)
	{
		m_vvpoNodes[i].clear();
	}
	Initialize();
}
void MainDataStructure::SetDomainData(const unsigned int& iDomainsCount,const unsigned int& iDomainID)
{
	m_iDomainsCount = iDomainsCount;
	m_iDomainID = iDomainID;
	m_vvpoNodes.resize(iDomainsCount);
}
void MainDataStructure::ReadInput(const string& sFileName)
{
	FILE* fpFile = fopen(sFileName.c_str(),"r");
	if(fpFile == NULL)
	{
		return;
	}
	string sRead;
	sRead = GetRealString(512,fpFile);
	sscanf(sRead.c_str(),"%d\t%d\t%d\n",&m_iXDomainsCount,&m_iYDomainsCount,&m_iZDomainsCount);
	double dTemp1 = 0.0;
	double dTemp2 = 0.0;
	double dTemp3 = 0.0;
	double dTemp4 = 0.0;
	double dTemp5 = 0.0;
	double dTemp6 = 0.0;
	sRead = GetRealString(512,fpFile);
	sscanf(sRead.c_str(),"%lf\t%lf\t%lf\n",&dTemp1,&dTemp2,&dTemp3);
	m_oProblemBox.SetXMin(dTemp1);
	m_oProblemBox.SetYMin(dTemp2);
	m_oProblemBox.SetZMin(dTemp3);
	sRead = GetRealString(512,fpFile);
	sscanf(sRead.c_str(),"%lf\t%lf\t%lf\n",&dTemp1,&dTemp2,&dTemp3);
	m_oProblemBox.SetXMax(dTemp1);
	m_oProblemBox.SetYMax(dTemp2);
	m_oProblemBox.SetZMax(dTemp3);
	SetDomainBoxes();
	unsigned int iNodesCount = 0;
	sRead = GetRealString(512,fpFile);
	sscanf(sRead.c_str(),"%d\t%d\n",&iNodesCount,&m_iOutputFrequency);
	unsigned int i = 0;
	unsigned int j = 0;
	unsigned int iDomainID = 0;
	unsigned int iNodeID = 0;
	unsigned int iNeighboursCount = 0;
	unsigned int iCategory = 0;
	unsigned int iMaxDomainID = 0;
	DislocationNode* poNode = NULL;
	DislocationSegment* poArm = NULL;
	for(i = 0 ; i < iNodesCount ; i++)
	{
		sRead = GetRealString(512,fpFile);
		sscanf(sRead.c_str(),"%d,%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\n",&iDomainID,&iNodeID,&dTemp1,&dTemp2,&dTemp3,&dTemp4,&dTemp5,&dTemp6,&iNeighboursCount,&iCategory);
		poNode = new DislocationNode(dTemp1,dTemp2,dTemp3);
		poNode->SetID(iNodeID);
		poNode->SetDomainID(iDomainID);
		poNode->SetSurfaceNormal(Vector(dTemp4,dTemp5,dTemp6));
		poNode->SetCategory(iCategory);
		if(iMaxDomainID < iDomainID)
		{
			iMaxDomainID = iDomainID;
		}
		for(j = 0 ; j < iNeighboursCount ; j++)
		{
			sRead = GetRealString(512,fpFile);
			sscanf(sRead.c_str(),"\t%d,%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",&iDomainID,&iNodeID,&dTemp1,&dTemp2,&dTemp3,&dTemp4,&dTemp5,&dTemp6);
			poArm = new DislocationSegment;
			poArm->SetEndDomainID(iDomainID);
			poArm->SetEndNodeID(iNodeID);
			poArm->SetBurgersVector(Vector(dTemp1,dTemp2,dTemp3));
			poArm->SetSlipPlaneNormal(Vector(dTemp4,dTemp5,dTemp6));
			poNode->AddArm(poArm);
		}
		m_lpoLocalNodes.push_back(poNode);
	}
	SetLocalNodes(iMaxDomainID);
	fclose(fpFile);
}
void MainDataStructure::WriteOutput(const string& sFileName)
{
	// get total nodes count
	unsigned int iLocalNodesCount = (unsigned int)m_lpoLocalNodes.size();
	unsigned int iGlobalNodesCount = 0;
	MPI_Reduce(&iLocalNodesCount,&iGlobalNodesCount,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	FILE* fpFile = NULL;
	if(m_iDomainID == 0)
	{
		fpFile = fopen(sFileName.c_str(),"w");
		fprintf(fpFile,"%d\t%d\t%d\n",m_iXDomainsCount,m_iYDomainsCount,m_iZDomainsCount);
		fprintf(fpFile,"%lf\t%lf\t%lf\n",m_oProblemBox.GetXMin(),m_oProblemBox.GetYMin(),m_oProblemBox.GetZMin());
		fprintf(fpFile,"%lf\t%lf\t%lf\n",m_oProblemBox.GetXMax(),m_oProblemBox.GetYMax(),m_oProblemBox.GetZMax());
		fprintf(fpFile,"%d\t%d\n",iGlobalNodesCount,m_iOutputFrequency);
		fprintf(fpFile,"data_begin\n");
		fclose(fpFile);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	unsigned int i = 0;
	list<DislocationNode*>::iterator liNodes;
	for(i = 0 ; i < m_iDomainsCount ; i++)
	{
		if(m_iDomainID == i)
		{
			fpFile = fopen(sFileName.c_str(),"a");
			for(liNodes = m_lpoLocalNodes.begin() ; liNodes != m_lpoLocalNodes.end() ; liNodes++)
			{
				(*liNodes)->Write(fpFile);
			}
			fclose(fpFile);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}
void MainDataStructure::CommunicateGhosts()
{
	ClearGhostNodes();
	vector< list<double> > vldData;
	
	PackGhostNodes(&vldData);
	MPI_Barrier(MPI_COMM_WORLD);
	
	// send the data buffers lengths
	// issue length receives
	unsigned int* piBufferSizes = new unsigned int[m_iDomainsCount];
	MPI_Request* poIncomingRequests = new MPI_Request[m_iDomainsCount];
	unsigned int i = 0;
	for(i = 0 ; i < m_iDomainsCount ; i++)
	{
		if(i == m_iDomainID)
		{
			piBufferSizes[i] = 0;
			continue;
		}
		MPI_Irecv(&piBufferSizes[i],1,MPI_INT,i,GhostMessageLengthTag,MPI_COMM_WORLD,&poIncomingRequests[i]);
	}
	// issue length sends
	unsigned int iTemp = 0;
	MPI_Request* poOutgoingRequests = new MPI_Request[m_iDomainsCount];
	for(i = 0 ; i < m_iDomainsCount ; i++)
	{
		if(i == m_iDomainID)
		{
			continue;
		}
		iTemp = (unsigned int)vldData[i].size();
		MPI_Isend(&iTemp,1,MPI_INT,i,GhostMessageLengthTag,MPI_COMM_WORLD,&poOutgoingRequests[i]);
	}
	// wait for the communications to end
	MPI_Status* poIncomingStatuses = new MPI_Status[m_iDomainsCount];
	MPI_Status* poOutgoingStatuses = new MPI_Status[m_iDomainsCount];
	for(i = 0 ; i < m_iDomainsCount ; i++)
	{
		if(i == m_iDomainID)
		{
			continue;
		}
		MPI_Wait(&poIncomingRequests[i],&poIncomingStatuses[i]);
		MPI_Wait(&poOutgoingRequests[i],&poOutgoingStatuses[i]);
	}
    MPI_Barrier(MPI_COMM_WORLD);
    // send the data buffers
 	// issue data receives
 	vector<double*> vpdIncomingBuffers;
 	vpdIncomingBuffers.resize(m_iDomainsCount);
	for(i = 0 ; i < m_iDomainsCount ; i++)
	{
		if(i == m_iDomainID)
		{
			vpdIncomingBuffers[i] = NULL;
			continue;
		}
		vpdIncomingBuffers[i] = new double[piBufferSizes[i]];
		MPI_Irecv(vpdIncomingBuffers[i],piBufferSizes[i],MPI_DOUBLE,i,GhostMessageDataTag,MPI_COMM_WORLD,&poIncomingRequests[i]);
	}
	// issue data sends
	vector<double*> vpdOutgoingBuffers;
 	vpdOutgoingBuffers.resize(m_iDomainsCount);
 	list<double>::iterator liData;
 	unsigned int j = 0;
	for(i = 0 ; i < m_iDomainsCount ; i++)
	{
		if(i == m_iDomainID)
		{
			vldData[i].clear();
			vpdOutgoingBuffers[i] = NULL;
			continue;
		}
		iTemp = (unsigned int)vldData[i].size();
		vpdOutgoingBuffers[i] = new double[iTemp];
		// copy from list
		j = 0;
		for(liData = vldData[i].begin() ; liData != vldData[i].end() ; liData++)
		{
			vpdOutgoingBuffers[i][j] = (*liData);
			j = j + 1;
		}
		vldData[i].clear();
		MPI_Isend(vpdOutgoingBuffers[i],iTemp,MPI_DOUBLE,i,GhostMessageDataTag,MPI_COMM_WORLD,&poOutgoingRequests[i]);
	}
	// wait for the communications to end and copy the buffers contents to the list
	for(i = 0 ; i < m_iDomainsCount ; i++)
	{
		if(i == m_iDomainID)
		{
			continue;
		}
		MPI_Wait(&poIncomingRequests[i],&poIncomingStatuses[i]);
		MPI_Wait(&poOutgoingRequests[i],&poOutgoingStatuses[i]);
		delete vpdOutgoingBuffers[i];
		// copy the buffer contents to the corresponding list
		for(j = 0 ; j < piBufferSizes[i] ; j++)
		{
			vldData[i].push_back(vpdIncomingBuffers[i][j]);
		}
		delete vpdIncomingBuffers[i];
	}
	vpdIncomingBuffers.clear();
	vpdOutgoingBuffers.clear();
	// free all communication allocated buffers
	delete [] poIncomingRequests;
	delete [] poOutgoingRequests;
	delete [] poIncomingStatuses;
	delete [] poOutgoingStatuses;
	delete [] piBufferSizes;
	MPI_Barrier(MPI_COMM_WORLD);
	
	UnpackGhostNodes(&vldData);
	RegisterNodes();
	MPI_Barrier(MPI_COMM_WORLD);
}
bool MainDataStructure::IsTagLower(const unsigned int& iDomainID1,const unsigned int& iNodeID1,const unsigned int& iDomainID2,const unsigned int& iNodeID2)
{
	if(iDomainID1 < iDomainID2)
	{
		return true;
	}
	else if(iDomainID1 > iDomainID2)
	{
		return false;
	}
	else
	{
		if(iNodeID1 < iNodeID2)
		{
			return true;
		}
		else if(iNodeID1 > iNodeID2)
		{
			return false;
		}
	}
	return true;
}
list<DislocationNode*>* MainDataStructure::GetLocalNodes()
{
	return &m_lpoLocalNodes;
}
list<DislocationNode*>* MainDataStructure::GetRemoteNodes()
{
	return &m_lpoRemoteNodes;
}
DislocationNode* MainDataStructure::GetNode(DislocationNodeTag* poTag) const
{
	return m_vvpoNodes[poTag->GetDomainID()][poTag->GetNodeID()];
}
unsigned int MainDataStructure::GetDomainID() const
{
	return m_iDomainID;
}
unsigned int MainDataStructure::GetLocalNodeID()
{
	unsigned int iID = 0;
	if(m_liFreeNodeTags.empty())
	{
		unsigned int iNewIDsPackSize = 50;
		iID = m_iFirstFreeNodeIndex;
		unsigned int i = 0;
		for(i = 0 ; i < iNewIDsPackSize ; i++)
		{
			m_liFreeNodeTags.push_back(iID + i);
		}
		m_iFirstFreeNodeIndex = m_liFreeNodeTags.back();
	}
	iID = m_liFreeNodeTags.front();
	m_liFreeNodeTags.pop_front();
	// if the first free node index was used, update it for the next call
	if(iID >= m_iFirstFreeNodeIndex)
	{
		m_iFirstFreeNodeIndex = iID + 1;
	}
	return iID;
}
unsigned int MainDataStructure::GetDomainsCount() const
{
	return m_iDomainsCount;
}
void MainDataStructure::RecycleNodeID(const unsigned int& iNodeID)
{
	m_liFreeNodeTags.push_back(iNodeID);
}
void MainDataStructure::AddNode(DislocationNode* poNode)
{
	if(poNode == NULL)
	{
		return;
	}
	m_lpoLocalNodes.push_back(poNode);
}
void MainDataStructure::AddAndRegisterNode(DislocationNode* poNode)
{
	if(poNode == NULL)
	{
		return;
	}
	// make sure that the node does not exist
	unsigned int iDomainID = poNode->GetDomainID();
	unsigned int iNodeID = poNode->GetID();
	unsigned int iSize = (unsigned int)m_vvpoNodes[iDomainID].size();
	if((iSize > iNodeID) && (m_vvpoNodes[iDomainID][iNodeID] != NULL))
	{
		printf("@ %d: error: adding a node (%d,%d) that already exists\n",m_iDomainID,iDomainID,iNodeID);
		fflush(NULL);
		exit(1);
	}
	m_lpoLocalNodes.push_back(poNode);
	if(iSize <= iNodeID)
	{
		unsigned int iNewSize = iNodeID + 1;
		m_vvpoNodes[iDomainID].resize(iNewSize);
		unsigned int i = 0;
		for(i = iSize ; i < iNewSize ; i++)
		{
			m_vvpoNodes[iDomainID][i] = NULL;
		}
		m_iFirstFreeNodeIndex = iNodeID + 1;
	}
	m_vvpoNodes[iDomainID][iNodeID] = poNode;
	// update the first free node index
	if(iNodeID >= m_iFirstFreeNodeIndex)
	{
		m_iFirstFreeNodeIndex = iNodeID + 1;
	}
}
void MainDataStructure::SetTime(const double& dValue)
{
	m_dTime = dValue;
}
void MainDataStructure::SetTimeStep(const double& dValue)
{
	m_dTimeStep = dValue;
}
double MainDataStructure::GetTime() const
{
	return m_dTime;
}
double MainDataStructure::GetTimeStep() const
{
	return m_dTimeStep;
}
Matrix* MainDataStructure::GetAppliedStress()
{
	return &m_oAppliedStress;
}
Matrix* MainDataStructure::GetTotalStrain()
{
	return &m_oTotalStrain;
}
Matrix* MainDataStructure::GetPlasticStrain()
{
	return &m_oPlasticStrain;
}
Matrix* MainDataStructure::GetPlasticStrainIncrement()
{
	return &m_oPlasticStrainIncrement;
}
bool MainDataStructure::Is111Vector(Vector oVector)
{
	oVector.Normalize();
	double dX = oVector.GetX();
	double dY = oVector.GetY();
	double dZ = oVector.GetZ();
	double dTolerance = 1.0E-6;
	if(fabs(fabs(dX) - fabs(dY)) < dTolerance)
	{
		if(fabs(fabs(dY) - fabs(dZ)) < dTolerance)
		{
			return true;
		}
	}
	return false;
}
bool MainDataStructure::Is110Vector(Vector oVector)
{
	oVector.Normalize();
	double dX = oVector.GetX();
	double dY = oVector.GetY();
	double dZ = oVector.GetZ();
	double dTolerance = 1.0E-6;
	if(fabs(fabs(dX) - fabs(dY)) < dTolerance)
	{
		if(fabs(dZ) < dTolerance)
		{
			return true;
		}
	}
	
	if(fabs(fabs(dY) - fabs(dZ)) < dTolerance)
	{
		if(fabs(dX) < dTolerance)
		{
			return true;
		}
	}
	
	if(fabs(fabs(dZ) - fabs(dX)) < dTolerance)
	{
		if(fabs(dY) < dTolerance)
		{
			return true;
		}
	}
	return false;
}
bool MainDataStructure::Is112Vector(Vector oVector)
{
	oVector.Normalize();
	double dX = oVector.GetX();
	double dY = oVector.GetY();
	double dZ = oVector.GetZ();
	double dTolerance = 1.0E-6;
	double dOne = 1.0/sqrt(6.0);
	double dTwo = 2.0*dOne;
	if(fabs(fabs(dX) - dTwo) < dTolerance)
	{
		if(fabs(fabs(dY) - dOne) < dTolerance)
		{
			if(fabs(fabs(dZ) - dOne) < dTolerance)
			{
				return true;
			}
		}
	}
	
	if(fabs(fabs(dY) - dTwo) < dTolerance)
	{
		if(fabs(fabs(dZ) - dOne) < dTolerance)
		{
			if(fabs(fabs(dX) - dOne) < dTolerance)
			{
				return true;
			}
		}
	}
	
	if(fabs(fabs(dZ) - dTwo) < dTolerance)
	{
		if(fabs(fabs(dX) - dOne) < dTolerance)
		{
			if(fabs(fabs(dY) - dOne) < dTolerance)
			{
				return true;
			}
		}
	}
	return false;
}
bool MainDataStructure::Is100Vector(Vector oVector)
{
	oVector.Normalize();
	double dX = oVector.GetX();
	double dY = oVector.GetY();
	double dZ = oVector.GetZ();
	double dTolerance = 1.0E-6;
	if(fabs(fabs(dX) - 1.0) < dTolerance)
	{
		if(fabs(dY) < dTolerance)
		{
			if(fabs(dZ) < dTolerance)
			{
				return true;
			}
		}
	}
	
	if(fabs(fabs(dY) - 1.0) < dTolerance)
	{
		if(fabs(dZ) < dTolerance)
		{
			if(fabs(dX) < dTolerance)
			{
				return true;
			}
		}
	}
	
	if(fabs(fabs(dZ) - 1.0) < dTolerance)
	{
		if(fabs(dX) < dTolerance)
		{
			if(fabs(dY) < dTolerance)
			{
				return true;
			}
		}
	}
	return false;
}
AxisAlignedBoundingBox* MainDataStructure::GetSimulationBox()
{
	return &m_oProblemBox;
}
AxisAlignedBoundingBox* MainDataStructure::GetDomainBox() const
{
	return m_poDomainBox;
}
vector<AxisAlignedBoundingBox*>* MainDataStructure::GetDomainBoxes()
{
	return &m_vpoDomains;
}
void MainDataStructure::Initialize()
{
	m_iDomainsCount = 0;
	m_iDomainID = 0;
	m_iXDomainsCount = 0;
	m_iYDomainsCount = 0;
	m_iZDomainsCount = 0;
	m_iOutputFrequency = 0;
	m_oProblemBox.Reset();
	m_poDomainBox = NULL;
	m_lpoLocalNodes.clear();
	m_lpoRemoteNodes.clear();
	m_vpoDomains.clear();
	m_vvpoNodes.clear();
	m_dTime = 0.0;
	m_dTimeStep = 0.0;
	m_oAppliedStress.SetSize(3,3);
	m_oTotalStrain.SetSize(3,3);
	m_oPlasticStrain.SetSize(3,3);
	m_oPlasticStrainIncrement.SetSize(3,3);
	m_liFreeNodeTags.clear();
	m_iFirstFreeNodeIndex = 0;
}
void MainDataStructure::SetDomainBoxes()
{
	list<AxisAlignedBoundingBox*> lpoDomains = m_oProblemBox.UniformPartition(m_iXDomainsCount,m_iYDomainsCount,m_iZDomainsCount);
	list<AxisAlignedBoundingBox*>::iterator liBoxes;
	m_vpoDomains.resize(lpoDomains.size());
	unsigned int i = 0;
	for(liBoxes = lpoDomains.begin() ; liBoxes != lpoDomains.end() ; liBoxes++)
	{
		m_vpoDomains[i] = (*liBoxes);
		i = i + 1;
	}
	lpoDomains.clear();
	m_poDomainBox = m_vpoDomains[m_iDomainID];
}
void MainDataStructure::SetLocalNodes(const unsigned int& iMaxDomainID)
{
	unsigned int iNodesCount = (unsigned int)m_lpoLocalNodes.size();
	unsigned int* piDomainIndices = new unsigned int[iNodesCount];
	unsigned int* piNodeIndices = new unsigned int[iNodesCount];
	unsigned int* piDomainCounts = new unsigned int[m_iDomainsCount];
	unsigned int i = 0;
	for(i = 0 ; i < iNodesCount ; i++)
	{
		piDomainIndices[i] = 0;
		piNodeIndices[i] = 0;
	}
	for(i = 0 ; i < m_iDomainsCount ; i++)
	{
		piDomainCounts[i] = 0;
	}
	// processor 0 loops over all the nodes and creates tags for them
	list<DislocationNode*>::iterator liNodes;
	list<AxisAlignedBoundingBox*>::iterator liDomains;
	unsigned int iNodeID = 0;
	unsigned int iSize = (unsigned int)m_vpoDomains.size();
	unsigned int j = 0;
	i = 0;
	if(m_iDomainID == 0)
	{
		for(liNodes = m_lpoLocalNodes.begin() ; liNodes != m_lpoLocalNodes.end() ; liNodes++)
		{
			for(j = 0 ; j < iSize ; j++)
			{
				if(m_vpoDomains[j]->IsPointInside(*(*liNodes),1.0E-6))
				{
					piDomainIndices[i] = j;
					piNodeIndices[i] = piDomainCounts[j];
					piDomainCounts[j] = piDomainCounts[j] + 1;
					break;
				}
			}
			i = i + 1;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// now the tags are generated, distribute them to all domains
	MPI_Bcast(piDomainIndices,iNodesCount,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(piNodeIndices,iNodesCount,MPI_INT,0,MPI_COMM_WORLD);
	
	// update node and arm end tags
	vector< map<unsigned int,DislocationNodeTag*>* > vpmoTagMappings;
	vpmoTagMappings.resize(iMaxDomainID + 1);
	for(i = 0 ; i <= iMaxDomainID ; i++)
	{
		vpmoTagMappings[i] = new map<unsigned int,DislocationNodeTag*>;
		vpmoTagMappings[i]->clear();
	}
		
	DislocationNodeTag* poTag = NULL;
	i = 0;
	for(liNodes = m_lpoLocalNodes.begin() ; liNodes != m_lpoLocalNodes.end() ; liNodes++)
	{
		poTag = new DislocationNodeTag;
		poTag->SetDomainID(piDomainIndices[i]);
		poTag->SetNodeID(piNodeIndices[i]);
 		vpmoTagMappings[(*liNodes)->GetDomainID()]->insert(pair<unsigned int,DislocationNodeTag*>((*liNodes)->GetID(),poTag));
		(*liNodes)->SetDomainID(piDomainIndices[i]);
		(*liNodes)->SetID(piNodeIndices[i]);
		i = i + 1;
	}
		
	list<DislocationSegment*>* plpoArms;
	list<DislocationSegment*>::iterator liArms;
	map<unsigned int,DislocationNodeTag*>::iterator liMappings;
	for(liNodes = m_lpoLocalNodes.begin() ; liNodes != m_lpoLocalNodes.end() ; liNodes++)
	{
		plpoArms = (*liNodes)->GetArms();
		for(liArms = plpoArms->begin() ; liArms != plpoArms->end() ; liArms++)
		{
			liMappings = vpmoTagMappings[(*liArms)->GetEndDomainID()]->find((*liArms)->GetEndNodeID());
			poTag = (*liMappings).second;
			(*liArms)->SetEndDomainID(poTag->GetDomainID());
			(*liArms)->SetEndNodeID(poTag->GetNodeID());
		}
	}
		
	// free the communication memory
	delete [] piDomainIndices;
	delete [] piNodeIndices;
	delete [] piDomainCounts;
	
	// free the tag mapping
	for(i = 0 ; i <= iMaxDomainID ; i++)
	{
		for(liMappings = vpmoTagMappings[i]->begin() ; liMappings != vpmoTagMappings[i]->end() ; liMappings++)
		{
			if((*liMappings).second != NULL)
			{
				delete (*liMappings).second;
			}
		}
		vpmoTagMappings[i]->clear();
		delete vpmoTagMappings[i];
	}
	vpmoTagMappings.clear();
		
	// clear all remote nodes
	liNodes = m_lpoLocalNodes.begin();
	while(liNodes != m_lpoLocalNodes.end())
	{
		if((*liNodes)->GetDomainID() != m_iDomainID)
		{
			delete (*liNodes);
			liNodes = m_lpoLocalNodes.erase(liNodes);
		}
		else
		{
			liNodes++;
		}
	}	
	MPI_Barrier(MPI_COMM_WORLD);
}
void MainDataStructure::ClearGhostNodes()
{
	list<DislocationNode*>::iterator liNodes;
	for(liNodes = m_lpoRemoteNodes.begin() ; liNodes != m_lpoRemoteNodes.end() ; liNodes++)
	{
		if((*liNodes) != NULL)
		{
			delete (*liNodes);
		}
	}
	m_lpoRemoteNodes.clear();
}
void MainDataStructure::PackGhostNodes(vector< list<double> >* pvldData)
{
	pvldData->clear();
	pvldData->resize(m_iDomainsCount);
	vector<bool> vbIsNodeAdded;
	vbIsNodeAdded.resize(m_iDomainsCount);
	// initialize the lists
	unsigned int i = 0;
	for(i = 0 ; i < m_iDomainsCount ; i++)
	{
		pvldData->at(i).clear();
	}
	
	// pack the data
	list<DislocationNode*>::iterator liNodes;
	list<DislocationSegment*>::iterator liArms;
	list<DislocationSegment*>* plpoArms = NULL;
	unsigned int iTargetDomain = 0;
	for(liNodes = m_lpoLocalNodes.begin() ; liNodes != m_lpoLocalNodes.end() ; liNodes++)
	{
		for(i = 0 ; i < m_iDomainsCount ; i++)
		{
			vbIsNodeAdded[i] = false;
		}
		plpoArms = (*liNodes)->GetArms();
		for(liArms = plpoArms->begin() ; liArms != plpoArms->end() ; liArms++)
		{
			iTargetDomain = (*liArms)->GetEndDomainID();
			if(iTargetDomain != m_iDomainID)
			{
				if(!vbIsNodeAdded[iTargetDomain])
				{
					(*liNodes)->Pack(&(pvldData->at(iTargetDomain)));
					vbIsNodeAdded[iTargetDomain] = true;
				}
			}
		}
	}
	vbIsNodeAdded.clear();
}
void MainDataStructure::UnpackGhostNodes(vector< list<double> >* pvldData)
{
	unsigned int i = 0;
	list<double>* pldDomainData = NULL;
	DislocationNode* poNode = NULL;
	for(i = 0 ; i < m_iDomainsCount ; i++)
	{
		pldDomainData = &pvldData->at(i);
		while(!pldDomainData->empty())
		{
			poNode = new DislocationNode;
			poNode->Unpack(pldDomainData);
			m_lpoRemoteNodes.push_back(poNode);
		}
	}
}
void MainDataStructure::RegisterNodes()
{
	vector<unsigned int> viMaxNodeIndex;
	viMaxNodeIndex.resize(m_iDomainsCount);
	unsigned int i = 0;
	for(i = 0 ; i < m_iDomainsCount ; i++)
	{
		viMaxNodeIndex[i] = 0;
		m_vvpoNodes[i].clear();
	}
	
	// estimate the sizes
	list<DislocationNode*>::iterator liNodes;
	unsigned int iNodeID = 0;
	for(liNodes = m_lpoLocalNodes.begin() ; liNodes != m_lpoLocalNodes.end() ; liNodes++)
	{
		iNodeID = (*liNodes)->GetID();
		if(viMaxNodeIndex[m_iDomainID] < iNodeID)
		{
			viMaxNodeIndex[m_iDomainID] = iNodeID;
		}
	}
	
	unsigned int iDomainID = 0;
	for(liNodes = m_lpoRemoteNodes.begin() ; liNodes != m_lpoRemoteNodes.end() ; liNodes++)
	{
		iNodeID = (*liNodes)->GetID();
		iDomainID = (*liNodes)->GetDomainID();
		if(viMaxNodeIndex[iDomainID] < iNodeID)
		{
			viMaxNodeIndex[iDomainID] = iNodeID;
		}
	}
	
	unsigned int j = 0;
	for(i = 0 ; i < m_iDomainsCount ; i++)
	{
		viMaxNodeIndex[i] = viMaxNodeIndex[i] + 1;
		m_vvpoNodes[i].resize(viMaxNodeIndex[i]);
		for(j = 0 ; j < viMaxNodeIndex[i] ; j++)
		{
			m_vvpoNodes[i][j] = NULL;
		}
	}
	viMaxNodeIndex.clear();
	
	// register the nodes
	for(liNodes = m_lpoLocalNodes.begin() ; liNodes != m_lpoLocalNodes.end() ; liNodes++)
	{
		m_vvpoNodes[m_iDomainID][(*liNodes)->GetID()] = (*liNodes);
	}
	
	for(liNodes = m_lpoRemoteNodes.begin() ; liNodes != m_lpoRemoteNodes.end() ; liNodes++)
	{
 		m_vvpoNodes[(*liNodes)->GetDomainID()][(*liNodes)->GetID()] = (*liNodes);
	}
	// update the first free node index
	m_iFirstFreeNodeIndex = (unsigned int)m_vvpoNodes[m_iDomainID].size();
}
list<DislocationNode*>::iterator MainDataStructure::RemoveLocalNode(list<DislocationNode*>::iterator liNode)
{
	if((*liNode) != NULL)
	{
		RecycleNodeID((*liNode)->GetID());
		delete (*liNode);
	}
	return m_lpoLocalNodes.erase(liNode);
}
list<DislocationNode*>::iterator MainDataStructure::RemoveRemoteNode(list<DislocationNode*>::iterator liNode)
{
	if((*liNode) != NULL)
	{
		delete (*liNode);
	}
	return m_lpoRemoteNodes.erase(liNode);
}
void MainDataStructure::UpdateConnectivityTags(vector< map<unsigned int,DislocationNodeTag*>* >* pvmpoTagMaps)
{
	list<DislocationNode*>::iterator liNodes;
	list<DislocationSegment*>::iterator liArms;
	list<DislocationSegment*>* plpoArms = NULL;
	map<unsigned int,DislocationNodeTag*>* poDomainMap = NULL;
	map<unsigned int,DislocationNodeTag*>::iterator miItem;
	for(liNodes = m_lpoLocalNodes.begin() ; liNodes != m_lpoLocalNodes.end() ; liNodes++)
	{
		plpoArms = (*liNodes)->GetArms();
		for(liArms = plpoArms->begin() ; liArms != plpoArms->end() ; liArms++)
		{
			poDomainMap = pvmpoTagMaps->at((*liArms)->GetEndDomainID());
			miItem = poDomainMap->find((*liArms)->GetEndNodeID());
			if(miItem != poDomainMap->end())
			{
				(*liArms)->SetEndDomainID(miItem->second->GetDomainID());
				(*liArms)->SetEndNodeID(miItem->second->GetNodeID());
			}
		}
	}
}




