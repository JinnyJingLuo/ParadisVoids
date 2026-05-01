#include "BalancingServer.h"
#include "mpi.h"

namespace BalancingSystem
{
	BalancingServer* BalancingServer::m_poInstance = NULL;
	BalancingServer* BalancingServer::GetInstance()
	{
		if(m_poInstance == NULL)
		{
			m_poInstance = new BalancingServer;
		}
		return m_poInstance;
	}
	BalancingServer::~BalancingServer()
	{
		Reset();
	}
	void BalancingServer::SetDataStructure(MainDataStructure* poDataStructure)
	{
		m_poDataStructure = poDataStructure;
	}
	void BalancingServer::Migrate()
	{
		CommunicateMigratingNodes();
		// at this point, we are done with all of the node communications and tag map
		// generation. it is required now to broadcast the tag maps from all of the processors to
		// all of the processors
		CommunicateTagMaps();
		MPI_Barrier(MPI_COMM_WORLD);
	}
	BalancingServer::BalancingServer()
	{
		Initialize();
	}
	void BalancingServer::Reset()
	{
		m_poDataStructure = NULL;
		FreeTagMaps();
	}
	void BalancingServer::Initialize()
	{
		m_poDataStructure = NULL;
		m_lpoTagMaps.clear();
	}
	void BalancingServer::PackMigratingNodes(vector< list<double> >* pvldData)
	{
		list<DislocationNode*>* plpoNodes =  m_poDataStructure->GetLocalNodes();
		list<DislocationNode*>::iterator liNodes;
		vector<AxisAlignedBoundingBox*>* pvpoDomainBoxes = m_poDataStructure->GetDomainBoxes();
		list<AxisAlignedBoundingBox*>::iterator liBoxes;
		AxisAlignedBoundingBox* poLocalDomainBox = m_poDataStructure->GetDomainBox();
		double dTolerance = 1.0E-3;
		unsigned int i = 0;
		unsigned int iSize = (unsigned int)pvpoDomainBoxes->size();
		AxisAlignedBoundingBox* poTargetDomainBox = NULL;
		unsigned int iTargetDomainID = 0;
		liNodes = plpoNodes->begin();
		pvldData->clear();
		pvldData->resize(m_poDataStructure->GetDomainsCount());
		while(liNodes != plpoNodes->end())
		{
			if(poLocalDomainBox->IsPointInside(*(*liNodes),dTolerance))
			{
				liNodes++;
				continue;
			}
			// if the local node is not in the local domain box, then this node needs to be 
			// moved to the suitable domain
			poTargetDomainBox = NULL;
			iTargetDomainID = 0;
			for(i = 0 ; i < iSize ; i++)
			{
				if(pvpoDomainBoxes->at(i)->IsPointInside(*(*liNodes),dTolerance))
				{
					poTargetDomainBox = pvpoDomainBoxes->at(i);
					iTargetDomainID = i;
					break;
				}
			}
			// make sure that the target domain was found
			if(poTargetDomainBox == NULL)
			{
				printf("@ %d: error: target domain not found while migrating node (%d,%d) at (%lf,%lf,%lf)\n",m_poDataStructure->GetDomainID(),(*liNodes)->GetDomainID(),(*liNodes)->GetID(),(*liNodes)->GetX(),(*liNodes)->GetY(),(*liNodes)->GetZ());
				fflush(NULL);
				exit(1);
			}
			// now we have a node and we know its target domain, pack this node into the 
			// suitable buffer and remove the node from the current domain
			(*liNodes)->Pack(&pvldData->at(iTargetDomainID));
			liNodes = m_poDataStructure->RemoveLocalNode(liNodes);
		}
	}
	void BalancingServer::CommunicateMigratingNodes()
	{
		vector< list<double> > vldData;
		PackMigratingNodes(&vldData);
		unsigned int iDomainsCount = m_poDataStructure->GetDomainsCount();
		unsigned int iDomainID = m_poDataStructure->GetDomainID();
		MPI_Barrier(MPI_COMM_WORLD);
		// send the data buffers lengths
		// issue length receives
		unsigned int* piBufferSizes = new unsigned int[iDomainsCount];
		MPI_Request* poIncomingRequests = new MPI_Request[iDomainsCount];
		unsigned int i = 0;
		for(i = 0 ; i < iDomainsCount ; i++)
		{
			if(i == iDomainID)
			{
				piBufferSizes[i] = 0;
				continue;
			}
			MPI_Irecv(&piBufferSizes[i],1,MPI_INT,i,MigratingNodesMessageLengthTag,MPI_COMM_WORLD,&poIncomingRequests[i]);
		}
		// issue length sends
		unsigned int iTemp = 0;
		MPI_Request* poOutgoingRequests = new MPI_Request[iDomainsCount];
		for(i = 0 ; i < iDomainsCount ; i++)
		{
			if(i == iDomainID)
			{
				continue;
			}
			iTemp = (unsigned int)vldData[i].size();
			MPI_Isend(&iTemp,1,MPI_INT,i,MigratingNodesMessageLengthTag,MPI_COMM_WORLD,&poOutgoingRequests[i]);
		}
		// wait for the communications to end
		MPI_Status* poIncomingStatuses = new MPI_Status[iDomainsCount];
		MPI_Status* poOutgoingStatuses = new MPI_Status[iDomainsCount];
		for(i = 0 ; i < iDomainsCount ; i++)
		{
			if(i == iDomainID)
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
		vpdIncomingBuffers.resize(iDomainsCount);
		for(i = 0 ; i < iDomainsCount ; i++)
		{
			if(i == iDomainID)
			{
				vpdIncomingBuffers[i] = NULL;
				continue;
			}
			vpdIncomingBuffers[i] = new double[piBufferSizes[i]];
			MPI_Irecv(vpdIncomingBuffers[i],piBufferSizes[i],MPI_DOUBLE,i,MigratingNodesMessageDataTag,MPI_COMM_WORLD,&poIncomingRequests[i]);
		}
		// issue data sends
		vector<double*> vpdOutgoingBuffers;
		vpdOutgoingBuffers.resize(iDomainsCount);
		list<double>::iterator liData;
		unsigned int j = 0;
		for(i = 0 ; i < iDomainsCount ; i++)
		{
			if(i == iDomainID)
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
			MPI_Isend(vpdOutgoingBuffers[i],iTemp,MPI_DOUBLE,i,MigratingNodesMessageDataTag,MPI_COMM_WORLD,&poOutgoingRequests[i]);
		}
		// wait for the communications to end and copy the buffers contents to the list
		for(i = 0 ; i < iDomainsCount ; i++)
		{
			if(i == iDomainID)
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
		UnpackMigratingNodes(&vldData);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	void BalancingServer::UnpackMigratingNodes(vector< list<double> >* pvldData)
	{
		unsigned int i = 0;
		unsigned int iDomainsCount = m_poDataStructure->GetDomainsCount();
		// initialize the tag maps
		FreeTagMaps();
		// unpack the nodes
		list<double>* pldDomainData = NULL;
		DislocationNode* poNode = NULL;
		unsigned int iDomainID = m_poDataStructure->GetDomainID();
		unsigned int iNewNodeID = 0;
		DislocationNodeTagMap* poTagMap = NULL;
		for(i = 0 ; i < iDomainsCount ; i++)
		{
			pldDomainData = &pvldData->at(i);
			while(!pldDomainData->empty())
			{
				poNode = new DislocationNode;
				poNode->Unpack(pldDomainData);
				// get new ID
				iNewNodeID = m_poDataStructure->GetLocalNodeID();
				// create and set the tag map
				poTagMap = new DislocationNodeTagMap;
				poTagMap->SetSourceDomainID(poNode->GetDomainID());
				poTagMap->SetSourceNodeID(poNode->GetID());
				poTagMap->SetTargetDomainID(iDomainID);
				poTagMap->SetTargetNodeID(iNewNodeID);
				m_lpoTagMaps.push_back(poTagMap);
				// add the node to the local nodes list
				poNode->SetDomainID(iDomainID);
				poNode->SetID(iNewNodeID);
				m_poDataStructure->AddNode(poNode);
			}
		}
	}
	void BalancingServer::FreeTagMaps()
	{
		list<DislocationNodeTagMap*>::iterator liMaps;
		for(liMaps = m_lpoTagMaps.begin() ; liMaps != m_lpoTagMaps.end() ; liMaps++)
		{
			if((*liMaps) != NULL)
			{
				delete (*liMaps);
			}
		}
		m_lpoTagMaps.clear();
	}
	void BalancingServer::CommunicateTagMaps()
	{
		list<unsigned int> liMapData;
		PackTagMaps(&liMapData);
		unsigned int iDomainsCount = m_poDataStructure->GetDomainsCount();
		unsigned int iDomainID = m_poDataStructure->GetDomainID();
		MPI_Barrier(MPI_COMM_WORLD);
		// send the data buffers lengths
		// issue length receives
		unsigned int* piBufferSizes = new unsigned int[iDomainsCount];
		MPI_Request* poIncomingRequests = new MPI_Request[iDomainsCount];
		unsigned int i = 0;
		for(i = 0 ; i < iDomainsCount ; i++)
		{
			if(i == iDomainID)
			{
				piBufferSizes[i] = 0;
				continue;
			}
			MPI_Irecv(&piBufferSizes[i],1,MPI_INT,i,TagMapsMessageLengthTag,MPI_COMM_WORLD,&poIncomingRequests[i]);
		}
		// issue length sends
		unsigned int iTemp = 0;
		MPI_Request* poOutgoingRequests = new MPI_Request[iDomainsCount];
		iTemp = (unsigned int)liMapData.size();
		for(i = 0 ; i < iDomainsCount ; i++)
		{
			if(i == iDomainID)
			{
				continue;
			}
			MPI_Isend(&iTemp,1,MPI_INT,i,TagMapsMessageLengthTag,MPI_COMM_WORLD,&poOutgoingRequests[i]);
		}
		// wait for the communications to end
		MPI_Status* poIncomingStatuses = new MPI_Status[iDomainsCount];
		MPI_Status* poOutgoingStatuses = new MPI_Status[iDomainsCount];
		for(i = 0 ; i < iDomainsCount ; i++)
		{
			if(i == iDomainID)
			{
				continue;
			}
			MPI_Wait(&poIncomingRequests[i],&poIncomingStatuses[i]);
			MPI_Wait(&poOutgoingRequests[i],&poOutgoingStatuses[i]);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		// send the data buffers
		// issue data receives
		vector<unsigned int*> vpiIncomingBuffers;
		vpiIncomingBuffers.resize(iDomainsCount);
		for(i = 0 ; i < iDomainsCount ; i++)
		{
			if(i == iDomainID)
			{
				vpiIncomingBuffers[i] = NULL;
				continue;
			}
			vpiIncomingBuffers[i] = new unsigned int[piBufferSizes[i]];
			MPI_Irecv(vpiIncomingBuffers[i],piBufferSizes[i],MPI_INT,i,TagMapsMessageDataTag,MPI_COMM_WORLD,&poIncomingRequests[i]);
		}
		// issue data sends
		iTemp = (unsigned int)liMapData.size();
		unsigned int* piOutgoingBuffer = new unsigned int[iTemp];
		list<unsigned int>::iterator liData;
		i = 0;
		for(liData = liMapData.begin() ; liData != liMapData.end() ; liData++)
		{
			piOutgoingBuffer[i] = (*liData);
			i = i + 1;
		}
		liMapData.clear();
			
		for(i = 0 ; i < iDomainsCount ; i++)
		{
			if(i == iDomainID)
			{
				continue;
			}
			MPI_Isend(piOutgoingBuffer,iTemp,MPI_INT,i,TagMapsMessageDataTag,MPI_COMM_WORLD,&poOutgoingRequests[i]);
		}
		// wait for the communications to end and copy the buffers contents to the list
		unsigned int j = 0;
		for(i = 0 ; i < iDomainsCount ; i++)
		{
			if(i == iDomainID)
			{
				continue;
			}
			MPI_Wait(&poIncomingRequests[i],&poIncomingStatuses[i]);
			MPI_Wait(&poOutgoingRequests[i],&poOutgoingStatuses[i]);
			
			// copy the buffer contents to the corresponding list
			for(j = 0 ; j < piBufferSizes[i] ; j++)
			{
				liMapData.push_back(vpiIncomingBuffers[i][j]);
			}
			delete [] vpiIncomingBuffers[i];
		}
		delete [] piOutgoingBuffer;
		vpiIncomingBuffers.clear();
		// free all communication allocated buffers
		delete [] poIncomingRequests;
		delete [] poOutgoingRequests;
		delete [] poIncomingStatuses;
		delete [] poOutgoingStatuses;
		delete [] piBufferSizes;
		MPI_Barrier(MPI_COMM_WORLD);
		UnpackTagMaps(&liMapData);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	void BalancingServer::PackTagMaps(list<unsigned int>* pliMapData)
	{
		pliMapData->clear();
		list<DislocationNodeTagMap*>::iterator liMaps;
		for(liMaps = m_lpoTagMaps.begin() ; liMaps != m_lpoTagMaps.end() ; liMaps++)
		{
			(*liMaps)->Pack(pliMapData);
		}
	}
	void BalancingServer::UnpackTagMaps(list<unsigned int>* pliMapData)
	{			
		list<unsigned int>::iterator liData = pliMapData->begin();
		DislocationNodeTagMap* poTagMap = NULL;
		while(!pliMapData->empty())
		{
			poTagMap = new DislocationNodeTagMap;
			poTagMap->Unpack(pliMapData);
			m_lpoTagMaps.push_back(poTagMap);
		}			
		// now the tag map list has all the tag maps, place them in a vector of maps
		vector< map<unsigned int,DislocationNodeTag*>* > vmpoTagMaps;
		unsigned int iDomainsCount = m_poDataStructure->GetDomainsCount();
		vmpoTagMaps.resize(iDomainsCount);
		unsigned int i = 0;
		for(i = 0 ; i < iDomainsCount ; i++)
		{
			vmpoTagMaps[i] = new map<unsigned int,DislocationNodeTag*>;
		}			
		list<DislocationNodeTagMap*>::iterator liMaps;
		DislocationNodeTag* poTag = NULL;
		for(liMaps = m_lpoTagMaps.begin() ; liMaps != m_lpoTagMaps.end() ; liMaps++)
		{
			poTag = new DislocationNodeTag;
			poTag->SetDomainID((*liMaps)->GetTargetDomainID());
			poTag->SetNodeID((*liMaps)->GetTargetNodeID());
			vmpoTagMaps[(*liMaps)->GetSourceDomainID()]->insert(pair<unsigned int,DislocationNodeTag*>((*liMaps)->GetSourceNodeID(),poTag));
		}
		// now update all the tags in the local nodes arms			
		m_poDataStructure->UpdateConnectivityTags(&vmpoTagMaps);			
		// clear the tag map table
		map<unsigned int,DislocationNodeTag*>::iterator miMaps;
		for(i = 0 ; i < iDomainsCount ; i++)
		{
			for(miMaps = vmpoTagMaps[i]->begin() ; miMaps != vmpoTagMaps[i]->end() ; miMaps++)
			{
				if(miMaps->second != NULL)
				{
					delete miMaps->second;
				}
			}
			vmpoTagMaps[i]->clear();
			delete vmpoTagMaps[i];
		}			
		vmpoTagMaps.clear();
		// finally, clear the tag map list
		FreeTagMaps();
	}
}



