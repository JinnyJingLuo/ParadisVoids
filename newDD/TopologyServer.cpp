#include "TopologyServer.h"
#include "ChangeConnectionOperation.h"
#include "Line.h"
#include "mpi.h"
#include "cmath"

using namespace GeometrySystem;

namespace TopologySystem
{
	TopologyServer* TopologyServer::m_poInstance = NULL;
	TopologyServer* TopologyServer::GetInstance()
	{
		if(m_poInstance == NULL)
		{
			m_poInstance = new TopologyServer;
		}
		return m_poInstance;
	}
	TopologyServer::~TopologyServer()
	{
		Reset();
	}
	void TopologyServer::SetDataStructure(MainDataStructure* poDataStructure)
	{
		m_poDataStructure = poDataStructure;
	}
	void TopologyServer::Remesh()
	{
		if(m_poDataStructure == NULL)
		{
			return;
		}
		ClearOperations();
		list<DislocationNode*>* plpoNodes = m_poDataStructure->GetLocalNodes();
		list<DislocationNode*>::iterator liNodes;
		list<DislocationSegment*>* plpoArms = NULL;
		list<DislocationSegment*>::iterator liArms;
		DislocationNode* poNode = NULL;
		DislocationNode* poNeighbour = NULL;
		DislocationNode* poNewNode = NULL;
		Vector oLine;
		double dLength = 0.0;
		unsigned int iDomainID = m_poDataStructure->GetDomainID();
		unsigned int iNeighbourDomainID = 0;
		DislocationSegment* poArm1 = NULL;
		DislocationSegment* poArm2 = NULL;
		DislocationSegment* poArm1Reverse = NULL;
		DislocationSegment* poArm2Reverse = NULL;
		Vector oNormal;
		Vector oBurgers;
		double dMaxSegment = 200.0;
		ChangeConnectionOperation* poOperation = NULL;
		for(liNodes = plpoNodes->begin() ; liNodes != plpoNodes->end() ; liNodes++)
		{
			poNode = (*liNodes);
			plpoArms = (*liNodes)->GetArms();
			for(liArms = plpoArms->begin() ; liArms != plpoArms->end() ; liArms++)
			{
				poNeighbour = m_poDataStructure->GetNode((*liArms)->GetEndTag());
				if(poNeighbour == NULL)
				{
					printf("@ %d : warning: null neighbour (%d,%d) found in remesh",iDomainID,(*liArms)->GetEndTag()->GetDomainID(),(*liArms)->GetEndTag()->GetNodeID());
					continue;
				}
				if(!poNode->IsNodeSuperior(poNeighbour))
				{
					continue;
				}
				oLine.SetByPoints(*poNode,*poNeighbour);
				dLength = oLine.Length();
				if(dLength <= dMaxSegment)
				{
					continue;
				}
				SplitArm(poNode,poNeighbour,(*liArms));
			}
		}
		CommunicateOperations();
		ProcessOperations();
		ClearOperations();
	}
	void TopologyServer::HandleCollisions()
	{
		MPI_Barrier(MPI_COMM_WORLD);
		ClearOperations();
		GenerateNodeCollisionLists();
		HandleSegmentSegmentCollisions();
		HandleHingeJointCollisions();
		//HandleTriangularLoops();
		m_iCollisionCallsCount = m_iCollisionCallsCount + 1;
		CommunicateOperations();
		MPI_Barrier(MPI_COMM_WORLD);
		ProcessOperations();
		ClearOperations();
		MPI_Barrier(MPI_COMM_WORLD);
	}
	TopologyServer::TopologyServer()
	{
		Initialize();
	}
	void TopologyServer::Reset()
	{
		ClearOperations();
		Initialize();
	}
	void TopologyServer::Initialize()
	{
		m_poDataStructure = NULL;
		m_lpoOperations.clear();
		m_vlpoNodes.clear();
		m_iCollisionCallsCount = 0;
	}
	void TopologyServer::ClearOperations()
	{
		list<TopologicalOperation*>::iterator liOperations;
		for(liOperations = m_lpoOperations.begin() ; liOperations != m_lpoOperations.end() ; liOperations++)
		{
			if((*liOperations) != NULL)
			{
				delete (*liOperations);
			}
		}
		m_lpoOperations.clear();
	}
	void TopologyServer::CommunicateOperations()
	{
		vector< list<double> > vldData;
		PackOperations(&vldData);
		MPI_Barrier(MPI_COMM_WORLD);
		
		// send the data buffers lengths
		// issue length receives
		unsigned int iDomainsCount = m_poDataStructure->GetDomainsCount();
		unsigned int iDomainID = m_poDataStructure->GetDomainID();
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
			MPI_Irecv(&piBufferSizes[i],1,MPI_INT,i,OperationMessageLengthTag,MPI_COMM_WORLD,&poIncomingRequests[i]);
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
			MPI_Isend(&iTemp,1,MPI_INT,i,OperationMessageLengthTag,MPI_COMM_WORLD,&poOutgoingRequests[i]);
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
			MPI_Irecv(vpdIncomingBuffers[i],piBufferSizes[i],MPI_DOUBLE,i,OperationMessageDataTag,MPI_COMM_WORLD,&poIncomingRequests[i]);
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
			MPI_Isend(vpdOutgoingBuffers[i],iTemp,MPI_DOUBLE,i,OperationMessageDataTag,MPI_COMM_WORLD,&poOutgoingRequests[i]);
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
		
		UnpackOperations(&vldData);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	void TopologyServer::PackOperations(vector< list<double> >* pvldBuffers)
	{
		pvldBuffers->clear();
		unsigned int iDomainsCount = m_poDataStructure->GetDomainsCount();
		pvldBuffers->resize(iDomainsCount);
		list<TopologicalOperation*>::iterator liOperations;
		unsigned int iDomainID = 0;
		for(liOperations = m_lpoOperations.begin() ; liOperations != m_lpoOperations.end() ; liOperations++)
		{
			iDomainID = (*liOperations)->GetResponsibleDomainID();
			(*liOperations)->Pack(&pvldBuffers->at(iDomainID));
			delete (*liOperations);
		}
		m_lpoOperations.clear();
	}
	void TopologyServer::UnpackOperations(vector< list<double> >* pvldBuffers)
	{
		unsigned int i = 0;
		unsigned int iSize = (unsigned int)pvldBuffers->size();
		list<double>* pldDomainData = NULL;
		TopologicalOperation* poOperation = NULL;
		for(i = 0 ; i < iSize ; i++)
		{
			pldDomainData = &pvldBuffers->at(i);
			while(!pldDomainData->empty())
			{
				poOperation = TopologicalOperation::CreateAndUnpackOperation(pldDomainData);
				m_lpoOperations.push_back(poOperation);
			}
		}
	}
	void TopologyServer::ProcessOperations()
	{
		list<TopologicalOperation*>::iterator liOperations;
		for(liOperations = m_lpoOperations.begin() ; liOperations != m_lpoOperations.end() ; liOperations++)
		{
			if((*liOperations)->GetType() == ChangeConnection)
			{
				ProcessChangeConnectionOperation((ChangeConnectionOperation*)(*liOperations));
			}
		}
	}
	void TopologyServer::GenerateNodeCollisionLists()
	{
		ClearNodeCollisionLists();
		// all the segments that span domains will not be included in the colliding
		// segments list
		unsigned int iDomainID = m_poDataStructure->GetDomainID();
		unsigned int iCollisionResolution = 2;
		list<AxisAlignedBoundingBox*> lpoPartitions = m_poDataStructure->GetDomainBox()->UniformPartition(iCollisionResolution,iCollisionResolution,iCollisionResolution);
		m_vlpoNodes.resize(lpoPartitions.size());
		list<DislocationNode*>* plpoNodes = m_poDataStructure->GetLocalNodes();
		list<DislocationNode*>::iterator liNodes;
		double dTolerance = 1.0E-2;
		list<AxisAlignedBoundingBox*>::iterator liPartitions;
		unsigned int i = 0;
		for(liNodes = plpoNodes->begin() ; liNodes != plpoNodes->end() ; liNodes++)
		{
			// enable all node collisions
			(*liNodes)->EnableCollision();
			i = 0;
			for(liPartitions = lpoPartitions.begin() ; liPartitions != lpoPartitions.end() ; liPartitions++)
			{
				if((*liPartitions)->IsPointInside(*(*liNodes),dTolerance))
				{
					m_vlpoNodes[i].push_back((*liNodes));
					break;
				}
				i = i + 1;
			}
		}
		// there is no need for the partitions any more
		for(liPartitions = lpoPartitions.begin() ; liPartitions != lpoPartitions.end() ; liPartitions++)
		{
			if((*liPartitions) != NULL)
			{
				delete (*liPartitions);
			}
		}
		lpoPartitions.clear();
	}
	void TopologyServer::ClearNodeCollisionLists()
	{
		unsigned int iSize = m_vlpoNodes.size();
		unsigned int i = 0;
		for(i = 0 ; i < iSize ; i++)
		{
			m_vlpoNodes[i].clear();
		}
		m_vlpoNodes.clear();
	}
	bool TopologyServer::IsNodeCollisionValid(DislocationNode* poNode)
	{
		if(poNode == NULL)
		{
			return false;
		}
		if(!poNode->IsAllowedToCollide())
		{
			return false;
		}
		return true;
	}
	bool TopologyServer::AreSegmentsCloseAndApproaching(DislocationNode* poNode1,DislocationNode* poNode2,DislocationNode* poNode3,DislocationNode* poNode4,Point& oNearPoint1,Point& oNearPoint2)
	{
		oNearPoint1.Set(0.0,0.0,0.0);
		oNearPoint2.Set(0.0,0.0,0.0);
		double dLocalCoordinate1 = 0.0;
		double dLocalCoordinate2 = 0.0;
		double dToleranceSquared = 1.0E-6;
		Vector oM(*poNode1,*poNode2);
		Vector oN(*poNode3,*poNode4);
		// extremely short segments cannot collide
		double dA = oM.LengthSquared();
		if(dA < dToleranceSquared)
		{
			return false;
		}
		double dB = oN.LengthSquared();
		if(dB < dToleranceSquared)
		{
			return false;
		}
		double dC = oM*oN;
		Vector oV13(*poNode1,*poNode3);
		Vector oV14(*poNode1,*poNode4);
		Vector oV23(*poNode2,*poNode3);
		Vector oV24(*poNode2,*poNode4);
		double dD = oM*oV13;
		double dE = oN*oV13*(-1.0);
		double dF = dA*dB - dC*dC;
		
		// get the 2 nearest points
		if(fabs(dF) < dToleranceSquared)
		{	
			// the segments are parallel
			double dTempCoordinate1[4] = {0.0,1.0,dD/dA,oV14*oN/dA};
			double dTempCoordinate2[4] = {dE/dB,oV23*oN*(-1.0/dB),0.0,1.0};
			double dDistanceSquared[4] = {0.0,0.0,0.0,0.0};
			unsigned int i = 0;
			for(i = 0 ; i < 4 ; i++)
			{
				if(dTempCoordinate1[i] < 0.0)
				{
					dTempCoordinate1[i] = 0.0;
				}
				else if(dTempCoordinate1[i] > 1.0)
				{
					dTempCoordinate1[i] = 1.0;
				}
				if(dTempCoordinate2[i] < 0.0)
				{
					dTempCoordinate2[i] = 0.0;
				}
				else if(dTempCoordinate2[i] > 1.0)
				{
					dTempCoordinate2[i] = 1.0;
				}
				dDistanceSquared[i] = Vector(*poNode1 + oM*dTempCoordinate1[i],*poNode3 + oN*dTempCoordinate2[i]).LengthSquared();
			}
			
			double dMinimumDistanceSquared = dDistanceSquared[0];
			dLocalCoordinate1 = dTempCoordinate1[0];
			dLocalCoordinate2 = dTempCoordinate2[0];
			for(i = 1 ; i < 4 ; i++)
			{
				if(dMinimumDistanceSquared < dDistanceSquared[i])
				{
					dMinimumDistanceSquared = dDistanceSquared[i];
					dLocalCoordinate1 = dTempCoordinate1[i];
					dLocalCoordinate2 = dTempCoordinate2[i];
				}
			}
		}
		else
		{
			dLocalCoordinate1 = (dC*dE + dD*dB)/dF;
			dLocalCoordinate2 = (dE + dC*dLocalCoordinate1)/dB;
		}
		
		if(dLocalCoordinate1 < 0.0)
		{
			dLocalCoordinate1 = 0.0;
		}
		if(dLocalCoordinate1 > 1.0)
		{
			dLocalCoordinate1 = 1.0;
		}
		
		if(dLocalCoordinate2 < 0.0)
		{
			dLocalCoordinate2 = 0.0;
		}
		if(dLocalCoordinate2 > 1.0)
		{
			dLocalCoordinate2 = 1.0;
		}
		
		// now check if the segments are close
		oNearPoint1 = *poNode1 + oM*dLocalCoordinate1;
		oNearPoint2 = *poNode3 + oN*dLocalCoordinate2;
		Vector oNewNodesSpacing(oNearPoint2,oNearPoint1);
		Vector oVel12 = poNode2->GetVelocity() - poNode1->GetVelocity();
		Vector oVel34 = poNode4->GetVelocity() - poNode3->GetVelocity();
		Vector oVel13 = poNode1->GetVelocity() - poNode3->GetVelocity();
		Vector oRelativeVelocity = oVel13 + oVel12*dLocalCoordinate1 - oVel34*dLocalCoordinate2;
		dD = oNewNodesSpacing.LengthSquared();
		double dDDot = 2.0*(oNewNodesSpacing*oRelativeVelocity);
		double dAbsoluteMinimum = 1.0E-2;
		double dCollisionSpacingSquared = 25.0;
		if(dD < dCollisionSpacingSquared)
		{
			if(dDDot < 0.0)
			{
				return true;
			}
		}
		if(dD < dAbsoluteMinimum)
		{
			return true;
		}
		return false;
	}
	bool TopologyServer::AreNodesCloseAndApproaching(DislocationNode* poNode1,DislocationNode* poNode2)
	{
		double dDistance = poNode1->Distance(*poNode2);
		double dAbsoluteMinimum = 1.0E-10;
		double dCollisionSpacing = 5.0;
		if(dDistance < dAbsoluteMinimum)
		{
			return true;
		}
		if(dDistance > dCollisionSpacing)
		{
			return false;
		}
		Vector oVelocity1 = poNode1->GetVelocity();
		Vector oVelocity2 = poNode2->GetVelocity();
		Vector oRelativeVelocity = oVelocity2 - oVelocity1;
		Vector oLine(*poNode1,*poNode2);
		oLine.Normalize();
		if(oLine*oRelativeVelocity < 0.0)
		{
			return true;
		}
		return false;
	}
	bool TopologyServer::IsTriangleCollapsing(DislocationNode* poNode1,DislocationNode* poNode2,DislocationNode* poNode3)
	{
		Vector oV1(*poNode1,*poNode2);
		Vector oV2(*poNode1,*poNode3);
		double dL1 = oV1.Length();
		double dL2 = oV2.Length();
		oV1.Normalize();
		oV2.Normalize();			
		double dCosTheta = oV1*oV2;
		if(dCosTheta < 0.85)						// the angle is too large for merging
		{
			return false;
		}
		if(dCosTheta > 0.95)						// the angle is too small, merging is allowed
		{
			return true;
		}
		Vector oV1Dot = poNode2->GetVelocity() - poNode1->GetVelocity();
		Vector oV2Dot = poNode3->GetVelocity() - poNode1->GetVelocity();
		oV1Dot = oV1Dot*(1.0/dL1);
		oV2Dot = oV2Dot*(1.0/dL2);
		double dR = (oV1Dot*oV2 + oV1*oV2Dot) - dCosTheta*(oV1Dot*oV1 + oV2Dot*oV2);
		if(dR > 0)									// the angle is moderate but decreasing
		{
			return true;
		}
		return false;
	}
	bool TopologyServer::GetExactCollisionPoint(DislocationNode* poNode1,DislocationNode* poNode2,Point& oCollisionPoint)
	{
		oCollisionPoint.Set(0.0,0.0,0.0);
		if(poNode1 == NULL)
		{
			return false;
		}
		if(poNode2 == NULL)
		{
			return false;
		}
		if(poNode1->GetNeighboursCount() != 2)
		{
			return false;
		}
		if(poNode2->GetNeighboursCount() != 2)
		{
			return false;
		}
		unsigned int iConstraintType = 0;
		Vector oNormal1;
		poNode1->GetDynamicConstraint(iConstraintType,oNormal1);
		if(iConstraintType != 1)
		{
			return false;
		}
		Vector oNormal2;
		poNode2->GetDynamicConstraint(iConstraintType,oNormal2);
		if(iConstraintType != 1)
		{
			return false;
		}
		// now we have 2 nodes that move in a plane, if they lie on the same plane, the collision
		// point is simply the mid point of the line connecting them. the possibility of the nodes
		// lying on two parallel planes is ruled out because this case is filtered before calling
		// this function
		oNormal1.Normalize();
		oNormal2.Normalize();
		double dTolerance = 1.0E-6;
		oCollisionPoint = (*poNode1 + *poNode2)*(0.5);
		if((oNormal1^oNormal2).Length() > dTolerance)
		{
			Matrix oSystem(5,5);
			oSystem.Set(1,1,1.0);
			oSystem.Set(1,2,0.0);
			oSystem.Set(1,3,0.0);
			oSystem.Set(1,4,oNormal1.GetX());
			oSystem.Set(1,5,oNormal2.GetX());
			
			oSystem.Set(2,1,0.0);
			oSystem.Set(2,2,1.0);
			oSystem.Set(2,3,0.0);
			oSystem.Set(2,4,oNormal1.GetY());
			oSystem.Set(2,5,oNormal2.GetY());
			
			oSystem.Set(3,1,0.0);
			oSystem.Set(3,2,0.0);
			oSystem.Set(3,3,1.0);
			oSystem.Set(3,4,oNormal1.GetZ());
			oSystem.Set(3,5,oNormal2.GetZ());
			
			oSystem.Set(4,1,oNormal1.GetX());
			oSystem.Set(4,2,oNormal1.GetY());
			oSystem.Set(4,3,oNormal1.GetZ());
			oSystem.Set(4,4,0.0);
			oSystem.Set(4,5,0.0);
			
			oSystem.Set(5,1,oNormal2.GetX());
			oSystem.Set(5,2,oNormal2.GetY());
			oSystem.Set(5,3,oNormal2.GetZ());
			oSystem.Set(5,4,0.0);
			oSystem.Set(5,5,0.0);
			
			Matrix oRHS(5,1);
			oRHS.Set(1,1,oCollisionPoint.GetX());
			oRHS.Set(2,1,oCollisionPoint.GetY());
			oRHS.Set(3,1,oCollisionPoint.GetZ());
			oRHS.Set(4,1,oNormal1*Vector(*poNode1));
			oRHS.Set(5,1,oNormal2*Vector(*poNode2));
			
			Matrix oX = oSystem.Solve(oRHS);
			oCollisionPoint.SetX(oX.Get(1,1));
			oCollisionPoint.SetY(oX.Get(2,1));
			oCollisionPoint.SetZ(oX.Get(3,1));
		}
		return true;
	}
	Point TopologyServer::GetExactCollisionPoint(const Point& oPoint1,const Point& oPoint2,const Vector& oNormal1,const Vector& oNormal2)
	{
		// we have 2 points that move in a plane, if they lie on the same plane, the collision
		// point is simply the mid point of the line connecting them. the possibility of the nodes
		// lying on two parallel planes is ruled out because this case is filtered before calling
		// this function
		// the function assumes that the normals passed to it are normalized
		double dTolerance = 1.0E-6;
		Point oCollisionPoint = (oPoint1 + oPoint2)*(0.5);
		if((oNormal1^oNormal2).Length() > dTolerance)
		{
			Matrix oSystem(5,5);
			oSystem.Set(1,1,1.0);
			oSystem.Set(1,2,0.0);
			oSystem.Set(1,3,0.0);
			oSystem.Set(1,4,oNormal1.GetX());
			oSystem.Set(1,5,oNormal2.GetX());
			
			oSystem.Set(2,1,0.0);
			oSystem.Set(2,2,1.0);
			oSystem.Set(2,3,0.0);
			oSystem.Set(2,4,oNormal1.GetY());
			oSystem.Set(2,5,oNormal2.GetY());
			
			oSystem.Set(3,1,0.0);
			oSystem.Set(3,2,0.0);
			oSystem.Set(3,3,1.0);
			oSystem.Set(3,4,oNormal1.GetZ());
			oSystem.Set(3,5,oNormal2.GetZ());
			
			oSystem.Set(4,1,oNormal1.GetX());
			oSystem.Set(4,2,oNormal1.GetY());
			oSystem.Set(4,3,oNormal1.GetZ());
			oSystem.Set(4,4,0.0);
			oSystem.Set(4,5,0.0);
			
			oSystem.Set(5,1,oNormal2.GetX());
			oSystem.Set(5,2,oNormal2.GetY());
			oSystem.Set(5,3,oNormal2.GetZ());
			oSystem.Set(5,4,0.0);
			oSystem.Set(5,5,0.0);
			
			Matrix oRHS(5,1);
			oRHS.Set(1,1,oCollisionPoint.GetX());
			oRHS.Set(2,1,oCollisionPoint.GetY());
			oRHS.Set(3,1,oCollisionPoint.GetZ());
			oRHS.Set(4,1,oNormal1*Vector(oPoint1));
			oRHS.Set(5,1,oNormal2*Vector(oPoint2));
			
			Matrix oX = oSystem.Solve(oRHS);
			oCollisionPoint.SetX(oX.Get(1,1));
			oCollisionPoint.SetY(oX.Get(2,1));
			oCollisionPoint.SetZ(oX.Get(3,1));
		}
		return oCollisionPoint;
	}
	void TopologyServer::HandleSegmentSegmentCollisions()
	{
		list<DislocationNode*>* plpoNodes = m_poDataStructure->GetLocalNodes();
		list<DislocationSegment*>* plpoArms1 = NULL;
		list<DislocationSegment*>* plpoArms2 = NULL;
		list<DislocationNode*>::iterator liNode1;
		list<DislocationNode*>::iterator liNode2;
		list<DislocationSegment*>::iterator liArm1;
		list<DislocationSegment*>::iterator liArm2;
		DislocationNode* poNode1 = NULL;
		DislocationNode* poNode2 = NULL;
		DislocationNode* poNode3 = NULL;
		DislocationNode* poNode4 = NULL;
		bool bHasCollided = false;		// because every node is allowed one collision a time
		Point oNearPoint1;
		Point oNearPoint2;
		Vector oNormal1;
		Vector oNormal2;
		Vector oBurgers1;
		Vector oBurgers2;
		double dTolerance = 1.0E-6;
		double dTemp = 0.0;
		DislocationSegment* poArm1 = NULL;
		DislocationSegment* poArm2 = NULL;
		DislocationSegment* poArm3 = NULL;
		DislocationSegment* poArm4 = NULL;
		Point oCollisionPoint;
		DislocationNode* poNewNode = NULL;
		unsigned int iNeighbourDomainID = 0;
		ChangeConnectionOperation* poOperation = NULL;
		unsigned int iDomainID = m_poDataStructure->GetDomainID();
		double dMaximumPlaneSpacing = 1.0;
		for(liNode1 = plpoNodes->begin() ; liNode1 != plpoNodes->end() ; liNode1++)
		{
			poNode1 = (*liNode1);
			if(!IsNodeCollisionValid(poNode1))
			{
				continue;
			}
			if(poNode1->GetNeighboursCount() > 2)
			{
				continue;
			}
			bool bHasCollided = false;
			liNode2 = liNode1;
			liNode2++;
			while(liNode2 != plpoNodes->end())
			{
				if(bHasCollided)
				{
					break;
				}
				poNode3 = (*liNode2);
				if(!IsNodeCollisionValid(poNode3))
				{
					liNode2++;
					continue;
				}
				if(poNode3->GetNeighboursCount() > 2)
				{
					liNode2++;
					continue;
				}
				if(!poNode1->IsNodeSuperior(poNode3))
				{
					liNode2++;
					continue;
				}
				// now the distinct arms of nodes 1 and 3 may collide, check for that
				plpoArms1 = poNode1->GetArms();
				for(liArm1 = plpoArms1->begin() ; liArm1 != plpoArms1->end() ; liArm1++)
				{
					if(bHasCollided)
					{
						break;
					}
					poNode2 = m_poDataStructure->GetNode((*liArm1)->GetEndTag());
					if(!IsNodeCollisionValid(poNode2))
					{
						continue;
					}
					if(poNode2->GetNeighboursCount() > 2)
					{
						continue;
					}
					if(!IsNodeSuperior(poNode1,poNode2))
					{
						continue;
					}
					if(poNode2->IsSameNode(poNode3))
					{
						continue;
					}
					plpoArms2 = poNode3->GetArms();
					for(liArm2 = plpoArms2->begin() ; liArm2 != plpoArms2->end() ; liArm2++)
					{
						if(bHasCollided)
						{
							break;
						}
						poNode4 = m_poDataStructure->GetNode((*liArm2)->GetEndTag());
						if(!IsNodeCollisionValid(poNode4))
						{
							continue;
						}
						if(poNode4->GetNeighboursCount() > 2)
						{
							continue;
						}
						if(!IsNodeSuperior(poNode3,poNode4))
						{
							continue;
						}
						if(poNode4->IsSameNode(poNode1))
						{
							continue;
						}
						if(poNode4->IsSameNode(poNode2))
						{
							continue;
						}
						// are segments close enough for collision ?? are they approaching each 
						// other or moving away from each other ??
						if(!AreSegmentsCloseAndApproaching(poNode1,poNode2,poNode3,poNode4,oNearPoint1,oNearPoint2))
						{
							continue;
						}
						// now these 2 segments are supposed to collide, make sure that the extra
						// kinematic constraints will not cause problems
						oNormal1 = (*liArm1)->GetSlipPlaneNormal();
						oNormal2 = (*liArm2)->GetSlipPlaneNormal();
						oNormal1.Normalize();
						oNormal2.Normalize();
						if((oNormal1^oNormal2).Length() < dTolerance)
						{
							// if the two segments are on the same plane, collide them
							dTemp = oNormal1*Vector(*poNode1,*poNode3);
							if(fabs(dTemp) > dMaximumPlaneSpacing)
							{
								continue;
							}
						}
						// now the segments can definitely collide, but if the nearest point
						// on either of the two segments is too close to one of the segment
						// nodes, skip the collision this time, because we don't want to split
						// the nodes and place the new nodes at the same position
						if(oNearPoint1.Distance(*poNode1) < dTolerance)
						{
							continue;
						}
						if(oNearPoint1.Distance(*poNode2) < dTolerance)
						{
							continue;
						}
						if(oNearPoint2.Distance(*poNode3) < dTolerance)
						{
							continue;
						}
						if(oNearPoint2.Distance(*poNode4) < dTolerance)
						{
							continue;
						}
						// split the arms
						oCollisionPoint = GetExactCollisionPoint(oNearPoint1,oNearPoint2,oNormal1,oNormal2);
						// now we have the kinemtically correct collision location, create a node
						// and connect it to the 4 nodes of the 2 segments
						poNewNode = new DislocationNode;
						*poNewNode = oCollisionPoint;
						poNewNode->SetID(m_poDataStructure->GetLocalNodeID());
						poNewNode->SetCategory(UnconstrainedNode);
						poNewNode->SetDomainID(iDomainID);
						poNewNode->SetSurfaceNormal(Vector(0.0,0.0,0.0));
						// set its arms
						poArm1 = new DislocationSegment;
						poArm2 = new DislocationSegment;
						poArm3 = new DislocationSegment;
						poArm4 = new DislocationSegment;
						poArm1->SetEndDomainID(poNode1->GetDomainID());
						poArm1->SetEndNodeID(poNode1->GetID());
						poArm2->SetEndDomainID(poNode2->GetDomainID());
						poArm2->SetEndNodeID(poNode2->GetID());
						poArm3->SetEndDomainID(poNode3->GetDomainID());
						poArm3->SetEndNodeID(poNode3->GetID());
						poArm4->SetEndDomainID(poNode4->GetDomainID());
						poArm4->SetEndNodeID(poNode4->GetID());
						poArm1->SetSlipPlaneNormal(oNormal1);
						poArm2->SetSlipPlaneNormal(oNormal1);
						poArm3->SetSlipPlaneNormal(oNormal2);
						poArm4->SetSlipPlaneNormal(oNormal2);
						oBurgers1 = (*liArm1)->GetBurgersVector();
						oBurgers2 = (*liArm2)->GetBurgersVector();
						poArm1->SetBurgersVector(oBurgers1*(-1.0));
						poArm2->SetBurgersVector(oBurgers1);
						poArm3->SetBurgersVector(oBurgers2*(-1.0));
						poArm4->SetBurgersVector(oBurgers2);
						poNewNode->AddArm(poArm1);
						poNewNode->AddArm(poArm2);
						poNewNode->AddArm(poArm3);
						poNewNode->AddArm(poArm4);
						// change the connections of the arms
						poArm1 = poNode1->GetArmByTag(poNode2->GetDomainID(),poNode2->GetID());
						poArm2 = poNode2->GetArmByTag(poNode1->GetDomainID(),poNode1->GetID());
						poArm3 = poNode3->GetArmByTag(poNode4->GetDomainID(),poNode4->GetID());
						poArm4 = poNode4->GetArmByTag(poNode3->GetDomainID(),poNode3->GetID());
						poArm1->SetEndDomainID(poNewNode->GetDomainID());
						poArm1->SetEndNodeID(poNewNode->GetID());
						poArm2->SetEndDomainID(poNewNode->GetDomainID());
						poArm2->SetEndNodeID(poNewNode->GetID());
						poArm3->SetEndDomainID(poNewNode->GetDomainID());
						poArm3->SetEndNodeID(poNewNode->GetID());
						poArm4->SetEndDomainID(poNewNode->GetDomainID());
						poArm4->SetEndNodeID(poNewNode->GetID());
						// add the new node to the data structure
						m_poDataStructure->AddAndRegisterNode(poNewNode);
						// both nodes 1 and 3 are guaranteed to be in this domain, as well as
						// the target node. check for nodes 2 and 4 only
						iNeighbourDomainID = poNode2->GetDomainID();
						if(iNeighbourDomainID != iDomainID)
						{
							poOperation = new ChangeConnectionOperation;
							poOperation->SetResponsibleDomainID(iNeighbourDomainID);
							poOperation->SetSourceData(iNeighbourDomainID,poNode2->GetID());
							poOperation->SetOldTargetData(poNode1->GetDomainID(),poNode1->GetID());
							poOperation->SetNewTargetData(poNewNode->GetDomainID(),poNewNode->GetID());
							m_lpoOperations.push_back(poOperation);
						}
						iNeighbourDomainID = poNode4->GetDomainID();
						if(iNeighbourDomainID != iDomainID)
						{
							poOperation = new ChangeConnectionOperation;
							poOperation->SetResponsibleDomainID(iNeighbourDomainID);
							poOperation->SetSourceData(iNeighbourDomainID,poNode4->GetID());
							poOperation->SetOldTargetData(poNode3->GetDomainID(),poNode3->GetID());
							poOperation->SetNewTargetData(poNewNode->GetDomainID(),poNewNode->GetID());
							m_lpoOperations.push_back(poOperation);
						}
						poNewNode->DisableCollision();
						bHasCollided = true;
					}
				}
				liNode2++;
			}
		}
	}
	void TopologyServer::HandleHingeJointCollisions()
	{
		list<DislocationNode*>* plpoNodes = m_poDataStructure->GetLocalNodes();
		list<DislocationSegment*>* plpoArms = NULL;
		list<DislocationNode*>::iterator liNode;
		list<DislocationSegment*>::iterator liArm1;
		list<DislocationSegment*>::iterator liArm2;
		DislocationNode* poNode1 = NULL;
		DislocationNode* poNode2 = NULL;
		DislocationNode* poNode3 = NULL;
		bool bHasCollided = false;		// because every node is allowed one collision a time
		Point oNearPoint1;
		Point oNearPoint2;
		Vector oNormal1;
		Vector oNormal2;
		Vector oBurgers1;
		Vector oBurgers2;
		double dTolerance = 1.0E-6;
		double dTemp = 0.0;
		DislocationSegment* poArm1 = NULL;
		DislocationSegment* poArm2 = NULL;
		DislocationSegment* poArm3 = NULL;
		DislocationSegment* poArm4 = NULL;
		Point oCollisionPoint;
		DislocationNode* poNewNode = NULL;
		unsigned int iNeighbourDomainID = 0;
		ChangeConnectionOperation* poOperation = NULL;
		unsigned int iDomainID = m_poDataStructure->GetDomainID();
		Vector oNewBurgersVector;
		list<DislocationSegment*>* plpoModifiedArms;
		list<DislocationSegment*>::iterator liModifiedArms;
		for(liNode = plpoNodes->begin() ; liNode != plpoNodes->end() ; liNode++)
		{
			poNode1 = (*liNode);
			if(!IsNodeCollisionValid(poNode1))
			{
				continue;
			}
			plpoArms = poNode1->GetArms();
			for(liArm1 = plpoArms->begin() ; liArm1 != plpoArms->end() ; liArm1++)
			{
				if((*liArm1)->GetEndDomainID() != iDomainID)
				{
					continue;
				}
				poNode2 = m_poDataStructure->GetNode((*liArm1)->GetEndTag());
				if(!IsNodeCollisionValid(poNode2))
				{
					continue;
				}
				liArm2 = liArm1;
				liArm2++;
				while(liArm2 != plpoArms->end())
				{
					if((*liArm2)->GetEndDomainID() != iDomainID)
					{
						liArm2++;
						continue;
					}
					poNode3 = m_poDataStructure->GetNode((*liArm2)->GetEndTag());
					if(!IsNodeCollisionValid(poNode3))
					{
						liArm2++;
						continue;
					}
					if(!IsTriangleCollapsing(poNode1,poNode2,poNode3))
					{
						liArm2++;
						continue;
					}
					if(!GetExactCollisionPoint(poNode2,poNode3,oCollisionPoint))
					{
						liArm2++;
						continue;
					}
					// now we have the collision location, move node 2 to that location
					poNode2->SetX(oCollisionPoint.GetX());
					poNode2->SetY(oCollisionPoint.GetY());
					poNode2->SetZ(oCollisionPoint.GetZ());
					// change the burgers vector of arm 1 to be the total of the 2 arms burgers vectors
					oNewBurgersVector = (*liArm1)->GetBurgersVector() + (*liArm2)->GetBurgersVector();
					(*liArm1)->SetBurgersVector(oNewBurgersVector);
					// loop over all the arms of node 3, except the one of the node 1, and change their
					// connections  to node 2
					plpoModifiedArms = poNode3->GetArms();
					for(liModifiedArms = plpoModifiedArms->begin() ; liModifiedArms != plpoModifiedArms->end() ; liModifiedArms++)
					{
						if(poNode1->IsSameAsNode((*liModifiedArms)->GetEndTag()))
						{
							continue;
						}
						// change the connections
					}
					// if them sum of the arms burgers vectors is non zero, then we are done, otherwise,
					// remove the arm from node 1 to node 2

					// create the change connections operations and process them

// 					dNewPosition[0] = oCollisionPoint.GetX();
// 					dNewPosition[1] = oCollisionPoint.GetY();
// 					dNewPosition[2] = oCollisionPoint.GetZ();
// 					MergeNode(poHome,OPCLASS_COLLISION,poNode2,poNode3,dNewPosition,&poTargetNode,&iMergeStatus,1);
// 					if((iMergeStatus & MERGE_SUCCESS) == 0)
// 					{
// 						continue;
// 					}
// 					if(poTargetNode != NULL)
// 					{
// 						poTargetNode->flags |= NO_COLLISIONS;
// 					}
					liArm2++;
				}
			}
		}
	}
	void TopologyServer::HandleTriangularLoops()
	{
// 		unsigned int i = 0;
// 		unsigned int j = 0;
// 		unsigned int k = 0;
// 		unsigned int l = 0;
// 		Node_t* poNode1 = NULL;
// 		Node_t* poNode2 = NULL;
// 		Node_t* poNode3 = NULL;
// 		Node_t* poMergeNode1 = NULL;
// 		Node_t* poMergeNode2 = NULL;
// 		Node_t* poTargetNode = NULL;
// 		bool bAreNeighbours = false;
// 		bool bFoundMergeNodes = false;
// 		double dNewPosition[3];
// 		int iMergeStatus = 0;
// 		for(i = 0 ; i < poHome->newNodeKeyPtr ; i++)
// 		{
// 			poNode1 = poHome->nodeKeys[i];
// 			if(!IsNodeCollisionValid(poNode1))
// 			{
// 				continue;
// 			}
// 			for(j = 0 ; j < poNode1->numNbrs ; j++)
// 			{
// 				if(poNode1->myTag.domainID != poNode1->nbrTag[j].domainID)
// 				{
// 					continue;
// 				}
// 				for(k = j + 1 ; k < poNode1->numNbrs ; k++)
// 				{
// 					if(poNode1->myTag.domainID != poNode1->nbrTag[k].domainID)
// 					{
// 						continue;
// 					}
// 					poNode2 = GetNodeFromTag(poHome,poNode1->nbrTag[j]);
// 					if(!IsNodeCollisionValid(poNode2))
// 					{
// 						continue;
// 					}
// 					poNode3 = GetNodeFromTag(poHome,poNode1->nbrTag[k]);
// 					if(!IsNodeCollisionValid(poNode3))
// 					{
// 						continue;
// 					}
// 					// check to see if the nodes 2 and 3 are neighbours
// 					bAreNeighbours = false;
// 					for(l = 0 ; l < poNode2->numNbrs ; l++)
// 					{
// 						if(poNode2->nbrTag[l].index == poNode3->myTag.index)
// 						{
// 							if(poNode2->nbrTag[l].domainID == poNode3->myTag.domainID)
// 							{
// 								bAreNeighbours = true;
// 								break;
// 							}
// 						}
// 					}
// 					if(!bAreNeighbours)
// 					{
// 						continue;
// 					}
// 					// now we have a triangle, make sure that one of the 3 nodes has exactly 2 arms
// 					// this node is to be removed
// 					poMergeNode1 = NULL;
// 					poMergeNode2 = NULL;
// 					bFoundMergeNodes = false;
// 					if(poNode1->numNbrs == 2)
// 					{
// 						poMergeNode1 = poNode1;
// 						if((poNode2->flags & NO_MESH_COARSEN) == 0)
// 						{
// 							poMergeNode2 = poNode2;
// 							bFoundMergeNodes = true;
// 						}
// 						else if((poNode3->flags & NO_MESH_COARSEN) == 0)
// 						{
// 							poMergeNode2 = poNode3;
// 							bFoundMergeNodes = true;
// 						}
// 					}
// 					if(!bFoundMergeNodes)
// 					{
// 						if(poNode2->numNbrs == 2)
// 						{
// 							poMergeNode1 = poNode2;
// 							if((poNode1->flags & NO_MESH_COARSEN) == 0)
// 							{
// 								poMergeNode2 = poNode1;
// 								bFoundMergeNodes = true;
// 							}
// 							else if((poNode3->flags & NO_MESH_COARSEN) == 0)
// 							{
// 								poMergeNode2 = poNode3;
// 								bFoundMergeNodes = true;
// 							}
// 						}
// 					}
// 					if(!bFoundMergeNodes)
// 					{
// 						if(poNode3->numNbrs == 2)
// 						{
// 							poMergeNode1 = poNode3;
// 							if((poNode1->flags & NO_MESH_COARSEN) == 0)
// 							{
// 								poMergeNode2 = poNode1;
// 								bFoundMergeNodes = true;
// 							}
// 							else if((poNode2->flags & NO_MESH_COARSEN) == 0)
// 							{
// 								poMergeNode2 = poNode2;
// 								bFoundMergeNodes = true;
// 							}
// 						}
// 					}
// 					if(!bFoundMergeNodes)
// 					{
// 						continue;
// 					}
// 					// now everything is set, we have the two nodes to be merged, merge and place in 
// 					// the position of the second merge node
// 					dNewPosition[0] = poMergeNode2->x;
// 					dNewPosition[1] = poMergeNode2->y;
// 					dNewPosition[2] = poMergeNode2->z;
// 			
// 					MergeNode(poHome,OPCLASS_COLLISION,poMergeNode1,poMergeNode2,dNewPosition,&poTargetNode,&iMergeStatus,1);
// 					if((iMergeStatus & MERGE_SUCCESS) == 0)
// 					{
// 						continue;
// 					}
// 					if(poTargetNode != NULL)
// 					{
// 						poTargetNode->flags |= NO_COLLISIONS;
// 						poTargetNode->flags |= NO_MESH_COARSEN;
// 					}
// 				}
// 			}
// 		}
	}
	void TopologyServer::ProcessChangeConnectionOperation(ChangeConnectionOperation* poOperation)
	{
		DislocationNodeTag* poTag = new DislocationNodeTag;
		poTag->SetDomainID(poOperation->GetSourceDomainID());
		poTag->SetNodeID(poOperation->GetSourceNodeID());
		DislocationNode* poNode = m_poDataStructure->GetNode(poTag);
		if(poNode == NULL)
		{
			printf("@ %d: error: couldn't find local node (%d,%d) while processing a change connection operation\n",m_poDataStructure->GetDomainID(),poTag->GetDomainID(),poTag->GetNodeID());
			fflush(NULL);
			exit(1);
		}
		DislocationSegment* poArm = poNode->GetArmByTag(poOperation->GetOldTargetDomainID(),poOperation->GetOldTargetNodeID());
		poArm->SetEndDomainID(poOperation->GetNewTargetDomainID());
		poArm->SetEndNodeID(poOperation->GetNewTargetNodeID());
	}
	DislocationNode* TopologyServer::SplitArm(DislocationNode* poNode1,DislocationNode* poNode2,DislocationSegment* poArm)
	{
		Vector oLine(*poNode1,*poNode2);
		Point oSplittingPoint = *poNode1 + oLine*(0.5);
		return SplitArm(poNode1,poNode2,poArm,oSplittingPoint);
	}
	DislocationNode* TopologyServer::SplitArm(DislocationNode* poNode1,DislocationNode* poNode2,DislocationSegment* poArm,const Point& oSplittingPoint)
	{
		// this function only works if the first node is a local node
		unsigned int iDomainID = m_poDataStructure->GetDomainID();
		if(poNode1->GetDomainID() != iDomainID)
		{
			return NULL;
		}
		// create the new node
		DislocationNode* poNewNode = new DislocationNode;
		*poNewNode = oSplittingPoint;
		poNewNode->SetID(m_poDataStructure->GetLocalNodeID());
		poNewNode->SetCategory(UnconstrainedNode);
		poNewNode->SetDomainID(iDomainID);
		poNewNode->SetSurfaceNormal(Vector(0.0,0.0,0.0));
		// set its arms
 		DislocationSegment* poArm1 = new DislocationSegment;
 		DislocationSegment* poArm2 = new DislocationSegment;
		poArm1->SetEndDomainID(poNode1->GetDomainID());
		poArm1->SetEndNodeID(poNode1->GetID());
		poArm2->SetEndDomainID(poNode2->GetDomainID());
		poArm2->SetEndNodeID(poNode2->GetID());
		Vector oNormal = poArm->GetSlipPlaneNormal();
		poArm1->SetSlipPlaneNormal(oNormal);
		poArm2->SetSlipPlaneNormal(oNormal);
		Vector oBurgers = poArm->GetBurgersVector();
		poArm1->SetBurgersVector(oBurgers*(-1.0));
		poArm2->SetBurgersVector(oBurgers);
		poNewNode->AddArm(poArm1);
		poNewNode->AddArm(poArm2);
		// change the connections of the arms
		DislocationSegment* poArm1Reverse = poNode1->GetArmByTag(poNode2->GetDomainID(),poNode2->GetID());
		DislocationSegment* poArm2Reverse = poNode2->GetArmByTag(poNode1->GetDomainID(),poNode1->GetID());
		poArm1Reverse->SetEndDomainID(poNewNode->GetDomainID());
		poArm1Reverse->SetEndNodeID(poNewNode->GetID());
		poArm2Reverse->SetEndDomainID(poNewNode->GetDomainID());
		poArm2Reverse->SetEndNodeID(poNewNode->GetID());
		// add the new node to the data structure
		m_poDataStructure->AddAndRegisterNode(poNewNode);
		// if the segment spans domains, then notify the other domain with the operation
		unsigned int iNeighbourDomainID = poNode2->GetDomainID();
		if(iNeighbourDomainID != iDomainID)
		{
			ChangeConnectionOperation* poOperation = new ChangeConnectionOperation;
			poOperation->SetResponsibleDomainID(iNeighbourDomainID);
			poOperation->SetSourceData(iNeighbourDomainID,poNode2->GetID());
			poOperation->SetOldTargetData(poNode1->GetDomainID(),poNode1->GetID());
			poOperation->SetNewTargetData(poNewNode->GetDomainID(),poNewNode->GetID());
			m_lpoOperations.push_back(poOperation);
		}
		return poNewNode;
	}
	bool TopologyServer::IsNodeSuperior(DislocationNode* poNode1,DislocationNode* poNode2) const
	{
		bool bIsSuperior = poNode1->IsNodeSuperior(poNode2);
		if(m_iCollisionCallsCount%2 != 0)
		{
			bIsSuperior = !bIsSuperior;
		}
		return bIsSuperior;
	}
}




