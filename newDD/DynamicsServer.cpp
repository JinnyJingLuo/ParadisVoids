#include "DynamicsServer.h"
#include "cmath"
#include "mpi.h"

namespace DynamicsSystem
{
	DynamicsServer* DynamicsServer::m_poInstance = NULL;
	DynamicsServer* DynamicsServer::GetInstance()
	{
		if(m_poInstance == NULL)
		{
			m_poInstance = new DynamicsServer;
		}
		return m_poInstance;
	}
	DynamicsServer::~DynamicsServer()
	{
		Reset();
	}
	void DynamicsServer::SetDataStructure(MainDataStructure* poDataStructure)
	{
		m_poDataStructure = poDataStructure;
	}
	void DynamicsServer::ApplyLoad()
	{
		// for now, constant strain rate loading, at a strain rate of 50 in the positive z direction
		// is the only available load type, also the stiffness is always 2.1E11 Pa;
		double dStrainRate = 50.0;
		double dStiffness = 2.1E11;
		Matrix* poPlasticStrainIncrement = m_poDataStructure->GetPlasticStrainIncrement();
		Matrix oTotalStrainIncrement(3,3);
		oTotalStrainIncrement.Set(3,3,1.0);
		oTotalStrainIncrement = oTotalStrainIncrement*dStrainRate*m_poDataStructure->GetTimeStep();
		Matrix oStressIncrement = (oTotalStrainIncrement - *poPlasticStrainIncrement)*dStiffness;
		Matrix* poAppliedStress = m_poDataStructure->GetAppliedStress();
		*poAppliedStress = *poAppliedStress + oStressIncrement;
		Matrix* poTotalStrain = m_poDataStructure->GetTotalStrain();
		*poTotalStrain = *poTotalStrain + oTotalStrainIncrement;
	}
	void DynamicsServer::ComputeForces()
	{
		ResetForces();
		ComputeAppliedForces();
		ComputeSelfForces();
		ComputeInteractionsForces();
	}
	void DynamicsServer::ComputeVelocities()
	{
		list<DislocationNode*>* plpoNodes = m_poDataStructure->GetLocalNodes();
		list<DislocationNode*>::iterator liNodes;
		for(liNodes = plpoNodes->begin() ; liNodes != plpoNodes->end() ; liNodes++)
		{
			(*liNodes)->SetVelocity(CalculateFCCVelocity((*liNodes)));
		}
		CommunicateVelocities();
	}
	void DynamicsServer::ComputeTimeStep()
	{
		// get the maximum velocity magnitude
		list<DislocationNode*>* plpoNodes = m_poDataStructure->GetLocalNodes();
		list<DislocationNode*>::iterator liNodes;
		double dMaximumVelocitySquared = 0.0;
		double dTemp = 0.0;
		for(liNodes = plpoNodes->begin() ; liNodes != plpoNodes->end() ; liNodes++)
		{
			dTemp = (*liNodes)->GetVelocity().LengthSquared();
			if(dTemp > dMaximumVelocitySquared)
			{
				dMaximumVelocitySquared = dTemp;
			}
		}
		plpoNodes = m_poDataStructure->GetRemoteNodes();
		for(liNodes = plpoNodes->begin() ; liNodes != plpoNodes->end() ; liNodes++)
		{
			dTemp = (*liNodes)->GetVelocity().LengthSquared();
			if(dTemp > dMaximumVelocitySquared)
			{
				dMaximumVelocitySquared = dTemp;
			}
		}
		// get the maximum of all processes
		MPI_Allreduce(&dMaximumVelocitySquared,&dTemp,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
		dMaximumVelocitySquared = dTemp;
		double dMaximumVelocity = sqrt(dMaximumVelocitySquared);
		double dTimeStep = 0.0;
		double dMaxFlightDistance = 25.0;
		double dMaximumTimeStep = 1.0E-7;
		double dTolerance = 1.0E-6;
		if(dMaximumVelocity <= dTolerance)
		{
			dTimeStep = dMaximumTimeStep;
		}
		else
		{
			dTimeStep = dMaxFlightDistance/dMaximumVelocity;
		}
		if(dTimeStep > dMaximumTimeStep)
		{
			dTimeStep = dMaximumTimeStep;
		}
		m_poDataStructure->SetTimeStep(dTimeStep);
		m_poDataStructure->SetTime(m_poDataStructure->GetTime() + dTimeStep);
	}
	void DynamicsServer::IntegrateMotion()
	{
		double dTimeStep = m_poDataStructure->GetTimeStep();
		list<DislocationNode*>* plpoNodes = m_poDataStructure->GetLocalNodes();
		list<DislocationNode*>::iterator liNodes;
		Vector oDisplacement;
		for(liNodes = plpoNodes->begin() ; liNodes != plpoNodes->end() ; liNodes++)
		{
			oDisplacement = (*liNodes)->GetVelocity()*dTimeStep;
			(*liNodes)->Shift(oDisplacement.GetX(),oDisplacement.GetY(),oDisplacement.GetZ());
		}
		plpoNodes = m_poDataStructure->GetRemoteNodes();
		for(liNodes = plpoNodes->begin() ; liNodes != plpoNodes->end() ; liNodes++)
		{
			oDisplacement = (*liNodes)->GetVelocity()*dTimeStep;
			(*liNodes)->Shift(oDisplacement.GetX(),oDisplacement.GetY(),oDisplacement.GetZ());
		}
	}
	DynamicsServer::DynamicsServer()
	{
		Initialize();
	}
	void DynamicsServer::Reset()
	{
		m_poDataStructure = NULL;
	}
	void DynamicsServer::Initialize()
	{
		m_poDataStructure = NULL;
	}
	void DynamicsServer::ResetForces()
	{
		list<DislocationNode*>* plpoNodes = m_poDataStructure->GetLocalNodes();
		list<DislocationNode*>::iterator liNodes;
		for(liNodes = plpoNodes->begin() ; liNodes != plpoNodes->end() ; liNodes++)
		{
			(*liNodes)->ResetForces();
		}
		plpoNodes = m_poDataStructure->GetRemoteNodes();
		for(liNodes = plpoNodes->begin() ; liNodes != plpoNodes->end() ; liNodes++)
		{
			(*liNodes)->ResetForces();
		}
	}
	void DynamicsServer::ComputeAppliedForces()
	{
		list<DislocationNode*>* plpoNodes = m_poDataStructure->GetLocalNodes();
		list<DislocationSegment*>* plpoArms = NULL;
		list<DislocationNode*>::iterator liNodes;
		list<DislocationSegment*>::iterator liArms;
		Vector oForce;
		DislocationNode* poNode = NULL;
		DislocationNode* poNeighbour = NULL;
		Matrix* poAppliedStress = m_poDataStructure->GetAppliedStress();
		for(liNodes = plpoNodes->begin() ; liNodes != plpoNodes->end() ; liNodes++)
		{
			poNode = (*liNodes);
			plpoArms = poNode->GetArms();
			for(liArms = plpoArms->begin() ; liArms != plpoArms->end() ; liArms++)
			{
				poNeighbour = m_poDataStructure->GetNode((*liArms)->GetEndTag());
				if(poNeighbour == NULL)
				{
					printf("@ %d : error : neighbour node (%d,%d) not found when computing applied forces\n",m_poDataStructure->GetDomainID(),(*liArms)->GetEndDomainID(),(*liArms)->GetEndNodeID());
					exit(1);
				}
				if(!poNode->IsNodeSuperior(poNeighbour))
				{
					continue;
				}
				oForce = CalculatePeachKoehlerForce(poNode,poNeighbour,(*liArms),poAppliedStress)*0.5;
				poNode->AddForce(oForce);
				poNeighbour->AddForce(oForce);
			}
		}
	}
	void DynamicsServer::ComputeSelfForces()
	{
// 		f1[0] = 0.0;
// 		f1[1] = 0.0;
// 		f1[2] = 0.0;
// 	
// 		f2[0] = 0.0;
// 		f2[1] = 0.0;
// 		f2[2] = 0.0;
// 		
// 		real8 tx, ty, tz, L, La, S;
// 		real8 bs, bs2, bex, bey, bez, be2, fL, ft;
// 		
// 		GetUnitVector(1, x1, y1, z1, x2, y2, z2, &tx, &ty, &tz, &L);
// 		if(L < 1.0E-6)
// 		{
// 			return;
// 		}
// 		
// 		bs = bx*tx + by*ty + bz*tz;
// 		bex = bx-bs*tx; bey = by-bs*ty; bez=bz-bs*tz;
// 		be2 = (bex*bex+bey*bey+bez*bez);
// 		bs2 = bs*bs;
// 	
// 		La=sqrt(L*L+a*a);
// 		S = (-(2*NU*La+(1-NU)*a*a/La-(1+NU)*a)/L + (NU*log((La+L)/a)-(1-NU)*0.5*L/La))*MU/4/M_PI/(1-NU)*bs;
// 	
// 	
// 		/* Ecore = MU/(4*pi) log(a/a0) */
// 		/* Account for the self force due to the core energy change when 
// 		   the core radius changes from a0-->a, M. Tang, 7/19/2004 */ 
// 		fL = -Ecore*(bs2+be2/(1-NU));
// 		ft =  Ecore*2*bs*NU/(1-NU);   
// 		
// 		f2[0] = bex*(S+ft) + fL*tx;
// 		f2[1] = bey*(S+ft) + fL*ty;
// 		f2[2] = bez*(S+ft) + fL*tz;
// 	
// 		f1[0] = -f2[0];
// 		f1[1] = -f2[1];
// 		f1[2] = -f2[2];
	}
	void DynamicsServer::ComputeInteractionsForces()
	{
	
	}
	Vector DynamicsServer::CalculatePeachKoehlerForce(DislocationNode* poNode,DislocationNode* poNeighbour,DislocationSegment* poArm,Matrix* poStress)
	{
		Vector oSigmaB;
		Vector oBurgersVector = poArm->GetBurgersVector();
		double dBx = oBurgersVector.GetX();
		double dBy = oBurgersVector.GetY();
		double dBz = oBurgersVector.GetZ();
		oSigmaB.SetX(poStress->Get(1,1)*dBx + poStress->Get(1,2)*dBy + poStress->Get(1,3)*dBz);
		oSigmaB.SetY(poStress->Get(2,1)*dBx + poStress->Get(2,2)*dBy + poStress->Get(2,3)*dBz);
		oSigmaB.SetZ(poStress->Get(3,1)*dBx + poStress->Get(3,2)*dBy + poStress->Get(3,3)*dBz);
		Vector oLine(*poNode,*poNeighbour);
		Vector oForce = oSigmaB^oLine;
		return oForce;
	}
	Vector DynamicsServer::CalculateFCCVelocity(DislocationNode* poNode)
	{
		// if node is pinned or has sessile bugers vector, set the velocity to zero and return
		Vector oVelocity(0.0,0.0,0.0);
		if(poNode->GetCategory() == PinnedNode)
		{
			return oVelocity;
		}
		
		unsigned int iConstraintType = 0;
		Vector oConstraintVector;
		poNode->GetDynamicConstraint(iConstraintType,oConstraintVector);
		// node is dynamically fixed, set the velocity to zero and return
		if(iConstraintType == 3)
		{
			return oVelocity;
		}
		
		// since this is FCC mobility, nodes that are not on {111} slip planes cannot move, 
		// to be more accurate, they can only move along their arms, but this is disabled for now
		if(iConstraintType == 1)
		{
			if(!MainDataStructure::Is111Vector(oConstraintVector))
			{
				return oVelocity;
			}
		}
		
		// calculate the node's mobility
		// find total dislocation length times drag coefficent (L x B)
		double dTolerance = 1.0E-6;
		double dLB = 0.0;
		unsigned int i = 0;
		DislocationNode* poNeighbour = NULL;
		Vector oLine;
		double dLength = 0.0;
		Vector oBurgersVector;
		double dProjection = 0.0;
		double dMobility = 0.0;
		// fix these parameters for now
		double dScrewSegmentMobility = 4000.0;
		double dEdgeSegmentMobility  = 4000.0;
		list<DislocationSegment*>* plpoArms = poNode->GetArms();
		list<DislocationSegment*>::iterator liArms;
		for(liArms = plpoArms->begin() ; liArms != plpoArms->end() ; liArms++)
		{
			poNeighbour = m_poDataStructure->GetNode((*liArms)->GetEndTag());
			if(poNeighbour == NULL)
			{
				printf("@ %d: couldn't find neighbour node (%d,%d) in velocity calculation\n",m_poDataStructure->GetDomainID(),(*liArms)->GetEndTag()->GetDomainID(),(*liArms)->GetEndTag()->GetNodeID());
				exit(1);
			}
			if(poNode->GetCategory() == SurfaceNode)
			{
				if(poNeighbour->GetCategory() == SurfaceNode)
				{
					continue;
				}
			}
			oLine.SetByPoints(*poNode,*poNeighbour);
			dLength = oLine.Length();
			if(dLength < dTolerance)
			{
				// skip zero length segments
				continue;
			}
			oLine.Normalize();
			oBurgersVector = (*liArms)->GetBurgersVector();
			oBurgersVector.Normalize();
			dProjection = fabs(oBurgersVector*oLine);
			dMobility = dEdgeSegmentMobility + (dScrewSegmentMobility - dEdgeSegmentMobility)*dProjection;
			dLB = dLB + (dLength/dMobility);
		}
		
		// if dLB is zero, zero the node's velocity
		if(dLB < dTolerance)
		{
			return oVelocity;
		}
		dLB = 0.5*dLB;
	
		// velocity is proportional to total force per unit length
		oVelocity = poNode->GetForce()*(1.0/dLB);
		
		// now we have the velocity, project it according to the dynamic constraint, in case of
		// no dynamic constraints (constraint type 0), nothing needs to be done
		if(iConstraintType == 1)
		{
			// node moves normal to the constraint vector, subtract the component of the
			// velocity vector in the direction of the constraint vector
			oVelocity = oVelocity - oConstraintVector*(oVelocity*oConstraintVector);
		}
		else if(iConstraintType == 2)
		{
			// node moves along the constraint vector, project the obtained velocity in the 
			// direction of the constraint vector
			oVelocity = oConstraintVector*(oVelocity*oConstraintVector);
		}
		
		// if this node is a surface node, force it to move on the surface
		if(poNode->GetCategory() == SurfaceNode)
		{		
			Vector oSurfaceNormal = poNode->GetSurfaceNormal();
			if(iConstraintType == 1)
			{
				// node moves on a plane
				Vector oTargetVelocityDirection = oConstraintVector^oSurfaceNormal;
				oTargetVelocityDirection.Normalize();
				oVelocity = oTargetVelocityDirection*(oVelocity*oTargetVelocityDirection);
			}
			else if(iConstraintType == 2)
			{
				if(fabs(oSurfaceNormal*oConstraintVector) > dTolerance)
				{
					// velocity direction has a component in normal to the surface, which is 
					// inconsistent with the surface motion constraint, the node should not move at all
					oVelocity.Set(0.0,0.0,0.0);
				}
			}
			else
			{
				// node is free to move, just subtract the suface normal component of the velocity
				// vector from the velocity vector
				oVelocity = oVelocity - oSurfaceNormal*(oVelocity*oSurfaceNormal);
			}
		}
		return oVelocity;
	}
	void DynamicsServer::CommunicateVelocities()
	{
 		vector< list<double> > vldData;
 		unsigned int iDomainsCount = m_poDataStructure->GetDomainsCount();
 		unsigned int iDomainID = m_poDataStructure->GetDomainID();
 		PackVelocities(&vldData);
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
 			MPI_Irecv(&piBufferSizes[i],1,MPI_INT,i,VelocityMessageLengthTag,MPI_COMM_WORLD,&poIncomingRequests[i]);
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
			MPI_Isend(&iTemp,1,MPI_INT,i,VelocityMessageLengthTag,MPI_COMM_WORLD,&poOutgoingRequests[i]);
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
			MPI_Irecv(vpdIncomingBuffers[i],piBufferSizes[i],MPI_DOUBLE,i,VelocityMessageDataTag,MPI_COMM_WORLD,&poIncomingRequests[i]);
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
			MPI_Isend(vpdOutgoingBuffers[i],iTemp,MPI_DOUBLE,i,VelocityMessageDataTag,MPI_COMM_WORLD,&poOutgoingRequests[i]);
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
		
		UnpackVelocities(&vldData);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	void DynamicsServer::PackVelocities(vector< list<double> >* pvldData)
	{
		unsigned int iDomainsCount = m_poDataStructure->GetDomainsCount();
		unsigned int iDomainID = m_poDataStructure->GetDomainID();
		pvldData->clear();
		pvldData->resize(iDomainsCount);
		vector<bool> vbIsNodeAdded;
		vbIsNodeAdded.resize(iDomainsCount);
		// initialize the lists
		unsigned int i = 0;
		for(i = 0 ; i < iDomainsCount ; i++)
		{
			pvldData->at(i).clear();
		}
		
		// pack the data
		list<DislocationNode*>* plpoNodes = m_poDataStructure->GetLocalNodes();
		list<DislocationNode*>::iterator liNodes;
		list<DislocationSegment*>::iterator liArms;
		list<DislocationSegment*>* plpoArms = NULL;
		unsigned int iTargetDomain = 0;
		for(liNodes = plpoNodes->begin() ; liNodes != plpoNodes->end() ; liNodes++)
		{
			for(i = 0 ; i < iDomainsCount ; i++)
			{
				vbIsNodeAdded[i] = false;
			}
			plpoArms = (*liNodes)->GetArms();
			for(liArms = plpoArms->begin() ; liArms != plpoArms->end() ; liArms++)
			{
				iTargetDomain = (*liArms)->GetEndDomainID();
				if(iTargetDomain != iDomainID)
				{
					if(!vbIsNodeAdded[iTargetDomain])
					{
						(*liNodes)->PackVelocity(&(pvldData->at(iTargetDomain)));
						vbIsNodeAdded[iTargetDomain] = true;
					}
				}
			}
		}
		vbIsNodeAdded.clear();
	}
	void DynamicsServer::UnpackVelocities(vector< list<double> >* pvldData)
	{
		unsigned int iDomainsCount = m_poDataStructure->GetDomainsCount();
		unsigned int i = 0;
		list<double>* pldDomainData = NULL;
		DislocationNode* poNode = NULL;
		Vector oVelocity;
		DislocationNodeTag oTag;
		for(i = 0 ; i < iDomainsCount ; i++)
		{
			pldDomainData = &pvldData->at(i);
			while(!pldDomainData->empty())
			{
				oTag = DislocationNode::UnpackVelocity(pldDomainData,oVelocity);
				poNode = m_poDataStructure->GetNode(&oTag);
				if(poNode == NULL)
				{
					printf("@ %d : error : remote node (%d,%d) not found when unpacking velocities\n",m_poDataStructure->GetDomainID(),oTag.GetDomainID(),oTag.GetNodeID());
					exit(1);
				}
				poNode->SetVelocity(oVelocity);
			}
		}
	}
}


