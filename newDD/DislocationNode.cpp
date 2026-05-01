// Ahmed M. Hussein

#include "DislocationNode.h"
#include "DislocationSegment.h"
#include "MainDataStructure.h"

namespace DislocationSystem
{
	DislocationNode::DislocationNode()
	{
		Initialize();
	}
	DislocationNode::DislocationNode(const double& dX,const double& dY,const double& dZ)
	{
		Initialize();
		Set(dX,dY,dZ);
	}
	DislocationNode::DislocationNode(const DislocationNode& oNode)
	{
		*this = oNode;
	}
	DislocationNode::DislocationNode(const CategorizedGenericNode& oNode)
	{
		*this = oNode;
	}
	DislocationNode::DislocationNode(const GenericNode& oNode)
	{
		*this = oNode;
	}
	DislocationNode::DislocationNode(const Point& oNode)
	{
		*this = oNode;
	}
	DislocationNode::~DislocationNode()
	{
		Reset();
	}
	DislocationNode& DislocationNode::operator=(const DislocationNode& oNode)
	{
		Reset();
		m_dX = oNode.m_dX;
		m_dY = oNode.m_dY;
		m_dZ = oNode.m_dZ;
		m_iID = oNode.m_iID;
		m_iCategory = oNode.m_iCategory;
		m_oSurfaceNormal = oNode.m_oSurfaceNormal;
		m_oForce = oNode.m_oForce;
		m_oVelocity = oNode.m_oVelocity;
		list<DislocationSegment*>* plpoArms = ((DislocationNode&)oNode).GetArms();
		list<DislocationSegment*>::iterator liArms;
		for(liArms = plpoArms->begin() ; liArms != plpoArms->end() ; liArms++)
		{
			m_lpoArms.push_back(new DislocationSegment(*(*liArms)));
		}
		m_bAllowedToCollide = oNode.m_bAllowedToCollide;
		return *this;
	}
	DislocationNode& DislocationNode::operator=(const CategorizedGenericNode& oNode)
	{
		Initialize();
		m_dX = oNode.GetX();
		m_dY = oNode.GetY();
		m_dZ = oNode.GetZ();
		m_iID = oNode.GetID();
		m_iCategory = oNode.GetCategory();
		return *this;
	}
	DislocationNode& DislocationNode::operator=(const GenericNode& oNode)
	{
		Initialize();
		m_dX = oNode.GetX();
		m_dY = oNode.GetY();
		m_dZ = oNode.GetZ();
		m_iID = oNode.GetID();
		return *this;
	}
	DislocationNode& DislocationNode::operator=(const Point& oNode)
	{
		Initialize();
		m_dX = oNode.GetX();
		m_dY = oNode.GetY();
		m_dZ = oNode.GetZ();
		return *this;
	}
	void DislocationNode::Reset()
	{
		m_dX = 0.0;
		m_dY = 0.0;
		m_dZ = 0.0;
		m_iID = 0;
		m_iCategory = 0;
		m_iDomainID = 0;
		m_oSurfaceNormal.Set(0.0,0.0,0.0);
		m_oForce.Set(0.0,0.0,0.0);
		m_oVelocity.Set(0.0,0.0,0.0);
		list<DislocationSegment*>::iterator liArms;
		for(liArms = m_lpoArms.begin() ; liArms != m_lpoArms.end() ; liArms++)
		{
			if((*liArms) != NULL)
			{
				delete (*liArms);
			}
		}
		m_lpoArms.clear();
		m_bAllowedToCollide = true;
	}
	void DislocationNode::SetDomainID(const unsigned int& iID)
	{
		m_iDomainID = iID;
	}
    unsigned int DislocationNode::GetDomainID() const
    {
    	return m_iDomainID;
    }
	void DislocationNode::SetSurfaceNormal(const Vector& oNormal)
	{
		m_oSurfaceNormal = oNormal;
	}
	Vector DislocationNode::GetSurfaceNormal() const
	{
		return m_oSurfaceNormal;
	}
	void DislocationNode::SetForce(const Vector& oForce)
	{
		m_oForce = oForce;
	}
    Vector DislocationNode::GetForce() const
    {
    	return m_oForce;
    }
    void DislocationNode::SetVelocity(const Vector& oVelocity)
    {
    	m_oVelocity = oVelocity;
    }
    Vector DislocationNode::GetVelocity() const
    {
    	return m_oVelocity;
    }
	list<DislocationSegment*>* DislocationNode::GetArms()
	{
		return &m_lpoArms;
	}
	void DislocationNode::AddArm(DislocationSegment* poArm)
	{
		if(poArm != NULL)
		{
			m_lpoArms.push_back(poArm);
		}
	}
	void DislocationNode::RemoveArmByTag(DislocationNodeTag* poTag)
	{
		RemoveArmByTag(poTag->GetDomainID(),poTag->GetNodeID());
	}
	void DislocationNode::RemoveArmByTag(const unsigned int& iDomainID,const unsigned int& iNodeID)
	{
		list<DislocationSegment*>::iterator liArms = m_lpoArms.begin();
		while(liArms != m_lpoArms.end())
		{
			if(((*liArms)->GetEndNodeID() == iNodeID) && ((*liArms)->GetEndDomainID() == iDomainID))
			{
				delete (*liArms);
				liArms = m_lpoArms.erase(liArms);
			}
			else
			{
				liArms++;
			}
		}
	}
	void DislocationNode::RemoveArmByEndNode(DislocationNode* poNode)
	{
		if(poNode == NULL)
		{
			return;
		}
		RemoveArmByTag(poNode->m_iDomainID,poNode->m_iID);
	}
	DislocationSegment* DislocationNode::GetArmByTag(const unsigned int& iDomainID,const unsigned int& iNodeID)
	{
		list<DislocationSegment*>::iterator liArms;
		for(liArms = m_lpoArms.begin() ; liArms != m_lpoArms.end() ; liArms++)
		{
			if(((*liArms)->GetEndNodeID() == iNodeID) && ((*liArms)->GetEndDomainID() == iDomainID))
			{
				return (*liArms);
			}
		}
		return NULL;
	}
	void DislocationNode::Write(FILE* fpFile) const
	{
		fprintf(fpFile,"%d,%d\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%d\t%d\n",m_iDomainID,m_iID,m_dX,m_dY,m_dZ,m_oSurfaceNormal.GetX(),m_oSurfaceNormal.GetY(),m_oSurfaceNormal.GetZ(),(unsigned int)m_lpoArms.size(),m_iCategory);
		list<DislocationSegment*>::const_iterator liArms;
		for(liArms = m_lpoArms.begin() ; liArms != m_lpoArms.end() ; liArms++)
		{
			(*liArms)->Write(fpFile);
		}
	}
	void DislocationNode::Pack(list<double>* pldData) const
	{
		pldData->push_back((double)m_iDomainID);
		pldData->push_back((double)m_iID);
		pldData->push_back((double)m_iCategory);
		pldData->push_back(m_dX);
		pldData->push_back(m_dY);
		pldData->push_back(m_dZ);
		pldData->push_back(m_oSurfaceNormal.GetX());
		pldData->push_back(m_oSurfaceNormal.GetY());
		pldData->push_back(m_oSurfaceNormal.GetZ());
		// no need to pach the force
		pldData->push_back(m_oVelocity.GetX());
		pldData->push_back(m_oVelocity.GetY());
		pldData->push_back(m_oVelocity.GetZ());
		pldData->push_back((double)((unsigned int)m_lpoArms.size()));
		list<DislocationSegment*>::const_iterator liArms;
		for(liArms = m_lpoArms.begin() ; liArms != m_lpoArms.end() ; liArms++)
		{
			(*liArms)->Pack(pldData);
		}
		if(m_bAllowedToCollide)
		{
			pldData->push_back(1.0);
		}
		else
		{
			pldData->push_back(-1.0);
		}
	}
	void DislocationNode::Unpack(list<double>* pldData)
	{
		m_iDomainID = (unsigned int)pldData->front();
		pldData->pop_front();
		m_iID = (unsigned int)pldData->front();
		pldData->pop_front();
		m_iCategory = (unsigned int)pldData->front();
		pldData->pop_front();
		m_dX = pldData->front();
		pldData->pop_front();
		m_dY = pldData->front();
		pldData->pop_front();
		m_dZ = pldData->front();
		pldData->pop_front();
		m_oSurfaceNormal.SetX(pldData->front());
		pldData->pop_front();
		m_oSurfaceNormal.SetY(pldData->front());
		pldData->pop_front();
		m_oSurfaceNormal.SetZ(pldData->front());
		pldData->pop_front();
		m_oVelocity.SetX(pldData->front());
		pldData->pop_front();
		m_oVelocity.SetY(pldData->front());
		pldData->pop_front();
		m_oVelocity.SetZ(pldData->front());
		pldData->pop_front();
		unsigned int iArmsCount = (unsigned int)pldData->front();
		pldData->pop_front();
		unsigned int i = 0;
		DislocationSegment* poArm = NULL;
		for(i = 0 ; i < iArmsCount ; i++)
		{
			poArm = new DislocationSegment;
			poArm->Unpack(pldData);
			AddArm(poArm);
		}
		double dTemp = pldData->front();
		pldData->pop_front();
		if(dTemp > 0.0)
		{
			m_bAllowedToCollide = true;
		}
		else
		{
			m_bAllowedToCollide = false;
		}
		// the force wasn't packed, so reset it
		m_oForce.Set(0.0,0.0,0.0);
	}
	void DislocationNode::PackVelocity(list<double>* pldData) const
	{
		pldData->push_back((double)m_iDomainID);
		pldData->push_back((double)m_iID);
		pldData->push_back(m_oVelocity.GetX());
		pldData->push_back(m_oVelocity.GetY());
		pldData->push_back(m_oVelocity.GetZ());
	}
	DislocationNodeTag DislocationNode::UnpackVelocity(list<double>* pldData,Vector& oVelocity)
	{
		unsigned int iDomainID = (unsigned int)pldData->front();
		pldData->pop_front();
		unsigned int iNodeID = (unsigned int)pldData->front();
		pldData->pop_front();
		oVelocity.SetX(pldData->front());
		pldData->pop_front();
		oVelocity.SetY(pldData->front());
		pldData->pop_front();
		oVelocity.SetZ(pldData->front());
		pldData->pop_front();
		DislocationNodeTag oTag;
		oTag.SetDomainID(iDomainID);
		oTag.SetNodeID(iNodeID);
		return oTag;
	}
	void DislocationNode::ResetForces()
	{
		m_oForce.Set(0.0,0.0,0.0);
	}
	void DislocationNode::AddForce(const Vector& oForceIncrement)
	{
		m_oForce = m_oForce + oForceIncrement;
	}
	bool DislocationNode::IsNodeSuperior(DislocationNode* poNode) const
	{
		if(MainDataStructure::IsTagLower(m_iDomainID,m_iID,poNode->m_iDomainID,poNode->m_iID))
		{
			return true;
		}
		return false;
	}
	void DislocationNode::GetDynamicConstraint(unsigned int& iConstraintType,Vector& oConstraintVector) const
	{
		oConstraintVector.Set(0.0,0.0,0.0);
		iConstraintType = 0;
		double dTolerance = 1.0E-6;
		unsigned int i = 0;
		// get the node dynamic constraint
		unsigned int iArmsCount = (unsigned int)m_lpoArms.size();
		if(iArmsCount == 0)
		{
			oConstraintVector.Set(0.0,0.0,0.0);
			iConstraintType = 0;
		}
		else if(iArmsCount == 1)
		{
			oConstraintVector = m_lpoArms.front()->GetSlipPlaneNormal();
			iConstraintType = 1;
		}
		else if(iArmsCount == 2)
		{
			Vector oNormal1 = m_lpoArms.front()->GetSlipPlaneNormal();
			Vector oNormal2 = m_lpoArms.back()->GetSlipPlaneNormal();
			oNormal1.Normalize();
			oNormal2.Normalize();
			Vector oCrossProduct = oNormal1^oNormal2;
			if(oCrossProduct.Length() < dTolerance)
			{
				iConstraintType = 1;
				oConstraintVector = oNormal1;
			}
			else
			{
				iConstraintType = 2;
				oConstraintVector = oCrossProduct;
			}
		}
		else
		{
			list<Vector> loNormals;
			list<Vector>::iterator liNormals;
			loNormals.push_back(m_lpoArms.front()->GetSlipPlaneNormal());
			loNormals.front().Normalize();
			Vector oTempNormal;
			Vector oCrossProduct;
			bool bAddNormal = false;
			list<DislocationSegment*>::const_iterator liArms;
			for(liArms = m_lpoArms.begin() ; liArms != m_lpoArms.end() ; liArms++)
			{
				oTempNormal = (*liArms)->GetSlipPlaneNormal();
				oTempNormal.Normalize();
				bAddNormal = false;
				for(liNormals = loNormals.begin() ; liNormals != loNormals.end() ; liNormals++)
				{
					oCrossProduct = oTempNormal^(*liNormals);
					if(oCrossProduct.Length() > dTolerance)
					{
						bAddNormal = true;
						break;
					}
				}
				if(bAddNormal)
				{
					loNormals.push_back(oTempNormal);
				}
			}
			unsigned int iDynamicConstraintsCount = (unsigned int)loNormals.size();
			if(iDynamicConstraintsCount == 1)
			{
				iConstraintType = 1;
				oConstraintVector = loNormals.front();
			}
			else if(iDynamicConstraintsCount == 2)
			{
				iConstraintType = 2;
				oConstraintVector = loNormals.front()^loNormals.back();
			}
			else
			{
				iConstraintType = 3;
				oConstraintVector.Set(0.0,0.0,0.0);
			}
		}
		oConstraintVector.Normalize();
	}
	bool DislocationNode::IsAllowedToCollide() const
	{
		return m_bAllowedToCollide;
	}
	void DislocationNode::EnableCollision()
	{
		m_bAllowedToCollide = true;
	}
	void DislocationNode::DisableCollision()
	{
		m_bAllowedToCollide = false;
	}
	unsigned int DislocationNode::GetNeighboursCount() const
	{
		return (unsigned int)m_lpoArms.size();
	}
	bool DislocationNode::IsSameNode(DislocationNode* poNode) const
	{
		if(m_iID == poNode->m_iID)
		{
			if(m_iDomainID == poNode->m_iDomainID)
			{
				return true;
			}
		}
		return false;
	}
	bool DislocationNode::IsSameAsNode(DislocationNodeTag* poTag) const
	{
		if(m_iID == poTag->GetNodeID())
		{
			if(m_iDomainID == poTag->GetDomainID())
			{
				return true;
			}
		}
		return false;
	}
	bool DislocationNode::IsMasterNode() const
	{
		// a master node is a node that is superior to all of its neighbours
		list<DislocationSegment*>::const_iterator liArms;
		for(liArms = m_lpoArms.begin() ; liArms != m_lpoArms.end() ; liArms++)
		{
			if(!MainDataStructure::IsTagLower(m_iDomainID,m_iID,(*liArms)->GetEndDomainID(),(*liArms)->GetEndNodeID()))
			{
				return false;
			}
		}
		return true;
	}
	void DislocationNode::DislocationNode::Initialize()
	{
		m_dX = 0.0;
		m_dY = 0.0;
		m_dZ = 0.0;
		m_iID = 0;
		m_iCategory = 0;
		m_iDomainID = 0;
		m_oSurfaceNormal.Set(0.0,0.0,0.0);
		m_oForce.Set(0.0,0.0,0.0);
		m_oVelocity.Set(0.0,0.0,0.0);
		m_lpoArms.clear();
		m_bAllowedToCollide = true;
	}
}



