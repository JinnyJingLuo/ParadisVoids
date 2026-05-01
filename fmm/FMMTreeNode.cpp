#include "FMMTreeNode.h"


namespace FMM
{
	FMMTreeNode::FMMTreeNode()
	{
		Initialize();
	}
	FMMTreeNode::FMMTreeNode(const FMMTreeNode& oNode)
	{
		*this = oNode;
	}
	FMMTreeNode::~FMMTreeNode()
	{
		Reset();
	}
	FMMTreeNode& FMMTreeNode::operator=(const FMMTreeNode& oNode)
	{
		Reset();
		m_oBox = oNode.m_oBox;
		m_iLevel = oNode.m_iLevel;
		m_iID = oNode.m_iID;
		m_poParent = oNode.m_poParent;
		
		list< FMMTreeNode* >::const_iterator liChildren;
		for(liChildren = oNode.m_lpoChildren.begin() ; liChildren != oNode.m_lpoChildren.end() ; liChildren++)
		{
			m_lpoChildren.push_back(new FMMTreeNode(*(*liChildren)));
		}
		
		m_lpoSources = oNode.m_lpoSources;
		m_lpoTargets = oNode.m_lpoTargets;
		return *this;
	}
	void FMMTreeNode::Reset(const bool& bResetBox)
	{
		list< FMMTreeNode* >::iterator liChildren; 
		for(liChildren = m_lpoChildren.begin() ; liChildren != m_lpoChildren.end() ; liChildren++)
		{
			if((*liChildren) != NULL)
			{
				delete (*liChildren);
			}
		}
		m_iLevel = 0;
		m_iID = 0;
		m_poParent = NULL;
		m_lpoChildren.clear();
		m_lpoSources.clear();
		m_lpoTargets.clear();
		if(bResetBox)
		{
			m_oBox.Reset();
		}
	}
	void FMMTreeNode::SetBox(const double& dXRange,const double& dYRange,const double& dZRange)
	{
		SetBox(-0.5*dXRange,0.5*dXRange,-0.5*dYRange,0.5*dYRange,-0.5*dZRange,0.5*dZRange);
	}
	void FMMTreeNode::SetBox(const double& dXMin,const double& dXMax,const double& dYMin,const double& dYMax,const double& dZMin,const double& dZMax)
	{
		m_oBox.SetXMin(dXMin);
		m_oBox.SetXMax(dXMax);
		m_oBox.SetYMin(dYMin);
		m_oBox.SetYMax(dYMax);
		m_oBox.SetZMin(dZMin);
		m_oBox.SetZMax(dZMax);
	}
	void FMMTreeNode::Build(const unsigned int& iDepth,const unsigned int& iStartLevel)
	{
		Reset(false);
		// Recursively build till the required depth is reached. The level is zero based.
		// 1. set the level of this node
		m_iLevel = iStartLevel;
		if(m_iLevel == iDepth)
		{
			return;
		}
		double dXMin = m_oBox.GetXMin();
		double dXMax = m_oBox.GetXMax();
		double dYMin = m_oBox.GetYMin();
		double dYMax = m_oBox.GetYMax();
		double dZMin = m_oBox.GetZMin();
		double dZMax = m_oBox.GetZMax();
		double dMidX = 0.5*(dXMin + dXMax);
		double dMidY = 0.5*(dYMin + dYMax);
		double dMidZ = 0.5*(dZMin + dZMax);
		// 2. create a list of 8 children and set their boxes
		FMMTreeNode* poNode = NULL;
		// NNN
		poNode = new FMMTreeNode;
		poNode->SetBox(dXMin,dMidX,dYMin,dMidY,dZMin,dMidZ);
		m_lpoChildren.push_back(poNode);
		// PNN
		poNode = new FMMTreeNode;
		poNode->SetBox(dMidX,dXMax,dYMin,dMidY,dZMin,dMidZ);
		m_lpoChildren.push_back(poNode);
		// NPN
		poNode = new FMMTreeNode;
		poNode->SetBox(dXMin,dMidX,dMidY,dYMax,dZMin,dMidZ);
		m_lpoChildren.push_back(poNode);
		// PPN
		poNode = new FMMTreeNode;
		poNode->SetBox(dMidX,dXMax,dMidY,dYMax,dZMin,dMidZ);
		m_lpoChildren.push_back(poNode);
		// NNP
		poNode = new FMMTreeNode;
		poNode->SetBox(dXMin,dMidX,dYMin,dMidY,dMidZ,dZMax);
		m_lpoChildren.push_back(poNode);
		// PNP
		poNode = new FMMTreeNode;
		poNode->SetBox(dMidX,dXMax,dYMin,dMidY,dMidZ,dZMax);
		m_lpoChildren.push_back(poNode);
		// NPP
		poNode = new FMMTreeNode;
		poNode->SetBox(dXMin,dMidX,dMidY,dYMax,dMidZ,dZMax);
		m_lpoChildren.push_back(poNode);
		// PPP
		poNode = new FMMTreeNode;
		poNode->SetBox(dMidX,dXMax,dMidY,dYMax,dMidZ,dZMax);
		m_lpoChildren.push_back(poNode);
		
		// 3. set their parents and build their subtrees
		list< FMMTreeNode* >::iterator liChildren; 
		for(liChildren = m_lpoChildren.begin() ; liChildren != m_lpoChildren.end() ; liChildren++)
		{
			(*liChildren)->SetParent(this);
			(*liChildren)->Build(iDepth,m_iLevel + 1);
		}
		
		// only the root node sets the IDs
		unsigned int iStartID = 1;
		if(m_iLevel == 0)
		{
			SetIDs(iStartID);
			printf("last id is : %d\n",iStartID);
		}
	}
	void FMMTreeNode::SetIDs(unsigned int& iStartID)
	{
		m_iID = iStartID;
		iStartID++;
		list< FMMTreeNode* >::iterator liChildren; 
		for(liChildren = m_lpoChildren.begin() ; liChildren != m_lpoChildren.end() ; liChildren++)
		{
			(*liChildren)->SetIDs(iStartID);
		}
	}
	void FMMTreeNode::SetParent(FMMTreeNode* poParent)
	{
		m_poParent = poParent;
	}
	const AxisAlignedBoundingBox* FMMTreeNode::GetBox() const
	{
		return &m_oBox;
	}
	void FMMTreeNode::ClaimSources(list< Point* >* plpoSources)
	{
		// clear the current sources
		m_lpoSources.clear();
		// loop over all the given sources
		list< Point* >::iterator liSources = plpoSources->begin();
		while(liSources != plpoSources->end())
		{
			// see if the specific source is inside the node's box
			if(m_oBox.IsPointInside(*(*liSources)))
			{
				// if so, and this node is a leaf node, add it to the node's sources list
				m_lpoSources.push_back((*liSources));
				// and remove it from the input sources list
				liSources = plpoSources->erase(liSources);
			}
			else
			{
				// if not, move on
				liSources++;
			}
		}
		// if this node is a leaf node, we are done, otherwise, pass the sources 
		// list to the children to be claimed by them. At return, the current sources list 
		// should be empty
		if(!m_lpoChildren.empty())
		{
			list< FMMTreeNode* >::iterator liChildren; 
			for(liChildren = m_lpoChildren.begin() ; liChildren != m_lpoChildren.end() ; liChildren++)
			{
				(*liChildren)->ClaimSources(&m_lpoSources);
			}
		}
	}
	void FMMTreeNode::ClaimTargets(list< Point* >* plpoTargets)
	{
		// clear the current targets
		m_lpoTargets.clear();
		// loop over all the given targets
		list< Point* >::iterator liTargets = plpoTargets->begin();
		while(liTargets != plpoTargets->end())
		{
			// see if the specific target is inside the node's box
			if(m_oBox.IsPointInside(*(*liTargets)))
			{
				// if so, and this node is a leaf node, add it to the node's targets list
				m_lpoTargets.push_back((*liTargets));
				// and remove it from the input sources list
				liTargets = plpoTargets->erase(liTargets);
			}
			else
			{
				// if not, move on
				liTargets++;
			}
		}
		// if this node is a leaf node, we are done, otherwise, pass the targets 
		// list to the children to be claimed by them. At return, the current targets list 
		// should be empty
		if(!m_lpoChildren.empty())
		{
			list< FMMTreeNode* >::iterator liChildren; 
			for(liChildren = m_lpoChildren.begin() ; liChildren != m_lpoChildren.end() ; liChildren++)
			{
				(*liChildren)->ClaimTargets(&m_lpoTargets);
			}
		}
	}
	void FMMTreeNode::Initialize()
	{
		m_oBox.Reset();
		m_iLevel = 0;
		m_iID = 0;
		m_poParent = NULL;
		m_lpoChildren.clear();
		m_lpoSources.clear();
		m_lpoTargets.clear();
	}
}		
		



