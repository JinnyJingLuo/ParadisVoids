// Ahmed M. Hussein

#include "FEMNode.h"
#include "Tools.h"
#include "FEMSolidNode.h"
#include "FEMPotentialNode.h"
#include "FEMThermoMechanicalNode.h"


using namespace SupportSystem;

namespace FEMSystem
{
	FEMNode::FEMNode()
	{
		Initialize();
	}
	FEMNode::~FEMNode()
	{
		Reset();
	}
	FEMNode::FEMNode(const FEMNode& oNode)
	{
		*this = oNode;
	}
	FEMNode::FEMNode(const double& dX,const double& dY,const double& dZ)
	{
		Initialize();
		Set(dX,dY,dZ);
	}
	FEMNode& FEMNode::operator=(const FEMNode& oNode)
	{
		if(oNode.GetType() != GetType())
		{
			return *this;
		}
		GenericNode::operator=(oNode);
		m_bOnSurface = oNode.m_bOnSurface;
		return *this;
	}
	void FEMNode::Reset()
	{
		GenericNode::Reset();
		m_bOnSurface = false;
	}
	void FEMNode::SetOnSurface(const bool& bOnSurface)
	{
		m_bOnSurface = bOnSurface;
	}
	bool FEMNode::IsOnSurface() const
	{
		return m_bOnSurface;
	}
	void FEMNode::Initialize()
	{
		GenericNode::Initialize();
		m_bOnSurface = false;
	}
	FEMNode* FEMNode::CreateNodeByType(FEMNodeType eType)
	{
		if(eType == NULLFEMNode)
		{
			return NULL;
		}
		else if(eType == PotentialFEMNode)
		{
			return new FEMPotentialNode;
		}
		else if(eType == SolidFEMNode)
		{
			return new FEMSolidNode;
		}
		else if(eType == ThermoMechanicalFEMNode)
		{
		  return new FEMThermoMechanicalNode;
		}
		return NULL;
	}
	FEMNode* FEMNode::CreateNodeByTypeIndex(const unsigned int& iIndex)
	{
		if(iIndex == 1)
		{
			return CreateNodeByType(PotentialFEMNode);
		}
		if(iIndex == 2)
		{
			return CreateNodeByType(SolidFEMNode);
		}
		if(iIndex == 3)
		{
			return CreateNodeByType(ThermoMechanicalFEMNode);
		}
		return NULL;
	}
}



