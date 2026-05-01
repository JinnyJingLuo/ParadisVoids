// Ahmed M. Hussein

#include "FEMElement.h"
#include "FEMSolidElement.h"
#include "FEMPotentialElement.h"
#include "FEMThermoMechanicalElement.h"

namespace FEMSystem
{
	FEMElement::~FEMElement()
	{
		Reset();
	}
	void FEMElement::Reset()
	{
		if(m_poElementGeometry != NULL)
		{
			delete m_poElementGeometry;
		}
		m_voGaussPoints.clear();
		Initialize();
	}
	FEMElement& FEMElement::operator=(const FEMElement& oElement)
	{
		Reset();
		m_poElementGeometry = oElement.m_poElementGeometry->Clone();
		m_poMaterial = oElement.m_poMaterial;
		m_voGaussPoints = oElement.m_voGaussPoints;
		return *this;
	}
	vector<FEMNode*>* FEMElement::GetNodes()
	{
		return m_poElementGeometry->GetNodes();
	}
	FEMNode* FEMElement::GetNode(const unsigned int& iNodeIndex) const
	{
		return m_poElementGeometry->GetNode(iNodeIndex);
	}
	unsigned int FEMElement::GetNodesCount() const
	{
		return m_poElementGeometry->GetNodesCount();
	}
	double FEMElement::GetVolume() const
	{
		return m_poElementGeometry->GetVolume();
	}
	FEMElementGeometry* FEMElement::GetGeometry() const
	{
		return m_poElementGeometry;
	}
	void FEMElement::SetGeometry(FEMElementGeometry* poElementGeometry)
	{
		m_poElementGeometry = poElementGeometry;
		SetFacesLoadsSize();
	}
	FEMElement* FEMElement::CreateElementByType(FEMElementType eType)
	{
		if(eType == NullFEMElement)
		{
			return NULL;
		}
		else if(eType == PotentialFEMElement)
		{
			return new FEMPotentialElement;
		}
		else if(eType == SolidFEMElement)
		{
			return new FEMSolidElement;
		}
		else if(eType == ThermoMechanicalFEMElement)
		{
			return new FEMThermoMechanicalElement;
		}
		return NULL;
	}
	FEMElement* FEMElement::CreateElementByTypeIndex(const unsigned int& iIndex)
	{
		if(iIndex == 1)
		{
			return CreateElementByType(PotentialFEMElement);
		}
		if(iIndex == 2)
		{
			return CreateElementByType(SolidFEMElement);
		}
		if(iIndex == 3)
		{
			return CreateElementByType(ThermoMechanicalFEMElement);
		}
		return NULL;
	}
	void FEMElement::Initialize()
	{
		m_poElementGeometry = NULL;
		m_poMaterial = NULL;
		m_voGaussPoints.clear();
	}
	void FEMElement::DeleteGeometry()
	{
		if(m_poElementGeometry != NULL)
		{
			delete m_poElementGeometry;
		}
		m_poElementGeometry = NULL;
	}
	void FEMElement::SetMaterial(FEMMaterial* poMaterial)
	{
		m_poMaterial = poMaterial;
	}
	FEMMaterial* FEMElement::GetMaterial() const
	{
		return m_poMaterial;
	}
	void FEMElement::InitializeGaussPoints()
	{
 		vector< vector<double> > vvdGaussPointsCoordinates = m_poElementGeometry->GetBodyGaussPointsCoordinates();
 		vector<double> vdGaussPointsWeights = m_poElementGeometry->GetBodyGaussPointsWeights();
 		unsigned int i = 0;
 		unsigned int iGaussPointsCount = (unsigned int)vvdGaussPointsCoordinates.size();
 		m_voGaussPoints.resize(iGaussPointsCount);
 		for(i = 0; i < iGaussPointsCount ; i++)
 		{
			m_voGaussPoints[i].SetCoordinates(vvdGaussPointsCoordinates[i][0],vvdGaussPointsCoordinates[i][1],vvdGaussPointsCoordinates[i][2]);
			m_voGaussPoints[i].SetWeight(vdGaussPointsWeights[i]);
			m_voGaussPoints[i].SetElement(this);
			m_voGaussPoints[i].Update();
 		}
	}
	void FEMElement::UpdateGaussPoints()
	{
 		unsigned int i = 0;
 		unsigned int iGaussPointsCount = (unsigned int)m_voGaussPoints.size();
 		for(i = 0; i < iGaussPointsCount ; i++)
 		{
			m_voGaussPoints[i].Update();
 		}
	}
}



