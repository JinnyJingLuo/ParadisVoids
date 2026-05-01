// Ahmed M. Hussein

#include "FEMElementGeometry.h"
#include "FEMHexahedralElement.h"
#include "FEMTetrahedralElement.h"
#include "float.h"

namespace FEMSystem
{
	FEMElementGeometry::~FEMElementGeometry()
	{
		Reset();
	}
	FEMElementGeometry& FEMElementGeometry::operator=(const FEMElementGeometry& oElementGeometry)
	{
		if(oElementGeometry.GetType() != GetType())
		{
			return *this;
		}
		Set(oElementGeometry.m_vpoNodes);
		return *this;
	}
	void FEMElementGeometry::Reset()
	{
		Initialize();
	}
	void FEMElementGeometry::Set(const vector<FEMNode*>& vpoNodes)
	{
		unsigned int iNodesCount = GetNodesCount();
		if(vpoNodes.size() != iNodesCount)
		{
			return;
		}
		unsigned int i = 0;
		m_vpoNodes.resize(iNodesCount);
		for(i = 0; i < iNodesCount ; i++)
		{
			m_vpoNodes[i] = vpoNodes[i];
		}
		UpdateAxisAlignedBoundingBox();
	}
	void FEMElementGeometry::UpdateAxisAlignedBoundingBox()
	{
		unsigned int i = 0;
		double dMinX = DBL_MAX;
		double dMaxX = -DBL_MAX;
		double dMinY = DBL_MAX;
		double dMaxY = -DBL_MAX;
		double dMinZ = DBL_MAX;
		double dMaxZ = -DBL_MAX;
		double dTemp = 0.0;
		unsigned int iNodesCount = GetNodesCount();
		for(i = 0; i < iNodesCount ; i++)
		{
			dTemp = m_vpoNodes[i]->GetX();
			if(dTemp < dMinX)
			{
				dMinX = dTemp;
			}
			if(dTemp > dMaxX)
			{
				dMaxX = dTemp;
			}

			dTemp = m_vpoNodes[i]->GetY();
			if(dTemp < dMinY)
			{
				dMinY = dTemp;
			}
			if(dTemp > dMaxY)
			{
				dMaxY = dTemp;
			}

			dTemp = m_vpoNodes[i]->GetZ();
			if(dTemp < dMinZ)
			{
				dMinZ = dTemp;
			}
			if(dTemp > dMaxZ)
			{
				dMaxZ = dTemp;
			}
		}
		double dTolerance = 1E-5*min((dMaxX - dMinX),min((dMaxY - dMinY),(dMaxZ - dMinZ)));

		m_oBoundingBox.SetXMin(dMinX - dTolerance);
		m_oBoundingBox.SetXMax(dMaxX + dTolerance);
		m_oBoundingBox.SetYMin(dMinY - dTolerance);
		m_oBoundingBox.SetYMax(dMaxY + dTolerance);
		m_oBoundingBox.SetZMin(dMinZ - dTolerance);
		m_oBoundingBox.SetZMax(dMaxZ + dTolerance);
	}
	FEMNode* FEMElementGeometry::GetNode(const unsigned int& iNodeIndex) const
	{
		if(iNodeIndex >= GetNodesCount())
		{
			return NULL;
		}
		return m_vpoNodes[iNodeIndex];
	}
	vector<FEMNode*>* FEMElementGeometry::GetNodes()
	{
		return &m_vpoNodes;
	}
	bool FEMElementGeometry::IsInAxisAlignedBoundingBox(Point* poPoint) const
	{
		return m_oBoundingBox.IsPointInside(*poPoint);
	}
	void FEMElementGeometry::Initialize()
	{
		m_vpoNodes.clear();
		m_oBoundingBox.Reset();
	}
	FEMElementGeometry* FEMElementGeometry::CreateElementGeometryByType(FEMElementGeometryType eType)
	{
		if(eType == NullFEMElementGeometry)
		{
			return NULL;
		}
		if(eType == HexahedralFEMElement)
		{
			return new FEMHexahedralElement;
		}
		if(eType == TetrahedralFEMElement)
		{
			//return new FEMTetrahedralElement;
		}
		return NULL;
	}
	FEMElementGeometry* FEMElementGeometry::CreateElementGeometryByTypeIndex(const unsigned int& iIndex)
	{
		if(iIndex == 1)
		{
			return CreateElementGeometryByType(TetrahedralFEMElement);
		}
		if(iIndex == 2)
		{
			return CreateElementGeometryByType(HexahedralFEMElement);
		}
		return NULL;
	}
	AxisAlignedBoundingBox* FEMElementGeometry::GetBox()
	{
		return &m_oBoundingBox;
	}
}




