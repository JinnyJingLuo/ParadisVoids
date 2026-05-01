#include "FEMBoundaryElementFace.h"
#include "FEMSolidNode.h"
#include "MathServices.h"

using namespace EZ;

namespace FEMSystem
{
	FEMBoundaryElementFace::FEMBoundaryElementFace()
	{
		Initialize();
	}
	FEMBoundaryElementFace::FEMBoundaryElementFace(const FEMBoundaryElementFace& oFace)
	{
		*this = oFace;
	}
	FEMBoundaryElementFace::~FEMBoundaryElementFace()
	{
		Reset();
	}
	FEMBoundaryElementFace& FEMBoundaryElementFace::operator=(const FEMBoundaryElementFace& oFace)
	{
		m_vpoNodes = oFace.m_vpoNodes;
		m_voGaussPoints = oFace.m_voGaussPoints;
		m_voGaussPointsNormals = oFace.m_voGaussPointsNormals;
		m_vdGaussPointsAreas = oFace.m_vdGaussPointsAreas;
		m_oShapeFunctions = oFace.m_oShapeFunctions;
		m_voGaussPointsStresses = oFace.m_voGaussPointsStresses;
		return *this;
	}
	void FEMBoundaryElementFace::Reset()
	{
		Initialize();
	}
	void FEMBoundaryElementFace::Set(FEMElementGeometry* poGeometry,const unsigned int& iFaceIndex,const unsigned int& iPointsCount)
	{
		Reset();
		// generate the face Gauss points coordinates, normals, weights and shape function values
		double dXi = 0.0;
		double dEta = 0.0;
		vector<double> vdGaussPointsLocations;
		vector<double> vdGaussPointsWeights;
		MathServices::GenerateGaussPoints(vdGaussPointsLocations,vdGaussPointsWeights,iPointsCount);
		Matrix oNaturalCoordinates;
		m_voGaussPoints.resize(iPointsCount*iPointsCount);
		m_voGaussPointsNormals.resize(iPointsCount*iPointsCount);
		m_vdGaussPointsAreas.resize(iPointsCount*iPointsCount);
		m_voGaussPointsStresses.resize(iPointsCount*iPointsCount);
		unsigned int iNodesCount = poGeometry->GetNodesCount();
		m_oShapeFunctions.SetSize(iPointsCount*iPointsCount,iNodesCount);
		Matrix oTempShapeFunctions;
		unsigned int i = 0;
		unsigned int j = 0;
		unsigned int k = 0;
		int iIndex = 0;
		Vector oNormal;
		for(i = 0 ; i < iPointsCount ; i++)
		{
			for(j = 0 ; j < iPointsCount ; j++)
			{
				iIndex = i*iPointsCount + j;
				oNaturalCoordinates = poGeometry->GetFaceNaturalCoordinates(iFaceIndex,vdGaussPointsLocations[i],vdGaussPointsLocations[j]);
				m_voGaussPoints[iIndex] = poGeometry->GetPointOnFace(iFaceIndex,oNaturalCoordinates);
				oNormal = poGeometry->GetFaceOuterNormal(iFaceIndex,oNaturalCoordinates);
				m_vdGaussPointsAreas[iIndex] = vdGaussPointsWeights[i]*vdGaussPointsWeights[j]*oNormal.Length();
				oNormal.Normalize();
				m_voGaussPointsNormals[iIndex] = oNormal;
				oTempShapeFunctions = poGeometry->GetShapeFunctions(oNaturalCoordinates.Get(1,1),oNaturalCoordinates.Get(2,1),oNaturalCoordinates.Get(3,1));
				for(k = 1 ; k <= iNodesCount ; k++)
				{
					m_oShapeFunctions.Set(iIndex + 1,k,oTempShapeFunctions.Get(1,k));
				}
				m_voGaussPointsStresses[iIndex].SetSize(3,3);
			}
		}
		// get the element's nodes
		vector<FEMNode*>* pvpoNodes = poGeometry->GetNodes();
		m_vpoNodes.resize(iNodesCount);
		for(i = 0 ; i < iNodesCount ; i++)
		{
			m_vpoNodes[i] = pvpoNodes->at(i);
		}
	}
	void FEMBoundaryElementFace::Initialize()
	{
		m_vpoNodes.clear();
		m_voGaussPoints.clear();
		m_voGaussPointsNormals.clear();
		m_vdGaussPointsAreas.clear();
		m_oShapeFunctions.Reset();
		m_voGaussPointsStresses.clear();
	}
	void FEMBoundaryElementFace::Print() const
	{
		unsigned int iNodesCount = (unsigned int)m_vpoNodes.size();
		unsigned int i = 0;
		printf("%d nodes : ",iNodesCount);
		for(i = 0 ; i < iNodesCount ; i++)
		{
			printf("%d,",m_vpoNodes[i]->GetID());
		}
		printf("\n");
		
		unsigned int iPointsCount = (unsigned int)m_voGaussPoints.size();
		printf("number of points : %d\n",iPointsCount);
		unsigned int j = 0;
		double dShapeFunctionSum = 0.0;
		double dArea = 0.0;
		for(i = 0 ; i < iPointsCount ; i++)
		{
			dShapeFunctionSum = 0.0;
			for(j = 1 ; j <= iNodesCount ; j++)
			{
				dShapeFunctionSum = dShapeFunctionSum + m_oShapeFunctions.Get(i + 1,j);
			}
			printf("sum of shape function values : %lf\n",dShapeFunctionSum);
			dArea = dArea + m_vdGaussPointsAreas[i];
		}
		printf("area : %lf\n",dArea);
	}
	vector<Point>* FEMBoundaryElementFace::GetGaussPoints()
	{
		return &m_voGaussPoints;
	}
	void FEMBoundaryElementFace::SetGaussPointStress(const unsigned int& iIndex,const Matrix& oStress)
	{
		if(iIndex >= m_voGaussPointsStresses.size())
		{
			return;
		}
		m_voGaussPointsStresses[iIndex] = oStress;
	}
	unsigned int FEMBoundaryElementFace::GetGaussPointsCount() const
	{
		return m_voGaussPoints.size();
	}
	void FEMBoundaryElementFace::ApplyNodalForcesFromStresses()
	{
		unsigned int iGaussPointsCount = (unsigned int)m_voGaussPoints.size();
		unsigned int iNodesCount = (unsigned int)m_vpoNodes.size();
		unsigned int i = 0;
		unsigned int j = 0;
		double dTx = 0.0;
		double dTy = 0.0;
		double dTz = 0.0;
		double dNx = 0.0;
		double dNy = 0.0;
		double dNz = 0.0;
		Vector oForce;
		for(i = 0 ; i < iGaussPointsCount ; i++)
		{
			// get the normal at this Gauss point
			dNx = m_voGaussPointsNormals[i].GetX();
			dNy = m_voGaussPointsNormals[i].GetY();
			dNz = m_voGaussPointsNormals[i].GetZ();
			// get the traction at this Gauss point
			dTx = m_voGaussPointsStresses[i].Get(1,1)*dNx + m_voGaussPointsStresses[i].Get(1,2)*dNy + m_voGaussPointsStresses[i].Get(1,3)*dNz;
			dTy = m_voGaussPointsStresses[i].Get(2,1)*dNx + m_voGaussPointsStresses[i].Get(2,2)*dNy + m_voGaussPointsStresses[i].Get(2,3)*dNz;
			dTz = m_voGaussPointsStresses[i].Get(3,1)*dNx + m_voGaussPointsStresses[i].Get(3,2)*dNy + m_voGaussPointsStresses[i].Get(3,3)*dNz;
			// multiply each traction by its allocated area
			oForce.Set(dTx,dTy,dTz);
			oForce = oForce*m_vdGaussPointsAreas[i];
			// loop over the nodes and distribute this Gauss point force over them, the negative sign is
			// for surface traction cancellation
			for(j = 1 ; j <= iNodesCount ; j++)
			{
				((FEMSolidNode*)m_vpoNodes[i])->AddForce(oForce*(-m_oShapeFunctions.Get(i + 1,j)));
			}
		}
	}
}



