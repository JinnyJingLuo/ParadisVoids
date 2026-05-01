// Ahmed M. Hussein

#include "FEMHexahedralElement.h"
#include "cmath"
#include "float.h"
#include "MathServices.h"
#include "Tools.h"

using namespace SupportSystem;

namespace FEMSystem
{
	unsigned int FEMHexahedralElement::NodesPerElement = 20;
	unsigned int FEMHexahedralElement::FacesCount = 6;
	unsigned int FEMHexahedralElement::GaussPointsCount = 27;
	unsigned int FEMHexahedralElement::FaceGaussPointsCount = 9;
	vector< vector<double> > FEMHexahedralElement::m_vvdGaussPointsCoordinates = GenerateBodyGaussPointsCoordinates();
	vector<double> FEMHexahedralElement::m_vdGaussPointsWeights = GenerateBodyGaussPointsWeights();
	vector< vector< vector<double> > > FEMHexahedralElement::m_vvdFaceGaussPointsCoordinates = GenerateFaceGaussPointsCoordinates();
	vector<double> FEMHexahedralElement::m_vdFaceGaussPointsWeights = GenerateFaceGaussPointsWeights();
	Matrix FEMHexahedralElement::m_oShapeFunctionsExtrapolations = GetShapeFunctionsExtrapolations();

	FEMHexahedralElement::FEMHexahedralElement()
	{
		Initialize();
	}
	FEMHexahedralElement::FEMHexahedralElement(const FEMHexahedralElement& oElement)
	{
		*this = oElement;
	}
	FEMHexahedralElement::~FEMHexahedralElement()
	{
		Reset();
	}
	FEMHexahedralElement& FEMHexahedralElement::operator=(const FEMHexahedralElement& oElement)
	{
		FEMElementGeometry::operator=(oElement);
		return *this;
	}
	void FEMHexahedralElement::Reset()
	{
		FEMElementGeometry::Reset();
	}
	void FEMHexahedralElement::Initialize()
	{
		FEMElementGeometry::Initialize();
	}
	FEMElementGeometryType FEMHexahedralElement::GetType() const
	{
		return HexahedralFEMElement;
	}
	double FEMHexahedralElement::GetVolume() const
 	{
 		double dVolume = 0.0;
 		unsigned int i = 0;
 		double dJacobian = 0.0;
 		for(i = 0; i < GaussPointsCount ; i++)
 		{
			dJacobian = GetJacobianMatrixDeterminant(m_vvdGaussPointsCoordinates[i][0],m_vvdGaussPointsCoordinates[i][1],m_vvdGaussPointsCoordinates[i][2]);
			dVolume = dVolume + dJacobian*m_vdGaussPointsWeights[i];
 		}
 		return dVolume;
 	}
	Point FEMHexahedralElement::GetPoint(const Matrix& oNaturalCoordinates) const
	{
		if(oNaturalCoordinates.GetRowsCount() != 3 || oNaturalCoordinates.GetColumnsCount() != 1)
		{
			return Point();
		}
		unsigned int i = 0;
		Matrix oShapeFunctions = GetShapeFunctions(oNaturalCoordinates.Get(1,1),oNaturalCoordinates.Get(2,1),oNaturalCoordinates.Get(3,1));
		double dX = 0.0;
		double dY = 0.0;
		double dZ = 0.0;
		double dShapeFunctionValue = 0.0;
		for(i = 0; i < NodesPerElement ; i++)
		{
			dShapeFunctionValue = oShapeFunctions.Get(1,i + 1);
			dX = dX + dShapeFunctionValue*m_vpoNodes[i]->GetX();
			dY = dY + dShapeFunctionValue*m_vpoNodes[i]->GetY();
			dZ = dZ + dShapeFunctionValue*m_vpoNodes[i]->GetZ();
		}
		return Point(dX,dY,dZ);
	}
	vector<Point> FEMHexahedralElement::GetPoints(const Matrix& oNaturalCoordinates) const
	{
 		vector<Point> voResults;
 		unsigned int iSize = oNaturalCoordinates.GetRowsCount();
 		voResults.resize(iSize);
 		if(oNaturalCoordinates.GetColumnsCount() != 3)
 		{
 			return voResults;
 		}
 		unsigned int i = 0;
 		Matrix oNodalCoordinates(NodesPerElement,3);
 		Matrix oPoint;
 		for(i = 0; i < NodesPerElement ; i++)
 		{
 			oNodalCoordinates.Set(i + 1,1,m_vpoNodes[i]->GetX());
 			oNodalCoordinates.Set(i + 1,2,m_vpoNodes[i]->GetY());
 			oNodalCoordinates.Set(i + 1,3,m_vpoNodes[i]->GetZ());
 		}
 		Matrix oShapeFunctions;
 		for(i = 0; i < iSize ; i++)
 		{
 			oShapeFunctions = GetShapeFunctions(oNaturalCoordinates.Get(i + 1,1),oNaturalCoordinates.Get(i + 1,2),oNaturalCoordinates.Get(i + 1,3));
 			oPoint = oShapeFunctions*oNodalCoordinates;
 			voResults[i].Set(oPoint.Get(1,1),oPoint.Get(1,2),oPoint.Get(1,3));
 		}
 		return voResults;
	}
	Vector FEMHexahedralElement::GetFaceOuterNormal(const unsigned int& iFaceID,const Matrix& oNaturalCoordinates) const
	{
		Vector oNormal;
		if(oNaturalCoordinates.GetRowsCount() != 3 || oNaturalCoordinates.GetColumnsCount() != 1)
		{
			return oNormal;
		}
		
		Point oElementCenter = GetCenterPoint();
		double dXi = oNaturalCoordinates.Get(1,1);
		double dEta = oNaturalCoordinates.Get(2,1);
		double dZeta = oNaturalCoordinates.Get(3,1);
		Point oFaceCenter;
		Matrix oJacobian;
		Vector oV1;
		Vector oV2;
		
		oFaceCenter = GetFaceCenter(iFaceID);
		if(iFaceID == 1)
		{
			oJacobian = GetJacobianMatrix(-1.0,dEta,dZeta);
			oV1.Set(oJacobian.Get(2,1),oJacobian.Get(2,2),oJacobian.Get(2,3));
			oV2.Set(oJacobian.Get(3,1),oJacobian.Get(3,2),oJacobian.Get(3,3));
		}
		else if(iFaceID == 2)
		{
			oJacobian = GetJacobianMatrix(1.0,dEta,dZeta);
			oV1.Set(oJacobian.Get(2,1),oJacobian.Get(2,2),oJacobian.Get(2,3));
			oV2.Set(oJacobian.Get(3,1),oJacobian.Get(3,2),oJacobian.Get(3,3));
		}
		else if(iFaceID == 3)
		{
			oJacobian = GetJacobianMatrix(dXi,-1.0,dZeta);
			oV1.Set(oJacobian.Get(1,1),oJacobian.Get(1,2),oJacobian.Get(1,3));
			oV2.Set(oJacobian.Get(3,1),oJacobian.Get(3,2),oJacobian.Get(3,3));
		}
		else if(iFaceID == 4)
		{
			oJacobian = GetJacobianMatrix(dXi,1.0,dZeta);
			oV1.Set(oJacobian.Get(1,1),oJacobian.Get(1,2),oJacobian.Get(1,3));
			oV2.Set(oJacobian.Get(3,1),oJacobian.Get(3,2),oJacobian.Get(3,3));
		}
		else if(iFaceID == 5)
		{
			oJacobian = GetJacobianMatrix(dXi,dEta,-1.0);
			oV1.Set(oJacobian.Get(1,1),oJacobian.Get(1,2),oJacobian.Get(1,3));
			oV2.Set(oJacobian.Get(2,1),oJacobian.Get(2,2),oJacobian.Get(2,3));
		}
		else if(iFaceID == 6)
		{
			oJacobian = GetJacobianMatrix(dXi,dEta,1.0);
			oV1.Set(oJacobian.Get(1,1),oJacobian.Get(1,2),oJacobian.Get(1,3));
			oV2.Set(oJacobian.Get(2,1),oJacobian.Get(2,2),oJacobian.Get(2,3));
		}
		Vector oTest(oElementCenter,oFaceCenter);
		oNormal = oV1^oV2;
		if(oNormal*oTest < 0.0)
		{
			oNormal = oNormal*(-1.0);
		}
		return oNormal;
	}
	unsigned int FEMHexahedralElement::GetNodesCount() const
	{
		return NodesPerElement;
	}
	Matrix FEMHexahedralElement::GetShapeFunctions(const double& dXi,const double& dEta,const double& dZeta) const
	{
		return EvaluateShapeFunctions(dXi,dEta,dZeta);
	}
	Matrix FEMHexahedralElement::GetShapeFunctionsDerivatives(const double& dXi,const double& dEta,const double& dZeta) const
	{
		Matrix oResult(3,NodesPerElement);
		
		vector<double> vdShapeFunctions;
		vdShapeFunctions.resize(NodesPerElement);
		unsigned int i = 0;
		double dFactor1 = 1.0/8.0;
		double dFactor2 = 1.0/4.0;

		// Xi derivatives
		vdShapeFunctions[0] = -dFactor1*(1-dXi)*(1-dEta)*(1-dZeta)-dFactor1*(-dXi-dEta-dZeta-2)*(1-dEta)*(1-dZeta);
		vdShapeFunctions[1] = dFactor1*(1+dXi)*(1-dEta)*(1-dZeta)+dFactor1*(dXi-dEta-dZeta-2)*(1-dEta)*(1-dZeta);
		vdShapeFunctions[2] = dFactor1*(1+dXi)*(1+dEta)*(1-dZeta)+dFactor1*(dXi+dEta-dZeta-2)*(1+dEta)*(1-dZeta);
		vdShapeFunctions[3] = -dFactor1*(1-dXi)*(1+dEta)*(1-dZeta)-dFactor1*(-dXi+dEta-dZeta-2)*(1+dEta)*(1-dZeta);
		vdShapeFunctions[4] = -dFactor1*(1-dXi)*(1-dEta)*(1+dZeta)-dFactor1*(-dXi-dEta+dZeta-2)*(1-dEta)*(1+dZeta);
		vdShapeFunctions[5] = dFactor1*(1+dXi)*(1-dEta)*(1+dZeta)+dFactor1*(dXi-dEta+dZeta-2)*(1-dEta)*(1+dZeta);
		vdShapeFunctions[6] = dFactor1*(1+dXi)*(1+dEta)*(1+dZeta)+dFactor1*(dXi+dEta+dZeta-2)*(1+dEta)*(1+dZeta);
		vdShapeFunctions[7] = -dFactor1*(1-dXi)*(1+dEta)*(1+dZeta)-dFactor1*(-dXi+dEta+dZeta-2)*(1+dEta)*(1+dZeta);
		vdShapeFunctions[8] = -2.0*dFactor2*dXi*(1-dEta)*(1-dZeta);
		vdShapeFunctions[9] = dFactor2*(1-dEta*dEta)*(1-dZeta);
		vdShapeFunctions[10] = -2.0*dFactor2*dXi*(1+dEta)*(1-dZeta);
		vdShapeFunctions[11] = -dFactor2*(1-dEta*dEta)*(1-dZeta);
		vdShapeFunctions[12] = -2.0*dFactor2*dXi*(1-dEta)*(1+dZeta);
		vdShapeFunctions[13] = dFactor2*(1-dEta*dEta)*(1+dZeta);
		vdShapeFunctions[14] = -2.0*dFactor2*dXi*(1+dEta)*(1+dZeta);
		vdShapeFunctions[15] = -dFactor2*(1-dEta*dEta)*(1+dZeta);
		vdShapeFunctions[16] = -dFactor2*(1-dZeta*dZeta)*(1-dEta);
		vdShapeFunctions[17] = dFactor2*(1-dZeta*dZeta)*(1-dEta);
		vdShapeFunctions[18] = dFactor2*(1-dZeta*dZeta)*(1+dEta);
		vdShapeFunctions[19] = -dFactor2*(1-dZeta*dZeta)*(1+dEta);

		for(i = 0; i < NodesPerElement ; i++)
		{
			oResult.Set(1,i + 1,vdShapeFunctions[i]);
		}
		// Eta derivatives
		vdShapeFunctions[0] = -dFactor1*(1-dXi)*(1-dEta)*(1-dZeta)-dFactor1*(-dXi-dEta-dZeta-2)*(1-dXi)*(1-dZeta);
		vdShapeFunctions[1] = -dFactor1*(1+dXi)*(1-dEta)*(1-dZeta)-dFactor1*(dXi-dEta-dZeta-2)*(1+dXi)*(1-dZeta);
		vdShapeFunctions[2] = dFactor1*(1+dXi)*(1+dEta)*(1-dZeta)+dFactor1*(dXi+dEta-dZeta-2)*(1+dXi)*(1-dZeta);
		vdShapeFunctions[3] = dFactor1*(1-dXi)*(1+dEta)*(1-dZeta)+dFactor1*(-dXi+dEta-dZeta-2)*(1-dXi)*(1-dZeta);
		vdShapeFunctions[4] = -dFactor1*(1-dXi)*(1-dEta)*(1+dZeta)-dFactor1*(-dXi-dEta+dZeta-2)*(1-dXi)*(1+dZeta);
		vdShapeFunctions[5] = -dFactor1*(1+dXi)*(1-dEta)*(1+dZeta)-dFactor1*(dXi-dEta+dZeta-2)*(1+dXi)*(1+dZeta);
		vdShapeFunctions[6] = dFactor1*(1+dXi)*(1+dEta)*(1+dZeta)+dFactor1*(dXi+dEta+dZeta-2)*(1+dXi)*(1+dZeta);
		vdShapeFunctions[7] = dFactor1*(1-dXi)*(1+dEta)*(1+dZeta)+dFactor1*(-dXi+dEta+dZeta-2)*(1-dXi)*(1+dZeta);
		vdShapeFunctions[8] = -dFactor2*(1-dXi*dXi)*(1-dZeta);
		vdShapeFunctions[9] = -2.0*dFactor2*dEta*(1+dXi)*(1-dZeta);
		vdShapeFunctions[10] =dFactor2*(1-dXi*dXi)*(1-dZeta);
		vdShapeFunctions[11] = -2.0*dFactor2*dEta*(1-dXi)*(1-dZeta);
		vdShapeFunctions[12] = -dFactor2*(1-dXi*dXi)*(1+dZeta);
		vdShapeFunctions[13] = -2.0*dFactor2*dEta*(1+dXi)*(1+dZeta);
		vdShapeFunctions[14] = dFactor2*(1-dXi*dXi)*(1+dZeta);
		vdShapeFunctions[15] = -2.0*dFactor2*dEta*(1-dXi)*(1+dZeta);
		vdShapeFunctions[16] = -dFactor2*(1-dZeta*dZeta)*(1-dXi);
		vdShapeFunctions[17] = -dFactor2*(1-dZeta*dZeta)*(1+dXi);
		vdShapeFunctions[18] = dFactor2*(1-dZeta*dZeta)*(1+dXi);
		vdShapeFunctions[19] = dFactor2*(1-dZeta*dZeta)*(1-dXi);

		for(i = 0; i < NodesPerElement ; i++)
		{
			oResult.Set(2,i + 1,vdShapeFunctions[i]);
		}
		// Zeta derivatives
		vdShapeFunctions[0] = -dFactor1*(1-dXi)*(1-dEta)*(1-dZeta)-dFactor1*(-dXi-dEta-dZeta-2)*(1-dXi)*(1-dEta);
		vdShapeFunctions[1] = -dFactor1*(1+dXi)*(1-dEta)*(1-dZeta)-dFactor1*(dXi-dEta-dZeta-2)*(1+dXi)*(1-dEta);
		vdShapeFunctions[2] = -dFactor1*(1+dXi)*(1+dEta)*(1-dZeta)-dFactor1*(dXi+dEta-dZeta-2)*(1+dXi)*(1+dEta);
		vdShapeFunctions[3] = -dFactor1*(1-dXi)*(1+dEta)*(1-dZeta)-dFactor1*(-dXi+dEta-dZeta-2)*(1-dXi)*(1+dEta);
		vdShapeFunctions[4] = dFactor1*(1-dXi)*(1-dEta)*(1+dZeta)+dFactor1*(-dXi-dEta+dZeta-2)*(1-dXi)*(1-dEta);
		vdShapeFunctions[5] = dFactor1*(1+dXi)*(1-dEta)*(1+dZeta)+dFactor1*(dXi-dEta+dZeta-2)*(1+dXi)*(1-dEta);
		vdShapeFunctions[6] = dFactor1*(1+dXi)*(1+dEta)*(1+dZeta)+dFactor1*(dXi+dEta+dZeta-2)*(1+dXi)*(1+dEta);
		vdShapeFunctions[7] = dFactor1*(1-dXi)*(1+dEta)*(1+dZeta)+dFactor1*(-dXi+dEta+dZeta-2)*(1-dXi)*(1+dEta);
		vdShapeFunctions[8] = -dFactor2*(1-dXi*dXi)*(1-dEta);
		vdShapeFunctions[9] = -dFactor2*(1-dEta*dEta)*(1+dXi);
		vdShapeFunctions[10] = -dFactor2*(1-dXi*dXi)*(1+dEta);
		vdShapeFunctions[11] = -dFactor2*(1-dEta*dEta)*(1-dXi);
		vdShapeFunctions[12] = dFactor2*(1-dXi*dXi)*(1-dEta);
		vdShapeFunctions[13] = dFactor2*(1-dEta*dEta)*(1+dXi);
		vdShapeFunctions[14] = dFactor2*(1-dXi*dXi)*(1+dEta);
		vdShapeFunctions[15] = dFactor2*(1-dEta*dEta)*(1-dXi);
		vdShapeFunctions[16] = -2.0*dFactor2*dZeta*(1-dXi)*(1-dEta);
		vdShapeFunctions[17] = -2.0*dFactor2*dZeta*(1+dXi)*(1-dEta);
		vdShapeFunctions[18] = -2.0*dFactor2*dZeta*(1+dXi)*(1+dEta);
		vdShapeFunctions[19] = -2.0*dFactor2*dZeta*(1-dXi)*(1+dEta);

		for(i = 0; i < NodesPerElement ; i++)
		{
			oResult.Set(3,i + 1,vdShapeFunctions[i]);
		}

		return oResult;
	}
	Matrix FEMHexahedralElement::GetJacobianMatrix(const double& dXi,const double& dEta,const double& dZeta) const
	{
		Matrix oJacobian;
		unsigned int i = 0;
		Matrix oNodalCoordinates(NodesPerElement,3);
		for(i = 0; i < NodesPerElement ; i++)
		{
			oNodalCoordinates.Set(i + 1,1,m_vpoNodes[i]->GetX());
			oNodalCoordinates.Set(i + 1,2,m_vpoNodes[i]->GetY());
			oNodalCoordinates.Set(i + 1,3,m_vpoNodes[i]->GetZ());
		}
		Matrix oDerivatives = GetShapeFunctionsDerivatives(dXi,dEta,dZeta);
		oJacobian = oDerivatives*oNodalCoordinates;
		return oJacobian;
	}
	double FEMHexahedralElement::GetJacobianMatrixDeterminant(const double& dXi,const double& dEta,const double& dZeta) const
	{
		Matrix oJacobian = GetJacobianMatrix(dXi,dEta,dZeta);
		double dDeterminant = 0.0;
		double dDet1 = oJacobian.Get(2,2)*oJacobian.Get(3,3) - oJacobian.Get(2,3)*oJacobian.Get(3,2);
		double dDet2 = oJacobian.Get(2,1)*oJacobian.Get(3,3) - oJacobian.Get(2,3)*oJacobian.Get(3,1);
		double dDet3 = oJacobian.Get(2,1)*oJacobian.Get(3,2) - oJacobian.Get(2,2)*oJacobian.Get(3,1);
		dDeterminant = oJacobian.Get(1,1)*dDet1 - oJacobian.Get(1,2)*dDet2 + oJacobian.Get(1,3)*dDet3;
		return dDeterminant;
	}
	Matrix FEMHexahedralElement::GetShapeFunctionsXYZDerivatives(const double& dXi,const double& dEta,const double& dZeta,double& dJacobian) const
	{
		Matrix oXYZDerivatives;
		Matrix oDerivatives = GetShapeFunctionsDerivatives(dXi,dEta,dZeta);
		Matrix oJacobian;
		unsigned int i = 0;
		Matrix oNodalCoordinates(NodesPerElement,3);
		for(i = 0; i < NodesPerElement ; i++)
		{
			oNodalCoordinates.Set(i + 1,1,m_vpoNodes[i]->GetX());
			oNodalCoordinates.Set(i + 1,2,m_vpoNodes[i]->GetY());
			oNodalCoordinates.Set(i + 1,3,m_vpoNodes[i]->GetZ());
		}
		oJacobian = oDerivatives*oNodalCoordinates;

		Matrix oInverseJacobian = Matrix::Invert3x3Matrix(oJacobian,dJacobian);

		oXYZDerivatives = oInverseJacobian*oDerivatives;
		return oXYZDerivatives;
	}	
	vector< vector<double> > FEMHexahedralElement::GetBodyGaussPointsCoordinates() const
	{
		return m_vvdGaussPointsCoordinates;
	}
	vector<double> FEMHexahedralElement::GetBodyGaussPointsWeights() const
	{
		return m_vdGaussPointsWeights;
	}
	vector< vector<double> > FEMHexahedralElement::GetFaceGaussPointsCoordinates(const unsigned int& iFaceIndex) const
	{
		vector< vector<double> > vvdCoordinates;
		vvdCoordinates.clear();
		if((iFaceIndex <= 0) || (iFaceIndex > FacesCount))
		{
			return vvdCoordinates;
		}
		vvdCoordinates = m_vvdFaceGaussPointsCoordinates[iFaceIndex - 1];
		return vvdCoordinates;
	}
	vector<double> FEMHexahedralElement::GetFaceGaussPointsWeights() const
	{
		return m_vdFaceGaussPointsWeights;
	}
	vector<unsigned int> FEMHexahedralElement::GetSurfaceFacesIndices() const
	{
		vector<unsigned int> viFacesIndices;
		viFacesIndices.reserve(FacesCount);
		unsigned int i = 0;
		for(i = 1 ; i <= FacesCount ; i++)
		{
			if(IsFaceOnSurface(i))
			{
				viFacesIndices.push_back(i);
			}
		}
		return viFacesIndices;
	}
	double FEMHexahedralElement::GetDistanceToElementCenterPoint(Point* poPoint) const
	{
		return poPoint->Distance(GetCenterPoint());
	}
	Matrix FEMHexahedralElement::GetShapeFunctionsExtrapolationsMatrix() const
	{
		return m_oShapeFunctionsExtrapolations;
	}
	Matrix FEMHexahedralElement::EvaluateShapeFunctions(const double& dXi,const double& dEta,const double& dZeta)
	{
		Matrix oResult(1,NodesPerElement);
		
		vector<double> vdShapeFunctions;
		vdShapeFunctions.resize(NodesPerElement);

		double dFactor1 = 1.0/8.0;
		double dFactor2 = 1.0/4.0;

		vdShapeFunctions[0] = dFactor1*(-dXi - dEta - dZeta - 2.0)*(1.0 - dXi)*(1.0 - dEta)*(1.0 - dZeta);
		vdShapeFunctions[1] = dFactor1*(dXi - dEta - dZeta - 2.0)*(1.0 + dXi)*(1.0 - dEta)*(1.0 - dZeta);
		vdShapeFunctions[2] = dFactor1*(dXi + dEta - dZeta - 2.0)*(1.0 + dXi)*(1.0 + dEta)*(1.0 - dZeta);
		vdShapeFunctions[3] = dFactor1*(-dXi + dEta - dZeta - 2.0)*(1.0 - dXi)*(1.0 + dEta)*(1.0 - dZeta);
		vdShapeFunctions[4] = dFactor1*(-dXi - dEta + dZeta - 2.0)*(1.0 - dXi)*(1.0 - dEta)*(1.0 + dZeta);
		vdShapeFunctions[5] = dFactor1*(dXi - dEta + dZeta - 2.0)*(1.0 + dXi)*(1.0 - dEta)*(1.0 + dZeta);
		vdShapeFunctions[6] = dFactor1*(dXi + dEta + dZeta - 2.0)*(1.0 + dXi)*(1.0 + dEta)*(1.0 + dZeta);
		vdShapeFunctions[7] = dFactor1*(-dXi + dEta + dZeta - 2.0)*(1.0 - dXi)*(1.0 + dEta)*(1.0 + dZeta);
		vdShapeFunctions[8] = dFactor2*(1.0 - dXi*dXi)*(1.0 - dEta)*(1.0 - dZeta);
		vdShapeFunctions[9] = dFactor2*(1.0 - dEta*dEta)*(1.0 + dXi)*(1.0 - dZeta);
		vdShapeFunctions[10] = dFactor2*(1.0 - dXi*dXi)*(1.0 + dEta)*(1.0 - dZeta);
		vdShapeFunctions[11] = dFactor2*(1.0 - dEta*dEta)*(1.0 - dXi)*(1.0 - dZeta);
		vdShapeFunctions[12] = dFactor2*(1.0 - dXi*dXi)*(1.0 - dEta)*(1.0 + dZeta);
		vdShapeFunctions[13] = dFactor2*(1.0 - dEta*dEta)*(1.0 + dXi)*(1.0 + dZeta);
		vdShapeFunctions[14] = dFactor2*(1.0 - dXi*dXi)*(1.0 + dEta)*(1.0 + dZeta);
		vdShapeFunctions[15] = dFactor2*(1.0 - dEta*dEta)*(1.0 - dXi)*(1.0 + dZeta);
		vdShapeFunctions[16] = dFactor2*(1.0 - dZeta*dZeta)*(1.0 - dXi)*(1.0 - dEta);
		vdShapeFunctions[17] = dFactor2*(1.0 - dZeta*dZeta)*(1.0 + dXi)*(1.0 - dEta);
		vdShapeFunctions[18] = dFactor2*(1.0 - dZeta*dZeta)*(1.0 + dXi)*(1.0 + dEta);
		vdShapeFunctions[19] = dFactor2*(1.0 - dZeta*dZeta)*(1.0 - dXi)*(1.0 + dEta);

		unsigned int i = 0;
		for(i = 0; i < NodesPerElement ; i++)
		{
			oResult.Set(1,i + 1,vdShapeFunctions[i]);
		}

		return oResult;
	}
	Matrix FEMHexahedralElement::GetShapeFunctionsExtrapolations()
	{
		// this function is tailored for the 20 nodes hex element with 27 Gauss points
		// based on the least squares solution for the mapping equations
		Matrix N(GaussPointsCount,NodesPerElement);
 		unsigned int i = 0;
 		unsigned int j = 0;
 		Matrix oShapeFunctions;
 		for(i = 0; i < GaussPointsCount ; i++)
 		{
 			oShapeFunctions = EvaluateShapeFunctions(m_vvdGaussPointsCoordinates[i][0],m_vvdGaussPointsCoordinates[i][1],m_vvdGaussPointsCoordinates[i][2]);
 			for(j = 1 ; j <= NodesPerElement ; j++)
 			{
 				N.Set(i + 1,j,oShapeFunctions.Get(1,j));
 			}
 		}
		Matrix A = N.GetTranspose()*N;
		Matrix B = A.GetInverse()*N.GetTranspose();
 		return B;
	}
	Matrix FEMHexahedralElement::GetShapeFunctionsExtrapolationsOld()
	{
		Matrix oExtrapolations(NodesPerElement,27);
		unsigned int i = 0;
		unsigned int j = 0;
		Matrix oTemp;
		double daCoordinates[60] = {-1.0,-1.0,-1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,-1.0,-1.0,-1.0,1.0,1.0,-1.0,1.0,1.0,1.0,1.0,-1.0,1.0,1.0,0.0,-1.0,-1.0,1.0,0.0,-1.0,0.0,1.0,-1.0,-1.0,0.0,-1.0,0.0,-1.0,1.0,1.0,0.0,1.0,0.0,1.0,1.0,-1.0,0.0,1.0,-1.0,-1.0,0.0,1.0,-1.0,0.0,1.0,1.0,0.0,-1.0,1.0,0.0};

		for(i = 0; i < NodesPerElement ; i++)
		{
			oTemp = EvaluateExtrapolatedShapeFunctions(daCoordinates[3*i],daCoordinates[3*i + 1],daCoordinates[3*i + 2]);
			for(j = 1; j <= 27 ; j++)
			{
				oExtrapolations.Set(i + 1,j,oTemp.Get(j,1));
			}
		}
		return oExtrapolations;
	}
	Matrix FEMHexahedralElement::EvaluateExtrapolatedShapeFunctions(const double& dXi,const double& dEta,const double& dZeta)
	{
		Matrix oN(27,1);
		double dK = sqrt(0.6);
		double dFactor1 = 1.0/8.0;
		double dFactor2 = 1.0/4.0;
		double dFactor3 = 1.0/2.0;
		

		oN.Set(1,1,-dFactor1*(1.0 - dXi/dK)*(1.0 - dEta/dK)*(1.0 - dZeta/dK)*dXi/dK*dEta/dK*dZeta/dK);
		oN.Set(2,1,dFactor2*(1.0 - dXi/dK)*(1.0 - dEta/dK)*(1.0 - dZeta/dK*dZeta/dK)*dXi/dK*dEta/dK);
		oN.Set(3,1,dFactor1*(1.0 - dXi/dK)*(1.0 - dEta/dK)*(1.0 + dZeta/dK)*dXi/dK*dEta/dK*dZeta/dK);
		oN.Set(4,1,dFactor2*(1.0 - dXi/dK)*(1.0 - dEta/dK*dEta/dK)*(1.0 - dZeta/dK)*dXi/dK*dZeta/dK);
		oN.Set(5,1,-dFactor3*(1.0 - dXi/dK)*(1.0 - dEta/dK*dEta/dK)*(1.0 - dZeta/dK*dZeta/dK)*dXi/dK);
		oN.Set(6,1,-dFactor2*(1.0 - dXi/dK)*(1.0 - dEta/dK*dEta/dK)*(1.0 + dZeta/dK)*dXi/dK*dZeta/dK);
		oN.Set(7,1,dFactor1*(1.0 - dXi/dK)*(1.0 + dEta/dK)*(1.0 - dZeta/dK)*dXi/dK*dEta/dK*dZeta/dK);
		oN.Set(8,1,-dFactor2*(1.0 - dXi/dK)*(1.0 + dEta/dK)*(1.0 - dZeta/dK*dZeta/dK)*dXi/dK*dEta/dK);
		oN.Set(9,1,-dFactor1*(1.0 - dXi/dK)*(1.0 + dEta/dK)*(1.0 + dZeta/dK)*dXi/dK*dEta/dK*dZeta/dK);
		oN.Set(10,1,dFactor2*(1.0 - dXi/dK*dXi/dK)*(1.0 - dEta/dK)*(1.0 - dZeta/dK)*dEta/dK*dZeta/dK);
		oN.Set(11,1,-dFactor3*(1.0 - dXi/dK*dXi/dK)*(1.0 - dEta/dK)*(1.0 - dZeta/dK*dZeta/dK)*dEta/dK);
		oN.Set(12,1,-dFactor2*(1.0 - dXi/dK*dXi/dK)*(1.0 - dEta/dK)*(1.0 + dZeta/dK)*dEta/dK*dZeta/dK);
		oN.Set(13,1,-dFactor3*(1.0 - dXi/dK*dXi/dK)*(1.0 - dEta/dK*dEta/dK)*(1.0 - dZeta/dK)*dZeta/dK);
		oN.Set(14,1,(1.0 - dXi/dK*dXi/dK)*(1.0 - dEta/dK*dEta/dK)*(1.0 - dZeta/dK*dZeta/dK));
		oN.Set(15,1,dFactor3*(1.0 - dXi/dK*dXi/dK)*(1.0 - dEta/dK*dEta/dK)*(1.0 + dZeta/dK)*dZeta/dK);
		oN.Set(16,1,-dFactor2*(1.0 - dXi/dK*dXi/dK)*(1.0 + dEta/dK)*(1.0 - dZeta/dK)*dEta/dK*dZeta/dK);
		oN.Set(17,1,dFactor3*(1.0 - dXi/dK*dXi/dK)*(1.0 + dEta/dK)*(1.0 - dZeta/dK*dZeta/dK)*dEta/dK);
		oN.Set(18,1,dFactor2*(1.0 - dXi/dK*dXi/dK)*(1.0 + dEta/dK)*(1.0 + dZeta/dK)*dEta/dK*dZeta/dK);
		oN.Set(19,1,dFactor1*(1.0 + dXi/dK)*(1.0 - dEta/dK)*(1.0 - dZeta/dK)*dXi/dK*dEta/dK*dZeta/dK);
		oN.Set(20,1,-dFactor2*(1.0 + dXi/dK)*(1.0 - dEta/dK)*(1.0 - dZeta/dK*dZeta/dK)*dXi/dK*dEta/dK);
		oN.Set(21,1,-dFactor1*(1.0 + dXi/dK)*(1.0 - dEta/dK)*(1.0 + dZeta/dK)*dXi/dK*dEta/dK*dZeta/dK);
		oN.Set(22,1,-dFactor2*(1.0 + dXi/dK)*(1.0 - dEta/dK*dEta/dK)*(1.0 - dZeta/dK)*dXi/dK*dZeta/dK);
		oN.Set(23,1,dFactor3*(1.0 + dXi/dK)*(1.0 - dEta/dK*dEta/dK)*(1.0 - dZeta/dK*dZeta/dK)*dXi/dK);
		oN.Set(24,1,dFactor2*(1.0 + dXi/dK)*(1.0 - dEta/dK*dEta/dK)*(1.0 + dZeta/dK)*dXi/dK*dZeta/dK);
		oN.Set(25,1,-dFactor1*(1.0 + dXi/dK)*(1.0 + dEta/dK)*(1.0 - dZeta/dK)*dXi/dK*dEta/dK*dZeta/dK);
		oN.Set(26,1,dFactor2*(1.0 + dXi/dK)*(1.0 + dEta/dK)*(1.0 - dZeta/dK*dZeta/dK)*dXi/dK*dEta/dK);
		oN.Set(27,1,dFactor1*(1.0 + dXi/dK)*(1.0 + dEta/dK)*(1.0 + dZeta/dK)*dXi/dK*dEta/dK*dZeta/dK);

		return oN;
	}
	vector< vector<double> > FEMHexahedralElement::GenerateBodyGaussPointsCoordinates()
	{
		vector< vector<double> > vvdCoordinates;
		vvdCoordinates.resize(GaussPointsCount);
		unsigned int i = 0;
		for(i = 0 ; i < GaussPointsCount ; i++)
		{
			vvdCoordinates[i].resize(3);
		}
		unsigned int j = 0;
		unsigned int k = 0;
		double dTemp = sqrt(3.0/5.0);
		double daCoordinates[3] = {-dTemp,0.0,dTemp};
		unsigned int iIndex = 0;
		for(i = 0 ; i < 3 ; i++)
		{
			for(j = 0 ; j < 3 ; j++)
			{
				for(k = 0 ; k < 3 ; k++)
				{
					vvdCoordinates[iIndex][0] = daCoordinates[i];
					vvdCoordinates[iIndex][1] = daCoordinates[j];
					vvdCoordinates[iIndex][2] = daCoordinates[k];
					iIndex = iIndex + 1;
				}
			}
		}
		return vvdCoordinates;
	}
	vector<double> FEMHexahedralElement::GenerateBodyGaussPointsWeights()
	{
		vector<double> vdWeights;
		vdWeights.resize(GaussPointsCount);
		unsigned int i = 0;
		unsigned int j = 0;
		unsigned int k = 0;
		double dTemp = 5.0/9.0;
		double daWeights[3] = {dTemp,8.0/9.0,dTemp};
		unsigned int iIndex = 0;
		for(i = 0 ; i < 3 ; i++)
		{
			for(j = 0 ; j < 3 ; j++)
			{
				for(k = 0 ; k < 3 ; k++)
				{
					vdWeights[iIndex] = daWeights[i]*daWeights[j]*daWeights[k];
					iIndex = iIndex + 1;
				}
			}
		}
		return vdWeights;
	}
	vector< vector< vector<double> > > FEMHexahedralElement::GenerateFaceGaussPointsCoordinates()
	{
		vector< vector< vector<double> > > vvvdCoordinates;
		vvvdCoordinates.resize(FacesCount);
		unsigned int i = 0;
		unsigned int j = 0;
		for(i = 0 ; i < FacesCount ; i++)
		{
			vvvdCoordinates[i].resize(FaceGaussPointsCount);
			for(j = 0 ; j < FaceGaussPointsCount ; j++)
			{
				vvvdCoordinates[i][j].resize(3);
			}
		}
		
		unsigned int k = 0;
		double dTemp = sqrt(3.0/5.0);
		double daCoordinates[3] = {-dTemp,0.0,dTemp};
		double dFixedCoordinate[6] = {-1.0,1.0,-1.0,1.0,-1.0,1.0};
		unsigned int iaFixedIndex[6] = {0,0,1,1,2,2};
		unsigned int iaJIndex[6] = {1,1,0,0,0,0};
		unsigned int iaKIndex[6] = {2,2,2,2,1,1};
		unsigned int iIndex = 0;
		for(i = 0 ; i < FacesCount ; i++)		// loop over the faces
		{
			iIndex = 0;
			for(j = 0 ; j < 3 ; j++)			// loop over the face Gauss points rows
			{
				for(k = 0 ; k < 3 ; k++)		// loop over the face Gauss points columns
				{
					vvvdCoordinates[i][iIndex][iaJIndex[i]] = daCoordinates[j];
					vvvdCoordinates[i][iIndex][iaKIndex[i]] = daCoordinates[k];
					vvvdCoordinates[i][iIndex][iaFixedIndex[i]] = dFixedCoordinate[i];
					iIndex = iIndex + 1;
				}
			}
		}
		return vvvdCoordinates;
	}
	vector<double> FEMHexahedralElement::GenerateFaceGaussPointsWeights()
	{
		vector<double> vdWeights;
		vdWeights.resize(FaceGaussPointsCount);
		unsigned int i = 0;		
		unsigned int j = 0;
		double dTemp = 5.0/9.0;
		double daWeights[3] = {dTemp,8.0/9.0,dTemp};
		unsigned int iIndex = 0;
		for(i = 0 ; i < 3 ; i++)			// loop over the face Gauss points rows
		{
			for(j = 0 ; j < 3 ; j++)		// loop over the face Gauss points columns
			{
				vdWeights[iIndex] = daWeights[i]*daWeights[j];
				iIndex = iIndex + 1;
			}
		}
		return vdWeights;
	}
	void FEMHexahedralElement::GenerateSurfacePatches(vector< vector<unsigned int> >& vviNodesIndices,vector<GenericNode*>& vpoMidPoints) const
	{
		vector<unsigned int> viFacesIndices = GetSurfaceFacesIndices();
		unsigned int iSize = viFacesIndices.size();
		vector<unsigned int> viNodesIndices;
		vviNodesIndices.resize(iSize);
		vpoMidPoints.resize(iSize);
		unsigned int i = 0;
		for(i = 0; i < iSize ; i++)
		{
			vpoMidPoints[i] = GetFaceQuadPatch(viFacesIndices[i],vviNodesIndices[i]);
		}
	}
	GenericNode* FEMHexahedralElement::GetFaceQuadPatch(const unsigned int& iFaceIndex,vector<unsigned int>& viOriginalNodesIndices) const
	{
		viOriginalNodesIndices.resize(8);
        Point oTempPoint = GetFaceCenter(iFaceIndex);
		if(iFaceIndex == 1)
		{
			viOriginalNodesIndices[0] = m_vpoNodes[3]->GetID();
			viOriginalNodesIndices[1] = m_vpoNodes[11]->GetID();
			viOriginalNodesIndices[2] = m_vpoNodes[0]->GetID();
			viOriginalNodesIndices[3] = m_vpoNodes[16]->GetID();
			viOriginalNodesIndices[4] = m_vpoNodes[4]->GetID();
			viOriginalNodesIndices[5] = m_vpoNodes[15]->GetID();
			viOriginalNodesIndices[6] = m_vpoNodes[7]->GetID();
			viOriginalNodesIndices[7] = m_vpoNodes[19]->GetID();
		}
		else if(iFaceIndex == 2)
		{
			viOriginalNodesIndices[0] = m_vpoNodes[1]->GetID();
			viOriginalNodesIndices[1] = m_vpoNodes[9]->GetID();
			viOriginalNodesIndices[2] = m_vpoNodes[2]->GetID();
			viOriginalNodesIndices[3] = m_vpoNodes[18]->GetID();
			viOriginalNodesIndices[4] = m_vpoNodes[6]->GetID();
			viOriginalNodesIndices[5] = m_vpoNodes[13]->GetID();
			viOriginalNodesIndices[6] = m_vpoNodes[5]->GetID();
			viOriginalNodesIndices[7] = m_vpoNodes[17]->GetID();
		}
		else if(iFaceIndex == 3)
		{
			viOriginalNodesIndices[0] = m_vpoNodes[0]->GetID();
			viOriginalNodesIndices[1] = m_vpoNodes[8]->GetID();
			viOriginalNodesIndices[2] = m_vpoNodes[1]->GetID();
			viOriginalNodesIndices[3] = m_vpoNodes[17]->GetID();
			viOriginalNodesIndices[4] = m_vpoNodes[5]->GetID();
			viOriginalNodesIndices[5] = m_vpoNodes[12]->GetID();
			viOriginalNodesIndices[6] = m_vpoNodes[4]->GetID();
			viOriginalNodesIndices[7] = m_vpoNodes[16]->GetID();
		}
		else if(iFaceIndex == 4)
		{
			viOriginalNodesIndices[0] = m_vpoNodes[2]->GetID();
			viOriginalNodesIndices[1] = m_vpoNodes[10]->GetID();
			viOriginalNodesIndices[2] = m_vpoNodes[3]->GetID();
			viOriginalNodesIndices[3] = m_vpoNodes[19]->GetID();
			viOriginalNodesIndices[4] = m_vpoNodes[7]->GetID();
			viOriginalNodesIndices[5] = m_vpoNodes[14]->GetID();
			viOriginalNodesIndices[6] = m_vpoNodes[6]->GetID();
			viOriginalNodesIndices[7] = m_vpoNodes[18]->GetID();
		}
		else if(iFaceIndex == 5)
		{
			viOriginalNodesIndices[0] = m_vpoNodes[0]->GetID();
			viOriginalNodesIndices[1] = m_vpoNodes[11]->GetID();
			viOriginalNodesIndices[2] = m_vpoNodes[3]->GetID();
			viOriginalNodesIndices[3] = m_vpoNodes[10]->GetID();
			viOriginalNodesIndices[4] = m_vpoNodes[2]->GetID();
			viOriginalNodesIndices[5] = m_vpoNodes[9]->GetID();
			viOriginalNodesIndices[6] = m_vpoNodes[1]->GetID();
			viOriginalNodesIndices[7] = m_vpoNodes[8]->GetID();
		}
		else if(iFaceIndex == 6)
		{
			viOriginalNodesIndices[0] = m_vpoNodes[5]->GetID();
			viOriginalNodesIndices[1] = m_vpoNodes[13]->GetID();
			viOriginalNodesIndices[2] = m_vpoNodes[6]->GetID();
			viOriginalNodesIndices[3] = m_vpoNodes[14]->GetID();
			viOriginalNodesIndices[4] = m_vpoNodes[7]->GetID();
			viOriginalNodesIndices[5] = m_vpoNodes[15]->GetID();
			viOriginalNodesIndices[6] = m_vpoNodes[4]->GetID();
			viOriginalNodesIndices[7] = m_vpoNodes[12]->GetID();
		}
		else
		{
			viOriginalNodesIndices.clear();
			return NULL;
		}
   	 	GenericNode* poMidPoint = new GenericNode;
        poMidPoint->Set(oTempPoint.GetX(),oTempPoint.GetY(),oTempPoint.GetZ());
		return poMidPoint;
	}
	bool FEMHexahedralElement::IsFaceOnSurface(const unsigned int& iFaceID) const
	{
		unsigned int viNodeIndices[48] = {1,12,4,20,8,16,5,17,2,10,3,19,7,14,6,18,1,9,2,18,6,13,5,17,4,11,3,19,7,15,8,20,1,9,2,10,3,11,4,12,5,13,6,14,7,15,8,16};
		unsigned int i = 0;
		bool bOnSurface = true;
		unsigned int iNodesPerFace = 8;
		for(i = 0; i < iNodesPerFace ; i++)
		{
			if(!m_vpoNodes[viNodeIndices[(iFaceID - 1)*iNodesPerFace + i] - 1]->IsOnSurface())
			{
				bOnSurface = false;
				break;
			}
		}
		return bOnSurface;
	}
	bool FEMHexahedralElement::IsPointInside(Point* poPoint,vector<double>& vdNaturalCoordinates) const
 	{
 		vdNaturalCoordinates.clear();
 		if(!IsInAxisAlignedBoundingBox(poPoint))
 		{
 			return false;
 		}

 		vector<double> vdCoordinates = GetNaturalCoordinates(poPoint);
 		double dLimit = 1.0 + 1E-4*m_oBoundingBox.GetMinimumDimension();
 		if(fabs(vdCoordinates[0]) <= dLimit)
 		{
 			if(fabs(vdCoordinates[1]) <= dLimit)
 			{
 				if(fabs(vdCoordinates[2]) <= dLimit)
 				{
 					vdNaturalCoordinates.resize(3);
 					vdNaturalCoordinates[0] = vdCoordinates[0];
 					vdNaturalCoordinates[1] = vdCoordinates[1];
 					vdNaturalCoordinates[2] = vdCoordinates[2];
 					return true;
 				}
 			}
 		}
 		return false;
 	}
	vector<double> FEMHexahedralElement::GetNearestNodeNaturalCoordinates(Point* poPoint) const
 	{
 		double dTemp = 0.0;
 		double dMinimum = DBL_MAX;
 		unsigned int iNearestNodeIndex = 0;
 		unsigned int i = 0;
 		for(i = 0 ; i < NodesPerElement ; i++)
 		{
 			dTemp = poPoint->Distance(*m_vpoNodes[i]);
 			if(dTemp < dMinimum)
 			{
 				dMinimum = dTemp;
 				iNearestNodeIndex = i;
 			}
 		}
 		double daNodeNaturalCoordinates[60] = {-1.0,-1.0,-1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,-1.0,-1.0,-1.0,1.0,1.0,-1.0,1.0,1.0,1.0,1.0,-1.0,1.0,1.0,0.0,-1.0,-1.0,1.0,0.0,-1.0,0.0,1.0,-1.0,-1.0,0.0,-1.0,0.0,-1.0,1.0,1.0,0.0,1.0,0.0,1.0,1.0,-1.0,0.0,1.0,-1.0,-1.0,0.0,1.0,-1.0,0.0,1.0,1.0,0.0,-1.0,1.0,0.0};
 		vector<double> vdCoordinates;
 		vdCoordinates.resize(3);
 		vdCoordinates[0] = daNodeNaturalCoordinates[3*iNearestNodeIndex];
 		vdCoordinates[1] = daNodeNaturalCoordinates[3*iNearestNodeIndex + 1];
 		vdCoordinates[2] = daNodeNaturalCoordinates[3*iNearestNodeIndex + 2];
 		return vdCoordinates;
 	}
	Point FEMHexahedralElement::GetFaceCenter(const unsigned int& iFaceID) const
 	{
		Matrix oNaturalCoordinates(3,1);
		oNaturalCoordinates.Set(1,1,0.0);
		oNaturalCoordinates.Set(2,1,0.0);
		oNaturalCoordinates.Set(3,1,0.0);
 		return GetPointOnFace(iFaceID,oNaturalCoordinates);
 	}
	vector<double> FEMHexahedralElement::GetNaturalCoordinates(Point* poPoint,const unsigned int& iMaxIterationsCount) const
 	{
 		double dTolerance = 1E-5*(m_oBoundingBox.GetMinimumDimension());
 		double dError = 100.0;
 		unsigned int iIterationsCount = 0;
 		vector<double> vdNearestNodeCoordinates = GetNearestNodeNaturalCoordinates(poPoint);
 		double daCurrentSolution[3] = {0.5*vdNearestNodeCoordinates[0],0.5*vdNearestNodeCoordinates[1],0.5*vdNearestNodeCoordinates[2]};
 		Matrix oSystem;
 		Matrix oRHS(3,1);
 		Point oCurrentPoint;
 		Point oNewPoint;
 		double dX = poPoint->GetX();
 		double dY = poPoint->GetY();
 		double dZ = poPoint->GetZ();
		Matrix oTemp(3,1);
 		while((dError > dTolerance) && (iIterationsCount < iMaxIterationsCount))
 		{
 			oSystem = GetJacobianMatrix(daCurrentSolution[0],daCurrentSolution[1],daCurrentSolution[2]);
 			oSystem = oSystem.GetTranspose();
			oTemp.Set(1,1,daCurrentSolution[0]);
			oTemp.Set(2,1,daCurrentSolution[1]);
			oTemp.Set(3,1,daCurrentSolution[2]);
 			oCurrentPoint = GetPoint(oTemp);
 			oRHS.Set(1,1,dX - oCurrentPoint.GetX());
 			oRHS.Set(2,1,dY - oCurrentPoint.GetY());
 			oRHS.Set(3,1,dZ - oCurrentPoint.GetZ());
 			oRHS = Matrix::Solve3x3System(oSystem,oRHS);
 			daCurrentSolution[0] = daCurrentSolution[0] + oRHS.Get(1,1);
 			daCurrentSolution[1] = daCurrentSolution[1] + oRHS.Get(2,1);
 			daCurrentSolution[2] = daCurrentSolution[2] + oRHS.Get(3,1);
			oTemp.Set(1,1,daCurrentSolution[0]);
			oTemp.Set(2,1,daCurrentSolution[1]);
			oTemp.Set(3,1,daCurrentSolution[2]);
 			oNewPoint = GetPoint(oTemp);
 			dError = (oCurrentPoint - oNewPoint).Distance();
 			iIterationsCount = iIterationsCount + 1;
 		}
 		vector<double> vdCoordinates;
 		vdCoordinates.resize(3);
 		vdCoordinates[0] = daCurrentSolution[0];
 		vdCoordinates[1] = daCurrentSolution[1];
 		vdCoordinates[2] = daCurrentSolution[2];
 		return vdCoordinates;
 	}
	bool FEMHexahedralElement::IsSurfaceElement() const
 	{
 		unsigned int i = 0;
 		unsigned int iCount = 0;
 		for(i = 0; i < NodesPerElement ; i++)
 		{
 			if(m_vpoNodes[i]->IsOnSurface())
 			{
 				iCount = iCount + 1;
 			}
 		}
 		if(iCount >= 8)
 		{
 			return true;
 		}
 		return false;
 	}
	bool FEMHexahedralElement::IsValid() const
 	{
 		unsigned int i = 0;
 		double dJacobian = 0.0;
 		for(i = 0; i < GaussPointsCount ; i++)
 		{
			dJacobian = GetJacobianMatrixDeterminant(m_vvdGaussPointsCoordinates[i][0],m_vvdGaussPointsCoordinates[i][1],m_vvdGaussPointsCoordinates[i][2]);
			if(dJacobian < 0.0)
			{
				return false;
			}
 		}
 		return true;
 	}
	Point FEMHexahedralElement::GetCenterPoint() const
 	{
		Matrix oNaturalCoordinates(3,1);
		oNaturalCoordinates.Set(1,1,0.0);
		oNaturalCoordinates.Set(2,1,0.0);
		oNaturalCoordinates.Set(3,1,0.0);
 		return GetPoint(oNaturalCoordinates);
 	}
	Point FEMHexahedralElement::GetPointOnFace(const unsigned int& iFaceID,const Matrix& oNaturalCoordinates) const
 	{
 		Point oPoint(0.0,0.0,0.0);
 		Matrix oTempNaturalCoordinates = oNaturalCoordinates;
 		if(iFaceID == 1)
 		{	
			oTempNaturalCoordinates.Set(1,1,-1.0);
 		}
 		else if(iFaceID == 2)
 		{
 			oTempNaturalCoordinates.Set(1,1,1.0);
 		}
 		else if(iFaceID == 3)
 		{
 			oTempNaturalCoordinates.Set(2,1,-1.0);
 		}
 		else if(iFaceID == 4)
 		{
 			oTempNaturalCoordinates.Set(2,1,1.0);
 		}
 		else if(iFaceID == 5)
 		{
 			oTempNaturalCoordinates.Set(3,1,-1.0);
 		}
 		else if(iFaceID == 6)
 		{
 			oTempNaturalCoordinates.Set(3,1,1.0);
 		}
		oPoint = GetPoint(oTempNaturalCoordinates);
 		return oPoint;
 	}
	Matrix FEMHexahedralElement::GetFaceNaturalCoordinates(const unsigned int& iFaceID,const double& dXi,const double& dEta) const
	{
 		Matrix oNaturalCoordinates(3,1);
 		if(iFaceID == 1)
 		{	
			oNaturalCoordinates.Set(1,1,-1.0);
			oNaturalCoordinates.Set(2,1,dXi);
			oNaturalCoordinates.Set(3,1,dEta);
 		}
 		else if(iFaceID == 2)
 		{
			oNaturalCoordinates.Set(1,1,1.0);
			oNaturalCoordinates.Set(2,1,dXi);
			oNaturalCoordinates.Set(3,1,dEta);
 		}
 		else if(iFaceID == 3)
 		{
			oNaturalCoordinates.Set(1,1,dXi);
			oNaturalCoordinates.Set(2,1,-1.0);
			oNaturalCoordinates.Set(3,1,dEta);
 		}
 		else if(iFaceID == 4)
 		{
			oNaturalCoordinates.Set(1,1,dXi);
			oNaturalCoordinates.Set(2,1,1.0);
			oNaturalCoordinates.Set(3,1,dEta);
 		}
 		else if(iFaceID == 5)
 		{
			oNaturalCoordinates.Set(1,1,dXi);
			oNaturalCoordinates.Set(2,1,dEta);
			oNaturalCoordinates.Set(3,1,-1.0);
 		}
 		else if(iFaceID == 6)
 		{
			oNaturalCoordinates.Set(1,1,dXi);
			oNaturalCoordinates.Set(2,1,dEta);
			oNaturalCoordinates.Set(3,1,1.0);
 		}
 		return oNaturalCoordinates;
	}
	unsigned int FEMHexahedralElement::GetFacesCount() const
	{
		return FacesCount;
	}
	void FEMHexahedralElement::Read(FILE* fpFile,const vector<FEMNode*>* pvpoNodes)
	{
 		unsigned int iNode1Index = 0;
 		unsigned int iNode2Index = 0;
 		unsigned int iNode3Index = 0;
 		unsigned int iNode4Index = 0;
 		unsigned int iNode5Index = 0;
 		unsigned int iNode6Index = 0;
 		unsigned int iNode7Index = 0;
 		unsigned int iNode8Index = 0;
 		unsigned int iNode9Index = 0;
 		unsigned int iNode10Index = 0;
 		unsigned int iNode11Index = 0;
 		unsigned int iNode12Index = 0;
 		unsigned int iNode13Index = 0;
 		unsigned int iNode14Index = 0;
 		unsigned int iNode15Index = 0;
 		unsigned int iNode16Index = 0;
 		unsigned int iNode17Index = 0;
 		unsigned int iNode18Index = 0;
 		unsigned int iNode19Index = 0;
 		unsigned int iNode20Index = 0;
 		vector<FEMNode*> vpoTempNodes;
 		vpoTempNodes.resize(NodesPerElement);
 		string sRead = SupportSystem::GetRealString(500,fpFile);
 		sscanf(sRead.c_str(),"%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",&iNode1Index,&iNode2Index,&iNode3Index,&iNode4Index,&iNode5Index,&iNode6Index,&iNode7Index,&iNode8Index,&iNode9Index,&iNode10Index,&iNode11Index,&iNode12Index,&iNode13Index,&iNode14Index,&iNode15Index,&iNode16Index,&iNode17Index,&iNode18Index,&iNode19Index,&iNode20Index);		
		vpoTempNodes[0] = pvpoNodes->at(iNode1Index - 1);
 		vpoTempNodes[1] = pvpoNodes->at(iNode2Index - 1);
 		vpoTempNodes[2] = pvpoNodes->at(iNode3Index - 1);
 		vpoTempNodes[3] = pvpoNodes->at(iNode4Index - 1);
 		vpoTempNodes[4] = pvpoNodes->at(iNode5Index - 1);
 		vpoTempNodes[5] = pvpoNodes->at(iNode6Index - 1);
 		vpoTempNodes[6] = pvpoNodes->at(iNode7Index - 1);
 		vpoTempNodes[7] = pvpoNodes->at(iNode8Index - 1);
 		vpoTempNodes[8] = pvpoNodes->at(iNode9Index - 1);
 		vpoTempNodes[9] = pvpoNodes->at(iNode10Index - 1);
 		vpoTempNodes[10] = pvpoNodes->at(iNode11Index - 1);
 		vpoTempNodes[11] = pvpoNodes->at(iNode12Index - 1);
 		vpoTempNodes[12] = pvpoNodes->at(iNode13Index - 1);
 		vpoTempNodes[13] = pvpoNodes->at(iNode14Index - 1);
 		vpoTempNodes[14] = pvpoNodes->at(iNode15Index - 1);
 		vpoTempNodes[15] = pvpoNodes->at(iNode16Index - 1);
 		vpoTempNodes[16] = pvpoNodes->at(iNode17Index - 1);
 		vpoTempNodes[17] = pvpoNodes->at(iNode18Index - 1);
 		vpoTempNodes[18] = pvpoNodes->at(iNode19Index - 1);
 		vpoTempNodes[19] = pvpoNodes->at(iNode20Index - 1);
		Set(vpoTempNodes);
	}
	void FEMHexahedralElement::Write(FILE* fpFile) const
	{
		unsigned int iNodesPerElement = (unsigned int)m_vpoNodes.size();
		unsigned int i = 0;
		fprintf(fpFile,"%d",m_vpoNodes[0]->GetID());
		for(i = 1 ; i < iNodesPerElement ; i++)
		{
			fprintf(fpFile,",%d",m_vpoNodes[i]->GetID());
		}
		fprintf(fpFile,"\n");
	}
	bool FEMHexahedralElement::IsMidSideNode(const unsigned int& iIndex) const
	{
		if(iIndex <= 8)
		{
			return false;
		}
		else if(iIndex <= 20)
		{
			return true;
		}
		return false;
	}
	double FEMHexahedralElement::GetMassLumpingFactor() const
	{
		// this is the ratio of the mass concentrated at the MIDSIDE nodes to the ratio of the mass
		// concentrated at the VERTEX nodes, this is the reciprocal of the definition given by Dhondt
		// and used in CalculiX
		return 3.428179636612958;
	}
	FEMElementGeometry* FEMHexahedralElement::Clone() const
	{
		return new FEMHexahedralElement(*this);
	}
}




