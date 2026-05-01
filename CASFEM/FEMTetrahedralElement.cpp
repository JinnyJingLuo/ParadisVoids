// Ahmed M. Hussein

#include "FEMTetrahedralElement.h"
#include "cmath"
#include "float.h"
#include "MathServices.h"
#include "Tools.h"

using namespace SupportSystem;

namespace FEMSystem
{
	// Matrix FEMTetrahedralElement::m_oShapeFunctionsExtrapolations = GetShapeFunctionsExtrapolations();
// 	unsigned int FEMTetrahedralElement::NodesPerElement = 10;
// 	unsigned int FEMTetrahedralElement::FacesCount = 4;
// 	unsigned int FEMTetrahedralElement::GaussPointsCount = 14;
// 	unsigned int FEMTetrahedralElement::FaceGaussPointsCount = 6;
// 	vector< vector<double> > FEMTetrahedralElement::m_vvdGaussPointsCoordinates = GenerateBodyGaussPointsCoordinates();
// 	vector<double> FEMTetrahedralElement::m_vdGaussPointsWeights = GenerateBodyGaussPointsWeights();
// 	vector< vector< vector<double> > > FEMTetrahedralElement::m_vvdFaceGaussPointsCoordinates = GenerateFaceGaussPointsCoordinates();
// 	vector<double> FEMTetrahedralElement::m_vdFaceGaussPointsWeights = GenerateFaceGaussPointsWeights();
// 		
// 	FEMTetrahedralElement::FEMTetrahedralElement()
// 	{
// 	
// 	}
// 	FEMTetrahedralElement::FEMTetrahedralElement(const FEMTetrahedralElement& oElement)
// 	{
// 		*this = oElement;
// 	}
// 	FEMTetrahedralElement::~FEMTetrahedralElement()
// 	{
// 	
// 	}
// 	FEMTetrahedralElement& FEMTetrahedralElement::operator=(const FEMTetrahedralElement& oElement)
// 	{
// 		FEMElementGeometry::operator=(oElement);
// 		return *this;
// 	}
// 	FEMElementGeometryType FEMTetrahedralElement::GetType() const
// 	{
// 		return TetrahedralFEMElement;
// 	}
// 	double FEMTetrahedralElement::GetVolume() const
//  	{
//  		double dVolume = 0.0;
//  		unsigned int i = 0;
//  		double dJacobian = 0.0;
//  		for(i = 0; i < GaussPointsCount ; i++)
//  		{
// 			dJacobian = GetJacobianMatrixDeterminant(m_vvdGaussPointsCoordinates[i][0],m_vvdGaussPointsCoordinates[i][1],m_vvdGaussPointsCoordinates[i][2]);
// 			dVolume = dVolume + dJacobian*m_vdGaussPointsWeights[i];
//  		}
//  		return dVolume;
//  	}
// 	Point FEMTetrahedralElement::GetPoint(const Matrix& oNaturalCoordinates) const
// 	{
// 		if(oNaturalCoordinates.GetRowsCount() != 3 || oNaturalCoordinates.GetColumnsCount() != 1)
// 		{
// 			return Point();
// 		}
// 		unsigned int i = 0;
// 		Matrix oShapeFunctions = GetShapeFunctions(oNaturalCoordinates.Get(1,1),oNaturalCoordinates.Get(2,1),oNaturalCoordinates.Get(3,1));
// 		double dX = 0.0;
// 		double dY = 0.0;
// 		double dZ = 0.0;
// 		double dShapeFunctionValue = 0.0;
// 		for(i = 0; i < NodesPerElement ; i++)
// 		{
// 			dShapeFunctionValue = oShapeFunctions.Get(1,i + 1);
// 			dX = dX + dShapeFunctionValue*m_vpoNodes[i]->GetX();
// 			dY = dY + dShapeFunctionValue*m_vpoNodes[i]->GetY();
// 			dZ = dZ + dShapeFunctionValue*m_vpoNodes[i]->GetZ();
// 		}
// 		return Point(dX,dY,dZ);
// 	}
// 	vector<Point> FEMTetrahedralElement::GetPoints(const Matrix& oNaturalCoordinates) const
// 	{
//  		vector<Point> voResults;
//  		unsigned int iSize = oNaturalCoordinates.GetRowsCount();
//  		voResults.resize(iSize);
//  		if(oNaturalCoordinates.GetColumnsCount() != 3)
//  		{
//  			return voResults;
//  		}
//  		unsigned int i = 0;
//  		Matrix oNodalCoordinates(NodesPerElement,3);
//  		Matrix oPoint;
//  		for(i = 0; i < NodesPerElement ; i++)
//  		{
//  			oNodalCoordinates.Set(i + 1,1,m_vpoNodes[i]->GetX());
//  			oNodalCoordinates.Set(i + 1,2,m_vpoNodes[i]->GetY());
//  			oNodalCoordinates.Set(i + 1,3,m_vpoNodes[i]->GetZ());
//  		}
//  		Matrix oShapeFunctions;
//  		for(i = 0; i < iSize ; i++)
//  		{
//  			oShapeFunctions = GetShapeFunctions(oNaturalCoordinates.Get(i + 1,1),oNaturalCoordinates.Get(i + 1,2),oNaturalCoordinates.Get(i + 1,3));
//  			oPoint = oShapeFunctions*oNodalCoordinates;
//  			voResults[i].Set(oPoint.Get(1,1),oPoint.Get(1,2),oPoint.Get(1,3));
//  		}
//  		return voResults;
// 	}
// 	Vector FEMTetrahedralElement::GetFaceOuterNormal(const unsigned int& iFaceID,const Matrix& oNaturalCoordinates) const
// 	{
// 		Vector oNormal;
// 		if(oNaturalCoordinates.GetRowsCount() != 3 || oNaturalCoordinates.GetColumnsCount() != 1)
// 		{
// 			return oNormal;
// 		}
// 		
// 		Point oElementCenter = GetCenterPoint();
// 		double dXi = oNaturalCoordinates.Get(1,1);
// 		double dEta = oNaturalCoordinates.Get(2,1);
// 		double dZeta = oNaturalCoordinates.Get(3,1);
// 		Point oFaceCenter;
// 		Matrix oJacobian;
// 		Vector oV1;
// 		Vector oV2;
// 		
// 		oFaceCenter = GetFaceCenter(iFaceID);
// 		if(iFaceID == 1)
// 		{
// 			oJacobian = GetJacobianMatrix(0.0,dEta,dZeta);
// 			oV1.Set(oJacobian.Get(2,1),oJacobian.Get(2,2),oJacobian.Get(2,3));
// 			oV2.Set(oJacobian.Get(3,1),oJacobian.Get(3,2),oJacobian.Get(3,3));
// 		}
// 		else if(iFaceID == 2)
// 		{
// 			oJacobian = GetJacobianMatrix(1.0,dEta,dZeta);
// 			oV1.Set(oJacobian.Get(2,1),oJacobian.Get(2,2),oJacobian.Get(2,3));
// 			oV2.Set(oJacobian.Get(3,1),oJacobian.Get(3,2),oJacobian.Get(3,3));
// 		}
// 		else if(iFaceID == 3)
// 		{
// 			oJacobian = GetJacobianMatrix(dXi,-1.0,dZeta);
// 			oV1.Set(oJacobian.Get(1,1),oJacobian.Get(1,2),oJacobian.Get(1,3));
// 			oV2.Set(oJacobian.Get(3,1),oJacobian.Get(3,2),oJacobian.Get(3,3));
// 		}
// 		else if(iFaceID == 4)
// 		{
// 			oJacobian = GetJacobianMatrix(dXi,1.0,dZeta);
// 			oV1.Set(oJacobian.Get(1,1),oJacobian.Get(1,2),oJacobian.Get(1,3));
// 			oV2.Set(oJacobian.Get(3,1),oJacobian.Get(3,2),oJacobian.Get(3,3));
// 		}
// 
// 		Vector oTest(oElementCenter,oFaceCenter);
// 		oNormal = oV1^oV2;
// 		if(oNormal*oTest < 0.0)
// 		{
// 			oNormal = oNormal*(-1.0);
// 		}
// 		return oNormal;
// 	}
// 	unsigned int FEMTetrahedralElement::GetNodesCount() const
// 	{
// 		return NodesPerElement;
// 	}
// 	Matrix FEMTetrahedralElement::GetShapeFunctions(const double& dXi,const double& dEta,const double& dZeta) const
// 	{
// 		Matrix oResult(1,NodesPerElement);
// 		
// 		vector<double> vdShapeFunctions;
// 		vdShapeFunctions.resize(NodesPerElement);
// 
// 		double dTheta = 1.0 - dXi - dEta - dZeta;
// 
// 		vdShapeFunctions[0] = dXi*(2.0*dXi - 1.0);
// 		vdShapeFunctions[1] = dEta*(2.0*dEta - 1.0);
// 		vdShapeFunctions[2] = dZeta*(2.0*dZeta - 1.0);
// 		vdShapeFunctions[3] = dTheta*(2.0*dTheta - 1.0);
// 		vdShapeFunctions[4] = 4.0*dXi*dEta;
// 		vdShapeFunctions[5] = 4.0*dEta*dZeta;
// 		vdShapeFunctions[6] = 4.0*dZeta*dXi;
// 		vdShapeFunctions[7] = 4.0*dXi*dTheta;
// 		vdShapeFunctions[8] = 4.0*dEta*dTheta;
// 		vdShapeFunctions[9] = 4.0*dZeta*dTheta;
// 
// 		unsigned int i = 0;
// 		for(i = 0; i < NodesPerElement ; i++)
// 		{
// 			oResult.Set(1,i + 1,vdShapeFunctions[i]);
// 		}
// 
// 		return oResult;
// 	}
// 	Matrix FEMTetrahedralElement::GetShapeFunctionsDerivatives(const double& dXi,const double& dEta,const double& dZeta) const
// 	{
// 		Matrix oResult(3,NodesPerElement);
// 		
// 		vector<double> vdShapeFunctions;
// 		vdShapeFunctions.resize(NodesPerElement);
// 		double dTheta = 1.0 - dXi - dEta - dZeta;
// 		unsigned int i = 0;
// 
// 		// Xi derivatives
// 		vdShapeFunctions[0] = 4.0*dXi - 1.0;
// 		vdShapeFunctions[1] = 0.0;
// 		vdShapeFunctions[2] = 0.0;
// 		vdShapeFunctions[3] = -4.0*dTheta + 1.0;
// 		vdShapeFunctions[4] = 4.0*dEta;
// 		vdShapeFunctions[5] = 0.0;
// 		vdShapeFunctions[6] = 4.0*dZeta;
// 		vdShapeFunctions[7] = 4.0*(dTheta - dXi);
// 		vdShapeFunctions[8] = -4.0*dEta;
// 		vdShapeFunctions[9] = -4.0*dZeta;
// 
// 		for(i = 0; i < NodesPerElement ; i++)
// 		{
// 			oResult.Set(1,i + 1,vdShapeFunctions[i]);
// 		}
// 		// Eta derivatives
// 		vdShapeFunctions[0] = 0.0;
// 		vdShapeFunctions[1] = 4.0*dEta - 1.0;
// 		vdShapeFunctions[2] = 0.0;
// 		vdShapeFunctions[3] = -4.0*dTheta + 1.0;
// 		vdShapeFunctions[4] = 4.0*dXi;
// 		vdShapeFunctions[5] = 4.0*dZeta;
// 		vdShapeFunctions[6] = 0.0;
// 		vdShapeFunctions[7] = -4.0*dXi;
// 		vdShapeFunctions[8] = 4.0*(dTheta - dEta);
// 		vdShapeFunctions[9] = -4.0*dZeta;
// 
// 		for(i = 0; i < NodesPerElement ; i++)
// 		{
// 			oResult.Set(2,i + 1,vdShapeFunctions[i]);
// 		}
// 		// Zeta derivatives
// 		vdShapeFunctions[0] = 0.0;
// 		vdShapeFunctions[1] = 0.0;
// 		vdShapeFunctions[2] = 4.0*dZeta - 1.0;
// 		vdShapeFunctions[3] = -4.0*dTheta + 1.0;
// 		vdShapeFunctions[4] = 0.0;
// 		vdShapeFunctions[5] = 4.0*dEta;
// 		vdShapeFunctions[6] = 4.0*dXi;
// 		vdShapeFunctions[7] = -4.0*dXi;
// 		vdShapeFunctions[8] = -4.0*dEta;
// 		vdShapeFunctions[9] = 4.0*(dTheta - dZeta);
// 
// 		for(i = 0; i < NodesPerElement ; i++)
// 		{
// 			oResult.Set(3,i + 1,vdShapeFunctions[i]);
// 		}
// 
// 		return oResult;
// 	}
// 	Matrix FEMTetrahedralElement::GetJacobianMatrix(const double& dXi,const double& dEta,const double& dZeta) const
// 	{
// 		Matrix oJacobian;
// 		unsigned int i = 0;
// 		Matrix oNodalCoordinates(NodesPerElement,3);
// 		for(i = 0; i < NodesPerElement ; i++)
// 		{
// 			oNodalCoordinates.Set(i + 1,1,m_vpoNodes[i]->GetX());
// 			oNodalCoordinates.Set(i + 1,2,m_vpoNodes[i]->GetY());
// 			oNodalCoordinates.Set(i + 1,3,m_vpoNodes[i]->GetZ());
// 		}
// 		Matrix oDerivatives = GetShapeFunctionsDerivatives(dXi,dEta,dZeta);
// 		oJacobian = oDerivatives*oNodalCoordinates;
// 		return oJacobian;
// 	}
// 	double FEMTetrahedralElement::GetJacobianMatrixDeterminant(const double& dXi,const double& dEta,const double& dZeta) const
// 	{
// 		Matrix oJacobian = GetJacobianMatrix(dXi,dEta,dZeta);
// 		double dDeterminant = 0.0;
// 		double dDet1 = oJacobian.Get(2,2)*oJacobian.Get(3,3) - oJacobian.Get(2,3)*oJacobian.Get(3,2);
// 		double dDet2 = oJacobian.Get(2,1)*oJacobian.Get(3,3) - oJacobian.Get(2,3)*oJacobian.Get(3,1);
// 		double dDet3 = oJacobian.Get(2,1)*oJacobian.Get(3,2) - oJacobian.Get(2,2)*oJacobian.Get(3,1);
// 		dDeterminant = oJacobian.Get(1,1)*dDet1 - oJacobian.Get(1,2)*dDet2 + oJacobian.Get(1,3)*dDet3;
// 		return dDeterminant;
// 	}
// 	Matrix FEMTetrahedralElement::GetShapeFunctionsXYZDerivatives(const double& dXi,const double& dEta,const double& dZeta,double& dJacobian) const
// 	{
// 		Matrix oXYZDerivatives;
// 		Matrix oDerivatives = GetShapeFunctionsDerivatives(dXi,dEta,dZeta);
// 		Matrix oJacobian;
// 		unsigned int i = 0;
// 		Matrix oNodalCoordinates(NodesPerElement,3);
// 		for(i = 0; i < NodesPerElement ; i++)
// 		{
// 			oNodalCoordinates.Set(i + 1,1,m_vpoNodes[i]->GetX());
// 			oNodalCoordinates.Set(i + 1,2,m_vpoNodes[i]->GetY());
// 			oNodalCoordinates.Set(i + 1,3,m_vpoNodes[i]->GetZ());
// 		}
// 		oJacobian = oDerivatives*oNodalCoordinates;
// 
// 		Matrix oInverseJacobian = Matrix::Invert3x3Matrix(oJacobian,dJacobian);
// 
// 		oXYZDerivatives = oInverseJacobian*oDerivatives;
// 		return oXYZDerivatives;
// 	}	
// 	vector< vector<double> > FEMTetrahedralElement::GetBodyGaussPointsCoordinates() const
// 	{
// 		return m_vvdGaussPointsCoordinates;
// 	}
// 	vector<double> FEMTetrahedralElement::GetBodyGaussPointsWeights() const
// 	{
// 		return m_vdGaussPointsWeights;
// 	}
// 	vector< vector<double> > FEMTetrahedralElement::GetFaceGaussPointsCoordinates(const unsigned int& iFaceIndex) const
// 	{
// 		vector< vector<double> > vvdCoordinates;
// 		vvdCoordinates.clear();
// 		if((iFaceIndex <= 0) || (iFaceIndex > FacesCount))
// 		{
// 			return vvdCoordinates;
// 		}
// 		vvdCoordinates = m_vvdFaceGaussPointsCoordinates[iFaceIndex - 1];
// 		return vvdCoordinates;
// 	}
// 	vector<double> FEMTetrahedralElement::GetFaceGaussPointsWeights() const
// 	{
// 		return m_vdFaceGaussPointsWeights;
// 	}
// 	vector<unsigned int> FEMTetrahedralElement::GetSurfaceFacesIndices() const
// 	{
// 		vector<unsigned int> viFacesIndices;
// 		viFacesIndices.reserve(FacesCount);
// 		unsigned int i = 0;
// 		for(i = 1 ; i <= FacesCount ; i++)
// 		{
// 			if(IsFaceOnSurface(i))
// 			{
// 				viFacesIndices.push_back(i);
// 			}
// 		}
// 		return viFacesIndices;
// 	}
// 	double FEMTetrahedralElement::GetDistanceToElementCenterPoint(Point* poPoint) const
// 	{
// 		return poPoint->Distance(GetCenterPoint());
// 	}
// 	Matrix FEMTetrahedralElement::GetShapeFunctionsExtrapolationsMatrix() const
// 	{
// 		return m_oShapeFunctionsExtrapolations;
// 	}
// 	Matrix FEMTetrahedralElement::GetShapeFunctionsExtrapolations()
// 	{
// 		Matrix oExtrapolations(NodesPerElement,27);
// 		unsigned int i = 0;
// 		unsigned int j = 0;
// 		Matrix oTemp;
// 		double daCoordinates[60] = {-1.0,-1.0,-1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,-1.0,-1.0,-1.0,1.0,1.0,-1.0,1.0,1.0,1.0,1.0,-1.0,1.0,1.0,0.0,-1.0,-1.0,1.0,0.0,-1.0,0.0,1.0,-1.0,-1.0,0.0,-1.0,0.0,-1.0,1.0,1.0,0.0,1.0,0.0,1.0,1.0,-1.0,0.0,1.0,-1.0,-1.0,0.0,1.0,-1.0,0.0,1.0,1.0,0.0,-1.0,1.0,0.0};
// 
// 		for(i = 0; i < NodesPerElement ; i++)
// 		{
// 			oTemp = EvaluateExtrapolatedShapeFunctions(daCoordinates[3*i],daCoordinates[3*i + 1],daCoordinates[3*i + 2]);
// 			for(j = 1; j <= 27 ; j++)
// 			{
// 				oExtrapolations.Set(i + 1,j,oTemp.Get(j,1));
// 			}
// 		}
// 		return oExtrapolations;
// 	}
// 	Matrix FEMTetrahedralElement::EvaluateExtrapolatedShapeFunctions(const double& dXi,const double& dEta,const double& dZeta)
// 	{
// 		Matrix oN(27,1);
// 		double dK = sqrt(0.6);
// 		double dFactor1 = 1.0/8.0;
// 		double dFactor2 = 1.0/4.0;
// 		double dFactor3 = 1.0/2.0;
// 		
// 
// 		oN.Set(1,1,-dFactor1*(1.0 - dXi/dK)*(1.0 - dEta/dK)*(1.0 - dZeta/dK)*dXi/dK*dEta/dK*dZeta/dK);
// 		oN.Set(2,1,dFactor2*(1.0 - dXi/dK)*(1.0 - dEta/dK)*(1.0 - dZeta/dK*dZeta/dK)*dXi/dK*dEta/dK);
// 		oN.Set(3,1,dFactor1*(1.0 - dXi/dK)*(1.0 - dEta/dK)*(1.0 + dZeta/dK)*dXi/dK*dEta/dK*dZeta/dK);
// 		oN.Set(4,1,dFactor2*(1.0 - dXi/dK)*(1.0 - dEta/dK*dEta/dK)*(1.0 - dZeta/dK)*dXi/dK*dZeta/dK);
// 		oN.Set(5,1,-dFactor3*(1.0 - dXi/dK)*(1.0 - dEta/dK*dEta/dK)*(1.0 - dZeta/dK*dZeta/dK)*dXi/dK);
// 		oN.Set(6,1,-dFactor2*(1.0 - dXi/dK)*(1.0 - dEta/dK*dEta/dK)*(1.0 + dZeta/dK)*dXi/dK*dZeta/dK);
// 		oN.Set(7,1,dFactor1*(1.0 - dXi/dK)*(1.0 + dEta/dK)*(1.0 - dZeta/dK)*dXi/dK*dEta/dK*dZeta/dK);
// 		oN.Set(8,1,-dFactor2*(1.0 - dXi/dK)*(1.0 + dEta/dK)*(1.0 - dZeta/dK*dZeta/dK)*dXi/dK*dEta/dK);
// 		oN.Set(9,1,-dFactor1*(1.0 - dXi/dK)*(1.0 + dEta/dK)*(1.0 + dZeta/dK)*dXi/dK*dEta/dK*dZeta/dK);
// 		oN.Set(10,1,dFactor2*(1.0 - dXi/dK*dXi/dK)*(1.0 - dEta/dK)*(1.0 - dZeta/dK)*dEta/dK*dZeta/dK);
// 		oN.Set(11,1,-dFactor3*(1.0 - dXi/dK*dXi/dK)*(1.0 - dEta/dK)*(1.0 - dZeta/dK*dZeta/dK)*dEta/dK);
// 		oN.Set(12,1,-dFactor2*(1.0 - dXi/dK*dXi/dK)*(1.0 - dEta/dK)*(1.0 + dZeta/dK)*dEta/dK*dZeta/dK);
// 		oN.Set(13,1,-dFactor3*(1.0 - dXi/dK*dXi/dK)*(1.0 - dEta/dK*dEta/dK)*(1.0 - dZeta/dK)*dZeta/dK);
// 		oN.Set(14,1,(1.0 - dXi/dK*dXi/dK)*(1.0 - dEta/dK*dEta/dK)*(1.0 - dZeta/dK*dZeta/dK));
// 		oN.Set(15,1,dFactor3*(1.0 - dXi/dK*dXi/dK)*(1.0 - dEta/dK*dEta/dK)*(1.0 + dZeta/dK)*dZeta/dK);
// 		oN.Set(16,1,-dFactor2*(1.0 - dXi/dK*dXi/dK)*(1.0 + dEta/dK)*(1.0 - dZeta/dK)*dEta/dK*dZeta/dK);
// 		oN.Set(17,1,dFactor3*(1.0 - dXi/dK*dXi/dK)*(1.0 + dEta/dK)*(1.0 - dZeta/dK*dZeta/dK)*dEta/dK);
// 		oN.Set(18,1,dFactor2*(1.0 - dXi/dK*dXi/dK)*(1.0 + dEta/dK)*(1.0 + dZeta/dK)*dEta/dK*dZeta/dK);
// 		oN.Set(19,1,dFactor1*(1.0 + dXi/dK)*(1.0 - dEta/dK)*(1.0 - dZeta/dK)*dXi/dK*dEta/dK*dZeta/dK);
// 		oN.Set(20,1,-dFactor2*(1.0 + dXi/dK)*(1.0 - dEta/dK)*(1.0 - dZeta/dK*dZeta/dK)*dXi/dK*dEta/dK);
// 		oN.Set(21,1,-dFactor1*(1.0 + dXi/dK)*(1.0 - dEta/dK)*(1.0 + dZeta/dK)*dXi/dK*dEta/dK*dZeta/dK);
// 		oN.Set(22,1,-dFactor2*(1.0 + dXi/dK)*(1.0 - dEta/dK*dEta/dK)*(1.0 - dZeta/dK)*dXi/dK*dZeta/dK);
// 		oN.Set(23,1,dFactor3*(1.0 + dXi/dK)*(1.0 - dEta/dK*dEta/dK)*(1.0 - dZeta/dK*dZeta/dK)*dXi/dK);
// 		oN.Set(24,1,dFactor2*(1.0 + dXi/dK)*(1.0 - dEta/dK*dEta/dK)*(1.0 + dZeta/dK)*dXi/dK*dZeta/dK);
// 		oN.Set(25,1,-dFactor1*(1.0 + dXi/dK)*(1.0 + dEta/dK)*(1.0 - dZeta/dK)*dXi/dK*dEta/dK*dZeta/dK);
// 		oN.Set(26,1,dFactor2*(1.0 + dXi/dK)*(1.0 + dEta/dK)*(1.0 - dZeta/dK*dZeta/dK)*dXi/dK*dEta/dK);
// 		oN.Set(27,1,dFactor1*(1.0 + dXi/dK)*(1.0 + dEta/dK)*(1.0 + dZeta/dK)*dXi/dK*dEta/dK*dZeta/dK);
// 
// 		return oN;
// 	}
// 	vector< vector<double> > FEMTetrahedralElement::GenerateBodyGaussPointsCoordinates()
// 	{
// 		vector< vector<double> > vvdCoordinates;
// 		vvdCoordinates.resize(GaussPointsCount);
// 		unsigned int i = 0;
// 		for(i = 0 ; i < GaussPointsCount ; i++)
// 		{
// 			vvdCoordinates[i].resize(3);
// 		}
// 		
// 		double dAlpha = 0.09273525031089122628655892066032137;
// 		double dBeta = .31088591926330060975814749494040332;
// 		double dGamma = 0.04550370412564965000000000000000000;
// 
// 		vvdCoordinates[0][0] = dAlpha;
// 		vvdCoordinates[0][1] = dAlpha;
// 		vvdCoordinates[0][2] = dAlpha;
// 		
// 		vvdCoordinates[1][0] = dAlpha;
// 		vvdCoordinates[1][1] = dAlpha;
// 		vvdCoordinates[1][2] = 1.0 - 3.0*dAlpha;
// 		
// 		vvdCoordinates[2][0] = dAlpha;
// 		vvdCoordinates[2][1] = 1.0 - 3.0*dAlpha;
// 		vvdCoordinates[2][2] = dAlpha;
// 		
// 		vvdCoordinates[3][0] = 1.0 - 3.0*dAlpha;
// 		vvdCoordinates[3][1] = dAlpha;
// 		vvdCoordinates[3][2] = dAlpha;
// 		
// 		vvdCoordinates[4][0] = dBeta;
// 		vvdCoordinates[4][1] = dBeta;
// 		vvdCoordinates[4][2] = dBeta;
// 		
// 		vvdCoordinates[5][0] = dBeta;
// 		vvdCoordinates[5][1] = dBeta;
// 		vvdCoordinates[5][2] = 1.0 - 3.0*dBeta;
// 		
// 		vvdCoordinates[6][0] = dBeta;
// 		vvdCoordinates[6][1] = 1.0 - 3.0*dBeta;
// 		vvdCoordinates[6][2] = dBeta;
// 		
// 		vvdCoordinates[7][0] = 1.0 - 3.0*dBeta;
// 		vvdCoordinates[7][1] = dBeta;
// 		vvdCoordinates[7][2] = dBeta;
// 		
// 		vvdCoordinates[8][0] = dGamma;
// 		vvdCoordinates[8][1] = dGamma;
// 		vvdCoordinates[8][2] = 0.5 - dGamma;
// 		
// 		vvdCoordinates[9][0] = dGamma;
// 		vvdCoordinates[9][1] = 0.5 - dGamma;
// 		vvdCoordinates[9][2] = dGamma;
// 		
// 		vvdCoordinates[10][0] = dGamma;
// 		vvdCoordinates[10][1] = 0.5 - dGamma;
// 		vvdCoordinates[10][2] = 0.5 - dGamma;
// 		
// 		vvdCoordinates[11][0] = 0.5 - dGamma;
// 		vvdCoordinates[11][1] = dGamma;
// 		vvdCoordinates[11][2] = dGamma;
// 		
// 		vvdCoordinates[12][0] = 0.5 - dGamma;
// 		vvdCoordinates[12][1] = dGamma;
// 		vvdCoordinates[12][2] = 0.5 - dGamma;
// 		
// 		vvdCoordinates[13][0] = 0.5 - dGamma;
// 		vvdCoordinates[13][1] = 0.5 - dGamma;
// 		vvdCoordinates[13][2] = dGamma;
// 
// 		return vvdCoordinates;
// 	}
// 	vector<double> FEMTetrahedralElement::GenerateBodyGaussPointsWeights()
// 	{
// 		vector<double> vdWeights;
// 		vdWeights.resize(GaussPointsCount);
// 
// 		double dAlpha = 0.07349304311636194934358694586367885;
// 		double dBeta = 0.11268792571801585036501492847638892;
// 		double dGamma = 0.04254602077708146686093208377328816;
// 		
// 		vdWeights[0] = dAlpha;
// 		vdWeights[1] = dAlpha;
// 		vdWeights[2] = dAlpha;
// 		vdWeights[3] = dAlpha;
// 		
// 		vdWeights[4] = dBeta;
// 		vdWeights[5] = dBeta;
// 		vdWeights[6] = dBeta;
// 		vdWeights[7] = dBeta;
// 		
// 		vdWeights[8] = dGamma;
// 		vdWeights[9] = dGamma;
// 		vdWeights[10] = dGamma;
// 		vdWeights[11] = dGamma;
// 		vdWeights[12] = dGamma;
// 		vdWeights[13] = dGamma;
// 
// 		return vdWeights;
// 	}
// 	vector< vector< vector<double> > > FEMTetrahedralElement::GenerateFaceGaussPointsCoordinates()
// 	{
// 		vector< vector< vector<double> > > vvvdCoordinates;
// 		vvvdCoordinates.resize(FacesCount);
// 		unsigned int i = 0;
// 		unsigned int j = 0;
// 		for(i = 0 ; i < FacesCount ; i++)
// 		{
// 			vvvdCoordinates[i].resize(FaceGaussPointsCount);
// 			for(j = 0 ; j < FaceGaussPointsCount ; j++)
// 			{
// 				vvvdCoordinates[i][j].resize(3);
// 			}
// 		}
// 		
// 		unsigned int k = 0;
// 		double dTemp = sqrt(3.0/5.0);
// 		double daCoordinates[3] = {-dTemp,0.0,dTemp};
// 		double dFixedCoordinate[6] = {-1.0,1.0,-1.0,1.0,-1.0,1.0};
// 		unsigned int iaFixedIndex[6] = {0,0,1,1,2,2};
// 		unsigned int iaJIndex[6] = {1,1,0,0,0,0};
// 		unsigned int iaKIndex[6] = {2,2,2,2,1,1};
// 		unsigned int iIndex = 0;
// 		for(i = 0 ; i < FacesCount ; i++)		// loop over the faces
// 		{
// 			iIndex = 0;
// 			for(j = 0 ; j < 3 ; j++)			// loop over the face Gauss points rows
// 			{
// 				for(k = 0 ; k < 3 ; k++)		// loop over the face Gauss points columns
// 				{
// 					vvvdCoordinates[i][iIndex][iaJIndex[i]] = daCoordinates[j];
// 					vvvdCoordinates[i][iIndex][iaKIndex[i]] = daCoordinates[k];
// 					vvvdCoordinates[i][iIndex][iaFixedIndex[i]] = dFixedCoordinate[i];
// 					iIndex = iIndex + 1;
// 				}
// 			}
// 		}
// 		return vvvdCoordinates;
// 	}
// 	vector<double> FEMTetrahedralElement::GenerateFaceGaussPointsWeights()
// 	{
// 		vector<double> vdWeights;
// 		vdWeights.resize(FaceGaussPointsCount);
// 		unsigned int i = 0;		
// 		unsigned int j = 0;
// 		double dTemp = 5.0/9.0;
// 		double daWeights[3] = {dTemp,8.0/9.0,dTemp};
// 		unsigned int iIndex = 0;
// 		for(i = 0 ; i < 3 ; i++)			// loop over the face Gauss points rows
// 		{
// 			for(j = 0 ; j < 3 ; j++)		// loop over the face Gauss points columns
// 			{
// 				vdWeights[iIndex] = daWeights[i]*daWeights[j];
// 				iIndex = iIndex + 1;
// 			}
// 		}
// 		return vdWeights;
// 	}
// 	void FEMTetrahedralElement::GenerateSurfacePatches(vector< vector<unsigned int> >& vviNodesIndices,vector<GenericNode*>& vpoMidPoints) const
// 	{
// 		vector<unsigned int> viFacesIndices = GetSurfaceFacesIndices();
// 		unsigned int iSize = viFacesIndices.size();
// 		vector<unsigned int> viNodesIndices;
// 		vviNodesIndices.resize(iSize);
// 		vpoMidPoints.resize(iSize);
// 		unsigned int i = 0;
// 		for(i = 0; i < iSize ; i++)
// 		{
// 			vpoMidPoints[i] = GetFaceQuadPatch(viFacesIndices[i],vviNodesIndices[i]);
// 		}
// 	}
// 	GenericNode* FEMTetrahedralElement::GetFaceQuadPatch(const unsigned int& iFaceIndex,vector<unsigned int>& viOriginalNodesIndices) const
// 	{
// 		viOriginalNodesIndices.resize(8);
//         Point oTempPoint = GetFaceCenter(iFaceIndex);
// 		if(iFaceIndex == 1)
// 		{
// 			viOriginalNodesIndices[0] = m_vpoNodes[3]->GetID();
// 			viOriginalNodesIndices[1] = m_vpoNodes[11]->GetID();
// 			viOriginalNodesIndices[2] = m_vpoNodes[0]->GetID();
// 			viOriginalNodesIndices[3] = m_vpoNodes[16]->GetID();
// 			viOriginalNodesIndices[4] = m_vpoNodes[4]->GetID();
// 			viOriginalNodesIndices[5] = m_vpoNodes[15]->GetID();
// 			viOriginalNodesIndices[6] = m_vpoNodes[7]->GetID();
// 			viOriginalNodesIndices[7] = m_vpoNodes[19]->GetID();
// 		}
// 		else if(iFaceIndex == 2)
// 		{
// 			viOriginalNodesIndices[0] = m_vpoNodes[1]->GetID();
// 			viOriginalNodesIndices[1] = m_vpoNodes[9]->GetID();
// 			viOriginalNodesIndices[2] = m_vpoNodes[2]->GetID();
// 			viOriginalNodesIndices[3] = m_vpoNodes[18]->GetID();
// 			viOriginalNodesIndices[4] = m_vpoNodes[6]->GetID();
// 			viOriginalNodesIndices[5] = m_vpoNodes[13]->GetID();
// 			viOriginalNodesIndices[6] = m_vpoNodes[5]->GetID();
// 			viOriginalNodesIndices[7] = m_vpoNodes[17]->GetID();
// 		}
// 		else if(iFaceIndex == 3)
// 		{
// 			viOriginalNodesIndices[0] = m_vpoNodes[0]->GetID();
// 			viOriginalNodesIndices[1] = m_vpoNodes[8]->GetID();
// 			viOriginalNodesIndices[2] = m_vpoNodes[1]->GetID();
// 			viOriginalNodesIndices[3] = m_vpoNodes[17]->GetID();
// 			viOriginalNodesIndices[4] = m_vpoNodes[5]->GetID();
// 			viOriginalNodesIndices[5] = m_vpoNodes[12]->GetID();
// 			viOriginalNodesIndices[6] = m_vpoNodes[4]->GetID();
// 			viOriginalNodesIndices[7] = m_vpoNodes[16]->GetID();
// 		}
// 		else if(iFaceIndex == 4)
// 		{
// 			viOriginalNodesIndices[0] = m_vpoNodes[2]->GetID();
// 			viOriginalNodesIndices[1] = m_vpoNodes[10]->GetID();
// 			viOriginalNodesIndices[2] = m_vpoNodes[3]->GetID();
// 			viOriginalNodesIndices[3] = m_vpoNodes[19]->GetID();
// 			viOriginalNodesIndices[4] = m_vpoNodes[7]->GetID();
// 			viOriginalNodesIndices[5] = m_vpoNodes[14]->GetID();
// 			viOriginalNodesIndices[6] = m_vpoNodes[6]->GetID();
// 			viOriginalNodesIndices[7] = m_vpoNodes[18]->GetID();
// 		}
// 		else if(iFaceIndex == 5)
// 		{
// 			viOriginalNodesIndices[0] = m_vpoNodes[0]->GetID();
// 			viOriginalNodesIndices[1] = m_vpoNodes[11]->GetID();
// 			viOriginalNodesIndices[2] = m_vpoNodes[3]->GetID();
// 			viOriginalNodesIndices[3] = m_vpoNodes[10]->GetID();
// 			viOriginalNodesIndices[4] = m_vpoNodes[2]->GetID();
// 			viOriginalNodesIndices[5] = m_vpoNodes[9]->GetID();
// 			viOriginalNodesIndices[6] = m_vpoNodes[1]->GetID();
// 			viOriginalNodesIndices[7] = m_vpoNodes[8]->GetID();
// 		}
// 		else if(iFaceIndex == 6)
// 		{
// 			viOriginalNodesIndices[0] = m_vpoNodes[5]->GetID();
// 			viOriginalNodesIndices[1] = m_vpoNodes[13]->GetID();
// 			viOriginalNodesIndices[2] = m_vpoNodes[6]->GetID();
// 			viOriginalNodesIndices[3] = m_vpoNodes[14]->GetID();
// 			viOriginalNodesIndices[4] = m_vpoNodes[7]->GetID();
// 			viOriginalNodesIndices[5] = m_vpoNodes[15]->GetID();
// 			viOriginalNodesIndices[6] = m_vpoNodes[4]->GetID();
// 			viOriginalNodesIndices[7] = m_vpoNodes[12]->GetID();
// 		}
// 		else
// 		{
// 			viOriginalNodesIndices.clear();
// 			return NULL;
// 		}
//    	 	GenericNode* poMidPoint = new GenericNode;
//         poMidPoint->Set(oTempPoint.GetX(),oTempPoint.GetY(),oTempPoint.GetZ());
// 		return poMidPoint;
// 	}
// 	bool FEMTetrahedralElement::IsFaceOnSurface(const unsigned int& iFaceID) const
// 	{
// 		unsigned int viNodeIndices[48] = {1,12,4,20,8,16,5,17,2,10,3,19,7,14,6,18,1,9,2,18,6,13,5,17,4,11,3,19,7,15,8,20,1,9,2,10,3,11,4,12,5,13,6,14,7,15,8,16};
// 		unsigned int i = 0;
// 		bool bOnSurface = true;
// 		unsigned int iNodesPerFace = 8;
// 		for(i = 0; i < iNodesPerFace ; i++)
// 		{
// 			if(!m_vpoNodes[viNodeIndices[(iFaceID - 1)*iNodesPerFace + i] - 1]->IsOnSurface())
// 			{
// 				bOnSurface = false;
// 				break;
// 			}
// 		}
// 		return bOnSurface;
// 	}
// 	bool FEMTetrahedralElement::IsPointInside(Point* poPoint,vector<double>& vdNaturalCoordinates) const
//  	{
//  		vdNaturalCoordinates.clear();
//  		if(!IsInAxisAlignedBoundingBox(poPoint))
//  		{
//  			return false;
//  		}
// 
//  		vector<double> vdCoordinates = GetNaturalCoordinates(poPoint);
//  		double dLimit = 1.0 + 1E-4*m_oBoundingBox.GetMinimumDimension();
//  		if(fabs(vdCoordinates[0]) <= dLimit)
//  		{
//  			if(fabs(vdCoordinates[1]) <= dLimit)
//  			{
//  				if(fabs(vdCoordinates[2]) <= dLimit)
//  				{
//  					vdNaturalCoordinates.resize(3);
//  					vdNaturalCoordinates[0] = vdCoordinates[0];
//  					vdNaturalCoordinates[1] = vdCoordinates[1];
//  					vdNaturalCoordinates[2] = vdCoordinates[2];
//  					return true;
//  				}
//  			}
//  		}
//  		return false;
//  	}
// 	vector<double> FEMTetrahedralElement::GetNearestNodeNaturalCoordinates(Point* poPoint) const
//  	{
//  		double dTemp = 0.0;
//  		double dMinimum = DBL_MAX;
//  		unsigned int iNearestNodeIndex = 0;
//  		unsigned int i = 0;
//  		for(i = 0 ; i < NodesPerElement ; i++)
//  		{
//  			dTemp = poPoint->Distance(*m_vpoNodes[i]);
//  			if(dTemp < dMinimum)
//  			{
//  				dMinimum = dTemp;
//  				iNearestNodeIndex = i;
//  			}
//  		}
//  		double daNodeNaturalCoordinates[60] = {-1.0,-1.0,-1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,-1.0,-1.0,-1.0,1.0,1.0,-1.0,1.0,1.0,1.0,1.0,-1.0,1.0,1.0,0.0,-1.0,-1.0,1.0,0.0,-1.0,0.0,1.0,-1.0,-1.0,0.0,-1.0,0.0,-1.0,1.0,1.0,0.0,1.0,0.0,1.0,1.0,-1.0,0.0,1.0,-1.0,-1.0,0.0,1.0,-1.0,0.0,1.0,1.0,0.0,-1.0,1.0,0.0};
//  		vector<double> vdCoordinates;
//  		vdCoordinates.resize(3);
//  		vdCoordinates[0] = daNodeNaturalCoordinates[3*iNearestNodeIndex];
//  		vdCoordinates[1] = daNodeNaturalCoordinates[3*iNearestNodeIndex + 1];
//  		vdCoordinates[2] = daNodeNaturalCoordinates[3*iNearestNodeIndex + 2];
//  		return vdCoordinates;
//  	}
// 	Point FEMTetrahedralElement::GetFaceCenter(const unsigned int& iFaceID) const
//  	{
// 		Matrix oNaturalCoordinates(3,1);
// 		oNaturalCoordinates.Set(1,1,0.0);
// 		oNaturalCoordinates.Set(2,1,0.0);
// 		oNaturalCoordinates.Set(3,1,0.0);
//  		return GetPointOnFace(iFaceID,oNaturalCoordinates);
//  	}
// 	vector<double> FEMTetrahedralElement::GetNaturalCoordinates(Point* poPoint,const unsigned int& iMaxIterationsCount) const
//  	{
//  		double dTolerance = 1E-5*(m_oBoundingBox.GetMinimumDimension());
//  		double dError = 100.0;
//  		unsigned int iIterationsCount = 0;
//  		vector<double> vdNearestNodeCoordinates = GetNearestNodeNaturalCoordinates(poPoint);
//  		double daCurrentSolution[3] = {0.5*vdNearestNodeCoordinates[0],0.5*vdNearestNodeCoordinates[1],0.5*vdNearestNodeCoordinates[2]};
//  		Matrix oSystem;
//  		Matrix oRHS(3,1);
//  		Point oCurrentPoint;
//  		Point oNewPoint;
//  		double dX = poPoint->GetX();
//  		double dY = poPoint->GetY();
//  		double dZ = poPoint->GetZ();
// 		Matrix oTemp(3,1);
//  		while((dError > dTolerance) && (iIterationsCount < iMaxIterationsCount))
//  		{
//  			oSystem = GetJacobianMatrix(daCurrentSolution[0],daCurrentSolution[1],daCurrentSolution[2]);
//  			oSystem = oSystem.GetTranspose();
// 			oTemp.Set(1,1,daCurrentSolution[0]);
// 			oTemp.Set(2,1,daCurrentSolution[1]);
// 			oTemp.Set(3,1,daCurrentSolution[2]);
//  			oCurrentPoint = GetPoint(oTemp);
//  			oRHS.Set(1,1,dX - oCurrentPoint.GetX());
//  			oRHS.Set(2,1,dY - oCurrentPoint.GetY());
//  			oRHS.Set(3,1,dZ - oCurrentPoint.GetZ());
//  			oRHS = Matrix::Solve3x3System(oSystem,oRHS);
//  			daCurrentSolution[0] = daCurrentSolution[0] + oRHS.Get(1,1);
//  			daCurrentSolution[1] = daCurrentSolution[1] + oRHS.Get(2,1);
//  			daCurrentSolution[2] = daCurrentSolution[2] + oRHS.Get(3,1);
// 			oTemp.Set(1,1,daCurrentSolution[0]);
// 			oTemp.Set(2,1,daCurrentSolution[1]);
// 			oTemp.Set(3,1,daCurrentSolution[2]);
//  			oNewPoint = GetPoint(oTemp);
//  			dError = (oCurrentPoint - oNewPoint).Distance();
//  			iIterationsCount = iIterationsCount + 1;
//  		}
//  		vector<double> vdCoordinates;
//  		vdCoordinates.resize(3);
//  		vdCoordinates[0] = daCurrentSolution[0];
//  		vdCoordinates[1] = daCurrentSolution[1];
//  		vdCoordinates[2] = daCurrentSolution[2];
//  		return vdCoordinates;
//  	}
// 	bool FEMTetrahedralElement::IsSurfaceElement() const
//  	{
//  		unsigned int i = 0;
//  		unsigned int iCount = 0;
//  		for(i = 0; i < NodesPerElement ; i++)
//  		{
//  			if(m_vpoNodes[i]->IsOnSurface())
//  			{
//  				iCount = iCount + 1;
//  			}
//  		}
//  		if(iCount >= 6)
//  		{
//  			return true;
//  		}
//  		return false;
//  	}
// 	bool FEMTetrahedralElement::IsValid() const
//  	{
//  		unsigned int i = 0;
//  		double dJacobian = 0.0;
//  		for(i = 0; i < GaussPointsCount ; i++)
//  		{
// 			dJacobian = GetJacobianMatrixDeterminant(m_vvdGaussPointsCoordinates[i][0],m_vvdGaussPointsCoordinates[i][1],m_vvdGaussPointsCoordinates[i][2]);
// 			if(dJacobian < 0.0)
// 			{
// 				return false;
// 			}
//  		}
//  		return true;
//  	}
// 	Point FEMTetrahedralElement::GetCenterPoint() const
//  	{
// 		Matrix oNaturalCoordinates(3,1);
// 		oNaturalCoordinates.Set(1,1,0.0);
// 		oNaturalCoordinates.Set(2,1,0.0);
// 		oNaturalCoordinates.Set(3,1,0.0);
//  		return GetPoint(oNaturalCoordinates);
//  	}
// 	Point FEMTetrahedralElement::GetPointOnFace(const unsigned int& iFaceID,const Matrix& oNaturalCoordinates) const
//  	{
//  		Point oPoint(0.0,0.0,0.0);
//  		Matrix oTempNaturalCoordinates = oNaturalCoordinates;
//  		if(iFaceID == 1)
//  		{	
// 			oTempNaturalCoordinates.Set(1,1,-1.0);
//  		}
//  		else if(iFaceID == 2)
//  		{
//  			oTempNaturalCoordinates.Set(1,1,1.0);
//  		}
//  		else if(iFaceID == 3)
//  		{
//  			oTempNaturalCoordinates.Set(2,1,-1.0);
//  		}
//  		else if(iFaceID == 4)
//  		{
//  			oTempNaturalCoordinates.Set(2,1,1.0);
//  		}
// 		oPoint = GetPoint(oTempNaturalCoordinates);
//  		return oPoint;
//  	}
// 	unsigned int FEMTetrahedralElement::GetFacesCount() const
// 	{
// 		return FacesCount;
// 	}
// 	void FEMTetrahedralElement::Read(FILE* fpFile,const vector<FEMNode*>* pvpoNodes)
// 	{
//  		unsigned int iNode1Index = 0;
//  		unsigned int iNode2Index = 0;
//  		unsigned int iNode3Index = 0;
//  		unsigned int iNode4Index = 0;
//  		unsigned int iNode5Index = 0;
//  		unsigned int iNode6Index = 0;
//  		unsigned int iNode7Index = 0;
//  		unsigned int iNode8Index = 0;
//  		unsigned int iNode9Index = 0;
//  		unsigned int iNode10Index = 0;
//  		vector<FEMNode*> vpoTempNodes;
//  		vpoTempNodes.resize(NodesPerElement);
//  		string sRead = SupportSystem::GetRealString(500,fpFile);
//  		sscanf(sRead.c_str(),"%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",&iNode1Index,&iNode2Index,&iNode3Index,&iNode4Index,&iNode5Index,&iNode6Index,&iNode7Index,&iNode8Index,&iNode9Index,&iNode10Index);		
// 		vpoTempNodes[0] = pvpoNodes->at(iNode1Index - 1);
//  		vpoTempNodes[1] = pvpoNodes->at(iNode2Index - 1);
//  		vpoTempNodes[2] = pvpoNodes->at(iNode3Index - 1);
//  		vpoTempNodes[3] = pvpoNodes->at(iNode4Index - 1);
//  		vpoTempNodes[4] = pvpoNodes->at(iNode5Index - 1);
//  		vpoTempNodes[5] = pvpoNodes->at(iNode6Index - 1);
//  		vpoTempNodes[6] = pvpoNodes->at(iNode7Index - 1);
//  		vpoTempNodes[7] = pvpoNodes->at(iNode8Index - 1);
//  		vpoTempNodes[8] = pvpoNodes->at(iNode9Index - 1);
//  		vpoTempNodes[9] = pvpoNodes->at(iNode10Index - 1);
// 		Set(vpoTempNodes);
// 	}
}




