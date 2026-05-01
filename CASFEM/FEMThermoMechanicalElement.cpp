// Ahmed M. Hussein
#include "FEMThermoMechanicalElement.h"
#include "FEMThermoMechanicalNode.h"
#include "Tools.h"

using namespace EZ;
using namespace std;
using namespace SupportSystem;

namespace FEMSystem
{
	unsigned int FEMThermoMechanicalElement::ThermoMechanicalDOFPerNode = 4;
	FEMThermoMechanicalElement::FEMThermoMechanicalElement()
	{
		Initialize();
	}
	FEMThermoMechanicalElement::FEMThermoMechanicalElement(const FEMThermoMechanicalElement& oElement)
	{
		*this = oElement;
	}
	FEMThermoMechanicalElement::~FEMThermoMechanicalElement()
	{
		Reset();
	}
	FEMThermoMechanicalElement& FEMThermoMechanicalElement::operator=(const FEMThermoMechanicalElement& oElement)
	{
		FEMSolidElement::operator=(oElement);
		m_vpoFacesTLoads = oElement.m_vpoFacesTLoads;
		m_poBodyTLoad = oElement.m_poBodyTLoad;
		return *this;
	}
	void FEMThermoMechanicalElement::Reset()
	{
		FEMSolidElement::Reset();
		Initialize();
	}
	void FEMThermoMechanicalElement::SetFacesLoadsSize()
	{
		FEMSolidElement::SetFacesLoadsSize();
		unsigned int iFacesPerElement = m_poElementGeometry->GetFacesCount();
 		m_vpoFacesTLoads.resize(iFacesPerElement);
		unsigned int i = 0;
 		for(i = 0; i < iFacesPerElement ; i++)
 		{
 			m_vpoFacesTLoads[i] = NULL;
 		}
	}
	Matrix FEMThermoMechanicalElement::GetConductionMatrix() const
	{
		unsigned int iElementDOF = GetThermalDegreesOfFreedomCount();
 		Matrix K(iElementDOF,iElementDOF);
 		vector< vector<double> > vvdGaussPointsCoordinates = m_poElementGeometry->GetBodyGaussPointsCoordinates();
 		vector<double> vdGaussPointsWeights = m_poElementGeometry->GetBodyGaussPointsWeights();
 		unsigned int iGaussPointsCount = (unsigned int)vvdGaussPointsCoordinates.size();
 		Matrix oMaterialConductionMatrix = m_poMaterial->GetConductionMatrix();
 		unsigned int i = 0;
 		Matrix B;
 		double dXi = 0.0;
 		double dEta = 0.0;
 		double dZeta = 0.0;
 		double dJacobian = 0.0;
 		double dTemp = 0.0;
 		for(i = 0; i < iGaussPointsCount ; i++)
 		{
 			dXi = vvdGaussPointsCoordinates[i][0];
 			dEta = vvdGaussPointsCoordinates[i][1];
 			dZeta = vvdGaussPointsCoordinates[i][2];
 			B = GetThermalGradientTransformation(dXi,dEta,dZeta,dJacobian);
 			dTemp = dJacobian*vdGaussPointsWeights[i];
 			K = K + (B.GetTranspose()*oMaterialConductionMatrix*B)*(dTemp);
 		}
 		return K;
	}
	Matrix FEMThermoMechanicalElement::GetStressTemperatureCouplingMatrix() const
	{
		unsigned int iMechanicalDOF = GetSolidDegreesOfFreedomCount();
		unsigned int iThermalDOF = GetThermalDegreesOfFreedomCount();
 		Matrix K(iMechanicalDOF,iThermalDOF);
 		vector< vector<double> > vvdGaussPointsCoordinates = m_poElementGeometry->GetBodyGaussPointsCoordinates();
 		vector<double> vdGaussPointsWeights = m_poElementGeometry->GetBodyGaussPointsWeights();
 		unsigned int iGaussPointsCount = (unsigned int)vvdGaussPointsCoordinates.size();
 		Matrix oMaterialCouplingMatrix = m_poMaterial->GetStressTemperatureCoupling();
 		unsigned int i = 0;
 		Matrix N;
 		Matrix B;
 		double dXi = 0.0;
 		double dEta = 0.0;
 		double dZeta = 0.0;
 		double dJacobian = 0.0;
 		double dTemp = 0.0;
 		for(i = 0; i < iGaussPointsCount ; i++)
 		{
 			dXi = vvdGaussPointsCoordinates[i][0];
 			dEta = vvdGaussPointsCoordinates[i][1];
 			dZeta = vvdGaussPointsCoordinates[i][2];
			N = GetThermalShapeFunctions(dXi,dEta,dZeta);
			B = GetStrainTransformation(dXi,dEta,dZeta,dJacobian);
 			dTemp = dJacobian*vdGaussPointsWeights[i];
 			K = K + (B.GetTranspose()*oMaterialCouplingMatrix*N)*(dTemp);
 		}
 		return K;
	}
	Matrix FEMThermoMechanicalElement::GetTemperatureStrainRateCouplingMatrix() const
	{
		unsigned int iMechanicalDOF = GetSolidDegreesOfFreedomCount();
		unsigned int iThermalDOF = GetThermalDegreesOfFreedomCount();
 		Matrix C(iThermalDOF,iMechanicalDOF);
 		vector< vector<double> > vvdGaussPointsCoordinates = m_poElementGeometry->GetBodyGaussPointsCoordinates();
 		vector<double> vdGaussPointsWeights = m_poElementGeometry->GetBodyGaussPointsWeights();
 		unsigned int iGaussPointsCount = (unsigned int)vvdGaussPointsCoordinates.size();
 		Matrix oMaterialCouplingMatrix = m_poMaterial->GetStressTemperatureCoupling();
 		unsigned int i = 0;
 		Matrix N;
 		Matrix B;
 		double dXi = 0.0;
 		double dEta = 0.0;
 		double dZeta = 0.0;
 		double dJacobian = 0.0;
 		double dTemp = 0.0;
 		for(i = 0; i < iGaussPointsCount ; i++)
 		{
 			dXi = vvdGaussPointsCoordinates[i][0];
 			dEta = vvdGaussPointsCoordinates[i][1];
 			dZeta = vvdGaussPointsCoordinates[i][2];
			N = GetThermalShapeFunctions(dXi,dEta,dZeta);
			B = GetStrainTransformation(dXi,dEta,dZeta,dJacobian);
 			dTemp = dJacobian*vdGaussPointsWeights[i];
 			C = C + (N.GetTranspose()*(oMaterialCouplingMatrix.GetTranspose()*B))*(dTemp);
 		}
 		C = C*m_poMaterial->GetReferenceTemperature();
 		return C;
	}
	Matrix FEMThermoMechanicalElement::GetThermalDampingMatrix() const
	{
		unsigned int iThermalDOF = GetThermalDegreesOfFreedomCount();
 		Matrix C(iThermalDOF,iThermalDOF);
 		vector< vector<double> > vvdGaussPointsCoordinates = m_poElementGeometry->GetBodyGaussPointsCoordinates();
 		vector<double> vdGaussPointsWeights = m_poElementGeometry->GetBodyGaussPointsWeights();
 		unsigned int iGaussPointsCount = (unsigned int)vvdGaussPointsCoordinates.size();
 		unsigned int i = 0;
 		Matrix N;
 		double dXi = 0.0;
 		double dEta = 0.0;
 		double dZeta = 0.0;
 		double dJacobian = 0.0;
 		double dTemp = 0.0;
 		for(i = 0; i < iGaussPointsCount ; i++)
 		{
 			dXi = vvdGaussPointsCoordinates[i][0];
 			dEta = vvdGaussPointsCoordinates[i][1];
 			dZeta = vvdGaussPointsCoordinates[i][2];
			N = GetThermalShapeFunctions(dXi,dEta,dZeta);
			dJacobian = m_poElementGeometry->GetJacobianMatrixDeterminant(dXi,dEta,dZeta);
 			dTemp = dJacobian*vdGaussPointsWeights[i];
 			C = C + (N.GetTranspose()*N)*dTemp;
 		}
 		C = C*(m_poMaterial->GetSpecificHeatCapacity()*m_poMaterial->GetMassDensity());
 		return C;
	}
	Matrix FEMThermoMechanicalElement::GetNodalFluxes() const
	{
	 	unsigned int iFluxesCount = 3;
		unsigned int iElementDOF = GetThermalDegreesOfFreedomCount();
		unsigned int iNodesPerElement = m_poElementGeometry->GetNodesCount();
 		Matrix U(iElementDOF,1);
 		unsigned int i = 0;
		vector<FEMNode*>* pvpoNodes = m_poElementGeometry->GetNodes();
 		for(i = 1; i <= iNodesPerElement ; i++)
 		{
 			U.Set(i,1,((FEMThermoMechanicalNode*)pvpoNodes->at(i - 1))->GetTemperature());
 		}
		vector< vector<double> > vvdGaussPointsCoordinates = m_poElementGeometry->GetBodyGaussPointsCoordinates();
 		unsigned int iGaussPointsCount = (unsigned int)vvdGaussPointsCoordinates.size();
 		Matrix oTempFluxes;
 		unsigned int j = 0;
 		unsigned int iCount = 0;
 		double dXi = 0.0;
 		double dEta = 0.0;
 		double dZeta = 0.0;
 		Matrix oGaussPointsFluxes(iGaussPointsCount,iFluxesCount);
 		double dDummyJacobian = 0.0;
 		Matrix B;
 		Matrix oConduction = m_poMaterial->GetConductionMatrix();
 		for(i = 0 ; i < iGaussPointsCount ; i++)
 		{
 			dXi = vvdGaussPointsCoordinates[i][0];
 			dEta = vvdGaussPointsCoordinates[i][1];
 			dZeta = vvdGaussPointsCoordinates[i][2];
			B = GetThermalGradientTransformation(dXi,dEta,dZeta,dDummyJacobian);
			oTempFluxes = oConduction*B*U;
			iCount = iCount + 1;
			for(j = 1; j <= iFluxesCount; j++)
			{
				oGaussPointsFluxes.Set(iCount,j,oTempFluxes.Get(j,1));
			}
 		}
 		Matrix oFluxes = m_poElementGeometry->GetShapeFunctionsExtrapolationsMatrix()*oGaussPointsFluxes;
 		return oFluxes;
	}
	double FEMThermoMechanicalElement::GetTemperature(const Matrix& oNaturalCoordinates) const
	{
	 	double dTemperature = 0.0;
		vector<FEMNode*>* pvpoNodes = m_poElementGeometry->GetNodes();
		unsigned int iNodesPerElement = m_poElementGeometry->GetNodesCount();
 		if((oNaturalCoordinates.GetRowsCount() != 3) || (oNaturalCoordinates.GetColumnsCount() != 1))
 		{
 			return dTemperature;
 		}
 		Matrix oShapeFunctions = m_poElementGeometry->GetShapeFunctions(oNaturalCoordinates.Get(1,1),oNaturalCoordinates.Get(2,1),oNaturalCoordinates.Get(3,1));
 		unsigned int i = 0;
 		for(i = 0; i < iNodesPerElement ; i++)
 		{
 			dTemperature = dTemperature + ((FEMThermoMechanicalNode*)pvpoNodes->at(i))->GetTemperature()*oShapeFunctions.Get(1,i + 1);
 		}
 		return dTemperature;
	}
	Vector FEMThermoMechanicalElement::GetFlux(const Matrix& oNaturalCoordinates) const
	{
	 	Vector oFlux(0.0,0.0,0.0);
		vector<FEMNode*>* pvpoNodes = m_poElementGeometry->GetNodes();
		unsigned int iNodesPerElement = m_poElementGeometry->GetNodesCount();
 		if((oNaturalCoordinates.GetRowsCount() != 3) || (oNaturalCoordinates.GetColumnsCount() != 1))
 		{
 			return oFlux;
 		}
 		Matrix oShapeFunctions = m_poElementGeometry->GetShapeFunctions(oNaturalCoordinates.Get(1,1),oNaturalCoordinates.Get(2,1),oNaturalCoordinates.Get(3,1));
 		unsigned int i = 0;
 		for(i = 0; i < iNodesPerElement ; i++)
 		{
 			oFlux = oFlux + (((FEMThermoMechanicalNode*)pvpoNodes->at(i))->GetFluxes())*oShapeFunctions.Get(1,i + 1);
 		}
 		return oFlux;
	}
	void FEMThermoMechanicalElement::SetFacesLoads(vector<FEMLoad*>* pvpoXLoads,vector<FEMLoad*>* pvpoYLoads,vector<FEMLoad*>* pvpoZLoads,vector<FEMLoad*>* pvpoTLoads)
	{
		unsigned int iFacesCount = m_poElementGeometry->GetFacesCount();
		if((pvpoXLoads->size() != iFacesCount) || (pvpoYLoads->size() != iFacesCount) || (pvpoZLoads->size() != iFacesCount) || (pvpoTLoads->size() != iFacesCount))
		{
			return;
		}
		unsigned int i = 0;
		for(i = 0 ; i < iFacesCount ; i++)
		{
			m_vpoFacesXLoads[i] = pvpoXLoads->at(i);
			m_vpoFacesYLoads[i] = pvpoYLoads->at(i);
			m_vpoFacesZLoads[i] = pvpoZLoads->at(i);
			m_vpoFacesTLoads[i] = pvpoTLoads->at(i);
		}
	}
	void FEMThermoMechanicalElement::SetBodyLoads(FEMLoad* poXLoad,FEMLoad* poYLoad,FEMLoad* poZLoad,FEMLoad* poTLoad)
	{
		m_poBodyXLoad = poXLoad;
		m_poBodyYLoad = poYLoad;
		m_poBodyZLoad = poZLoad;
		m_poBodyTLoad = poTLoad;
	}
	FEMElementType FEMThermoMechanicalElement::GetType() const
	{
		return ThermoMechanicalFEMElement;
	}
	unsigned int FEMThermoMechanicalElement::GetDegreesOfFreedomCount() const
	{
		return (m_poElementGeometry->GetNodesCount()*ThermoMechanicalDOFPerNode);
	}
	unsigned int FEMThermoMechanicalElement::GetThermalDegreesOfFreedomCount() const
	{
		return (m_poElementGeometry->GetNodesCount());
	}
	void FEMThermoMechanicalElement::ApplyLoads(const double& dTime)
	{
		ApplySolidLoads(dTime);
		ApplyThermalLoads(dTime);
	}
	void FEMThermoMechanicalElement::ApplySolidLoads(const double& dTime)
	{
		Matrix oShapeFunctions;
  		unsigned int i = 0;
  		unsigned int j = 0;
  		unsigned int k = 0;
  		double dTemp = 0.0;
		unsigned int iElementDOF = GetSolidDegreesOfFreedomCount();
  		Matrix oForces(iElementDOF,1);
  		Matrix oTempForces(iElementDOF,1);
 		// face tractions
 		Vector oNormal;
  		vector< vector<double> > vvdGaussPointsCoordinates;
  		vector<double> vdGaussPointsWeights = m_poElementGeometry->GetFaceGaussPointsWeights();
  		Vector oTractions;
  		double dReducedJacobian = 0.0;
  		double dX = 0.0;
 		double dY = 0.0;
 		double dZ = 0.0;
  		Point oPoint;
  		unsigned int iGaussPointsCount = 0;
  		Matrix oNaturalCoordinates(3,1);
		unsigned int iFacesPerElement = m_poElementGeometry->GetFacesCount();
		unsigned int iNodesPerElement = m_poElementGeometry->GetNodesCount();
  		for(i = 1 ; i <= iFacesPerElement ; i++)
  		{
  			vvdGaussPointsCoordinates = m_poElementGeometry->GetFaceGaussPointsCoordinates(i);
  			iGaussPointsCount = (unsigned int)vvdGaussPointsCoordinates.size();
  			for(j = 0 ; j < iGaussPointsCount ; j++)
  			{
  				oNaturalCoordinates.Set(1,1,vvdGaussPointsCoordinates[j][0]);
  				oNaturalCoordinates.Set(2,1,vvdGaussPointsCoordinates[j][1]);
  				oNaturalCoordinates.Set(3,1,vvdGaussPointsCoordinates[j][2]);
   				oNormal = m_poElementGeometry->GetFaceOuterNormal(i,oNaturalCoordinates);
   				dReducedJacobian = oNormal.Length();
 				oPoint = m_poElementGeometry->GetPointOnFace(i,oNaturalCoordinates);
 				dX = m_vpoFacesXLoads[i - 1]->Get(oPoint,dTime);
  				dY = m_vpoFacesYLoads[i - 1]->Get(oPoint,dTime);
  				dZ = m_vpoFacesZLoads[i - 1]->Get(oPoint,dTime);
  				oTractions.Set(dX,dY,dZ);
  				oShapeFunctions = m_poElementGeometry->GetShapeFunctions(vvdGaussPointsCoordinates[j][0],vvdGaussPointsCoordinates[j][1],vvdGaussPointsCoordinates[j][2]);
  				for(k = 0; k < iNodesPerElement ; k++)
  				{
  					dTemp = oShapeFunctions.Get(1,k + 1);
  					oTempForces.Set(3*k + 1,1,dTemp*oTractions.GetX());
  					oTempForces.Set(3*k + 2,1,dTemp*oTractions.GetY());
  					oTempForces.Set(3*k + 3,1,dTemp*oTractions.GetZ());
  				}
  				oForces = oForces + oTempForces*dReducedJacobian*vdGaussPointsWeights[j];
  			}
  		}
 
  		Vector oNodalForce;
		vector<FEMNode*>* pvpoNodes = m_poElementGeometry->GetNodes();
  		for(i = 0; i < iNodesPerElement ; i++)
  		{
  			oNodalForce.SetX(oForces.Get(3*i + 1,1));
  			oNodalForce.SetY(oForces.Get(3*i + 2,1));
  			oNodalForce.SetZ(oForces.Get(3*i + 3,1));
  			((FEMThermoMechanicalNode*)pvpoNodes->at(i))->AddForce(oNodalForce);
  		}
 		
 		// body forces
 		oForces.SetSize(iElementDOF,1);
 		vvdGaussPointsCoordinates = m_poElementGeometry->GetBodyGaussPointsCoordinates();
 		vdGaussPointsWeights = m_poElementGeometry->GetBodyGaussPointsWeights();
 		iGaussPointsCount = (unsigned int)vvdGaussPointsCoordinates.size();
 		double dJacobian = 0.0;
 		Vector oBodyForce;
 		for(i = 0 ; i < iGaussPointsCount ; i++)
 		{
			oNaturalCoordinates.Set(1,1,vvdGaussPointsCoordinates[i][0]);
			oNaturalCoordinates.Set(2,1,vvdGaussPointsCoordinates[i][1]);
			oNaturalCoordinates.Set(3,1,vvdGaussPointsCoordinates[i][2]);
			oPoint = m_poElementGeometry->GetPoint(oNaturalCoordinates);
			dX = m_poBodyXLoad->Get(oPoint,dTime);
			dY = m_poBodyYLoad->Get(oPoint,dTime);
			dZ = m_poBodyZLoad->Get(oPoint,dTime);
			oBodyForce.Set(dX,dY,dZ);
			oShapeFunctions = m_poElementGeometry->GetShapeFunctions(vvdGaussPointsCoordinates[i][0],vvdGaussPointsCoordinates[i][1],vvdGaussPointsCoordinates[i][2]);
			dJacobian = m_poElementGeometry->GetJacobianMatrixDeterminant(vvdGaussPointsCoordinates[i][0],vvdGaussPointsCoordinates[i][1],vvdGaussPointsCoordinates[i][2]);
			for(j = 0; j < iNodesPerElement ; j++)
			{
				dTemp = oShapeFunctions.Get(1,j + 1);
				oTempForces.Set(3*j + 1,1,dTemp*oBodyForce.GetX());
				oTempForces.Set(3*j + 2,1,dTemp*oBodyForce.GetY());
				oTempForces.Set(3*j + 3,1,dTemp*oBodyForce.GetZ());
			}
			oForces = oForces + oTempForces*dJacobian*vdGaussPointsWeights[i];
 		}
 
 		for(i = 0; i < iNodesPerElement ; i++)
 		{
 			oNodalForce.SetX(oForces.Get(3*i + 1,1));
 			oNodalForce.SetY(oForces.Get(3*i + 2,1));
 			oNodalForce.SetZ(oForces.Get(3*i + 3,1));
 			((FEMThermoMechanicalNode*)pvpoNodes->at(i))->AddForce(oNodalForce);
 		}
	}
	void FEMThermoMechanicalElement::ApplyThermalLoads(const double& dTime)
	{	
	 	Matrix oShapeFunctions;
  		unsigned int i = 0;
  		unsigned int j = 0;
  		unsigned int k = 0;
		unsigned int iElementDOF = GetThermalDegreesOfFreedomCount();
  		Matrix oFluxes(iElementDOF,1);
  		Matrix oTempFluxes(iElementDOF,1);
 		// face fluxes
 		Vector oNormal;
  		vector< vector<double> > vvdGaussPointsCoordinates;
  		vector<double> vdGaussPointsWeights = m_poElementGeometry->GetFaceGaussPointsWeights();
  		double dReducedJacobian = 0.0;
  		double dFlux = 0.0;
  		Point oPoint;
  		unsigned int iGaussPointsCount = 0;
  		Matrix oNaturalCoordinates(3,1);
		unsigned int iFacesPerElement = m_poElementGeometry->GetFacesCount();
		unsigned int iNodesPerElement = m_poElementGeometry->GetNodesCount();
  		for(i = 1 ; i <= iFacesPerElement ; i++)
  		{
  			vvdGaussPointsCoordinates = m_poElementGeometry->GetFaceGaussPointsCoordinates(i);
  			iGaussPointsCount = (unsigned int)vvdGaussPointsCoordinates.size();
  			for(j = 0 ; j < iGaussPointsCount ; j++)
  			{
  				oNaturalCoordinates.Set(1,1,vvdGaussPointsCoordinates[j][0]);
  				oNaturalCoordinates.Set(2,1,vvdGaussPointsCoordinates[j][1]);
  				oNaturalCoordinates.Set(3,1,vvdGaussPointsCoordinates[j][2]);
   				oNormal = m_poElementGeometry->GetFaceOuterNormal(i,oNaturalCoordinates);
   				dReducedJacobian = oNormal.Length();
 				oPoint = m_poElementGeometry->GetPointOnFace(i,oNaturalCoordinates);
 				dFlux = m_vpoFacesTLoads[i - 1]->Get(oPoint,dTime);
  				oShapeFunctions = m_poElementGeometry->GetShapeFunctions(vvdGaussPointsCoordinates[j][0],vvdGaussPointsCoordinates[j][1],vvdGaussPointsCoordinates[j][2]);
  				for(k = 1 ; k <= iNodesPerElement ; k++)
  				{
  					oTempFluxes.Set(k,1,oShapeFunctions.Get(1,k)*dFlux);
  				}
  				oFluxes = oFluxes + oTempFluxes*dReducedJacobian*vdGaussPointsWeights[j];
  			}
  		}
		vector<FEMNode*>* pvpoNodes = m_poElementGeometry->GetNodes();
  		for(i = 1 ; i <= iNodesPerElement ; i++)
  		{
  			((FEMThermoMechanicalNode*)pvpoNodes->at(i - 1))->AddFlux(oFluxes.Get(i,1));
  		}
 		// body sources
 		oFluxes.SetSize(iElementDOF,1);
 		vvdGaussPointsCoordinates = m_poElementGeometry->GetBodyGaussPointsCoordinates();
 		vdGaussPointsWeights = m_poElementGeometry->GetBodyGaussPointsWeights();
 		iGaussPointsCount = (unsigned int)vvdGaussPointsCoordinates.size();
 		double dJacobian = 0.0;
 		for(i = 0 ; i < iGaussPointsCount ; i++)
 		{
			oNaturalCoordinates.Set(1,1,vvdGaussPointsCoordinates[i][0]);
			oNaturalCoordinates.Set(2,1,vvdGaussPointsCoordinates[i][1]);
			oNaturalCoordinates.Set(3,1,vvdGaussPointsCoordinates[i][2]);
			oPoint = m_poElementGeometry->GetPoint(oNaturalCoordinates);
			dFlux = m_poBodyTLoad->Get(oPoint,dTime);
			oShapeFunctions = m_poElementGeometry->GetShapeFunctions(vvdGaussPointsCoordinates[i][0],vvdGaussPointsCoordinates[i][1],vvdGaussPointsCoordinates[i][2]);
			dJacobian = m_poElementGeometry->GetJacobianMatrixDeterminant(vvdGaussPointsCoordinates[i][0],vvdGaussPointsCoordinates[i][1],vvdGaussPointsCoordinates[i][2]);
			for(j = 1 ; j <= iNodesPerElement ; j++)
			{
				oTempFluxes.Set(j,1,oShapeFunctions.Get(1,j)*dFlux);
			}
			oFluxes = oFluxes + oTempFluxes*dJacobian*vdGaussPointsWeights[i];
 		}
 		for(i = 1; i <= iNodesPerElement ; i++)
 		{
 			((FEMThermoMechanicalNode*)pvpoNodes->at(i - 1))->AddFlux(oFluxes.Get(i,1));
 		}
	}
	void FEMThermoMechanicalElement::ReadLoads(FILE* fpFile,vector<FEMLoad*>* pvpoLoads)
	{
	 	unsigned int iFacesPerElement = m_poElementGeometry->GetFacesCount();
 		vector<unsigned int> viFacesXLoads;
 		vector<unsigned int> viFacesYLoads;
 		vector<unsigned int> viFacesZLoads;
 		vector<unsigned int> viFacesTLoads;
 		viFacesXLoads.resize(iFacesPerElement);
 		viFacesYLoads.resize(iFacesPerElement);
 		viFacesZLoads.resize(iFacesPerElement);
 		viFacesTLoads.resize(iFacesPerElement);
 		
 		unsigned int iBodyXLoad = 0;
 		unsigned int iBodyYLoad = 0;
 		unsigned int iBodyZLoad = 0;
 		unsigned int iBodyTLoad = 0;
 		
 		string sRead = GetRealString(500,fpFile);
 		sscanf(sRead.c_str(),"%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",&viFacesXLoads[0],&viFacesYLoads[0],&viFacesZLoads[0],&viFacesTLoads[0],&viFacesXLoads[1],&viFacesYLoads[1],&viFacesZLoads[1],&viFacesTLoads[1],&viFacesXLoads[2],&viFacesYLoads[2],&viFacesZLoads[2],&viFacesTLoads[2],&viFacesXLoads[3],&viFacesYLoads[3],&viFacesZLoads[3],&viFacesTLoads[3],&viFacesXLoads[4],&viFacesYLoads[4],&viFacesZLoads[4],&viFacesTLoads[4],&viFacesXLoads[5],&viFacesYLoads[5],&viFacesZLoads[5],&viFacesTLoads[5],&iBodyXLoad,&iBodyYLoad,&iBodyZLoad,&iBodyTLoad);
 		
 		vector<FEMLoad*> vpoFaceXLoads;
 		vector<FEMLoad*> vpoFaceYLoads;
 		vector<FEMLoad*> vpoFaceZLoads;
 		vector<FEMLoad*> vpoFaceTLoads;
 		vpoFaceXLoads.resize(iFacesPerElement);
 		vpoFaceYLoads.resize(iFacesPerElement);
 		vpoFaceZLoads.resize(iFacesPerElement);
 		vpoFaceTLoads.resize(iFacesPerElement);
 		
 		unsigned int i = 0;
 		for(i = 0 ; i < iFacesPerElement ; i++)
 		{
 			vpoFaceXLoads[i] = pvpoLoads->at(viFacesXLoads[i] - 1);
 			vpoFaceYLoads[i] = pvpoLoads->at(viFacesYLoads[i] - 1);
 			vpoFaceZLoads[i] = pvpoLoads->at(viFacesZLoads[i] - 1);
 			vpoFaceTLoads[i] = pvpoLoads->at(viFacesTLoads[i] - 1);
 		}
 		SetFacesLoads(&vpoFaceXLoads,&vpoFaceYLoads,&vpoFaceZLoads,&vpoFaceTLoads);
 		SetBodyLoads(pvpoLoads->at(iBodyXLoad - 1),pvpoLoads->at(iBodyYLoad - 1),pvpoLoads->at(iBodyZLoad - 1),pvpoLoads->at(iBodyTLoad - 1));
	}
	void FEMThermoMechanicalElement::Write(FILE* fpFile) const
	{
		fprintf(fpFile,"%d,%d,%d\n",m_poElementGeometry->GetType(),GetType(),m_poMaterial->GetID());
		m_poElementGeometry->Write(fpFile);
		unsigned int i = 0;
		unsigned int iFacesCount = (unsigned int)m_vpoFacesXLoads.size();
		for(i = 0 ; i < iFacesCount ; i++)
		{
			fprintf(fpFile,"%d,%d,%d,%d,",m_vpoFacesXLoads[i]->GetID(),m_vpoFacesYLoads[i]->GetID(),m_vpoFacesZLoads[i]->GetID(),m_vpoFacesTLoads[i]->GetID());
		}
		fprintf(fpFile,"%d,%d,%d,%d\n",m_poBodyXLoad->GetID(),m_poBodyYLoad->GetID(),m_poBodyZLoad->GetID(),m_poBodyTLoad->GetID());
	}
	vector<unsigned int> FEMThermoMechanicalElement::GetDOFIndices() const
	{
	  	vector<unsigned int> viDOFIndices;
  		unsigned int i = 0;
  		vector<FEMNode*>* pvpoNodes = m_poElementGeometry->GetNodes();
  		unsigned int iNodesPerElement = (unsigned int)pvpoNodes->size();
  		viDOFIndices.reserve(ThermoMechanicalDOFPerNode*iNodesPerElement);
  		FEMNode* poNode = NULL;
  		for(i = 0 ; i < iNodesPerElement ; i++)
  		{
  			poNode = pvpoNodes->at(i);
  			viDOFIndices.push_back(((FEMThermoMechanicalNode*)poNode)->GetXDOF()->GetIndex());
  			viDOFIndices.push_back(((FEMThermoMechanicalNode*)poNode)->GetYDOF()->GetIndex());
  			viDOFIndices.push_back(((FEMThermoMechanicalNode*)poNode)->GetZDOF()->GetIndex());
  			viDOFIndices.push_back(((FEMThermoMechanicalNode*)poNode)->GetTDOF()->GetIndex());
  		}
  		return viDOFIndices;
	}
	void FEMThermoMechanicalElement::Initialize()
	{
		FEMSolidElement::Initialize();
		m_vpoFacesTLoads.clear();
		m_poBodyTLoad = NULL;
	}
	Matrix FEMThermoMechanicalElement::GetThermalGradientTransformation(const double& dXi,const double& dEta,const double& dZeta,double& dJacobian) const
	{
		unsigned int iElementDOF = GetThermalDegreesOfFreedomCount();
 		Matrix B(3,iElementDOF);
 		Matrix oDerivatives = m_poElementGeometry->GetShapeFunctionsXYZDerivatives(dXi,dEta,dZeta,dJacobian); 
 		unsigned int i = 0;
 		for(i = 1 ; i <= iElementDOF ; i++)
 		{
 			B.Set(1,i,oDerivatives.Get(1,i));
 			B.Set(2,i,oDerivatives.Get(2,i));
 			B.Set(3,i,oDerivatives.Get(3,i));
 		}
 		return B;
	}
	Matrix FEMThermoMechanicalElement::GetThermalShapeFunctions(const double& dXi,const double& dEta,const double& dZeta) const
	{
		unsigned int iElementDOF = GetThermalDegreesOfFreedomCount();
 		Matrix N(1,iElementDOF);
 		Matrix oShapeFunctions = m_poElementGeometry->GetShapeFunctions(dXi,dEta,dZeta);
 		unsigned int i = 0;
 		for(i = 1 ; i <= iElementDOF ; i++)
 		{
 			N.Set(1,i,oShapeFunctions.Get(1,i));
 		}
 		return N;
	}
	Matrix FEMThermoMechanicalElement::GetThermoMechanicalStiffnessMatrix()
	{
		Matrix oStiffness = GetStiffnessMatrix();
		Matrix oConductivity = GetConductionMatrix();
		Matrix oCoupling = GetStressTemperatureCouplingMatrix();
		unsigned int iElementDOF = GetDegreesOfFreedomCount();
		unsigned int iNodesCount = m_poElementGeometry->GetNodesCount();
		Matrix K(iElementDOF,iElementDOF);
		unsigned int i = 0;
		unsigned int j = 0;
		unsigned int iRowStartingIndex = 0;
		unsigned int iColumnStartingIndex = 0;
		unsigned int iStiffnessCurrentRowIndex = 0;
		unsigned int iStiffnessCurrentColumnIndex = 0;
		unsigned int iConductionCurrentRowIndex = 0;
		unsigned int iConductionCurrentColumnIndex = 0;
		for(i = 0 ; i < iNodesCount ; i++)
		{
			iRowStartingIndex = ThermoMechanicalDOFPerNode*i + 1;
			iStiffnessCurrentRowIndex = 3*i + 1;
			iConductionCurrentRowIndex = i + 1;
			for(j = 0 ; j < iNodesCount ; j++)
			{
				iColumnStartingIndex = ThermoMechanicalDOFPerNode*j + 1;
				iStiffnessCurrentColumnIndex = 3*j + 1;
				iConductionCurrentColumnIndex = j + 1;
				
				K.Set(iRowStartingIndex,iColumnStartingIndex,oStiffness.Get(iStiffnessCurrentRowIndex,iStiffnessCurrentColumnIndex));
				K.Set(iRowStartingIndex,iColumnStartingIndex + 1,oStiffness.Get(iStiffnessCurrentRowIndex,iStiffnessCurrentColumnIndex + 1));
				K.Set(iRowStartingIndex,iColumnStartingIndex + 2,oStiffness.Get(iStiffnessCurrentRowIndex,iStiffnessCurrentColumnIndex + 2));
				K.Set(iRowStartingIndex,iColumnStartingIndex + 3,-oCoupling.Get(iStiffnessCurrentRowIndex,iConductionCurrentColumnIndex));
				
				K.Set(iRowStartingIndex + 1,iColumnStartingIndex,oStiffness.Get(iStiffnessCurrentRowIndex + 1,iStiffnessCurrentColumnIndex));
				K.Set(iRowStartingIndex + 1,iColumnStartingIndex + 1,oStiffness.Get(iStiffnessCurrentRowIndex + 1,iStiffnessCurrentColumnIndex + 1));
				K.Set(iRowStartingIndex + 1,iColumnStartingIndex + 2,oStiffness.Get(iStiffnessCurrentRowIndex + 1,iStiffnessCurrentColumnIndex + 2));
				K.Set(iRowStartingIndex + 1,iColumnStartingIndex + 3,-oCoupling.Get(iStiffnessCurrentRowIndex + 1,iConductionCurrentColumnIndex));
				
				K.Set(iRowStartingIndex + 2,iColumnStartingIndex,oStiffness.Get(iStiffnessCurrentRowIndex + 2,iStiffnessCurrentColumnIndex));
				K.Set(iRowStartingIndex + 2,iColumnStartingIndex + 1,oStiffness.Get(iStiffnessCurrentRowIndex + 2,iStiffnessCurrentColumnIndex + 1));
				K.Set(iRowStartingIndex + 2,iColumnStartingIndex + 2,oStiffness.Get(iStiffnessCurrentRowIndex + 2,iStiffnessCurrentColumnIndex + 2));
				K.Set(iRowStartingIndex + 2,iColumnStartingIndex + 3,-oCoupling.Get(iStiffnessCurrentRowIndex + 2,iConductionCurrentColumnIndex));
				
				K.Set(iRowStartingIndex + 3,iColumnStartingIndex,0.0);
				K.Set(iRowStartingIndex + 3,iColumnStartingIndex + 1,0.0);
				K.Set(iRowStartingIndex + 3,iColumnStartingIndex + 2,0.0);
				K.Set(iRowStartingIndex + 3,iColumnStartingIndex + 3,oConductivity.Get(iConductionCurrentRowIndex,iConductionCurrentColumnIndex));
			}
		}
		return K;
	}
	Matrix FEMThermoMechanicalElement::GetThermoMechanicalDampingMatrix() const
	{
		unsigned int iElementDOF = GetDegreesOfFreedomCount();
		unsigned int iNodesCount = m_poElementGeometry->GetNodesCount();
		Matrix oCoupling = GetTemperatureStrainRateCouplingMatrix();
		Matrix oDamping = GetThermalDampingMatrix();
		Matrix C(iElementDOF,iElementDOF);
		unsigned int i = 0;
		unsigned int j = 0;
		unsigned int iRowStartingIndex = 0;
		unsigned int iColumnStartingIndex = 0;
		unsigned int iCouplingCurrentRowIndex = 0;
		unsigned int iCouplingCurrentColumnIndex = 0;
		unsigned int iDampingCurrentRowIndex = 0;
		unsigned int iDampingCurrentColumnIndex = 0;
		for(i = 0 ; i < iNodesCount ; i++)
		{
			iRowStartingIndex = ThermoMechanicalDOFPerNode*i + 1;
			iCouplingCurrentRowIndex = i + 1;
			iDampingCurrentRowIndex = i + 1;
			
			for(j = 0 ; j < iNodesCount ; j++)
			{
				iColumnStartingIndex = ThermoMechanicalDOFPerNode*j + 1;
				iCouplingCurrentColumnIndex = 3*j + 1;
				iDampingCurrentColumnIndex = j + 1;
				
				C.Set(iRowStartingIndex,iColumnStartingIndex,0.0);
				C.Set(iRowStartingIndex,iColumnStartingIndex + 1,0.0);
				C.Set(iRowStartingIndex,iColumnStartingIndex + 2,0.0);
				C.Set(iRowStartingIndex,iColumnStartingIndex + 3,0.0);
				
				C.Set(iRowStartingIndex + 1,iColumnStartingIndex,0.0);
				C.Set(iRowStartingIndex + 1,iColumnStartingIndex + 1,0.0);
				C.Set(iRowStartingIndex + 1,iColumnStartingIndex + 2,0.0);
				C.Set(iRowStartingIndex + 1,iColumnStartingIndex + 3,0.0);
				
				C.Set(iRowStartingIndex + 2,iColumnStartingIndex,0.0);
				C.Set(iRowStartingIndex + 2,iColumnStartingIndex + 1,0.0);
				C.Set(iRowStartingIndex + 2,iColumnStartingIndex + 2,0.0);
				C.Set(iRowStartingIndex + 2,iColumnStartingIndex + 3,0.0);
				
				C.Set(iRowStartingIndex + 3,iColumnStartingIndex,oCoupling.Get(iCouplingCurrentRowIndex,iCouplingCurrentColumnIndex));
				C.Set(iRowStartingIndex + 3,iColumnStartingIndex + 1,oCoupling.Get(iCouplingCurrentRowIndex,iCouplingCurrentColumnIndex + 1));
				C.Set(iRowStartingIndex + 3,iColumnStartingIndex + 2,oCoupling.Get(iCouplingCurrentRowIndex,iCouplingCurrentColumnIndex + 2));
				C.Set(iRowStartingIndex + 3,iColumnStartingIndex + 3,oDamping.Get(iDampingCurrentRowIndex,iDampingCurrentColumnIndex));
			}
		}
		return C;
	}
	Matrix FEMThermoMechanicalElement::GetThermoMechanicalMassMatrix() const
	{
		Matrix oMass = GetMassMatrix();
		unsigned int iElementDOF = GetDegreesOfFreedomCount();
		unsigned int iNodesCount = m_poElementGeometry->GetNodesCount();
		Matrix M(iElementDOF,iElementDOF);
		unsigned int i = 0;
		unsigned int j = 0;
		unsigned int iRowStartingIndex = 0;
		unsigned int iColumnStartingIndex = 0;
		unsigned int iMassCurrentRowIndex = 0;
		unsigned int iMassCurrentColumnIndex = 0;
		for(i = 0 ; i < iNodesCount ; i++)
		{
			iRowStartingIndex = ThermoMechanicalDOFPerNode*i + 1;
			iMassCurrentRowIndex = 3*i + 1;
			
			for(j = 0 ; j < iNodesCount ; j++)
			{
				iColumnStartingIndex = ThermoMechanicalDOFPerNode*j + 1;
				iMassCurrentColumnIndex = 3*j + 1;
				
				M.Set(iRowStartingIndex,iColumnStartingIndex,oMass.Get(iMassCurrentRowIndex,iMassCurrentColumnIndex));
				M.Set(iRowStartingIndex,iColumnStartingIndex + 1,oMass.Get(iMassCurrentRowIndex,iMassCurrentColumnIndex + 1));
				M.Set(iRowStartingIndex,iColumnStartingIndex + 2,oMass.Get(iMassCurrentRowIndex,iMassCurrentColumnIndex + 2));
				M.Set(iRowStartingIndex,iColumnStartingIndex + 3,0.0);
				
				M.Set(iRowStartingIndex + 1,iColumnStartingIndex,oMass.Get(iMassCurrentRowIndex + 1,iMassCurrentColumnIndex));
				M.Set(iRowStartingIndex + 1,iColumnStartingIndex + 1,oMass.Get(iMassCurrentRowIndex + 1,iMassCurrentColumnIndex + 1));
				M.Set(iRowStartingIndex + 1,iColumnStartingIndex + 2,oMass.Get(iMassCurrentRowIndex + 1,iMassCurrentColumnIndex + 2));
				M.Set(iRowStartingIndex + 1,iColumnStartingIndex + 3,0.0);
				
				M.Set(iRowStartingIndex + 2,iColumnStartingIndex,oMass.Get(iMassCurrentRowIndex + 2,iMassCurrentColumnIndex));
				M.Set(iRowStartingIndex + 2,iColumnStartingIndex + 1,oMass.Get(iMassCurrentRowIndex + 2,iMassCurrentColumnIndex + 1));
				M.Set(iRowStartingIndex + 2,iColumnStartingIndex + 2,oMass.Get(iMassCurrentRowIndex + 2,iMassCurrentColumnIndex + 2));
				M.Set(iRowStartingIndex + 2,iColumnStartingIndex + 3,0.0);
				
				M.Set(iRowStartingIndex + 3,iColumnStartingIndex,0.0);
				M.Set(iRowStartingIndex + 3,iColumnStartingIndex + 1,0.0);
				M.Set(iRowStartingIndex + 3,iColumnStartingIndex + 2,0.0);
				M.Set(iRowStartingIndex + 3,iColumnStartingIndex + 3,0.0);
			}
		}
		return M;
	}
}



