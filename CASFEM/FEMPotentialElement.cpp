// Ahmed M. Hussein

#include "FEMPotentialElement.h"
#include "FEMPotentialNode.h"
#include "MathServices.h"
#include "Tools.h"
#include "cmath"

using namespace SupportSystem;

namespace FEMSystem
{
	unsigned int FEMPotentialElement::PotentialDOFPerNode = 1;
	FEMPotentialElement::FEMPotentialElement()
	{
		Initialize();
	}
	FEMPotentialElement::FEMPotentialElement(const FEMPotentialElement& oElement)
	{
		*this = oElement;
	}
	FEMPotentialElement::~FEMPotentialElement()
	{
		Reset();
	}
	FEMPotentialElement& FEMPotentialElement::operator=(const FEMPotentialElement& oElement)
	{
		FEMElement::operator=(oElement);
		m_vpoFacesLoads = oElement.m_vpoFacesLoads;
		m_poBodyLoad = oElement.m_poBodyLoad;
		return *this;
	}
	void FEMPotentialElement::Reset()
	{
		FEMElement::Reset();
		Initialize();
	}
	void FEMPotentialElement::SetFacesLoads(vector<FEMLoad*>* pvpoLoads)
	{
		unsigned int iFacesCount = m_poElementGeometry->GetFacesCount();
		if(pvpoLoads->size() != iFacesCount)
		{
			return;
		}
		unsigned int i = 0;
		for(i = 0 ; i < iFacesCount ; i++)
		{
			m_vpoFacesLoads[i] = pvpoLoads->at(i);
		}
	}
	void FEMPotentialElement::SetBodyLoad(FEMLoad* poLoad)
	{
		m_poBodyLoad = poLoad;
	}
	void FEMPotentialElement::Initialize()
	{
		FEMElement::Initialize();
		m_vpoFacesLoads.clear();
		m_poBodyLoad = NULL;
	}
	void FEMPotentialElement::SetFacesLoadsSize()
	{
		unsigned int iFacesPerElement = m_poElementGeometry->GetFacesCount();
 		m_vpoFacesLoads.resize(iFacesPerElement);
		unsigned int i = 0;
 		for(i = 0; i < iFacesPerElement ; i++)
 		{
 			m_vpoFacesLoads[i] = NULL;
 		}
	}
	FEMElementType FEMPotentialElement::GetType() const
	{
		return PotentialFEMElement;
	}
	unsigned int FEMPotentialElement::GetDegreesOfFreedomCount() const
	{
		return (m_poElementGeometry->GetNodesCount()*PotentialDOFPerNode);
	}
	void FEMPotentialElement::ApplyLoads(const double& dTime)
	{
 		Matrix oShapeFunctions;
  		unsigned int i = 0;
  		unsigned int j = 0;
  		unsigned int k = 0;
		unsigned int iElementDOF = GetDegreesOfFreedomCount();
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
 				dFlux = m_vpoFacesLoads[i - 1]->Get(oPoint,dTime);
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
  			((FEMPotentialNode*)pvpoNodes->at(i - 1))->AddFlux(oFluxes.Get(i,1));
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
			dFlux = m_poBodyLoad->Get(oPoint,dTime);
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
 			((FEMPotentialNode*)pvpoNodes->at(i - 1))->AddFlux(oFluxes.Get(i,1));
 		}
	}
	void FEMPotentialElement::ReadLoads(FILE* fpFile,vector<FEMLoad*>* pvpoLoads)
	{
 		unsigned int iFacesPerElement = m_poElementGeometry->GetFacesCount();
 		vector<unsigned int> viFacesLoads;
 		viFacesLoads.resize(iFacesPerElement);
 		unsigned int iBodyLoad = 0;	
 		string sRead = GetRealString(500,fpFile);
 		sscanf(sRead.c_str(),"%d,%d,%d,%d,%d,%d,%d\n",&viFacesLoads[0],&viFacesLoads[1],&viFacesLoads[2],&viFacesLoads[3],&viFacesLoads[4],&viFacesLoads[5],&iBodyLoad);
 		
 		vector<FEMLoad*> vpoFaceLoads;
 		vpoFaceLoads.resize(iFacesPerElement);
 		unsigned int i = 0;
 		for(i = 0 ; i < iFacesPerElement ; i++)
 		{
 			vpoFaceLoads[i] = pvpoLoads->at(viFacesLoads[i] - 1);
 		}
 		SetFacesLoads(&vpoFaceLoads);
 		SetBodyLoad(pvpoLoads->at(iBodyLoad - 1));
	}
	void FEMPotentialElement::Write(FILE* fpFile) const
	{
		fprintf(fpFile,"%d,%d,%d\n",m_poElementGeometry->GetType(),GetType(),m_poMaterial->GetID());
		m_poElementGeometry->Write(fpFile);
		unsigned int i = 0;
		unsigned int iFacesCount = (unsigned int)m_vpoFacesLoads.size();
		for(i = 0 ; i < iFacesCount ; i++)
		{
			fprintf(fpFile,"%d,",m_vpoFacesLoads[i]->GetID());
		}
		fprintf(fpFile,"%d\n",m_poBodyLoad->GetID());
	}
	Matrix FEMPotentialElement::GetConductionMatrix() const
	{
		unsigned int iElementDOF = GetDegreesOfFreedomCount();
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
 			B = GetGradientTransformation(dXi,dEta,dZeta,dJacobian);
 			dTemp = dJacobian*vdGaussPointsWeights[i];
 			K = K + (B.GetTranspose()*oMaterialConductionMatrix*B)*(dTemp);
 		}
 		return K;
	}
	Matrix FEMPotentialElement::GetNodalFluxes() const
	{
 		unsigned int iFluxesCount = 3;
		unsigned int iElementDOF = GetDegreesOfFreedomCount();
		unsigned int iNodesPerElement = m_poElementGeometry->GetNodesCount();
 		Matrix U(iElementDOF,1);
 		unsigned int i = 0;
		vector<FEMNode*>* pvpoNodes = m_poElementGeometry->GetNodes();
 		for(i = 1; i <= iNodesPerElement ; i++)
 		{
 			U.Set(i,1,((FEMPotentialNode*)pvpoNodes->at(i - 1))->GetPotential());
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
			B = GetGradientTransformation(dXi,dEta,dZeta,dDummyJacobian);
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
	double FEMPotentialElement::GetPotential(const Matrix& oNaturalCoordinates) const
	{
 		double dPotential = 0.0;
		vector<FEMNode*>* pvpoNodes = m_poElementGeometry->GetNodes();
		unsigned int iNodesPerElement = m_poElementGeometry->GetNodesCount();
 		if((oNaturalCoordinates.GetRowsCount() != 3) || (oNaturalCoordinates.GetColumnsCount() != 1))
 		{
 			return dPotential;
 		}
 		Matrix oShapeFunctions = m_poElementGeometry->GetShapeFunctions(oNaturalCoordinates.Get(1,1),oNaturalCoordinates.Get(2,1),oNaturalCoordinates.Get(3,1));
 		unsigned int i = 0;
 		for(i = 0; i < iNodesPerElement ; i++)
 		{
 			dPotential = dPotential + ((FEMPotentialNode*)pvpoNodes->at(i))->GetPotential()*oShapeFunctions.Get(1,i + 1);
 		}
 		return dPotential;
	}
	Vector FEMPotentialElement::GetFlux(const Matrix& oNaturalCoordinates) const
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
 			oFlux = oFlux + (((FEMPotentialNode*)pvpoNodes->at(i))->GetFluxes())*oShapeFunctions.Get(1,i + 1);
 		}
 		return oFlux;
	}
 	Matrix FEMPotentialElement::GetGradientTransformation(const double& dXi,const double& dEta,const double& dZeta,double& dJacobian) const
 	{
		unsigned int iElementDOF = GetDegreesOfFreedomCount();
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
  	vector<unsigned int> FEMPotentialElement::GetDOFIndices() const
  	{
  		vector<unsigned int> viDOFIndices;
  		unsigned int i = 0;
  		vector<FEMNode*>* pvpoNodes = m_poElementGeometry->GetNodes();
  		unsigned int iNodesPerElement = (unsigned int)pvpoNodes->size();
  		viDOFIndices.reserve(PotentialDOFPerNode*iNodesPerElement);
  		FEMNode* poNode = NULL;
  		for(i = 0 ; i < iNodesPerElement ; i++)
  		{
  			poNode = pvpoNodes->at(i);
  			viDOFIndices.push_back(((FEMPotentialNode*)poNode)->GetPotentialDOF()->GetIndex());
  		}
  		return viDOFIndices;
  	}
}




