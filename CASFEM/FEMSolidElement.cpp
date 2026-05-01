// Ahmed M. Hussein

#include "FEMSolidElement.h"
#include "FEMSolidNode.h"
#include "MathServices.h"
#include "Tools.h"
#include "cmath"

using namespace SupportSystem;

namespace FEMSystem
{
	unsigned int FEMSolidElement::SolidDOFPerNode = 3;
	FEMSolidElement::FEMSolidElement()
	{
		Initialize();
	}
	FEMSolidElement::FEMSolidElement(const FEMSolidElement& oElement)
	{
		*this = oElement;
	}
	FEMSolidElement::~FEMSolidElement()
	{
		Reset();
	}
	FEMSolidElement& FEMSolidElement::operator=(const FEMSolidElement& oElement)
	{
		FEMElement::operator=(oElement);
		m_vpoFacesXLoads = oElement.m_vpoFacesXLoads;
		m_vpoFacesYLoads = oElement.m_vpoFacesYLoads;
		m_vpoFacesZLoads = oElement.m_vpoFacesZLoads;
		
		m_poBodyXLoad = oElement.m_poBodyXLoad;
		m_poBodyYLoad = oElement.m_poBodyYLoad;
		m_poBodyZLoad = oElement.m_poBodyZLoad;
		return *this;
	}
	void FEMSolidElement::Reset()
	{
		FEMElement::Reset();
		Initialize();
	}
	void FEMSolidElement::SetFacesLoads(vector<FEMLoad*>* pvpoXLoads,vector<FEMLoad*>* pvpoYLoads,vector<FEMLoad*>* pvpoZLoads)
	{
		unsigned int iFacesCount = m_poElementGeometry->GetFacesCount();
		if((pvpoXLoads->size() != iFacesCount) || (pvpoYLoads->size() != iFacesCount) || (pvpoZLoads->size() != iFacesCount))
		{
			return;
		}
		unsigned int i = 0;
		for(i = 0 ; i < iFacesCount ; i++)
		{
			m_vpoFacesXLoads[i] = pvpoXLoads->at(i);
			m_vpoFacesYLoads[i] = pvpoYLoads->at(i);
			m_vpoFacesZLoads[i] = pvpoZLoads->at(i);
		}
	}
	void FEMSolidElement::SetBodyLoads(FEMLoad* poXLoad,FEMLoad* poYLoad,FEMLoad* poZLoad)
	{
		m_poBodyXLoad = poXLoad;
		m_poBodyYLoad = poYLoad;
		m_poBodyZLoad = poZLoad;
	}
	void FEMSolidElement::Initialize()
	{
		FEMElement::Initialize();
		m_vpoFacesXLoads.clear();
		m_vpoFacesYLoads.clear();
		m_vpoFacesZLoads.clear();
		m_poBodyXLoad = NULL;
		m_poBodyYLoad = NULL;
		m_poBodyZLoad = NULL;
	}
	void FEMSolidElement::SetFacesLoadsSize()
	{
		unsigned int iFacesPerElement = m_poElementGeometry->GetFacesCount();
 		m_vpoFacesXLoads.resize(iFacesPerElement);
 		m_vpoFacesYLoads.resize(iFacesPerElement);
 		m_vpoFacesZLoads.resize(iFacesPerElement);
		unsigned int i = 0;
 		for(i = 0; i < iFacesPerElement ; i++)
 		{
 			m_vpoFacesXLoads[i] = NULL;
 			m_vpoFacesYLoads[i] = NULL;
 			m_vpoFacesZLoads[i] = NULL;
 		}
	}
	FEMElementType FEMSolidElement::GetType() const
	{
		return SolidFEMElement;
	}
	unsigned int FEMSolidElement::GetDegreesOfFreedomCount() const
	{
		return GetSolidDegreesOfFreedomCount();
	}
	unsigned int FEMSolidElement::GetSolidDegreesOfFreedomCount() const
	{
		return (m_poElementGeometry->GetNodesCount()*SolidDOFPerNode);
	}
	void FEMSolidElement::ApplyLoads(const double& dTime)
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
  			((FEMSolidNode*)pvpoNodes->at(i))->AddForce(oNodalForce);
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
 			((FEMSolidNode*)pvpoNodes->at(i))->AddForce(oNodalForce);
 		}
	}
	void FEMSolidElement::ReadLoads(FILE* fpFile,vector<FEMLoad*>* pvpoLoads)
	{
 		unsigned int iFacesPerElement = m_poElementGeometry->GetFacesCount();
 		vector<unsigned int> viFacesXLoads;
 		vector<unsigned int> viFacesYLoads;
 		vector<unsigned int> viFacesZLoads;
 		viFacesXLoads.resize(iFacesPerElement);
 		viFacesYLoads.resize(iFacesPerElement);
 		viFacesZLoads.resize(iFacesPerElement);
 		
 		unsigned int iBodyXLoad = 0;
 		unsigned int iBodyYLoad = 0;
 		unsigned int iBodyZLoad = 0;
 		
 		string sRead = GetRealString(500,fpFile);;
 		sscanf(sRead.c_str(),"%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",&viFacesXLoads[0],&viFacesYLoads[0],&viFacesZLoads[0],&viFacesXLoads[1],&viFacesYLoads[1],&viFacesZLoads[1],&viFacesXLoads[2],&viFacesYLoads[2],&viFacesZLoads[2],&viFacesXLoads[3],&viFacesYLoads[3],&viFacesZLoads[3],&viFacesXLoads[4],&viFacesYLoads[4],&viFacesZLoads[4],&viFacesXLoads[5],&viFacesYLoads[5],&viFacesZLoads[5],&iBodyXLoad,&iBodyYLoad,&iBodyZLoad);
 		
 		vector<FEMLoad*> vpoFaceXLoads;
 		vector<FEMLoad*> vpoFaceYLoads;
 		vector<FEMLoad*> vpoFaceZLoads;
 		vpoFaceXLoads.resize(iFacesPerElement);
 		vpoFaceYLoads.resize(iFacesPerElement);
 		vpoFaceZLoads.resize(iFacesPerElement);
 		
 		unsigned int i = 0;
 		for(i = 0 ; i < iFacesPerElement ; i++)
 		{
 			vpoFaceXLoads[i] = pvpoLoads->at(viFacesXLoads[i] - 1);
 			vpoFaceYLoads[i] = pvpoLoads->at(viFacesYLoads[i] - 1);
 			vpoFaceZLoads[i] = pvpoLoads->at(viFacesZLoads[i] - 1);
 		}
 		SetFacesLoads(&vpoFaceXLoads,&vpoFaceYLoads,&vpoFaceZLoads);
 		SetBodyLoads(pvpoLoads->at(iBodyXLoad - 1),pvpoLoads->at(iBodyYLoad - 1),pvpoLoads->at(iBodyZLoad - 1));
	}
	void FEMSolidElement::Write(FILE* fpFile) const
	{
		fprintf(fpFile,"%d,%d,%d\n",m_poElementGeometry->GetType(),GetType(),m_poMaterial->GetID());
		m_poElementGeometry->Write(fpFile);
		unsigned int i = 0;
		unsigned int iFacesCount = (unsigned int)m_vpoFacesXLoads.size();
		for(i = 0 ; i < iFacesCount ; i++)
		{
			fprintf(fpFile,"%d,%d,%d,",m_vpoFacesXLoads[i]->GetID(),m_vpoFacesYLoads[i]->GetID(),m_vpoFacesZLoads[i]->GetID());
		}
		fprintf(fpFile,"%d,%d,%d\n",m_poBodyXLoad->GetID(),m_poBodyYLoad->GetID(),m_poBodyZLoad->GetID());
	}
	Matrix FEMSolidElement::GetStiffnessMatrix()
	{
		unsigned int iElementDOF = GetSolidDegreesOfFreedomCount();
 		Matrix K(iElementDOF,iElementDOF);
 		unsigned int i = 0;
 		unsigned int iGaussPointsCount = (unsigned int)m_voGaussPoints.size();
 		for(i = 0; i < iGaussPointsCount ; i++)
 		{
			K = K + m_voGaussPoints[i].GetStiffnessMatrix();
 		}
 		return K;
	}
 	Matrix FEMSolidElement::GetMassMatrix() const
 	{
 		unsigned int iElementDOF = GetSolidDegreesOfFreedomCount();
  		Matrix M(iElementDOF,iElementDOF);
 		vector< vector<double> > vvdGaussPointsCoordinates = m_poElementGeometry->GetBodyGaussPointsCoordinates();
 		vector<double> vdGaussPointsWeights = m_poElementGeometry->GetBodyGaussPointsWeights();
  		unsigned int i = 0;
  		unsigned int iGaussPointsCount = (unsigned int)vvdGaussPointsCoordinates.size();
  		double dMaterialDensity = m_poMaterial->GetMassDensity();
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
			N = GetInertiaTransformation(dXi,dEta,dZeta,dJacobian);
			dTemp = dJacobian*vdGaussPointsWeights[i];
			M = M + (N.GetTranspose()*N)*(dMaterialDensity*dTemp);
  		}
  		return M;
 	}
 	Matrix FEMSolidElement::GetLumpedMassMatrix() const
 	{
 		Matrix oMC = GetMassMatrix();
 		unsigned int i = 0;
 		unsigned int iNodesCount = m_poElementGeometry->GetNodesCount();
 		double dVertexNodesTotalMass = 0.0;
 		double dMidSideNodesTotalMass = 0.0;
 		for(i = 1 ; i <= iNodesCount ; i++)
 		{
 			if(m_poElementGeometry->IsMidSideNode(i))
 			{
				dMidSideNodesTotalMass = dMidSideNodesTotalMass + oMC.Get(i,i);
 			}
 			else
 			{
 				dVertexNodesTotalMass = dVertexNodesTotalMass + oMC.Get(i,i);
 			}
 		}
 		double dMassMatrixSum = oMC.SumAllEntries();
 		// lump the mass matrix as described in Dhondt G. The finite element mathod for three 
 		// dimensional thermomechanical applications, page 141. The mass lumping factor is the
 		// reciprocal of the one introduced by Dhondt though. 
 		double dBeta = m_poElementGeometry->GetMassLumpingFactor();
 		double dVertexNodesScalingFactor = dMassMatrixSum/dVertexNodesTotalMass*(1.0/(1.0 + dBeta));
 		double dMidSideNodesScalingFactor = dMassMatrixSum/dMidSideNodesTotalMass*(dBeta/(1.0 + dBeta));
 		// in case of no midside nodes, the dMidSideNodesScalingFactor will be infinity, but this
 		// is not a problem because it won't be used anyway
 		unsigned int iRowsCount = iNodesCount*SolidDOFPerNode;
 		Matrix oML(iNodesCount*SolidDOFPerNode,1);
 		double dEntry = 0.0;
 		unsigned int j = 0;
 		for(i = 1 ; i <= iNodesCount ; i++)
 		{
 			if(m_poElementGeometry->IsMidSideNode(i))
 			{
				dEntry = oMC.Get(i,i)*dMidSideNodesScalingFactor;
 			}
 			else
 			{
 				dEntry = oMC.Get(i,i)*dVertexNodesScalingFactor;
 			}
 			for(j = 0 ; j < SolidDOFPerNode ; j++)
 			{
				oML.Set(SolidDOFPerNode*i - j,1,dEntry);
 			}
 		}
 		return oML;
 	}
	Matrix FEMSolidElement::GetNodalStresses() const
	{
		unsigned int iStressesCount = 6;
 		unsigned int i = 0;
 		unsigned int iGaussPointsCount = (unsigned int)m_voGaussPoints.size();
 		Matrix oTempStresses;
 		Matrix oStressVector(6,1);
 		unsigned int j = 0;
 		unsigned int iCount = 0;
 		Matrix oGaussPointsStresses(iGaussPointsCount,iStressesCount);
 		for(i = 0 ; i < iGaussPointsCount ; i++)
 		{
 			oTempStresses = m_voGaussPoints[i].GetStress();
 			oStressVector.Set(1,1,oTempStresses.Get(1,1));
 			oStressVector.Set(2,1,oTempStresses.Get(2,2));
 			oStressVector.Set(3,1,oTempStresses.Get(3,3));
 			oStressVector.Set(4,1,oTempStresses.Get(1,2));
 			oStressVector.Set(5,1,oTempStresses.Get(2,3));
 			oStressVector.Set(6,1,oTempStresses.Get(3,1));
 			iCount = iCount + 1;
			for(j = 1; j <= iStressesCount; j++)
			{
				oGaussPointsStresses.Set(iCount,j,oStressVector.Get(j,1));
			}
 		}
 		Matrix oStress = m_poElementGeometry->GetShapeFunctionsExtrapolationsMatrix()*oGaussPointsStresses;
 		return oStress;
	}
	Vector FEMSolidElement::GetDisplacement(const Matrix& oNaturalCoordinates) const
	{
 		Vector oDisplacement(0.0,0.0,0.0);
		vector<FEMNode*>* pvpoNodes = m_poElementGeometry->GetNodes();
		unsigned int iNodesPerElement = m_poElementGeometry->GetNodesCount();
 		if((oNaturalCoordinates.GetRowsCount() != 3) || (oNaturalCoordinates.GetColumnsCount() != 1))
 		{
 			return oDisplacement;
 		}
 		Matrix oShapeFunctions = m_poElementGeometry->GetShapeFunctions(oNaturalCoordinates.Get(1,1),oNaturalCoordinates.Get(2,1),oNaturalCoordinates.Get(3,1));
 		unsigned int i = 0;
 		for(i = 0; i < iNodesPerElement ; i++)
 		{
 			oDisplacement = oDisplacement + ((FEMSolidNode*)pvpoNodes->at(i))->GetDisplacement()*oShapeFunctions.Get(1,i + 1);
 		}
 		return oDisplacement;
	}
	Matrix FEMSolidElement::GetStress(const Matrix& oNaturalCoordinates) const
	{
 		Matrix oStress(3,3);
		vector<FEMNode*>* pvpoNodes = m_poElementGeometry->GetNodes();
		unsigned int iNodesPerElement = m_poElementGeometry->GetNodesCount();
 		if((oNaturalCoordinates.GetRowsCount() != 3) || (oNaturalCoordinates.GetColumnsCount() != 1))
 		{
 			return oStress;
 		}
 		Matrix oShapeFunctions = m_poElementGeometry->GetShapeFunctions(oNaturalCoordinates.Get(1,1),oNaturalCoordinates.Get(2,1),oNaturalCoordinates.Get(3,1));
 		unsigned int i = 0;
 		for(i = 0; i < iNodesPerElement ; i++)
 		{
 			oStress = oStress + (((FEMSolidNode*)pvpoNodes->at(i))->GetStresses())*oShapeFunctions.Get(1,i + 1);
 		}
 		return oStress;
	}
 	Matrix FEMSolidElement::GetStrainTransformation(const double& dXi,const double& dEta,const double& dZeta,double& dJacobian) const
 	{
		unsigned int iElementDOF = GetSolidDegreesOfFreedomCount();
 		Matrix B(6,iElementDOF);
 		Matrix oDerivatives = m_poElementGeometry->GetShapeFunctionsXYZDerivatives(dXi,dEta,dZeta,dJacobian);
 
 		unsigned int i = 0;
 		double dNx = 0.0;
 		double dNy = 0.0;
 		double dNz = 0.0;
		unsigned int iNodesPerElement = m_poElementGeometry->GetNodesCount();
 		for(i = 0 ; i < iNodesPerElement ; i++)
 		{
 			dNx = oDerivatives.Get(1,i + 1);
 			dNy = oDerivatives.Get(2,i + 1);
 			dNz = oDerivatives.Get(3,i + 1);
 
 			B.Set(1,3*i + 1,dNx);
 			B.Set(2,3*i + 2,dNy);
 			B.Set(3,3*i + 3,dNz);
 
 			B.Set(4,3*i + 1,dNy);
 			B.Set(4,3*i + 2,dNx);
 
 			B.Set(5,3*i + 2,dNz);
 			B.Set(5,3*i + 3,dNy);
 
 			B.Set(6,3*i + 1,dNz);
 			B.Set(6,3*i + 3,dNx);
 		}
 		return B;
 	}
 	Matrix FEMSolidElement::GetInertiaTransformation(const double& dXi,const double& dEta,const double& dZeta,double& dJacobian) const
 	{
		unsigned int iElementDOF = GetSolidDegreesOfFreedomCount();
 		Matrix N(3,iElementDOF);
 		Matrix oShapeFunctions = m_poElementGeometry->GetShapeFunctions(dXi,dEta,dZeta);
 
 		unsigned int i = 0;
 		double dN = 0.0;
		unsigned int iNodesPerElement = m_poElementGeometry->GetNodesCount();
 		for(i = 0 ; i < iNodesPerElement ; i++)
 		{
 			dN = oShapeFunctions.Get(1,i + 1);
 
 			N.Set(1,3*i + 1,dN);
 			N.Set(2,3*i + 2,dN);
 			N.Set(3,3*i + 3,dN);
 		}
		dJacobian = m_poElementGeometry->GetJacobianMatrixDeterminant(dXi,dEta,dZeta);
  		return N;
  	}
  	vector<unsigned int> FEMSolidElement::GetDOFIndices() const
  	{
  		vector<unsigned int> viDOFIndices;
  		unsigned int i = 0;
  		vector<FEMNode*>* pvpoNodes = m_poElementGeometry->GetNodes();
  		unsigned int iNodesPerElement = (unsigned int)pvpoNodes->size();
  		viDOFIndices.reserve(SolidDOFPerNode*iNodesPerElement);
  		FEMNode* poNode = NULL;
  		for(i = 0 ; i < iNodesPerElement ; i++)
  		{
  			poNode = pvpoNodes->at(i);
  			viDOFIndices.push_back(((FEMSolidNode*)poNode)->GetXDOF()->GetIndex());
  			viDOFIndices.push_back(((FEMSolidNode*)poNode)->GetYDOF()->GetIndex());
  			viDOFIndices.push_back(((FEMSolidNode*)poNode)->GetZDOF()->GetIndex());
  		}
  		return viDOFIndices;
  	}
}




