#include "FEMGaussPoint.h"
#include "FEMElement.h"
#include "FEMSolidElement.h"
#include "FEMSolidNode.h"

namespace FEMSystem
{
	FEMGaussPoint::FEMGaussPoint()
	{
		Initialize();
	}
	FEMGaussPoint::FEMGaussPoint(const FEMGaussPoint& oPoint)
	{
		*this = oPoint;
	}
	FEMGaussPoint::~FEMGaussPoint()
	{
		Reset();
	}
	FEMGaussPoint& FEMGaussPoint::operator=(const FEMGaussPoint& oPoint)
	{
		m_dXi = oPoint.m_dXi;
		m_dEta = oPoint.m_dEta;
		m_dZeta = oPoint.m_dZeta;
		m_oStress = oPoint.m_oStress;
		m_oStrainTransformation = oPoint.m_oStrainTransformation;
		m_oPlasticityState = oPoint.m_oPlasticityState;
		m_poElement = oPoint.m_poElement;
		m_dWeight = oPoint.m_dWeight;
		m_dJacobian = oPoint.m_dJacobian;
		return *this;
	}
	void FEMGaussPoint::Reset()
	{
		Initialize();
	}
	void FEMGaussPoint::SetCoordinates(const double& dXi,const double& dEta,const double& dZeta)
	{
		m_dXi = dXi;
		m_dEta = dEta;
		m_dZeta = dZeta;
	}
	void FEMGaussPoint::Update()
	{
	 	Matrix oShapeFunctions = m_poElement->GetMaterial()->UpdateGaussPoint(this);
	}
	void FEMGaussPoint::SetElement(FEMElement* poElement)
	{
		m_poElement = poElement;
		if((m_poElement->GetType() == SolidFEMElement) || (m_poElement->GetType() == ThermoMechanicalFEMElement))
		{
			m_oStrainTransformation = ((FEMSolidElement*)m_poElement)->GetStrainTransformation(m_dXi,m_dEta,m_dZeta,m_dJacobian);
		}
	}
	void FEMGaussPoint::SetWeight(const double& dWeight)
	{
		m_dWeight = dWeight;
	}
	Matrix FEMGaussPoint::GetMaterialStiffnessMatrix()
	{
		return m_poElement->GetMaterial()->UpdateGaussPoint(this);
	}
	Matrix FEMGaussPoint::GetStiffnessMatrix()
	{
		Matrix oC = GetMaterialStiffnessMatrix();
		Matrix oK = m_oStrainTransformation.GetTranspose()*oC*m_oStrainTransformation;
		oK = oK*m_dJacobian*m_dWeight;
		return oK;
	}
	Matrix FEMGaussPoint::GetStress() const
	{
		return m_oStress;
	}
	Matrix FEMGaussPoint::GetTotalStrain() const
	{
		unsigned int iStrainsCount = 6;
		unsigned int iElementDOF = ((FEMSolidElement*)m_poElement)->GetSolidDegreesOfFreedomCount();
 		Matrix oDisplacements(iElementDOF,1);
 		unsigned int i = 0;
		vector<FEMNode*>* pvpoNodes = m_poElement->GetGeometry()->GetNodes();
		unsigned int iNodesPerElement = (unsigned int)pvpoNodes->size();
 		for(i = 0; i < iNodesPerElement ; i++)
 		{
 			oDisplacements.Set(3*i + 1,1,((FEMSolidNode*)pvpoNodes->at(i))->GetXDisplacement());
 			oDisplacements.Set(3*i + 2,1,((FEMSolidNode*)pvpoNodes->at(i))->GetYDisplacement());
 			oDisplacements.Set(3*i + 3,1,((FEMSolidNode*)pvpoNodes->at(i))->GetZDisplacement());
 		}
 		return m_oStrainTransformation*oDisplacements;
	}
	void FEMGaussPoint::SetStress(const Matrix& oStress)
	{
		m_oStress = oStress;
	}
	FEMPlasticityState* FEMGaussPoint::GetPlasticityState()
	{
		return &m_oPlasticityState;
	}
	Matrix FEMGaussPoint::GetVectorizedStress() const
	{
		Matrix oResult(6,1);
		oResult.Set(1,1,m_oStress.Get(1,1));
		oResult.Set(2,1,m_oStress.Get(2,2));
		oResult.Set(3,1,m_oStress.Get(3,3));
		oResult.Set(4,1,m_oStress.Get(1,2));
		oResult.Set(5,1,m_oStress.Get(2,3));
		oResult.Set(6,1,m_oStress.Get(3,1));
		return oResult;
	}
	void FEMGaussPoint::Initialize()
	{
		m_dXi = 0.0;
		m_dEta = 0.0;
		m_dZeta = 0.0;
		m_oStress.SetSize(3,3);
		m_oStrainTransformation.Reset();
		m_oPlasticityState.Reset();
		m_poElement = NULL;
		m_dWeight = 0.0;
		m_dJacobian = 0.0;
	}
	void FEMGaussPoint::SetVectorizedStress(const Matrix& oStress)
	{
		m_oStress.Set(1,1,oStress.Get(1,1));
		m_oStress.Set(1,2,oStress.Get(4,1));
		m_oStress.Set(1,3,oStress.Get(6,1));
		
		m_oStress.Set(2,1,oStress.Get(4,1));
		m_oStress.Set(2,2,oStress.Get(2,1));
		m_oStress.Set(2,3,oStress.Get(5,1));
		
		m_oStress.Set(3,1,oStress.Get(6,1));
		m_oStress.Set(3,2,oStress.Get(5,1));
		m_oStress.Set(3,3,oStress.Get(3,1));
	}
}

