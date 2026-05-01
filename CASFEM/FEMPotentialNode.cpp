// Ahmed M. Hussein

#include "FEMPotentialNode.h"
#include "Tools.h"

using namespace SupportSystem;

namespace FEMSystem
{
	unsigned int FEMPotentialNode::PotentialDOFCount = 1;
	FEMPotentialNode::FEMPotentialNode()
	{
		Initialize();
	}
	FEMPotentialNode::~FEMPotentialNode()
	{
		Reset();
	}
	FEMPotentialNode::FEMPotentialNode(const FEMPotentialNode& oNode)
	{
		*this = oNode;
	}
	FEMPotentialNode::FEMPotentialNode(const double& dX,const double& dY,const double& dZ)
	{
		Initialize();
		Set(dX,dY,dZ);
	}
	FEMPotentialNode& FEMPotentialNode::operator=(const FEMPotentialNode& oNode)
	{
        FEMNode::operator =(oNode);
        m_oPotentialDOF = oNode.m_oPotentialDOF;
		m_poLoad = oNode.m_poLoad;
		m_oFluxes = oNode.m_oFluxes;
		return *this;
	}
	unsigned int FEMPotentialNode::GetDOFCount() const
	{
		return PotentialDOFCount;
	}
	void FEMPotentialNode::Reset()
	{
		FEMNode::Reset();
		m_oPotentialDOF.Reset();
		m_oFluxes.Set(0.0,0.0,0.0);
		m_poLoad = NULL;
	}
	FEMDegreeOfFreedom* FEMPotentialNode::GetPotentialDOF()
	{
		return &m_oPotentialDOF;
	}
	void FEMPotentialNode::ResetLoads()
	{
		m_oPotentialDOF.SetPrimaryValue(0.0);
		m_oPotentialDOF.SetSecondaryValue(0.0);
	}
	FEMNode* FEMPotentialNode::Clone() const
	{
		return new FEMPotentialNode(*this);
	}
	void FEMPotentialNode::ApplyLoads(const double& dTime)
	{
		ResetLoads();
		if(m_poLoad != NULL)
		{
			m_oPotentialDOF.AddToConstraint(m_poLoad->Get(*this,dTime));
		}
	}
	void FEMPotentialNode::AddFlux(const double& dFlux)
	{
		m_oPotentialDOF.AddToSecondaryValue(dFlux);
	}
	void FEMPotentialNode::Initialize()
	{
		FEMNode::Initialize();
		m_oPotentialDOF.Reset();
		ResetLoads();
		m_oFluxes.Set(0.0,0.0,0.0);
		m_poLoad = NULL;
	}
	void FEMPotentialNode::SetLoad(FEMLoad* poLoad)
	{
		m_poLoad = poLoad;
	}
	FEMLoad* FEMPotentialNode::GetLoad() const
	{
		return m_poLoad;
	}
	FEMNodeType FEMPotentialNode::GetType() const
	{
		return PotentialFEMNode;
	}
	void FEMPotentialNode::ReadNode(FILE* fpFile,vector<FEMLoad*>* pvpoLoads)
	{
		unsigned int iTemp1 = 0;
		unsigned int iTemp2 = 0;
		unsigned int iSurface = 0;
        double dTemp = 0.0;
        
        // read the node position, condition and whether it is on the surface or not
        string sRead = GetRealString(500,fpFile);
		sscanf(sRead.c_str(),"%lf,%lf,%lf,%d,%d,%d\n",&m_dX,&m_dY,&m_dZ,&iTemp1,&iTemp2,&iSurface);
               
	   // read the initial conditions
	   sRead = GetRealString(500,fpFile);
	   sscanf(sRead.c_str(),"%lf\n",&dTemp);
	   m_oPotentialDOF.SetPrimaryValue(dTemp);
	   
	   // read the fluxes
	   sRead = GetRealString(500,fpFile);
	   double dFluxX = 0.0;
	   double dFluxY = 0.0;
	   double dFluxZ = 0.0;
	   sscanf(sRead.c_str(),"%lf,%lf,%lf,%lf\n",&dTemp,&dFluxX,&dFluxY,&dFluxZ);
	   
	   m_oPotentialDOF.SetSecondaryValue(dTemp);
	   m_oFluxes.Set(dFluxX,dFluxY,dFluxZ);
		
		if(iTemp1 == 1)
 		{
 			m_oPotentialDOF.SetCondition(true);
 		}
 		else if(iTemp1 == 0)
 		{
 			m_oPotentialDOF.SetCondition(false);
 		}
 		m_poLoad = pvpoLoads->at(iTemp2 - 1);
 		
 		if(iSurface == 1)
 		{
 			SetOnSurface();
 		}
	}
	void FEMPotentialNode::WriteNode(FILE* fpFile) const
	{
		// write node type
		fprintf(fpFile,"%d\n",GetType());
		int iCondition = 0;
		int iSurface = 0;
		if(m_oPotentialDOF.GetCondition())
		{
			iCondition = 1;
		}
		if(m_bOnSurface)
		{
			iSurface = 1;
		}
		fprintf(fpFile,"%e,%e,%e,%d,%d,%d\n",m_dX,m_dY,m_dZ,iCondition,m_poLoad->GetID(),iSurface);
		fprintf(fpFile,"%e\n",m_oPotentialDOF.GetPrimaryValue());
		fprintf(fpFile,"%e,%e,%e,%e\n",m_oPotentialDOF.GetSecondaryValue(),m_oFluxes.GetX(),m_oFluxes.GetY(),m_oFluxes.GetZ());
	}
	bool FEMPotentialNode::IsConstrained() const
	{
		return m_oPotentialDOF.GetCondition();
	}
	double FEMPotentialNode::GetPotential() const
	{
		return m_oPotentialDOF.GetPrimaryValue();
	}
	unsigned int FEMPotentialNode::SetDOFIndices(const unsigned int& iCurrentDOFIndex)
	{
		m_oPotentialDOF.SetIndex(iCurrentDOFIndex);
		return (iCurrentDOFIndex + 1);
	}
	vector<FEMDegreeOfFreedom*> FEMPotentialNode::GetDOFs()
	{
		vector<FEMDegreeOfFreedom*> vpoDOFs;
		vpoDOFs.resize(1);
		vpoDOFs[0] = &m_oPotentialDOF;
		return vpoDOFs;
	}
	void FEMPotentialNode::SetFluxes(const Vector& oFluxes)
	{
		m_oFluxes = oFluxes;
	}
	Vector FEMPotentialNode::GetFluxes() const
	{
		return m_oFluxes;
	}
}



