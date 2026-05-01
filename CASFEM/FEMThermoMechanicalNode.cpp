#include "FEMThermoMechanicalNode.h"
#include "Tools.h"

using namespace SupportSystem;

namespace FEMSystem
{
	unsigned int FEMThermoMechanicalNode::ThermoMechanicalDOFCount = 4;
	FEMThermoMechanicalNode::FEMThermoMechanicalNode()
	{
	  Initialize();

	}
	FEMThermoMechanicalNode::~FEMThermoMechanicalNode()
	{
	  		Reset();
	}
	FEMThermoMechanicalNode::FEMThermoMechanicalNode(const double& dX,const double& dY,const double& dZ)
	{
		Initialize();
		Set(dX,dY,dZ);
	}
	FEMThermoMechanicalNode::FEMThermoMechanicalNode(const FEMThermoMechanicalNode& oNode)
	{
		*this = oNode;
	}
	FEMThermoMechanicalNode& FEMThermoMechanicalNode::operator=(const FEMThermoMechanicalNode& oNode)
	{
	    FEMSolidNode::operator =(oNode);
	    m_oTDOF = oNode.m_oTDOF;
		m_poTLoad = oNode.m_poTLoad;
		m_oFluxes = oNode.m_oFluxes;
		m_dHeatingRate = oNode.m_dHeatingRate;
		m_dHeatingAcceleration = oNode.m_dHeatingAcceleration;
		return *this;
	}
	void FEMThermoMechanicalNode::Reset()
	{
		FEMSolidNode::Reset();
		m_oTDOF.Reset();
		m_oFluxes.Set(0.0,0.0,0.0);
		m_poTLoad = NULL;
		m_dHeatingRate = 0.0;
		m_dHeatingAcceleration = 0.0;
	}
	unsigned int FEMThermoMechanicalNode::GetDOFCount() const
	{
		return ThermoMechanicalDOFCount;
	}
	FEMNodeType FEMThermoMechanicalNode::GetType() const
	{
		return ThermoMechanicalFEMNode;
	}
	FEMDegreeOfFreedom* FEMThermoMechanicalNode::GetTDOF()
	{
		return &m_oTDOF;
	}
	void FEMThermoMechanicalNode::ResetLoads()
	{
		FEMSolidNode::ResetLoads();
		m_oTDOF.SetPrimaryValue(0.0);
		m_oTDOF.SetSecondaryValue(0.0);
	}
	FEMNode* FEMThermoMechanicalNode::Clone() const
	{
		return new FEMThermoMechanicalNode(*this);
	}
	void FEMThermoMechanicalNode::ApplyLoads(const double& dTime)
	{
		ResetLoads();
		FEMSolidNode::ApplyLoads(dTime);
		if(m_poTLoad != NULL)
		{
			m_oTDOF.AddToConstraint(m_poTLoad->Get(*this,dTime));
		}
	}
	void FEMThermoMechanicalNode::AddFlux(const double& dFlux)
	{
		m_oTDOF.AddToSecondaryValue(dFlux);
	}
	void FEMThermoMechanicalNode::SetTLoad(FEMLoad* poLoad)
	{
		m_poTLoad = poLoad;
	}
	FEMLoad* FEMThermoMechanicalNode::GetTLoad() const
	{
		return m_poTLoad;
	}
	void FEMThermoMechanicalNode::ReadNode(FILE* fpFile,vector<FEMLoad*>* pvpoLoads)
	{
		unsigned int iTemp1 = 0;
		unsigned int iTemp2 = 0;
		unsigned int iTemp3 = 0;
		unsigned int iTemp4 = 0;
		unsigned int iTemp5 = 0;
		unsigned int iTemp6 = 0;
		unsigned int iTemp7 = 0;
		unsigned int iTemp8 = 0;
		unsigned int iSurface = 0;
		double dUx = 0.0;
		double dUy = 0.0;
		double dUz = 0.0;
		double dT = 0.0;
		double dVx = 0.0;
		double dVy = 0.0;
		double dVz = 0.0;
		double dVT = 0.0;
		double dAx = 0.0;
		double dAy = 0.0;
		double dAz = 0.0;
		double dAT = 0.0;
		
		// read the node position, conditions and whether it is on the surface or not
 		string sRead = GetRealString(500,fpFile);
 		sscanf(sRead.c_str(),"%lf,%lf,%lf,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",&m_dX,&m_dY,&m_dZ,&iTemp1,&iTemp2,&iTemp3,&iTemp4,&iTemp5,&iTemp6,&iTemp7,&iTemp8,&iSurface);
		
		// read the initial conditions
		sRead = GetRealString(500,fpFile);
		sscanf(sRead.c_str(),"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",&dUx,&dUy,&dUz,&dT,&dVx,&dVy,&dVz,&dVT,&dAx,&dAy,&dAz,&dAT);
		m_oXDOF.SetPrimaryValue(dUx);
		m_oYDOF.SetPrimaryValue(dUy);
		m_oZDOF.SetPrimaryValue(dUz);
		m_oTDOF.SetPrimaryValue(dT);
		SetVelocity(Vector(dVx,dVy,dVz));
		SetAcceleration(Vector(dAx,dAy,dAz));
		SetHeatingRate(dVT);
		SetHeatingAcceleration(dAT);
		
		// read the forces and stresses
		sRead = GetRealString(500,fpFile);
		double dFluxZ = 0.0;
		sscanf(sRead.c_str(),"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",&dUx,&dUy,&dUz,&dT,&dVx,&dVy,&dVz,&dAx,&dAy,&dAz,&dVT,&dAT,&dFluxZ);
		
		m_oXDOF.SetSecondaryValue(dUx);
		m_oYDOF.SetSecondaryValue(dUy);
		m_oZDOF.SetSecondaryValue(dUz);
		m_oTDOF.SetSecondaryValue(dT);
		
		m_oStresses.SetSize(3,3);
		m_oStresses.Set(1,1,dVx);
		m_oStresses.Set(1,2,dAy);
		m_oStresses.Set(1,3,dAz);
		m_oStresses.Set(2,1,dAy);
		m_oStresses.Set(2,2,dVy);
		m_oStresses.Set(2,3,dAx);
		m_oStresses.Set(3,1,dAz);
		m_oStresses.Set(3,2,dAx);
		m_oStresses.Set(3,3,dVz);
		
		m_oFluxes.Set(dVT,dAT,dFluxZ);
		
 		// X condition
 		if(iTemp1 == 1)
 		{
 			m_oXDOF.SetCondition(true);
 		}
 		else if(iTemp1 == 0)
 		{
 			m_oXDOF.SetCondition(false);
 		}
 		m_poXLoad = pvpoLoads->at(iTemp2 - 1);
 		
 		// Y condition
 		if(iTemp3 == 1)
 		{
 			m_oYDOF.SetCondition(true);
 		}
 		else if(iTemp3 == 0)
 		{
 			m_oYDOF.SetCondition(false);
 		}
 		m_poYLoad = pvpoLoads->at(iTemp4 - 1);
 		
 		// Z condition
 		if(iTemp5 == 1)
 		{
 			m_oZDOF.SetCondition(true);
 		}
 		else if(iTemp5 == 0)
 		{
 			m_oZDOF.SetCondition(false);
 		}
 		m_poZLoad = pvpoLoads->at(iTemp6 - 1);
 		
 		// T condition
 		if(iTemp7 == 1)
 		{
 			m_oTDOF.SetCondition(true);
 		}
 		else if(iTemp7 == 0)
 		{
 			m_oTDOF.SetCondition(false);
 		}
 		m_poTLoad = pvpoLoads->at(iTemp8 - 1);
 		 		
 		if(iSurface == 1)
 		{
 			SetOnSurface();
 		}
	}
	void FEMThermoMechanicalNode::WriteNode(FILE* fpFile) const
	{
		// write node type
		fprintf(fpFile,"%d\n",GetType());
		int iXCondition = 0;
		int iYCondition = 0;
		int iZCondition = 0;
		int iTCondition = 0;
		int iSurface = 0;
		if(m_oXDOF.GetCondition())
		{
			iXCondition = 1;
		}
		if(m_oYDOF.GetCondition())
		{
			iYCondition = 1;
		}
		if(m_oZDOF.GetCondition())
		{
			iZCondition = 1;
		}
		if(m_oTDOF.GetCondition())
		{
			iTCondition = 1;
		}
		if(m_bOnSurface)
		{
			iSurface = 1;
		}
		fprintf(fpFile,"%e,%e,%e,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",m_dX,m_dY,m_dZ,iXCondition,m_poXLoad->GetID(),iYCondition,m_poYLoad->GetID(),iZCondition,m_poZLoad->GetID(),iTCondition,m_poTLoad->GetID(),iSurface);
		fprintf(fpFile,"%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n",m_oXDOF.GetPrimaryValue(),m_oYDOF.GetPrimaryValue(),m_oZDOF.GetPrimaryValue(),m_oTDOF.GetPrimaryValue(),m_oVelocity.GetX(),m_oVelocity.GetY(),m_oVelocity.GetZ(),m_dHeatingRate,m_oAcceleration.GetX(),m_oAcceleration.GetY(),m_oAcceleration.GetZ(),m_dHeatingAcceleration);
		fprintf(fpFile,"%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n",m_oXDOF.GetSecondaryValue(),m_oYDOF.GetSecondaryValue(),m_oZDOF.GetSecondaryValue(),m_oTDOF.GetSecondaryValue(),m_oStresses.Get(1,1),m_oStresses.Get(2,2),m_oStresses.Get(3,3),m_oStresses.Get(3,2),m_oStresses.Get(2,1),m_oStresses.Get(1,3),m_oFluxes.GetX(),m_oFluxes.GetY(),m_oFluxes.GetZ());
	}
	double FEMThermoMechanicalNode::GetTemperature() const
	{
		return m_oTDOF.GetPrimaryValue();
	}
	unsigned int FEMThermoMechanicalNode::SetDOFIndices(const unsigned int& iCurrentDOFIndex)
	{
		m_oXDOF.SetIndex(iCurrentDOFIndex);
		m_oYDOF.SetIndex(iCurrentDOFIndex + 1);
		m_oZDOF.SetIndex(iCurrentDOFIndex + 2);
		m_oTDOF.SetIndex(iCurrentDOFIndex + 3);
		return (iCurrentDOFIndex + 4);
	}
	vector<FEMDegreeOfFreedom*> FEMThermoMechanicalNode::GetDOFs()
	{
		vector<FEMDegreeOfFreedom*> vpoDOFs;
		vpoDOFs.resize(4);
		vpoDOFs[0] = &m_oXDOF;
		vpoDOFs[1] = &m_oYDOF;
		vpoDOFs[2] = &m_oZDOF;
		vpoDOFs[3] = &m_oTDOF;
		return vpoDOFs;
	}
	void FEMThermoMechanicalNode::SetFluxes(const Vector& oFluxes)
	{
		m_oFluxes = oFluxes;
	}
	Vector FEMThermoMechanicalNode::GetFluxes() const
	{
		return m_oFluxes;
	}
	void FEMThermoMechanicalNode::SetHeatingRate(const double& dRate)
	{
		m_dHeatingRate = dRate;
	}
	void FEMThermoMechanicalNode::SetHeatingAcceleration(const double& dAcceleration)
	{
		m_dHeatingAcceleration = dAcceleration;
	}
	double FEMThermoMechanicalNode::GetHeatingRate()
	{
		return m_dHeatingRate;
	}
	double FEMThermoMechanicalNode::GetHeatingAcceleration()
	{
		return m_dHeatingAcceleration;
	}
	void FEMThermoMechanicalNode::Initialize()
	{

		FEMSolidNode::Initialize();
		m_oTDOF.Reset();
		ResetLoads();
		m_oFluxes.Set(0.0,0.0,0.0);
		m_poTLoad = NULL;
		m_dHeatingRate = 0.0;
		m_dHeatingAcceleration = 0.0;

	}
}



