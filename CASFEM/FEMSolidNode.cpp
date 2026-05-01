// Ahmed M. Hussein
#include "FEMSolidNode.h"
#include "Tools.h"

using namespace SupportSystem;

namespace FEMSystem
{
	unsigned int FEMSolidNode::SolidDOFCount = 3;
	FEMSolidNode::FEMSolidNode()
	{
		Initialize();
	}
	FEMSolidNode::~FEMSolidNode()
	{
		Reset();
	}
	FEMSolidNode::FEMSolidNode(const FEMSolidNode& oNode)
	{
		*this = oNode;
	}
	FEMSolidNode::FEMSolidNode(const double& dX,const double& dY,const double& dZ)
	{
		Initialize();
		Set(dX,dY,dZ);
	}
	FEMSolidNode& FEMSolidNode::operator=(const FEMSolidNode& oNode)
	{
        FEMNode::operator =(oNode);
        m_oXDOF = oNode.m_oXDOF;
        m_oYDOF = oNode.m_oYDOF;
        m_oZDOF = oNode.m_oZDOF;
		m_oStresses = oNode.m_oStresses;
		m_poXLoad = oNode.m_poXLoad;
		m_poYLoad = oNode.m_poYLoad;
		m_poZLoad = oNode.m_poZLoad;
		m_oVelocity = oNode.m_oVelocity;
		m_oAcceleration = oNode.m_oAcceleration;
		return *this;
	}
	unsigned int FEMSolidNode::GetDOFCount() const
	{
		return SolidDOFCount;
	}
	FEMDegreeOfFreedom* FEMSolidNode::GetXDOF()
	{
		return &m_oXDOF;
	}
	FEMDegreeOfFreedom* FEMSolidNode::GetYDOF()
	{
		return &m_oYDOF;
	}
	FEMDegreeOfFreedom* FEMSolidNode::GetZDOF()
	{
		return &m_oZDOF;
	}
	void FEMSolidNode::Reset()
	{
		FEMNode::Reset();
		m_oXDOF.Reset();
		m_oYDOF.Reset();
		m_oZDOF.Reset();
		m_oStresses.Reset();
		m_poXLoad = NULL;
		m_poYLoad = NULL;
		m_poZLoad = NULL;
		m_oVelocity.Set(0.0,0.0,0.0);
		m_oAcceleration.Set(0.0,0.0,0.0);
	}
	void FEMSolidNode::ResetLoads()
	{
		m_oXDOF.SetPrimaryValue(0.0);
		m_oYDOF.SetPrimaryValue(0.0);
		m_oZDOF.SetPrimaryValue(0.0);
		m_oXDOF.SetSecondaryValue(0.0);
		m_oYDOF.SetSecondaryValue(0.0);
		m_oZDOF.SetSecondaryValue(0.0);
	}
	FEMNode* FEMSolidNode::Clone() const
	{
		return new FEMSolidNode(*this);
	}
	void FEMSolidNode::ApplyLoads(const double& dTime)
	{
		ResetLoads();
		if(m_poXLoad != NULL)
		{
			m_oXDOF.AddToConstraint(m_poXLoad->Get(*this,dTime));
		}

		if(m_poYLoad != NULL)
		{
			m_oYDOF.AddToConstraint(m_poYLoad->Get(*this,dTime));
		}
		
		if(m_poZLoad != NULL)
		{
			m_oZDOF.AddToConstraint(m_poZLoad->Get(*this,dTime));
		}
	}
	void FEMSolidNode::AddForce(Vector oForce)
	{
		m_oXDOF.AddToSecondaryValue(oForce.GetX());
		m_oYDOF.AddToSecondaryValue(oForce.GetY());
		m_oZDOF.AddToSecondaryValue(oForce.GetZ());
	}
	void FEMSolidNode::SetStresses(const Matrix& oStresses)
	{
		m_oStresses = oStresses;
	}
	Matrix FEMSolidNode::GetStresses() const
	{
		return m_oStresses;
	}
	void FEMSolidNode::Initialize()
	{
	  //		FEMNode::Reset();
	  FEMNode :: Initialize();
		m_oXDOF.Reset();
		m_oYDOF.Reset();
		m_oZDOF.Reset();
		ResetLoads();
		m_oStresses.Reset();
		m_poXLoad = NULL;
		m_poYLoad = NULL;
		m_poZLoad = NULL;
		m_oVelocity.Set(0.0,0.0,0.0);
		m_oAcceleration.Set(0.0,0.0,0.0);
	}
	void FEMSolidNode::SetXLoad(FEMLoad* poLoad)
	{
		m_poXLoad = poLoad;
	}
	void FEMSolidNode::SetYLoad(FEMLoad* poLoad)
	{
		m_poYLoad = poLoad;
	}
	void FEMSolidNode::SetZLoad(FEMLoad* poLoad)
	{
		m_poZLoad = poLoad;
	}
	FEMLoad* FEMSolidNode::GetXLoad() const
	{
		return m_poXLoad;
	}
	FEMLoad* FEMSolidNode::GetYLoad() const
	{
		return m_poYLoad;
	}
	FEMLoad* FEMSolidNode::GetZLoad() const
	{
		return m_poZLoad;
	}
	FEMNodeType FEMSolidNode::GetType() const
	{
		return SolidFEMNode;
	}
	void FEMSolidNode::ReadNode(FILE* fpFile,vector<FEMLoad*>* pvpoLoads)
	{
		unsigned int iTemp1 = 0;
		unsigned int iTemp2 = 0;
		unsigned int iTemp3 = 0;
		unsigned int iTemp4 = 0;
		unsigned int iTemp5 = 0;
		unsigned int iTemp6 = 0;
		unsigned int iSurface = 0;
		double dUx = 0.0;
		double dUy = 0.0;
		double dUz = 0.0;
		double dVx = 0.0;
		double dVy = 0.0;
		double dVz = 0.0;
		double dAx = 0.0;
		double dAy = 0.0;
		double dAz = 0.0;
		
		// read the node position, conditions and whether it is on the surface or not
 		string sRead = GetRealString(500,fpFile);
 		sscanf(sRead.c_str(),"%lf,%lf,%lf,%d,%d,%d,%d,%d,%d,%d\n",&m_dX,&m_dY,&m_dZ,&iTemp1,&iTemp2,&iTemp3,&iTemp4,&iTemp5,&iTemp6,&iSurface);
		
		// read the initial conditions
		sRead = GetRealString(500,fpFile);
		sscanf(sRead.c_str(),"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",&dUx,&dUy,&dUz,&dVx,&dVy,&dVz,&dAx,&dAy,&dAz);
		m_oXDOF.SetPrimaryValue(dUx);
		m_oYDOF.SetPrimaryValue(dUy);
		m_oZDOF.SetPrimaryValue(dUz);
		SetVelocity(Vector(dVx,dVy,dVz));
		SetAcceleration(Vector(dAx,dAy,dAz));
		
		// read the forces and stresses
		sRead = GetRealString(500,fpFile);
		sscanf(sRead.c_str(),"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",&dUx,&dUy,&dUz,&dVx,&dVy,&dVz,&dAx,&dAy,&dAz);
		
		m_oXDOF.SetSecondaryValue(dUx);
		m_oYDOF.SetSecondaryValue(dUy);
		m_oZDOF.SetSecondaryValue(dUz);
		
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
 		
 		if(iSurface == 1)
 		{
 			SetOnSurface();
 		}
	}
	void FEMSolidNode::WriteNode(FILE* fpFile) const
	{
		// write node type
		fprintf(fpFile,"%d\n",GetType());
		int iXCondition = 0;
		int iYCondition = 0;
		int iZCondition = 0;
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
		if(m_bOnSurface)
		{
			iSurface = 1;
		}
		fprintf(fpFile,"%e,%e,%e,%d,%d,%d,%d,%d,%d,%d\n",m_dX,m_dY,m_dZ,iXCondition,m_poXLoad->GetID(),iYCondition,m_poYLoad->GetID(),iZCondition,m_poZLoad->GetID(),iSurface);
		fprintf(fpFile,"%e,%e,%e,%e,%e,%e,%e,%e,%e\n",m_oXDOF.GetPrimaryValue(),m_oYDOF.GetPrimaryValue(),m_oZDOF.GetPrimaryValue(),m_oVelocity.GetX(),m_oVelocity.GetY(),m_oVelocity.GetZ(),m_oAcceleration.GetX(),m_oAcceleration.GetY(),m_oAcceleration.GetZ());
		fprintf(fpFile,"%e,%e,%e,%e,%e,%e,%e,%e,%e\n",m_oXDOF.GetSecondaryValue(),m_oYDOF.GetSecondaryValue(),m_oZDOF.GetSecondaryValue(),m_oStresses.Get(1,1),m_oStresses.Get(2,2),m_oStresses.Get(3,3),m_oStresses.Get(3,2),m_oStresses.Get(2,1),m_oStresses.Get(1,3));
	}
	bool FEMSolidNode::IsConstrained() const
	{
		if(m_oXDOF.GetCondition() || m_oYDOF.GetCondition() || m_oZDOF.GetCondition())
		{
			return true;
		}
		return false;
	}
	Vector FEMSolidNode::GetDisplacement() const
	{
		Vector oDisplacement(m_oXDOF.GetPrimaryValue(),m_oYDOF.GetPrimaryValue(),m_oZDOF.GetPrimaryValue());
		return oDisplacement;
	}
	double FEMSolidNode::GetXDisplacement() const
	{
		return m_oXDOF.GetPrimaryValue();
	}
	double FEMSolidNode::GetYDisplacement() const
	{
		return m_oYDOF.GetPrimaryValue();
	}
	double FEMSolidNode::GetZDisplacement() const
	{
		return m_oZDOF.GetPrimaryValue();
	}
	Vector FEMSolidNode::GetForce() const
	{
		Vector oForce(m_oXDOF.GetSecondaryValue(),m_oYDOF.GetSecondaryValue(),m_oZDOF.GetSecondaryValue());
		return oForce;
	}
	double FEMSolidNode::GetXForce() const
	{
		return m_oXDOF.GetSecondaryValue();
	}
	double FEMSolidNode::GetYForce() const
	{
		return m_oYDOF.GetSecondaryValue();
	}
	double FEMSolidNode::GetZForce() const
	{
		return m_oZDOF.GetSecondaryValue();
	}
	Vector FEMSolidNode::GetVelocity()
	{
		return m_oVelocity;
	}
	Vector FEMSolidNode::GetAcceleration()
	{
		return m_oAcceleration;
	}
	void FEMSolidNode::SetVelocity(const Vector& oVelocity)
	{
		m_oVelocity = oVelocity;
	}
	void FEMSolidNode::SetAcceleration(const Vector& oAcceleration)
	{
		m_oAcceleration = oAcceleration;
	}
	unsigned int FEMSolidNode::SetDOFIndices(const unsigned int& iCurrentDOFIndex)
	{
		m_oXDOF.SetIndex(iCurrentDOFIndex);
		m_oYDOF.SetIndex(iCurrentDOFIndex + 1);
		m_oZDOF.SetIndex(iCurrentDOFIndex + 2);
		return (iCurrentDOFIndex + 3);
	}
	vector<FEMDegreeOfFreedom*> FEMSolidNode::GetDOFs()
	{
		vector<FEMDegreeOfFreedom*> vpoDOFs;
		vpoDOFs.resize(3);
		vpoDOFs[0] = &m_oXDOF;
		vpoDOFs[1] = &m_oYDOF;
		vpoDOFs[2] = &m_oZDOF;
		return vpoDOFs;
	}
}




