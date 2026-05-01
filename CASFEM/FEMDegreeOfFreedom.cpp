// Ahmed M. Hussein

#include "iostream"
#include "FEMDegreeOfFreedom.h"

namespace FEMSystem
{
	FEMDegreeOfFreedom::FEMDegreeOfFreedom()
	{
		Initialize();
	}
	FEMDegreeOfFreedom::FEMDegreeOfFreedom(const FEMDegreeOfFreedom& oDOF)
	{
		*this = oDOF;
	}
	FEMDegreeOfFreedom::~FEMDegreeOfFreedom()
	{
	
	}
	FEMDegreeOfFreedom FEMDegreeOfFreedom::operator=(const FEMDegreeOfFreedom& oDOF)
	{
		m_bCondition = oDOF.m_bCondition;
		m_iIndex = oDOF.m_iIndex;
		m_dPrimaryValue = oDOF.m_dPrimaryValue;
		m_dSecondaryValue = oDOF.m_dSecondaryValue;
		return *this;
	}
	void FEMDegreeOfFreedom::Reset()
	{
		Initialize();
	}
	void FEMDegreeOfFreedom::SetCondition(const bool& bCondition)
	{
		m_bCondition = bCondition;
	}
	bool FEMDegreeOfFreedom::GetCondition() const
	{
		return m_bCondition;
	}
	void FEMDegreeOfFreedom::SetIndex(const unsigned int& iIndex)
	{
		m_iIndex = iIndex;
	}
	unsigned int FEMDegreeOfFreedom::GetIndex() const
	{
		return m_iIndex;
	}
	void FEMDegreeOfFreedom::SetPrimaryValue(const double& dValue)
	{
		m_dPrimaryValue = dValue;
	}
	double FEMDegreeOfFreedom::GetPrimaryValue() const
	{
		return m_dPrimaryValue;
	}
	void FEMDegreeOfFreedom::SetSecondaryValue(const double& dValue)
	{
		m_dSecondaryValue = dValue;
	}
	double FEMDegreeOfFreedom::GetSecondaryValue() const
	{
		return m_dSecondaryValue;
	}
	void FEMDegreeOfFreedom::AddToPrimaryValue(const double& dIncrement)
	{
		m_dPrimaryValue = m_dPrimaryValue + dIncrement;
	}
	void FEMDegreeOfFreedom::AddToSecondaryValue(const double& dIncrement)
	{
		m_dSecondaryValue = m_dSecondaryValue + dIncrement;
	}
	void FEMDegreeOfFreedom::SetConstraint(const double& dValue)
	{
		if(m_bCondition)
		{
			m_dPrimaryValue = dValue;
		}
		else
		{
			m_dSecondaryValue = dValue;
		}
	}
	void FEMDegreeOfFreedom::AddToConstraint(const double& dIncrement)
	{
		if(m_bCondition)
		{
			m_dPrimaryValue = m_dPrimaryValue + dIncrement;
		}
		else
		{
			m_dSecondaryValue = m_dSecondaryValue + dIncrement;
		}
	}
	void FEMDegreeOfFreedom::SetSolution(const double& dValue)
	{
		if(m_bCondition)
		{
			m_dSecondaryValue = dValue;
		}
		else
		{
			m_dPrimaryValue = dValue;
		}
	}
	double FEMDegreeOfFreedom::GetConstraintValue() const
	{
		if(m_bCondition)
		{
			return m_dPrimaryValue;
		}
		return m_dSecondaryValue;
	}
	double FEMDegreeOfFreedom::GetSolutionValue() const
	{
		if(m_bCondition)
		{
			return m_dSecondaryValue;
		}
		return m_dPrimaryValue;
	}
	void FEMDegreeOfFreedom::Initialize()
	{
		m_bCondition = false;
		m_iIndex = 0;
		double m_dPrimaryValue = 0.0;
		double m_dSecondaryValue = 0.0;
	}
}


