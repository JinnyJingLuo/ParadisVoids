// Ahmed M. Hussein

#ifndef FEMDEGREEOFFREEDOM_H_
#define FEMDEGREEOFFREEDOM_H_

namespace FEMSystem
{
	class FEMDegreeOfFreedom
	{
	public:
		FEMDegreeOfFreedom();
		FEMDegreeOfFreedom(const FEMDegreeOfFreedom& oDOF);
		virtual ~FEMDegreeOfFreedom();
		FEMDegreeOfFreedom operator=(const FEMDegreeOfFreedom& oDOF);
		void Reset();
		void SetCondition(const bool& bCondition);
		bool GetCondition() const;
		void SetIndex(const unsigned int& iIndex);
		unsigned int GetIndex() const;
		void SetPrimaryValue(const double& dValue);
		double GetPrimaryValue() const;
		void SetSecondaryValue(const double& dValue);
		double GetSecondaryValue() const;
		void AddToPrimaryValue(const double& dIncrement);
		void AddToSecondaryValue(const double& dIncrement);
		void SetConstraint(const double& dValue);
		void AddToConstraint(const double& dIncrement);
		void SetSolution(const double& dValue);
		double GetSolutionValue() const;
		double GetConstraintValue() const;
		
	private:
	
	protected:
		virtual void Initialize();
		// conditions : true for prescribed primary values, false for prescribed secondary values
		bool m_bCondition;
		unsigned int m_iIndex;
		double m_dPrimaryValue;
		double m_dSecondaryValue;
	};
}

#endif


