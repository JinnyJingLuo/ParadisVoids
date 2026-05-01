#ifndef FEMSOLVER_H_
#define FEMSOLVER_H_

#include "FEMElement.h"
#include "SparseMatrix.h"
#include "MainDataStructure.h"
#include "string"

using namespace std;
 
namespace FEMSystem
{
	enum FEMSolverType
	{
		NullFEMSolver = 0,
		PotentialFEMSolver = 1,
		SolidFEMSolver = 2,
		ThermoMechanicalFEMSolver = 3
	};
	
	class FEMSolver
	{
	public:
		virtual ~FEMSolver();
		virtual FEMSolver& operator=(const FEMSolver& oSolver);
		virtual void Reset();
		void SetDataStructure(MainDataStructure* poStructure);
		virtual void InitializeMatrices(const double& dTimeStep = 1.0) = 0;
		virtual void Solve(const double& dTime) = 0;
		virtual void WriteFEMSolution(const unsigned int& iStep) const;
		virtual void WriteFEMSolutionToParaview(const unsigned int& iStep) const = 0;
		static FEMSolver* CreateSolverByPhysics(FEMSolverType eType,const ProblemType& eProblemType);
		static FEMSolver* CreateSolverByPhysicsIndex(const unsigned int& iIndex,const ProblemType& eProblemType = StaticProblem);
		int DeterminePointLocation(Point* poPoint,FEMElement*& poElement,vector<double>& vdNaturalCoordinates) const;
		double GetTimeStep() const;
		string GetOutputFileName(const unsigned int& iStep) const;
		
	private:
	
	protected:
		virtual void Initialize();
		int DeterminePointLocationNaively(Point* poPoint,FEMElement*& poElement,vector<double>& vdNaturalCoordinates) const;
		MainDataStructure* m_poDataStructure;
		unsigned int m_iDOFCount;
		double m_dTimeStep;
	};
}

#endif
 


