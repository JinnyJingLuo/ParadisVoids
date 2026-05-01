#ifndef FEMEXPLICITDYNAMICSSOLIDSOLVER_H_
#define FEMEXPLICITDYNAMICSSOLIDSOLVER_H_

#include "FEMSolver.h"

namespace FEMSystem
{
	class FEMExplicitDynamicsSolidSolver : public FEMSolver
	{
	public:
		FEMExplicitDynamicsSolidSolver();
		FEMExplicitDynamicsSolidSolver(const FEMExplicitDynamicsSolidSolver& oSolver);
		virtual ~FEMExplicitDynamicsSolidSolver();
		virtual FEMExplicitDynamicsSolidSolver& operator=(const FEMExplicitDynamicsSolidSolver& oSolver);
		virtual void Reset();
		
		void InitializeMatrices(const double& dTimeStep);
		void Solve(const double& dTime);
		
		Matrix GetStress(Point* poPoint,unsigned int& iStatus);
		Vector GetDisplacement(Point* poPoint,unsigned int& iStatus);
		virtual void WriteFEMSolutionToParaview(const unsigned int& iStep) const;
		
	private:
	
	protected:
		virtual void Initialize();
		
		void GenerateMatrices(const unsigned int& iDOFCount);
		void PartitionMatrices(const SparseMatrix& oTMatrix,const SparseMatrix& oSMatrix,const Matrix& oMassMatrix);
		void GenerateNodalStateVectors();
		void UpdateNodalStateVectors();
		void GenerateNodalStressesMatrix();
		void AdjustTimeStep(const SparseMatrix& oStiffnessMatrix,const Matrix& oMassMatrix);
		void SetMethodParameters(const double& dTimeStep);
		
		// In each step, we solve the equation 
		// D u(i + 1) = R(n) - T u(n) - S u(n - 1)
		// T = K - a2 M		K is the stiffness mateix, M is the mass matrix
		// S = a0 M - a1 C	C is the damping matrix, if it exists
		// D = a0 M + a1 C								dynamic matrix
		// v(n) = a1*[u(n + 1) - u(n - 1)]				velocity
		// a(n) = a0*[u(n + 1) - 2 u(n) + u(n - 1)]		acceleration
		// u(n - 1) = u(n) -  dT v(n) + a3 a(n)
		// a0 = 1/dT/dT
		// a1 = 0.5/dT
		// a2 = 2.0*a0
		// a3 = 1.0/a2
		Matrix m_oForces1;
		Matrix m_oForces2;
		Matrix m_oNextDisplacements1;
		Matrix m_oNextDisplacements2;
		Matrix m_oCurrentDisplacements1;
		Matrix m_oCurrentDisplacements2;
		Matrix m_oPreviousDisplacements1;
		Matrix m_oPreviousDisplacements2;
		Matrix m_oVelocities1;
		Matrix m_oVelocities2;
		Matrix m_oAccelerations1;
		Matrix m_oAccelerations2;
	
		vector<unsigned int> m_viUnknownDisplacementsIndices;
		vector<unsigned int> m_viKnownDisplacementsIndices;
		vector<unsigned int> m_viUnknownDisplacementsIndicesReverseMap;
		vector<unsigned int> m_viKnownDisplacementsIndicesReverseMap;
		
		Matrix m_oDynamicMatrix1;
		Matrix m_oDynamicMatrix2;
		
		SparseMatrix m_oT11Matrix;
		SparseMatrix m_oT12Matrix;
		SparseMatrix m_oT21Matrix;
		SparseMatrix m_oT22Matrix;
		SparseMatrix m_oS11Matrix;
		SparseMatrix m_oS12Matrix;
		SparseMatrix m_oS21Matrix;
		SparseMatrix m_oS22Matrix;
		
		// explicit method parameters
		double m_dA0;
		double m_dA1;
		double m_dA2;
		double m_dA3;
	};
}

#endif



