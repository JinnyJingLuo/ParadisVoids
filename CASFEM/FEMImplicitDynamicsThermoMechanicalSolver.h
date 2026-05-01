// Ahmed M. Hussein

#ifndef FEMIMPLICITDYNAMICSTHERMOMECHANICALSOVLER_H_
#define FEMIMPLICITDYNAMICSTHERMOMECHANICALSOVLER_H_


#include "FEMSolver.h"

namespace FEMSystem
{
	class FEMImplicitDynamicsThermoMechanicalSolver : public FEMSolver
	{
	public:
		FEMImplicitDynamicsThermoMechanicalSolver();
		FEMImplicitDynamicsThermoMechanicalSolver(const FEMImplicitDynamicsThermoMechanicalSolver& oSolver);
		virtual ~FEMImplicitDynamicsThermoMechanicalSolver();
		virtual FEMImplicitDynamicsThermoMechanicalSolver& operator=(const FEMImplicitDynamicsThermoMechanicalSolver& oSolver);
		virtual void Reset();
		
		void InitializeMatrices(const double& dTimeStep = 1.0);
		void Solve(const double& dTime);
		
		Matrix GetStress(Point* poPoint,unsigned int& iStatus);
		Vector GetDisplacement(Point* poPoint,unsigned int& iStatus);
		double GetTemperature(Point* poPoint,unsigned int& iStatus);
		Vector GetFlux(Point* poPoint,unsigned int& iStatus);
		virtual void WriteFEMSolutionToParaview(const unsigned int& iStep) const;
		
	private:
	
	protected:
		virtual void Initialize();
		
		void GenerateMatrices(const unsigned int& iDOFCount);
		void PartitionMatrices(const SparseMatrix& oStiffnessMatrix,const SparseMatrix& oDampingMatrix,const SparseMatrix& oMassMatrix);
		void GenerateNodalStateVectors();
		void UpdateNodalStateVectors();
		void GenerateNodalStressesMatrix();
		void GenerateNodalFluxesMatrix();
		
		// In each step, we solve the equation 
		// D u(i + 1) = R(i + 1) + S(i + 1) + T(i + 1)
		// D = K + a0 M + a1 C							K is the stiffness mateix
		// S(i + 1) = M(a0 u(i) + a2 v(i) + a3 a(i))	M is the mass matrix
		// T(i + 1) = C(a1 u(i) + a4 v(i) + a5 a(i))	C is the damping matrix, if it exists
		// a(i + 1) = a0 (u(i + 1) - u(i)) - a2 v(i) - a3 a(i)
		// v(i + 1) = v(i) + a6 a(i) + a7 a(i + 1)
		// a0 = 1/alpha/dT/dT
		// a1 = delta/alpha/dT
		// a2 = 1.0/alpha/dT
		// a3 = 1.0/2/alpha - 1
		// a4 = delta/alpha - 1
		// a5 = dT/2 (delta/alpha - 2)
		// a6 = dT(1 - delta)
		// a7 = delta dT
		// delta >= 0.5, alpha >= 0.25(0.5 + delta)^2
		// delta is taken to be 0.5 always, alpha is taken to be 0.25
		Matrix m_oForces1;
		Matrix m_oForces2;
		Matrix m_oCurrentDisplacements1;
		Matrix m_oCurrentDisplacements2;
		Matrix m_oPreviousDisplacements1;
		Matrix m_oPreviousDisplacements2;
		Matrix m_oCurrentVelocities1;
		Matrix m_oCurrentVelocities2;
		Matrix m_oPreviousVelocities1;
		Matrix m_oPreviousVelocities2;
		Matrix m_oCurrentAccelerations1;
		Matrix m_oCurrentAccelerations2;
		Matrix m_oPreviousAccelerations1;
		Matrix m_oPreviousAccelerations2;
	
		vector<unsigned int> m_viUnknownDisplacementsIndices;
		vector<unsigned int> m_viKnownDisplacementsIndices;
		vector<unsigned int> m_viUnknownDisplacementsIndicesReverseMap;
		vector<unsigned int> m_viKnownDisplacementsIndicesReverseMap;
		
		SparseMatrix m_oDynamicMatrix11;
		SparseMatrix m_oDynamicMatrix12;
		SparseMatrix m_oDynamicMatrix21;
		SparseMatrix m_oDynamicMatrix22;
		
		SparseMatrix m_oMassMatrix11;
		SparseMatrix m_oMassMatrix12;
		SparseMatrix m_oMassMatrix21;
		SparseMatrix m_oMassMatrix22;
		
		SparseMatrix m_oDampingMatrix11;
		SparseMatrix m_oDampingMatrix12;
		SparseMatrix m_oDampingMatrix21;
		SparseMatrix m_oDampingMatrix22;
		
		// Newmark method parameters
		double m_dA0;
		double m_dA1;
		double m_dA2;
		double m_dA3;
		double m_dA4;
		double m_dA5;
		double m_dA6;
		double m_dA7;
	};
}


#endif


