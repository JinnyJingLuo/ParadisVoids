#ifndef FEMPOTENTIALSOLVER_H_
#define FEMPOTENTIALSOLVER_H_

 
 #include "FEMSolver.h"
 
 namespace FEMSystem
 {
 	class FEMPotentialSolver : public FEMSolver
 	{
 	public:
 		FEMPotentialSolver();
 		FEMPotentialSolver(const FEMPotentialSolver& oSolver);
 		virtual ~FEMPotentialSolver();
 		virtual FEMPotentialSolver& operator=(const FEMPotentialSolver& oSolver);
 		virtual void Reset();
 		
 		void InitializeMatrices(const double& dTimeStep = 1.0);
 		void Solve(const double& dTime);
 		
 		Vector GetFlux(Point* poPoint,unsigned int& iStatus);
 		double GetPotential(Point* poPoint,unsigned int& iStatus);
		virtual void WriteFEMSolutionToParaview(const unsigned int& iStep) const;
 	
 	private:
 	
 	protected:
		virtual void Initialize();
		
 		void GenerateMatrices(const unsigned int& iDOFCount);
 		void PartitionMatrices(const SparseMatrix& oConductionMatrix);
 		void GenerateNodalStateVectors(const double& dTime);
 		void UpdateNodalStateVectors();
 		void GenerateNodalFluxesMatrix();
 		
 		Matrix m_oFluxes1;
 		Matrix m_oFluxes2;
 		Matrix m_oPotentials1;
 		Matrix m_oPotentials2;

 		vector<unsigned int> m_viUnknownPotentialsIndices;
 		vector<unsigned int> m_viKnownPotentialsIndices;
 		vector<unsigned int> m_viUnknownPotentialsIndicesReverseMap;
 		vector<unsigned int> m_viKnownPotentialsIndicesReverseMap;
 		
 		SparseMatrix m_oConductionMatrix11;
 		SparseMatrix m_oConductionMatrix12;
 		SparseMatrix m_oConductionMatrix21;
 		SparseMatrix m_oConductionMatrix22;
 	};
 }
 
 #endif
 
 
