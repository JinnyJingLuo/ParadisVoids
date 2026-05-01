 #ifndef FEMSOLIDSOLVER_H_
 #define FEMSOLIDSOLVER_H_
 
 #include "FEMSolver.h"
 
 namespace FEMSystem
 {
 	class FEMSolidSolver : public FEMSolver
 	{
 	public:
 		FEMSolidSolver();
 		FEMSolidSolver(const FEMSolidSolver& oSolver);
 		virtual ~FEMSolidSolver();
 		virtual FEMSolidSolver& operator=(const FEMSolidSolver& oSolver);
 		virtual void Reset();
 		
 		void InitializeMatrices(const double& dTimeStep = 1.0);
 		void Solve(const double& dTime);
 		virtual void WriteFEMSolutionToParaview(const unsigned int& iStep) const;
 		
 		Matrix GetStress(Point* poPoint,unsigned int& iStatus);
 		Vector GetDisplacement(Point* poPoint,unsigned int& iStatus);
 	
 	private:
 	
 	protected:		
		virtual void Initialize();
 		virtual void GenerateMatrices();
 		virtual void PartitionMatrices(const SparseMatrix& oStiffnessMatrix);
 		virtual void GenerateNodalStateVectors();
 		virtual void UpdateNodalStateVectors();
 		virtual void GenerateNodalStressesMatrix();
 		
 		Matrix m_oForces1;
 		Matrix m_oForces2;
 		Matrix m_oDisplacements1;
 		Matrix m_oDisplacements2;

 		vector<unsigned int> m_viUnknownDisplacementsIndices;
 		vector<unsigned int> m_viKnownDisplacementsIndices;
 		vector<unsigned int> m_viUnknownDisplacementsIndicesReverseMap;
 		vector<unsigned int> m_viKnownDisplacementsIndicesReverseMap;
 		
 		SparseMatrix m_oStiffnessMatrix11;
 		SparseMatrix m_oStiffnessMatrix12;
 		SparseMatrix m_oStiffnessMatrix21;
 		SparseMatrix m_oStiffnessMatrix22;	
 	};
 }
 
 #endif
 
 
