// Ahmed M. Hussein

#ifndef FEMTHERMOMECHANICALSOVLER_H_
#define FEMTHERMOMECHANICALSOVLER_H_


#include "FEMSolidSolver.h"

namespace FEMSystem
{
	class FEMThermoMechanicalSolver : public FEMSolidSolver
	{
	public:
 		FEMThermoMechanicalSolver();
 		FEMThermoMechanicalSolver(const FEMThermoMechanicalSolver& oSolver);
 		virtual ~FEMThermoMechanicalSolver();
 		virtual FEMThermoMechanicalSolver& operator=(const FEMThermoMechanicalSolver& oSolver);
 		virtual void Reset();
		
		void InitializeMatrices(const double& dTimeStep = 1.0);
		void Solve(const double& dTime);
		virtual void WriteFEMSolutionToParaview(const unsigned int& iStep) const;
		
		double GetTemperature(Point* poPoint,unsigned int& iStatus);
		Vector GetFlux(Point* poPoint,unsigned int& iStatus);
		
	private:
	
	protected:
		virtual void Initialize();
		virtual void GenerateMatrices();
		virtual void PartitionMatrices(const SparseMatrix& oStiffnessMatrix);
		virtual void GenerateNodalStateVectors();
		virtual void UpdateNodalStateVectors();
		virtual void GenerateNodalFluxesMatrix();
	};
}


#endif


