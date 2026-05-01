#ifndef DYNAMICSSERVER_H_
#define DYNAMICSSERVER_H_

#include "MainDataStructure.h"
#include "Matrix.h"

using namespace EZ;

namespace DynamicsSystem
{
	class DynamicsServer
	{
	public:
		static DynamicsServer* GetInstance();
		~DynamicsServer();
		void SetDataStructure(MainDataStructure* poDataStructure);
		void ApplyLoad();
		void ComputeForces();
		void ComputeVelocities();
		void ComputeTimeStep();
		void IntegrateMotion();
		
	private:
	
	protected:
		static DynamicsServer* m_poInstance;
		DynamicsServer();
		void Reset();
		void Initialize();
		void ResetForces();
		void ComputeAppliedForces();
		void ComputeSelfForces();
		void ComputeInteractionsForces();
		Vector CalculatePeachKoehlerForce(DislocationNode* poNode,DislocationNode* poNeighbour,DislocationSegment* poArm,Matrix* poStress);
		Vector CalculateFCCVelocity(DislocationNode* poNode);
		void CommunicateVelocities();
		void PackVelocities(vector< list<double> >* pvldData);
		void UnpackVelocities(vector< list<double> >* pvldData);
		MainDataStructure* m_poDataStructure;
	};
}

#endif


