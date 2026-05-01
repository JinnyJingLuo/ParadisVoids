#ifndef TOPOLOGYSERVER_H_
#define TOPOLOGYSERVER_H_

#include "MainDataStructure.h"
#include "TopologicalOperation.h"
#include "ChangeConnectionOperation.h"

namespace TopologySystem
{
	class TopologyServer
	{
	public:
		static TopologyServer* GetInstance();
		~TopologyServer();
		void SetDataStructure(MainDataStructure* poDataStructure);
		void Remesh();
		void HandleCollisions();
		
	private:
	
	protected:
		static TopologyServer* m_poInstance;
		TopologyServer();
		void Reset();
		void Initialize();
		void ClearOperations();
		void CommunicateOperations();
		void PackOperations(vector< list<double> >* pvldBuffers);
		void UnpackOperations(vector< list<double> >* pvldBuffers);
		void ProcessOperations();
		void GenerateNodeCollisionLists();
		void ClearNodeCollisionLists();
		void HandleSegmentSegmentCollisions();
		void HandleHingeJointCollisions();
		void HandleTriangularLoops();
		static bool IsNodeCollisionValid(DislocationNode* poNode);
		static bool AreSegmentsCloseAndApproaching(DislocationNode* poNode1,DislocationNode* poNode2,DislocationNode* poNode3,DislocationNode* poNode4,Point& oNearPoint1,Point& oNearPoint2);
		static bool AreNodesCloseAndApproaching(DislocationNode* poNode1,DislocationNode* poNode2);
		DislocationNode* SplitArm(DislocationNode* poNode1,DislocationNode* poNode2,DislocationSegment* poArm);
		DislocationNode* SplitArm(DislocationNode* poNode1,DislocationNode* poNode2,DislocationSegment* poArm,const Point& oSplittingPoint);
		void ProcessChangeConnectionOperation(ChangeConnectionOperation* poOperation);
		bool IsNodeSuperior(DislocationNode* poNode1,DislocationNode* poNode2) const;
		bool IsTriangleCollapsing(DislocationNode* poNode1,DislocationNode* poNode2,DislocationNode* poNode3);
		bool GetExactCollisionPoint(DislocationNode* poNode1,DislocationNode* poNode2,Point& oCollisionPoint);
		Point GetExactCollisionPoint(const Point& oPoint1,const Point& oPoint2,const Vector& oNormal1,const Vector& oNormal2);
		MainDataStructure* m_poDataStructure;
		list<TopologicalOperation*> m_lpoOperations;
		vector< list<DislocationNode*> > m_vlpoNodes;
		unsigned int m_iCollisionCallsCount;
	};
}

#endif


