#ifndef DISLOCATIONSIMULATOR_H_
#define DISLOCATIONSIMULATOR_H_

#include "MainDataStructure.h"
#include "TopologyServer.h"
#include "DynamicsServer.h"
#include "SurfaceServer.h"
#include "BalancingServer.h"

using namespace std;
using namespace TopologySystem;
using namespace DynamicsSystem;
using namespace BalancingSystem;

class DislocationSimulator
{
public:
	static DislocationSimulator* GetInstance();
	~DislocationSimulator();
	void Start(const string& sInputFileName,const unsigned int& iDomainsCount,const unsigned int& iDomainID);
	void End();
	void Load();
	void Move();
	void Write();
	void Remesh();
	void Rebalance();
	void Collide();
	
private:

protected:
	static DislocationSimulator* m_poInstance;
	void Reset();
	void Initialize();
	MainDataStructure* m_poDataStructure;
	TopologyServer* m_poTopologyServer;
	DynamicsServer* m_poDynamicsServer;
	SurfaceServer* m_poSurfaceServer;
	BalancingServer* m_poBalancingServer;
	unsigned int m_iOutputFileIndex;
};
	
#endif


