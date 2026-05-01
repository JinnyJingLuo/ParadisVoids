#include "DislocationSimulator.h"
#include "MainDataStructure.h"
#include "mpi.h"
#include "string"
#include "TopologyServer.h"
#include "DynamicsServer.h"
#include "sys/stat.h"

using namespace std;
using namespace TopologySystem;
using namespace DynamicsSystem;

DislocationSimulator* DislocationSimulator::m_poInstance = NULL;
DislocationSimulator* DislocationSimulator::GetInstance()
{
	if(m_poInstance == NULL)
	{
		m_poInstance = new DislocationSimulator;
	}
	return m_poInstance;
}
DislocationSimulator::~DislocationSimulator()
{
	Reset();
}
void DislocationSimulator::Start(const string& sInputFileName,const unsigned int& iDomainsCount,const unsigned int& iDomainID)
{
	// get server instances
	m_poDataStructure = MainDataStructure::GetInstance();
	m_poTopologyServer = TopologyServer::GetInstance();
	m_poDynamicsServer = DynamicsServer::GetInstance();
	m_poSurfaceServer = SurfaceServer::GetInstance();
	m_poBalancingServer = BalancingServer::GetInstance();
	// initialize data structure
	m_poDataStructure->SetDomainData(iDomainsCount,iDomainID);
	m_poDataStructure->ReadInput(sInputFileName);
	m_poDataStructure->CommunicateGhosts();
	m_poTopologyServer->SetDataStructure(m_poDataStructure);
	m_poDynamicsServer->SetDataStructure(m_poDataStructure);
	m_poSurfaceServer->SetDataStructure(m_poDataStructure);
	m_poBalancingServer->SetDataStructure(m_poDataStructure);
	// remesh
	m_poTopologyServer->Remesh();
	m_poDataStructure->CommunicateGhosts();
	// create the output directory
	mkdir("OUTPUT",S_IRWXU | S_IRWXG);
	mkdir("OUTPUT//restart",S_IRWXU | S_IRWXG);
	// write initial microstructure after remeshing
	char cOutputFileName[512];
	m_iOutputFileIndex = 0;
	sprintf(cOutputFileName,"OUTPUT//restart//restart_%06d.ddo",m_iOutputFileIndex);
	m_poDataStructure->WriteOutput(cOutputFileName);
}
void DislocationSimulator::End()
{
	Reset();
}
void DislocationSimulator::Load()
{
	m_poDynamicsServer->ApplyLoad();
	m_poDynamicsServer->ComputeForces();
}
void DislocationSimulator::Move()
{
	m_poDynamicsServer->ComputeVelocities();
	m_poDynamicsServer->ComputeTimeStep();
	m_poDynamicsServer->IntegrateMotion();
	m_poSurfaceServer->CheckNodes();
	m_poDataStructure->CommunicateGhosts();
}
void DislocationSimulator::Write()
{
	m_iOutputFileIndex = m_iOutputFileIndex + 1;
	char cOutputFileName[512];
	sprintf(cOutputFileName,"OUTPUT//restart//restart_%06d.ddo",m_iOutputFileIndex);
	m_poDataStructure->WriteOutput(cOutputFileName);
	if(m_poDataStructure->GetDomainID() == 0)
	{
		printf("wrote restart file %s\n",cOutputFileName);
		fflush(NULL);
	}
}
void DislocationSimulator::Remesh()
{
	m_poTopologyServer->Remesh();
	m_poDataStructure->CommunicateGhosts();
}
void DislocationSimulator::Rebalance()
{
	m_poBalancingServer->Migrate();
	m_poDataStructure->CommunicateGhosts();
}
void DislocationSimulator::Collide()
{
	m_poTopologyServer->HandleCollisions();
	m_poDataStructure->CommunicateGhosts();
}
void DislocationSimulator::Reset()
{
	if(m_poDynamicsServer != NULL)
	{
		delete m_poDynamicsServer;
	}
	if(m_poTopologyServer != NULL)
	{
		delete m_poTopologyServer;
	}
	if(m_poSurfaceServer != NULL)
	{
		delete m_poSurfaceServer;
	}
	if(m_poBalancingServer != NULL)
	{
		delete m_poBalancingServer;
	}
	if(m_poDataStructure != NULL)
	{
		delete m_poDataStructure;
	}
	m_poDynamicsServer = NULL;
	m_poTopologyServer = NULL;
	m_poSurfaceServer = NULL;
	m_poBalancingServer = NULL;
	m_poDataStructure = NULL;
}
void DislocationSimulator::Initialize()
{
	m_poDataStructure = NULL;
	m_poTopologyServer = NULL;
	m_poDynamicsServer = NULL;
	m_poSurfaceServer = NULL;
	m_poBalancingServer = NULL;
	m_iOutputFileIndex = 0;
}
	
int main(int argc,char** argv)
{
	if(argc != 2)
	{
		printf("error: incorrect number of arguments\n");
		return 1;
	}
	string sInputFileName = argv[1];
	MPI_Init(&argc,&argv);
	int iDomainsCount = 0;
	int iDomainID = 0;
	MPI_Comm_rank(MPI_COMM_WORLD,&iDomainID);
	MPI_Comm_size(MPI_COMM_WORLD,&iDomainsCount);
	DislocationSimulator* poSimulator = DislocationSimulator::GetInstance();
	poSimulator->Start(sInputFileName,iDomainsCount,iDomainID);
	MPI_Barrier(MPI_COMM_WORLD);
    for(int i = 0 ; i < 1000 ; i++)
    {
		poSimulator->Load();
		MPI_Barrier(MPI_COMM_WORLD);
		poSimulator->Move();
		MPI_Barrier(MPI_COMM_WORLD);
		// calculate properties
		MPI_Barrier(MPI_COMM_WORLD);
		// handle cross slip
		MPI_Barrier(MPI_COMM_WORLD);
		// generate output
		MPI_Barrier(MPI_COMM_WORLD);
		poSimulator->Write();
		MPI_Barrier(MPI_COMM_WORLD);
		poSimulator->Collide();
		MPI_Barrier(MPI_COMM_WORLD);
		poSimulator->Remesh();
		MPI_Barrier(MPI_COMM_WORLD);
		// rebalance
		MPI_Barrier(MPI_COMM_WORLD);
		poSimulator->Rebalance();
		MPI_Barrier(MPI_COMM_WORLD);
		// send ghosts
		MPI_Barrier(MPI_COMM_WORLD);
    }
// 		
// 		ParadisCrossSlipServer::Start();
//         TimerStart(poHome,INITIALIZE);
//         Initialize(poHome,sControlFileName,bIsFirst);
//         TimerStop(poHome,INITIALIZE);
//         poHome->cycle = poHome->param->cycleStart;
//         MPI_Barrier(MPI_COMM_WORLD);
//         
//     unsigned int iCycleEnd = home->param->cycleStart + home->param->maxstep;
//     home->iBulkCrossSlipEventsCount = 0;
//     home->iSurfaceCrossSlipEventsCount = 0;
//     home->iAttractiveCrossSlipEventsCount = 0;
//     home->iRepulsiveCrossSlipEventsCount = 0;
//     home->dTotalBulkCrossSlippedChainsLength = 0.0;
//     home->dTotalSurfaceCrossSlippedChainsLength = 0.0;
//     home->dTotalAttractiveCrossSlippedChainsLength = 0.0;
//     home->dTotalRepulsiveCrossSlippedChainsLength = 0.0;
//     
//     // set the surface server
//     ParadisSurfaceServer oSurfaceServer;
//     oSurfaceServer.Set(home);
// 	
// 	while(home->cycle < iCycleEnd)
// 	{      
// 		oSurfaceServer.CheckNodes(home);
// 		ParadisStep(home);
// 		//oSurfaceServer.DetectVirtualNodes(home);
// 		//if(home->cycle%home->param->savecnfreq == 0)
// 		//{
// 		//	oSurfaceServer.RemoveSurfaceArms(home);
// 		//}
// 		if(home->cycle%home->param->CrossSlipHandlingFrequency == 0)
// 		{
// 			printf("bulk cross slip events count : %d\n",home->iBulkCrossSlipEventsCount);
// 			printf("surface cross slip events count : %d\n",home->iSurfaceCrossSlipEventsCount);
// 			printf("attractive cross slip events count : %d\n",home->iAttractiveCrossSlipEventsCount);
// 			printf("repulsive cross slip events count : %d\n",home->iRepulsiveCrossSlipEventsCount);
// 			printf("bulk cross slip total length : %E\n",home->dTotalBulkCrossSlippedChainsLength);
// 			printf("surface cross slip total length : %E\n",home->dTotalSurfaceCrossSlippedChainsLength);
// 			printf("attractive cross slip total length : %E\n",home->dTotalAttractiveCrossSlippedChainsLength);
// 			printf("repulsive cross slip total length : %E\n",home->dTotalRepulsiveCrossSlippedChainsLength);
// 		}
//     }
//     ParadisFinish(home);
	delete poSimulator;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
    return 0;
}



