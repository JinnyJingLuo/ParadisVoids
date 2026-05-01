#ifndef BALANCINGSERVER_H_
#define BALANCINGSERVER_H_

#include "MainDataStructure.h"
#include "DislocationNodeTagMap.h"

using namespace std;

namespace BalancingSystem
{
	class BalancingServer
	{
	public:
		static BalancingServer* GetInstance();
		~BalancingServer();
		void SetDataStructure(MainDataStructure* poDataStructure);
		void Migrate();
		
	private:
	
	protected:
		static BalancingServer* m_poInstance;
		BalancingServer();
		void Reset();
		void Initialize();
		void PackMigratingNodes(vector< list<double> >* pvldData);
		void CommunicateMigratingNodes();
		void UnpackMigratingNodes(vector< list<double> >* pvldData);
		void FreeTagMaps();
		void CommunicateTagMaps();
		void PackTagMaps(list<unsigned int>* pliMapData);
		void UnpackTagMaps(list<unsigned int>* pliMapData);
		MainDataStructure* m_poDataStructure;
		list<DislocationNodeTagMap*> m_lpoTagMaps;
	};
}


#endif


