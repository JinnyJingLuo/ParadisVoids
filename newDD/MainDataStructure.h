#ifndef MAINDATASTRUCTURE_H_
#define MAINDATASTRUCTURE_H_

#include "vector"
#include "list"
#include "map"
#include "DislocationNode.h"
#include "DislocationSegment.h"
#include "AxisAlignedBoundingBox.h"
#include "Matrix.h"

using namespace std;
using namespace DislocationSystem;
using namespace EZ;
using namespace GeometrySystem;

enum CommunicationMessageTag
{
	NullCommunicationMessageTag = 0,
	GhostMessageLengthTag = 1,
	GhostMessageDataTag = 2,
	OperationMessageLengthTag = 3,
	OperationMessageDataTag = 4,
	VelocityMessageLengthTag = 5,
	VelocityMessageDataTag = 6,
	MigratingNodesMessageLengthTag = 7,
	MigratingNodesMessageDataTag = 8,
	TagMapsMessageLengthTag = 9,
	TagMapsMessageDataTag = 10
};

class MainDataStructure
{
public:
	static MainDataStructure* GetInstance();
	~MainDataStructure();
	void SetDomainData(const unsigned int& iDomainsCount,const unsigned int& iDomainID);
	void ReadInput(const string& sFileName);
	void WriteOutput(const string& sFileName);
	void CommunicateGhosts();
	static bool IsTagLower(const unsigned int& iDomainID1,const unsigned int& iNodeID1,const unsigned int& iDomainID2,const unsigned int& iNodeID2);
	list<DislocationNode*>* GetLocalNodes();
	list<DislocationNode*>* GetRemoteNodes();
	DislocationNode* GetNode(DislocationNodeTag* poTag) const;
	unsigned int GetDomainID() const;
	unsigned int GetLocalNodeID();
	unsigned int GetDomainsCount() const;
	void AddNode(DislocationNode* poNode);
	void AddAndRegisterNode(DislocationNode* poNode);
	void SetTime(const double& dValue);
	void SetTimeStep(const double& dValue);
	double GetTime() const;
	double GetTimeStep() const;
	Matrix* GetAppliedStress();
	Matrix* GetTotalStrain();
	Matrix* GetPlasticStrain();
	Matrix* GetPlasticStrainIncrement();
	static bool Is111Vector(Vector oVector);
	static bool Is110Vector(Vector oVector);
	static bool Is112Vector(Vector oVector);
	static bool Is100Vector(Vector oVector);
	AxisAlignedBoundingBox* GetSimulationBox();
	AxisAlignedBoundingBox* GetDomainBox() const;
	vector<AxisAlignedBoundingBox*>* GetDomainBoxes();
	list<DislocationNode*>::iterator RemoveLocalNode(list<DislocationNode*>::iterator liNode);
	list<DislocationNode*>::iterator RemoveRemoteNode(list<DislocationNode*>::iterator liNode);
	void RecycleNodeID(const unsigned int& iNodeID);
	void UpdateConnectivityTags(vector< map<unsigned int,DislocationNodeTag*>* >* pvmpoTagMaps);
	
private:

protected:
	static MainDataStructure* m_poInstance;
	MainDataStructure();
	void Reset();
	void Initialize();
	void SetDomainBoxes();
	void SetLocalNodes(const unsigned int& iMaxDomainID);
	void ClearGhostNodes();
	void PackGhostNodes(vector< list<double> >* pvldData);
	void UnpackGhostNodes(vector< list<double> >* pvldData);
	void RegisterNodes();
	unsigned int m_iDomainsCount;
	unsigned int m_iDomainID;
	unsigned int m_iXDomainsCount;
	unsigned int m_iYDomainsCount;
	unsigned int m_iZDomainsCount;
	unsigned int m_iOutputFrequency;
	double m_dTime;
	double m_dTimeStep;
	AxisAlignedBoundingBox m_oProblemBox;
	AxisAlignedBoundingBox* m_poDomainBox;
	vector<AxisAlignedBoundingBox*> m_vpoDomains;
	list<DislocationNode*> m_lpoLocalNodes;
	list<DislocationNode*> m_lpoRemoteNodes;
	vector< vector<DislocationNode*> > m_vvpoNodes;
	Matrix m_oAppliedStress;
	Matrix m_oTotalStrain;
	Matrix m_oPlasticStrain;
	Matrix m_oPlasticStrainIncrement;
	list<unsigned int> m_liFreeNodeTags;
	unsigned int m_iFirstFreeNodeIndex;
};


#endif


