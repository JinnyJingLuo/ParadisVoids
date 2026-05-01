//// Ahmed M. Hussein
//
//#ifndef LOOPFILTER_H_
//#define LOOPFILTER_H_
//
//#include "Graph.h"
//#include "DislocationChain.h"
//#include "AxisAlignedBoundingBox.h"
// using namespace DislocationSystem;
// using namespace GeometrySystem;
//
//#define DualGraphNode GraphNode<DislocationChain*,DislocationNetworkNode*>
//#define DualGraphEdge GraphEdge<DislocationChain*,DislocationNetworkNode*>
//
// class LoopFilter
//{
// public:
//	static LoopFilter* CreateInstance();
//	~LoopFilter();
//	void Reset();
//	void Set(const string& sInputFileName,const unsigned int&
// iProcessesCount); 	void ExtractLoops(const string& sOutputFileName);
// void GenerateModifiedRestartFile(const string& sOutputFileName);
//
//	void GenerateDislocationChains();
//	void SeparateInternalChains();
//	void BuildDualGraph();
//	void ExtractChains();
//
// private:
//
// protected:
//	static LoopFilter* m_poLoopFilterInstance;
//	LoopFilter();
//	void Initialize();
//	void WriteLoop(DislocationChain* poChain,FILE* fpFile);
//	bool IsSurfaceLoop(DislocationChain* poChain);
//	void SetBoundingBox(const string& sInputFileName);
//	static void SeekRestartFileToData(FILE* fpFile);
//	void FixLoopBurgersVector(DislocationChain* poChain);
//	void UpdateNodeIndices(const vector<bool>& vbHasNodeSurvived);
//	void ExtractPath(list<DualGraphNode*>* plpoPath);
//	int EncodeConnection(DislocationChain* poCurrentChain,DislocationChain*
// poNextChain) const; 	void WritePath(list<DualGraphNode*>* plpoPath) const;
//
//	void Fire();
//
//	Graph<DislocationNode,DislocationSegment> m_oNetwork;
//	Graph<DislocationChain*,DislocationNetworkNode*> m_oDualGraph;
//	list< DislocationChain* > m_lpoSurvivingChains;
//	AxisAlignedBoundingBox m_oBox;
//};
//
//
//#endif
//
//
