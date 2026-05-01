// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef GRAPH_H_
#define GRAPH_H_

#include "list"
#include "Randomizer.h"
#include "math.h"
#include "GraphNode.h"
#include "GraphEdge.h"
#include "GraphChain.h"
#include "queue"
#include "limits.h"


using namespace std;
using namespace EZ;

namespace GraphSystem
{
	template <typename NodeType,typename EdgeType> class Graph // to send any class needed 
	{
 	public:
 		Graph()
 		{
			Reset();
 		}
 		~Graph()
 		{
			Reset();
 		}
		Graph(const Graph& oGraph)
		{
			*this = oGraph;
		}
		Graph(list< GraphNode<NodeType,EdgeType>* > lpoNodes,list< GraphEdge<NodeType,EdgeType>* > lpoEdges)
		{
			Set(lpoNodes,lpoEdges);
		}
		Graph& operator=(const Graph& oGraph)
		{
			Set(oGraph.m_lpoNodes,oGraph.m_lpoEdges);
			return *this;
		}
		void Reset()
		{
		 	m_lpoNodes.clear();
 			m_lpoEdges.clear();
		}
		void Collapse()
		{
			ClearNodes();
			ClearEdges();
		}
		void ClearNodes()
		{
			typename list< GraphNode<NodeType,EdgeType>* >::iterator liNodes; // ask
			for(liNodes = m_lpoNodes.begin() ; liNodes != m_lpoNodes.end() ; liNodes++)
			{
				if((*liNodes) != NULL)
				{
					delete (*liNodes);
				}
			}
			m_lpoNodes.clear();
		}
		void ClearEdges()
		{
			typename list< GraphEdge<NodeType,EdgeType>* >::iterator liEdges;
			for(liEdges = m_lpoEdges.begin() ; liEdges != m_lpoEdges.end() ; liEdges++)
			{
				if((*liEdges) != NULL)
				{
					delete (*liEdges);
				}
			}
			m_lpoEdges.clear();
		}
 		void Set(list< GraphNode<NodeType,EdgeType>* > lpoNodes,list< GraphEdge<NodeType,EdgeType>* > lpoEdges)
 		{
 			m_lpoNodes.clear();
 			m_lpoEdges.clear();

			typename list< GraphNode<NodeType,EdgeType>* >::iterator liNodes;
			for(liNodes = lpoNodes.begin() ; liNodes != lpoNodes.end() ; liNodes++)
			{
 				m_lpoNodes.push_back((*liNodes));
			}
			
 			typename list< GraphEdge<NodeType,EdgeType>* >::iterator liEdges;
 			for(liEdges = lpoEdges.begin() ; liEdges != lpoEdges.end() ; liEdges++)
 			{
				m_lpoEdges.push_back((*liEdges));
 			}
 		}
		list< GraphNode<NodeType,EdgeType>* >* GetNodes()
		{
			return &m_lpoNodes;
		}
		list< GraphEdge<NodeType,EdgeType>* >* GetEdges()
		{
			return &m_lpoEdges;
		}
		unsigned int GetNodesCount()
		{
			return (unsigned int)m_lpoNodes.size();
		}
		void AddNode(GraphNode<NodeType,EdgeType>* poNode)
		{
			if(poNode != NULL)
			{
				m_lpoNodes.push_back(poNode);
			}
		}
		void AddEdge(GraphEdge<NodeType,EdgeType>* poEdge)
		{
			if(poEdge != NULL)
			{
				m_lpoEdges.push_back(poEdge);
			}
		}
		void RemoveNode(typename list< GraphNode<NodeType,EdgeType>* >::iterator liNode)
		{
			m_lpoNodes.erase(liNode);
		}
		void RemoveNode(GraphNode<NodeType,EdgeType>* plpoNode,bool bDeleteNode = false)
		{
			typename list< GraphNode<NodeType,EdgeType>* >::iterator liNodes = m_lpoNodes.begin();
			while(liNodes != m_lpoNodes.end())
			{
				if((*liNodes) == plpoNode)
				{
					if(bDeleteNode)
					{
						delete (*liNodes);
					}
					liNodes = m_lpoNodes.erase(liNodes);
				}
				else
				{
					liNodes++;
				}
			}
		}
		void RemoveEdge(typename list< GraphEdge<NodeType,EdgeType>* >::iterator liEdge)
		{
			m_lpoEdges.erase(liEdge);
		}
 		void GenerateRandomUndirectedGraph(const unsigned int& iNodesCount,const double& dDensity)
 		{
 			if(iNodesCount <= 0 || dDensity <= 0.0)
 			{
 				return;
 			}
 			unsigned int iEdgesCount = (unsigned int)floor(0.5*dDensity*(double)(iNodesCount*(iNodesCount - 1)) + 0.5);
 			GenerateRandomGraph(iNodesCount,iEdgesCount,false);
 		}
 		void GenerateRandomDirectedGraph(const unsigned int& iNodesCount,const double& dDensity)
 		{
 			if(iNodesCount <= 0 || dDensity <= 0.0)
 			{
 				return;
 			}
 			unsigned int iEdgesCount = (unsigned int)floor(dDensity*(double)(iNodesCount*(iNodesCount - 1)) + 0.5); //ask check on 0.5 (directed) 
 			GenerateRandomGraph(iNodesCount,iEdgesCount,true);
 		}
 		list< GraphNode<NodeType,EdgeType>* > GetShortestPath(GraphNode<NodeType,EdgeType>* poStartNode,GraphNode<NodeType,EdgeType>* poEndNode)
 		{
 			list< GraphNode<NodeType,EdgeType>* > lpoPath;
 			lpoPath.clear();
 			unsigned int iSize = m_lpoNodes.size();
 			unsigned int iSourceIndex = poStartNode->GetID() - 1;
 			unsigned int iTargetID = poEndNode->GetID();
 			vector<int> viStatus;
 			vector<unsigned int> viDepth;
 			vector< GraphNode<NodeType,EdgeType>* > vpoParent;
 			viStatus.resize(iSize);
 			viDepth.resize(iSize);
 			vpoParent.resize(iSize);
 			unsigned int i = 0;
 			for(i = 0; i < iSize ; i++)
			{
				viStatus[i] = 0;
				viDepth[i] = UINT_MAX;
				vpoParent[i] = NULL;
			}
			viStatus[iSourceIndex] = 1;
			viDepth[iSourceIndex] = 0;
			queue< GraphNode<NodeType,EdgeType>* > qpoLevelNodes;
			qpoLevelNodes.push(poStartNode);
			GraphNode<NodeType,EdgeType>* poParentNode = NULL;
			GraphNode<NodeType,EdgeType>* poChildNode = NULL;
			list< GraphNode<NodeType,EdgeType>* > lpoNeighbours;
			typename list< GraphNode<NodeType,EdgeType>* >::iterator liIterator;
			unsigned int iParentIndex = 0;
			unsigned int iChildIndex = 0;
			bool bIsFound = false;
			while(!qpoLevelNodes.empty())
			{
				lpoNeighbours.clear();
				poParentNode = qpoLevelNodes.front();
				qpoLevelNodes.pop();
				iParentIndex = poParentNode->GetID() - 1;
				lpoNeighbours = poParentNode->GetNeighbours();
				for(liIterator = lpoNeighbours.begin() ; liIterator != lpoNeighbours.end() ; liIterator++)
				{
					poChildNode = *liIterator;
					iChildIndex = poChildNode->GetID() - 1;
					if(viStatus[iChildIndex] == 0)
					{
						viStatus[iChildIndex] = 1;
						viDepth[iChildIndex] = viDepth[iParentIndex] + 1;
						vpoParent[iChildIndex] = poParentNode;
						qpoLevelNodes.push(poChildNode);
						if(poChildNode == poEndNode)
						{
							bIsFound = true;
							break;
						}
					}
				}
				viStatus[iParentIndex] = 2;
				if(bIsFound)
				{
					break;
				}
			}
			if(!bIsFound)
			{
				return lpoPath;
			}
			poChildNode = qpoLevelNodes.back();
			unsigned int iDepth = viDepth[iTargetID - 1] + 1;
			vector< GraphNode<NodeType,EdgeType>* > vpoPath;
			vpoPath.resize(iDepth);
			vpoPath[0] = poChildNode;
			iChildIndex = poChildNode->GetID() - 1;
			for(i = 1; i < iDepth ; i++)
			{
				poParentNode = vpoParent[iChildIndex];
				if(poParentNode != NULL)
				{
					iChildIndex = poParentNode->GetID() - 1;
				}
				vpoPath[i] = poParentNode;
			}

			for(i = 0; i < iDepth ; i++)
			{
				lpoPath.push_back(vpoPath[iDepth - i - 1]);
			}
  			return lpoPath;
 		}
 		Graph<NodeType,EdgeType> GenerateNodeTree(GraphNode<NodeType,EdgeType>* poRootNode)
 		{
 			Graph<NodeType,EdgeType> oTree;
 			unsigned int iSize = m_lpoNodes.size();
 			if(poRootNode == NULL)
 			{
 				return oTree;
 			}
			vector<int> viStatus;
			vector< GraphNode<NodeType,EdgeType>* > vpoParent;
			vector< list< GraphNode<NodeType,EdgeType>* > > vlpoChildren;
			viStatus.resize(iSize);
			vpoParent.resize(iSize);
			vlpoChildren.resize(iSize);
			unsigned int i = 0;
			for(i = 0; i < iSize ; i++)
			{
				viStatus[i] = 0;
				vpoParent[i] = NULL;
				vlpoChildren[i].clear();
			}
			unsigned int iRootIndex = poRootNode->GetID() - 1;
 			viStatus[iRootIndex] = 1;
 			queue< GraphNode<NodeType,EdgeType>* > qpoLevelNodes;
 			qpoLevelNodes.push(poRootNode);
			GraphNode<NodeType,EdgeType>* poParentNode = NULL;
			GraphNode<NodeType,EdgeType>* poChildNode = NULL;
			list< GraphNode<NodeType,EdgeType>* > lpoNeighbours;
			typename list< GraphNode<NodeType,EdgeType>* >::iterator liNodes;
			unsigned int iParentIndex = 0;
			unsigned int iChildIndex = 0;
			while(!qpoLevelNodes.empty())
			{
				lpoNeighbours.clear();
				poParentNode = qpoLevelNodes.front();
				qpoLevelNodes.pop();
				iParentIndex = poParentNode->GetID() - 1;
				lpoNeighbours = poParentNode->GetNeighbours();
				for(liNodes = lpoNeighbours.begin() ; liNodes != lpoNeighbours.end() ; liNodes++)
				{
					poChildNode = *liNodes;
					iChildIndex = poChildNode->GetID() - 1;
					if(viStatus[iChildIndex] == 0)
					{
						viStatus[iChildIndex] = 1;
						vpoParent[iChildIndex] = poParentNode;
						vlpoChildren[iParentIndex].push_back(poChildNode);
						qpoLevelNodes.push(poChildNode);
					}
				}
				viStatus[iParentIndex] = 2;
			}
			list< GraphNode<NodeType,EdgeType>* > lpoNodes;
			list< GraphEdge<NodeType,EdgeType>* > lpoEdges;
			GraphNode<NodeType,EdgeType>* poTempNode = NULL;
			vector< GraphNode<NodeType,EdgeType>* > vpoNodesMap;
			vpoNodesMap.resize(iSize);
			for(liNodes = m_lpoNodes.begin() ; liNodes != m_lpoNodes.end() ; liNodes++)
			{
				i = (*liNodes)->GetID() - 1;
				vpoNodesMap[i] = NULL;
				if(viStatus[i] == 2)
				{
					poTempNode = new GraphNode<NodeType,EdgeType>;
					poTempNode->EquateData((*liNodes));
					poTempNode->SetID((*liNodes)->GetID());
					lpoNodes.push_back(poTempNode);
					vpoNodesMap[i] = poTempNode;
				}
			}
 			GraphEdge<NodeType,EdgeType>* poNewEdge = NULL;
 			GraphEdge<NodeType,EdgeType>* poOldEdge = NULL;
 			typename list< GraphNode<NodeType,EdgeType>* >::iterator liChildNodes;
 			for(liNodes = m_lpoNodes.begin() ; liNodes != m_lpoNodes.end() ; liNodes++)
 			{
 				i = (*liNodes)->GetID() - 1;
 				if(viStatus[i] == 2)
 				{
 					for(liChildNodes = vlpoChildren[i].begin(); liChildNodes != vlpoChildren[i].end() ; liChildNodes++)
 					{
 						if(vpoNodesMap[i]->ConnectTo(vpoNodesMap[(*liChildNodes)->GetID() - 1],poNewEdge,true))
 						{
							if((*liNodes)->IsConnected((*liChildNodes),poOldEdge))
							{
								poNewEdge->EquateData(poOldEdge);
							}
 						}
 						lpoEdges.push_back(poNewEdge);
					}
 				}
 			}
 			oTree.Set(lpoNodes,lpoEdges);
 			unsigned int iNodesCount = m_lpoNodes.size();
 			for(i = 0; i < iNodesCount ; i++)
 			{
 				vlpoChildren[i].clear();
 			}
 			lpoNodes.clear();
 			lpoEdges.clear();
 			viStatus.clear();
 			vpoParent.clear();
 			vlpoChildren.clear();
 			return oTree;
 		}
		void GenerateSequentialNodeIDs()
		{
			unsigned int i = 0;
			typename list< GraphNode<NodeType,EdgeType>* >::iterator liNodes;
			for(liNodes = m_lpoNodes.begin() ; liNodes != m_lpoNodes.end() ; liNodes++)
			{
				i = i + 1;
				(*liNodes)->SetID(i);
			}
		}
		list< GraphChain<NodeType,EdgeType>* > GenerateGraphChains()
		{
			// this function does NOT work if the graph has directed edges
			vector<bool> vbHasBeenExamined;
			vbHasBeenExamined.resize(m_lpoNodes.size());
			typename list< GraphNode<NodeType,EdgeType>* >::iterator liNodes;
			for(liNodes = m_lpoNodes.begin() ; liNodes != m_lpoNodes.end() ; liNodes++)
			{
				vbHasBeenExamined[(*liNodes)->GetID() - 1] = false;
			}
			// start chains detection
			list< GraphChain<NodeType,EdgeType>* > lpoChains;
			GraphChain<NodeType,EdgeType> oPath;
			GraphChain<NodeType,EdgeType>* poTempPath;
			unsigned int iChainLength = 0;
			unsigned int i = 0;
			for(liNodes = m_lpoNodes.begin() ; liNodes != m_lpoNodes.end() ; liNodes++)
			{
				if(!vbHasBeenExamined[(*liNodes)->GetID() - 1])
				{
					oPath = GenerateNodeChain((*liNodes));
					iChainLength = oPath.GetSize();
					if(iChainLength < 2)
					{
						continue;
					}
					poTempPath = new GraphChain<NodeType,EdgeType>;
					oPath.ResetIterator();
					for(i = 0 ; i < iChainLength ; i++)
					{
						vbHasBeenExamined[oPath.GetCurrentNode()->GetID() - 1] = true;
						poTempPath->AddNodeAtEnd(oPath.GetCurrentNode());
						oPath.IncrementIterator();
					}
					lpoChains.push_back(poTempPath);
				}
			}
			// generate stub chains (chains with junction nodes at both ends)
			// reset the examined conditions
			for(liNodes = m_lpoNodes.begin() ; liNodes != m_lpoNodes.end() ; liNodes++)
			{
				vbHasBeenExamined[(*liNodes)->GetID() - 1] = false;
			}
			
			list< GraphNode<NodeType,EdgeType>* > lpoNeighbours;
			typename list< GraphNode<NodeType,EdgeType>* >::iterator liNeighbours;
			for(liNodes = m_lpoNodes.begin() ; liNodes != m_lpoNodes.end() ; liNodes++)
			{
				if(!vbHasBeenExamined[(*liNodes)->GetID() - 1])
				{
					if((*liNodes)->GetAllNeighboursCount() > 2)
					{
						lpoNeighbours = (*liNodes)->GetNeighbours();
						for(liNeighbours = lpoNeighbours.begin() ; liNeighbours != lpoNeighbours.end() ; liNeighbours++)
						{
							if((*liNeighbours)->GetAllNeighboursCount() > 2)
							{
								if(!vbHasBeenExamined[(*liNeighbours)->GetID() - 1])
								{
									// now we have a stub chain, create and add it
									poTempPath = new GraphChain<NodeType,EdgeType>;
									poTempPath->AddNodeAtEnd((*liNodes));
									poTempPath->AddNodeAtEnd((*liNeighbours));
									lpoChains.push_back(poTempPath);
								}
							}
						}
					}
				}
				vbHasBeenExamined[(*liNodes)->GetID() - 1] = true;
			}
			return lpoChains;
		}
		static GraphChain<NodeType,EdgeType> GenerateNodeChain(GraphNode<NodeType,EdgeType>* poNode)
		{
			// this function does NOT work if the graph has directed edges
			GraphChain<NodeType,EdgeType> oChain;
			if(poNode == NULL)
			{
				return oChain;
			}
			oChain.AddNodeAtEnd(poNode);
			unsigned int iRootNodeNeighboursCount = poNode->GetAllNeighboursCount();
			if(iRootNodeNeighboursCount < 1 || iRootNodeNeighboursCount > 2)
			{
				return oChain;
			}
			// now, a proper chain (length > 1) can be formed
			GraphNode<NodeType,EdgeType>* poNodeBeforeLast = poNode;
			list< GraphNode<NodeType,EdgeType>* > lpoNeighbours = poNodeBeforeLast->GetNeighbours();
			GraphNode<NodeType,EdgeType>* poLastNode = NULL;
			// get the nodes on the first side of the initial node
			poLastNode = lpoNeighbours.front();
			oChain.AddNodeAtEnd(poLastNode);
			GraphNode<NodeType,EdgeType>* poTempNode = NULL;
			bool bProceed = true;
			bool bIsLoop = false;
			while(bProceed)
			{
				poTempNode = poLastNode->GetOtherNeighbour(poNodeBeforeLast);
				if(poTempNode == NULL)
				{
					bProceed = false;
				}
				else
				{
					oChain.AddNodeAtEnd(poTempNode);
					poNodeBeforeLast = poLastNode;
					poLastNode = poTempNode;
					if(poTempNode == poNode)
					{
						bProceed = false;
						bIsLoop = true;
					}
				}
			}

			if(iRootNodeNeighboursCount == 2 && !bIsLoop)
			{
				// get the nodes on the second side of the initial node
				poNodeBeforeLast = poNode;
				poLastNode = lpoNeighbours.back();
				oChain.AddNodeAtBeginning(poLastNode);
				poTempNode = NULL;
				bProceed = true;	
				while(bProceed)
				{
					poTempNode = poLastNode->GetOtherNeighbour(poNodeBeforeLast);
					if(poTempNode == NULL || poTempNode == poNode)
					{
						bProceed = false;
					}
					else
					{
						oChain.AddNodeAtBeginning(poTempNode);
						poNodeBeforeLast = poLastNode;
						poLastNode = poTempNode;
						// the following condition should never be satisfied
						if(poTempNode == poNode)
						{
							bProceed = false;
							bIsLoop = true;
						}
					}
				}
			}
			return oChain;
		}
		void RemoveRepeatedEdges()
		{
			typename list< GraphEdge<NodeType,EdgeType>* >::iterator liOuterEdges;
			typename list< GraphEdge<NodeType,EdgeType>* >::iterator liInnerEdges;
			GraphEdge<NodeType,EdgeType>* poEdge = NULL;
			for(liOuterEdges = m_lpoEdges.begin() ; liOuterEdges != m_lpoEdges.end() ; liOuterEdges++)
			{
				poEdge = (*liOuterEdges);
				if(poEdge == NULL)
				{
					continue;
				}
				liInnerEdges = liOuterEdges;
				liInnerEdges++;
				while(liInnerEdges != m_lpoEdges.end())
				{
					if(poEdge->IsSimilar((*liInnerEdges)))
					{
						liInnerEdges = m_lpoEdges.erase(liInnerEdges);
					}
					else
					{
						liInnerEdges++;
					}
				}
			}
		}
		void AbstractChains(const bool& bDeleteObjects,bool bTakeRightNodeEdgeData = true)
		{
			// this function does NOT work if the graph has directed edges
			GraphNode<NodeType,EdgeType>* poNode = NULL;
			GraphNode<NodeType,EdgeType>* poRightNode = NULL;
			GraphNode<NodeType,EdgeType>* poLeftNode = NULL;
			GraphEdge<NodeType,EdgeType>* poRightEdge = NULL;
			GraphEdge<NodeType,EdgeType>* poLeftEdge = NULL;
			GraphEdge<NodeType,EdgeType>* poNewEdge = NULL;
			list< GraphNode<NodeType,EdgeType>* > lpoNeighbours;
			EdgeType oEdgeData;
			typename list< GraphNode<NodeType,EdgeType>* >::iterator liNodes;
			for(liNodes = m_lpoNodes.begin() ; liNodes != m_lpoNodes.end() ; liNodes++)
			{
				poNode = (*liNodes);
				if(poNode->GetAllNeighboursCount() != 2)
				{
					continue;
				}
				lpoNeighbours = poNode->GetNeighbours();
				poRightNode = lpoNeighbours.front();
				poLeftNode = lpoNeighbours.back();
				if((poRightNode->GetAllNeighboursCount() != 2) && (poLeftNode->GetAllNeighboursCount() != 2))
				{
					continue;
				}
				// the two terminal nodes shouldn't be connected
				if(poRightNode->IsConnected(poLeftNode,poNewEdge))
				{
					continue;
				}
				if(!(poNode->IsConnected(poRightNode,poRightEdge)))
				{
					continue;
				}
				if(!(poNode->IsConnected(poLeftNode,poLeftEdge)))
				{
					continue;
				}
				if(bTakeRightNodeEdgeData)
				{
					oEdgeData = poRightEdge->GetData();
				}
				else
				{
					oEdgeData = poLeftEdge->GetData();
				}
				// disconnect the edges and the nodes
				poRightEdge->Isolate();
				poLeftEdge->Isolate();
				// add the new edge
				if(poRightNode->ConnectTo(poLeftNode,poNewEdge))
				{
					poNewEdge->SetData(oEdgeData);
					m_lpoEdges.push_back(poNewEdge);
				}
			}
			// now we are done with isolating the non essential nodes and edges, now, remove them
			// from the abstract chain graph, notice that the abstract chain graph object does NOT
			// own these nodes or edges, hence, we cannot delete them
			liNodes = m_lpoNodes.begin();
			while(liNodes != m_lpoNodes.end())
			{
				if((*liNodes)->IsIsolated())
				{
					if(bDeleteObjects)
					{
						delete (*liNodes);
					}
					liNodes = m_lpoNodes.erase(liNodes);
				}
				else
				{
					liNodes++;
				}
			}
			typename list< GraphEdge<NodeType,EdgeType>* >::iterator liEdges;
			liEdges = m_lpoEdges.begin();
			while(liEdges != m_lpoEdges.end())
			{
				if((*liEdges)->IsIsolated())
				{
					if(bDeleteObjects)
					{
						delete (*liEdges);
					}
					liEdges = m_lpoEdges.erase(liEdges);
				}
				else
				{
					liEdges++;
				}
			}
		}
 	private:
 
 	protected:
 		void GenerateRandomGraph(const unsigned int& iNodesCount,const unsigned int& iEdgesCount,bool bIsDirected)
 		{
 			unsigned int i = 0;
 			vector< GraphNode<NodeType,EdgeType>* > vpoNodes;
 			vpoNodes.resize(iNodesCount);
 			for(i = 0; i < iNodesCount ; i++)
 			{
 				vpoNodes[i] = new GraphNode<NodeType,EdgeType>;
 			}
 			vector< GraphEdge<NodeType,EdgeType>* > vpoEdges;
 			vpoEdges.resize(iEdgesCount);
 			bool bAccepted = false;
 			GraphNode<NodeType,EdgeType>* poFirstNode = NULL;//ask
 			GraphNode<NodeType,EdgeType>* poLastNode = NULL;
 			GraphEdge<NodeType,EdgeType>* poTempEdge = NULL;
 			for(i = 0; i < iEdgesCount ; i++)
 			{
 				bAccepted = false;
 				while(!bAccepted)
 				{
 					poFirstNode = vpoNodes[EZ::Randomizer::RandomInteger(1,iNodesCount) -1];
 					poLastNode = vpoNodes[EZ::Randomizer::RandomInteger(1,iNodesCount) -1];
 					if(poFirstNode == poLastNode)
 					{
 						continue;
 					}
 					if(poFirstNode->ConnectTo(poLastNode,poTempEdge,bIsDirected))
 					{
 						bAccepted = true;
 					}
				}
 				vpoEdges[i] = poTempEdge;
 			}
		
 			unsigned int iFinalNodesCount = vpoNodes.size();
 			unsigned int iFinalEdgesCount = vpoEdges.size();
 			list< GraphNode<NodeType,EdgeType>* > lpoNodes;
 			for(i = 0; i < iFinalNodesCount ; i++)
 			{
 				lpoNodes.push_back(vpoNodes[i]);
 			}
 			list< GraphEdge<NodeType,EdgeType>* > lpoEdges;
 			for(i = 0; i < iFinalEdgesCount ; i++)
 			{
 				lpoEdges.push_back(vpoEdges[i]);
 			}
 			vpoNodes.clear();
 			vpoEdges.clear();
 			Set(lpoNodes,lpoEdges);
 			GenerateSequentialNodeIDs();
 		}
 		list< GraphNode<NodeType,EdgeType>* > m_lpoNodes;
 		list< GraphEdge<NodeType,EdgeType>* > m_lpoEdges;
 	};
}

#endif


