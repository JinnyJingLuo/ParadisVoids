// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef GRAPHNODE_H_
#define GRAPHNODE_H_

#include "list"
#include "GraphEdge.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

using namespace std;

namespace GraphSystem
{
	template <typename NodeType,typename EdgeType> class GraphNode
	{
	public:
		GraphNode()
		{
			m_lpoEdges.clear();
			m_iID = 0;
		}
		GraphNode(const GraphNode& oNode)
		{
			*this = oNode;
		}
		~GraphNode()
		{
			m_iID = 0;
			Isolate();
		}
		GraphNode(NodeType oData)
		{
			Set(oData);
		}
		GraphNode(NodeType* poData)
		{
			Set(poData);
		}
		GraphNode& operator=(const GraphNode& oNode)
		{
			m_oData = oNode.m_oData;
			m_lpoEdges = oNode.m_lpoEdges;
			m_iID = oNode.m_iID;
			return *this;
		}
		void Set(NodeType oData)
		{
			m_oData = oData;
		}
		void Set(NodeType* poData)
		{
			m_oData = *poData;
		}
		NodeType GetData() const
		{
			return m_oData;
		}
		NodeType* GetDataPointer()
		{
			return &m_oData;
		}
		void SetID(const unsigned int& iID)
		{
			m_iID = iID;
		}
		unsigned int GetID() const
		{
			return m_iID;
		}
		bool IsConnected(GraphNode<NodeType,EdgeType>* poNode,GraphEdge<NodeType,EdgeType>*& poEdge)
		{
			poEdge = NULL;
			typename list< GraphEdge<NodeType,EdgeType>* >::iterator liIterator;
			for(liIterator = m_lpoEdges.begin() ; liIterator != m_lpoEdges.end() ; liIterator++)
			{
				if((*liIterator)->GetOther(this) == poNode)
				{
					if(!(*liIterator)->IsDirected())
					{
						poEdge = (*liIterator);
						return true;
					}
					else
					{
						if(((*liIterator)->GetFirst() == this) && ((*liIterator)->GetLast() == poNode))
						{
							poEdge = (*liIterator);
							return true;
						}
						return false;
					}
				}
			}
			return false;
		}
		bool ConnectTo(GraphNode<NodeType,EdgeType>* poNode,GraphEdge<NodeType,EdgeType>*& poEdge,bool bIsDirected = false)
		{
			poEdge = NULL;
			if(!IsConnected(poNode,poEdge))
			{
				poEdge = new GraphEdge<NodeType,EdgeType>;
				poEdge->Set(this,poNode);
				if(bIsDirected)
				{
					poEdge->SetDirected();
				}
				poNode->ConnectFrom(poEdge);
				m_lpoEdges.push_back(poEdge);
				return true;
			}
			return false;
		}
		void DisconnectFrom(GraphNode<NodeType,EdgeType>* poNode)
		{
			GraphEdge<NodeType,EdgeType>* poEdge = NULL;
			if(IsConnected(poNode,poEdge))
			{
				typename list< GraphEdge<NodeType,EdgeType>* >::iterator liIterator;
				liIterator = m_lpoEdges.begin();
				while(liIterator != m_lpoEdges.end())
				{
					if(*liIterator == poEdge)
					{
						(*liIterator)->RemoveNode(this);
						liIterator = m_lpoEdges.erase(liIterator);
					}
					else
					{
						liIterator++;
					}
				}
			}
		}
		void Isolate()
		{
			typename list< GraphEdge<NodeType,EdgeType>* >::iterator liIterator;
			for(liIterator = m_lpoEdges.begin() ; liIterator != m_lpoEdges.end() ; liIterator++)
			{
				(*liIterator)->RemoveNode(this);
			}
			m_lpoEdges.clear();
		}
		void RemoveEdge(GraphEdge<NodeType,EdgeType>* poEdge)
		{
			typename list< GraphEdge<NodeType,EdgeType>* >::iterator liIterator;
			liIterator = m_lpoEdges.begin();
			while(liIterator != m_lpoEdges.end())
			{
				if((*liIterator) == poEdge)
				{
					(*liIterator)->RemoveNode(this);
					liIterator = m_lpoEdges.erase(liIterator);
				}
				else
				{
					liIterator++;
				}
			}
		}
		void RemoveEdgeSoftly(GraphEdge<NodeType,EdgeType>* poEdge)
		{
			typename list< GraphEdge<NodeType,EdgeType>* >::iterator liIterator;
			liIterator = m_lpoEdges.begin();
			while(liIterator != m_lpoEdges.end())
			{
				if((*liIterator) == poEdge)
				{
					liIterator = m_lpoEdges.erase(liIterator);
				}
				else
				{
					liIterator++;
				}
			}
		}
		list< GraphEdge<NodeType,EdgeType>* >* GetEdges()
		{
			return &m_lpoEdges;
		}
		list< GraphNode<NodeType,EdgeType>* > GetNeighbours()
		{
			typename list< GraphEdge<NodeType,EdgeType>* >::iterator liIterator;
 			list< GraphNode<NodeType,EdgeType>* > lpoNeighbours;
			GraphNode<NodeType,EdgeType>* poNode = NULL;
			lpoNeighbours.clear();
			for(liIterator = m_lpoEdges.begin() ; liIterator != m_lpoEdges.end() ; liIterator++)
			{
				poNode = (*liIterator)->GetEnd(this);
				if((poNode != NULL) && (poNode != this))
				{
					lpoNeighbours.push_back(poNode);
				}
			}
			lpoNeighbours.sort();
			lpoNeighbours.unique();
			return lpoNeighbours;
		}
		GraphNode<NodeType,EdgeType>* Clone() const
		{
			GraphNode<NodeType,EdgeType>* poClone = new GraphNode<NodeType,EdgeType>(*this);
			return poClone;
		}
		void EquateData(GraphNode<NodeType,EdgeType>* poNode)
		{
			if(poNode == NULL)
			{
				return;
			}
			m_oData = poNode->m_oData;
		}
		unsigned int GetAllNeighboursCount() const
		{
			return m_lpoEdges.size();
		}
		GraphNode<NodeType,EdgeType>* GetOtherNeighbour(GraphNode<NodeType,EdgeType>* poNode)
		{
			if(m_lpoEdges.size() != 2)
			{
				return NULL;
			}
			GraphNode<NodeType,EdgeType>* poFirstEdgeOtherNode = m_lpoEdges.front()->GetOther(this);
			GraphNode<NodeType,EdgeType>* poSecondEdgeOtherNode = m_lpoEdges.back()->GetOther(this);
			if(poFirstEdgeOtherNode == NULL || poSecondEdgeOtherNode == NULL)
			{
				return NULL;
			}
			if(poFirstEdgeOtherNode == poNode)
			{
				return poSecondEdgeOtherNode;
			}
			
			if(poSecondEdgeOtherNode == poNode)
			{
				return poFirstEdgeOtherNode;
			}
			return NULL;
		}
		bool IsIsolated() const
		{
			return m_lpoEdges.empty();
		}
		GraphNode<NodeType,EdgeType>* DuplicateNode()
		{
			// clone without copying the edges
			GraphNode<NodeType,EdgeType>* poNode = new GraphNode<NodeType,EdgeType>;
			poNode->m_oData = m_oData;
			poNode->m_iID = m_iID;
			return poNode;
		}
		void ConnectFrom(GraphEdge<NodeType,EdgeType>* poEdge)
		{
			// no need to check for connection existence because the only
			// place where this function is called from is  the ConnectTo function
			// which makes that check before calling anyway
			m_lpoEdges.push_back(poEdge);
		}
	private:

	protected:
 		NodeType m_oData;
 		list< GraphEdge<NodeType,EdgeType>* > m_lpoEdges;
 		unsigned int m_iID;
	};
}

#endif


