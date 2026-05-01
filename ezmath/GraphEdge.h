// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef	GRAPHEDGE_H_
#define GRAPHEDGE_H_

#include "GraphNode.h"
#include "stdio.h"

namespace GraphSystem
{
	template <typename NodeType,typename EdgeType> class GraphNode;
	template <typename NodeType,typename EdgeType> class GraphEdge
	{
	public:
		GraphEdge()
		{
			m_poFirst = NULL;
			m_poLast = NULL;
			m_bIsDirected = false;
		}
		~GraphEdge()
		{

		}
		GraphEdge(const GraphEdge& oEdge)
		{
			*this = oEdge;
		}
		GraphEdge(EdgeType oData,GraphNode<NodeType,EdgeType>* poFirst,GraphNode<NodeType,EdgeType>* poLast)
		{
			Set(oData,poFirst,poLast);
		}
		GraphEdge& operator=(const GraphEdge& oEdge)
		{
			m_poFirst = oEdge.m_poFirst;
			m_poLast = oEdge.m_poLast;
			m_bIsDirected = oEdge.m_bIsDirected;
			m_oData = oEdge.m_oData;
			return *this;
		}
		void Set(GraphNode<NodeType,EdgeType>* poFirst,GraphNode<NodeType,EdgeType>* poLast)
		{
			m_poFirst = poFirst;
			m_poLast = poLast;
		}
		void Set(EdgeType oData,GraphNode<NodeType,EdgeType>* poFirst,GraphNode<NodeType,EdgeType>* poLast)
		{
			m_poFirst = poFirst;
			m_poLast = poLast;
			m_oData = oData;
		}
		void SetFirst(GraphNode<NodeType,EdgeType>* poNode)
		{
			m_poFirst = poNode;
		}
		void SetLast(GraphNode<NodeType,EdgeType>* poNode)
		{
			m_poLast = poNode;
		}
		void SetData(EdgeType oData)
		{
			m_oData = oData;
		}
		void SetData(EdgeType* poData)
		{
			m_oData = *poData;
		}
		void EquateData(GraphEdge<NodeType,EdgeType>* poEdge)
		{
			m_oData = poEdge->m_oData;
		}
		void SetDirected()
		{
			m_bIsDirected = true;
		}
		void SetUndirected()
		{
			m_bIsDirected = false;
		}
		GraphNode<NodeType,EdgeType>* GetFirst() const
		{
			return m_poFirst;
		}
		GraphNode<NodeType,EdgeType>* GetLast() const
		{
			return m_poLast;
		}
		bool IsFirst(GraphNode<NodeType,EdgeType>* poNode) const
		{
			if(poNode == m_poFirst)
			{
				return true;
			}
			return false;
		}
		bool IsLast(GraphNode<NodeType,EdgeType>* poNode) const
		{
			if(poNode == m_poLast)
			{
				return true;
			}
			return false;
		}
		EdgeType GetData() const
		{
			return m_oData;
		}
		EdgeType* GetDataPointer()
		{
			return &m_oData;
		}
		bool IsDirected() const
		{
			return m_bIsDirected;
		}
		GraphNode<NodeType,EdgeType>* GetOther(GraphNode<NodeType,EdgeType>* poNode) const
		{
			if(m_poFirst == poNode)
			{
				return m_poLast;
			}
			if(m_poLast == poNode)
			{
				return m_poFirst;
			}
			return NULL;
		}
		GraphNode<NodeType,EdgeType>* GetEnd(GraphNode<NodeType,EdgeType>* poNode) const
		{
			if(m_bIsDirected)
			{
				if(poNode == m_poFirst)
				{
					return m_poLast;
				}
				else
				{
					return NULL;
				}
			}
			else
			{
				return GetOther(poNode);
			}
		}
		bool IsSimilar(GraphEdge<NodeType,EdgeType>* poEdge) const
		{
			if(poEdge == NULL)
			{
				return false;
			}
			if(m_bIsDirected == poEdge->m_bIsDirected)
			{
				if(m_poFirst == poEdge->m_poFirst)
				{
					if(m_poLast == poEdge->m_poLast)
					{
						return true;
					}
				}
				if(!m_bIsDirected)
				{
					if(m_poFirst == poEdge->m_poLast)
					{
						if(m_poLast == poEdge->m_poFirst)
						{
							return true;
						}
					}
				}
			}
			return false;
		}
		bool IsLoop() const
		{
			if(m_poFirst == m_poLast)
			{
				return true;
			}
			return false;
		}
		void RemoveNode(GraphNode<NodeType,EdgeType>* poNode)
		{
			if(m_poFirst == poNode)
			{
				m_poFirst = NULL;
			}
			if(m_poLast == poNode)
			{
				m_poLast = NULL;
			}
		}
		void Isolate()
		{
			m_poFirst->RemoveEdge(this);
			m_poLast->RemoveEdge(this);
			m_poFirst = NULL;
			m_poLast = NULL;
		}
		bool IsIsolated() const
		{
			if((m_poFirst == NULL) && (m_poLast == NULL))
			{
				return true;
			}
			return false;
		}
		bool SwitchConnection(GraphNode<NodeType,EdgeType>* poOldNode,GraphNode<NodeType,EdgeType>* poNewNode)
		{
			// this function switches the connection from the old node to the new one, the edge
			// direction does not matter, will return false if the old node is not connected to 
			// the edge
			if(m_poFirst == poOldNode)
			{
				m_poFirst = poNewNode;
				poOldNode->RemoveEdgeSoftly(this);
				return true;
			}
			if(m_poLast == poOldNode)
			{
				m_poLast = poNewNode;
				poOldNode->RemoveEdgeSoftly(this);
				return true;
			}
			return false;
		}
		bool InsertNode(GraphNode<NodeType,EdgeType>* poNewNode,GraphEdge<NodeType,EdgeType>*& poNewEdge)
		{
			// this function takes the new node and switches the connection from the first node
			// to it, and then connects it to the last node using a new edge that it creates
			// if the new node is one of the end nodes, or if it is already connected to any of them,
			// the function returns false and does nothing
			if(m_poFirst == poNewNode)
			{
				return false;
			}
			if(m_poLast == poNewNode)
			{
				return false;
			}
			GraphEdge<NodeType,EdgeType>* poDummyEdge = NULL;
			if(m_poFirst->IsConnected(poNewNode,poDummyEdge))
			{
				return false;
			}
			if(poNewNode->IsConnected(m_poLast,poDummyEdge))
			{
				return false;
			}

			// the last node must lost this edge
			m_poLast->RemoveEdgeSoftly(this);
			poNewNode->ConnectTo(m_poLast,poNewEdge,m_bIsDirected);
			m_poLast = poNewNode;
			poNewNode->ConnectFrom(this);
			return true;
		}
	private:

	protected:
		GraphNode<NodeType,EdgeType>* m_poFirst;
		GraphNode<NodeType,EdgeType>* m_poLast;
		EdgeType m_oData;
		bool m_bIsDirected;
	};
}

#endif


