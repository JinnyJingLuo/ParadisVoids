// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef GRAPHCHAIN_H_
#define GRAPHCHAIN_H_

#include "stdio.h"
#include "GraphNode.h"

namespace GraphSystem
{
	template <typename NodeType,typename EdgeType> class GraphChain //ask
	{
	public:
		GraphChain()
		{
			m_lpoChain.clear();
			ResetIterator();
		}
		GraphChain(const GraphChain& oChain)
		{
			*this = oChain; //ask
		}
		~GraphChain()
		{
			m_lpoChain.clear();
			ResetIterator();
		}
		GraphChain& operator=(const GraphChain& oChain)
		{
			m_lpoChain = oChain.m_lpoChain; //ask
			ResetIterator();
			return *this;
		}
		void AddNodeAtBeginning(GraphNode<NodeType,EdgeType>* poNode)
		{
			m_lpoChain.push_front(poNode);
		}
		void AddNodeAtEnd(GraphNode<NodeType,EdgeType>* poNode)
		{
			m_lpoChain.push_back(poNode);
		}
		void RemoveNodeFromBeginning()
		{
			m_lpoChain.pop_front();
		}
		void RemoveNodeFromEnd()
		{
			m_lpoChain.pop_back();
		}
		void Clear()
		{
			m_lpoChain.clear();
			ResetIterator();
		}
		unsigned int GetSize() const
		{
			return (unsigned int)m_lpoChain.size();
		}
		GraphNode<NodeType,EdgeType>* GetCurrentNode() const
		{
			return (*m_liIterator);
		}
		GraphNode<NodeType,EdgeType>* GetPreviousNode()
		{
			GraphNode<NodeType,EdgeType>* poNode = NULL;
			if(m_liIterator != m_lpoChain.begin())
			{
				DecrementIterator();
				poNode = (*m_liIterator);
				IncrementIterator();
			}
			return poNode;
		}
		GraphNode<NodeType,EdgeType>* GetNextNode()
		{
			GraphNode<NodeType,EdgeType>* poNode = NULL;
			m_liIterator++;
			if(m_liIterator != m_lpoChain.end())
			{
				poNode = (*m_liIterator);
			}
			m_liIterator--;
			return poNode;
		}
		void IncrementIterator()
		{
			m_liIterator++;
			if(m_liIterator == m_lpoChain.end())
			{
				m_liIterator--;
			}
		}
		void IncrementIterator(const unsigned int& iCount)
		{
			unsigned int i = 0;
			for(i = 0 ; i < iCount ; i++)
			{
				IncrementIterator();
			}
		}
		void DecrementIterator()
		{
			if(m_liIterator != m_lpoChain.begin())
			{
				m_liIterator--;
			}
		}
		void DecrementIterator(const unsigned int& iCount)
		{
			unsigned int i = 0;
			for(i = 0 ; i < iCount ; i++)
			{
				DecrementIterator();
			}
		}
		void ResetIterator()
		{
			m_liIterator = m_lpoChain.begin();
		}
		GraphNode<NodeType,EdgeType>* GetFirst() const
		{
			return m_lpoChain.front();
		}
		GraphNode<NodeType,EdgeType>* GetLast() const
		{
			return m_lpoChain.back();
		}
		bool IsIteratorAtFirstNode() const
		{
			if(m_liIterator == m_lpoChain.begin())
			{
				return true;
			}
			return false;
		}
		bool IsIteratorAtLastNode()
		{
			m_liIterator++;
			if(m_liIterator == m_lpoChain.end())
			{
				m_liIterator--;
				return true;
			}
			else
			{
				m_liIterator--;
				return false;
			}
		}
		GraphChain<NodeType,EdgeType>* Split()
		{
			GraphChain<NodeType,EdgeType>* poNewChain = NULL;
			if(IsIteratorAtFirstNode() || IsIteratorAtLastNode())
			{
				return poNewChain;
			}
			// add the split node (should be included in both chains)
			poNewChain = new GraphChain<NodeType,EdgeType>;
			poNewChain->AddNodeAtEnd((*m_liIterator));
			m_liIterator++;
			while(m_liIterator != m_lpoChain.end())
			{
				poNewChain->AddNodeAtEnd((*m_liIterator));
				m_liIterator = m_lpoChain.erase(m_liIterator);
			}
			ResetIterator();
			return poNewChain;
		}
		GraphChain<NodeType,EdgeType>* SplitEvenly()
		{
			unsigned int iSplitNodeIndex = GetMiddleNodeIndex();
			unsigned int iChainSize = GetSize();
			unsigned int i = 0;
			GraphChain<NodeType,EdgeType>* poNewChain = NULL;
			ResetIterator();
			for(i = 0 ; i < iChainSize ; i++)
			{
				if(i == iSplitNodeIndex)
				{
					poNewChain = Split();
					break;
				}
				IncrementIterator();
			}
			ResetIterator();
			return poNewChain;
		}
		void RemoveCurrentNode()
		{
			if(IsIteratorAtFirstNode())
			{
				RemoveNodeFromBeginning();
				ResetIterator();
				return;
			}
			if(IsIteratorAtLastNode())
			{
				RemoveNodeFromEnd();
				m_liIterator = m_lpoChain.end();
				m_liIterator--;
				return;
			}
			m_liIterator = m_lpoChain.erase(m_liIterator);
		}
		unsigned int GetMiddleNodeIndex() const
		{
			unsigned int iChainSize = GetSize();
			unsigned int iMiddleNodeIndex = 0;
			if(iChainSize%2 == 0)
			{
				iMiddleNodeIndex = iChainSize/2 + 1;
			}
			else
			{
				iMiddleNodeIndex = (iChainSize + 1)/2;
			}
			iMiddleNodeIndex = iMiddleNodeIndex - 1;
			return iMiddleNodeIndex;
		}
		GraphNode<NodeType,EdgeType>* GetMiddleNode()
		{
			unsigned int iMiddleNodeIndex = GetMiddleNodeIndex();
			unsigned int iChainSize = GetSize();
			unsigned int i = 0;
			typename list< GraphNode<NodeType,EdgeType>* >::iterator liTempIterator = m_lpoChain.begin();
			for(i = 0 ; i < iChainSize ; i++)
			{
				if(i == iMiddleNodeIndex)
				{
					return (*liTempIterator);
				}
				liTempIterator++;
			}
			return NULL;
		}
		bool IsLoop() const
		{
			if(GetFirst() == GetLast())
			{
				return true;
			}
			return false;
		}
		void MoveIteratorToEnd()
		{
			m_liIterator = m_lpoChain.end();
			m_liIterator--;
		}
		bool PreStitch(GraphChain<NodeType,EdgeType>* poChain)
		{
			if(IsStitchable(poChain,true))
			{
				return false;
			}
			poChain->MoveIteratorToEnd();
			unsigned int i = 0;
			unsigned int iChainSize = poChain->GetSize();
			for(i = 0 ; i < iChainSize ; i++)
			{
				AddNodeAtBeginning(poChain->GetCurrentNode());
				poChain->DecrementIterator();
			}
		}
		bool PostStitch(GraphChain<NodeType,EdgeType>* poChain)
		{
			if(IsStitchable(poChain,false))
			{
				return false;
			}
			poChain->ResetIterator();
			unsigned int i = 0;
			unsigned int iChainSize = poChain->GetSize();
			for(i = 0 ; i < iChainSize ; i++)
			{
				AddNodeAtEnd(poChain->GetCurrentNode());
				poChain->IncrementIterator();
			}
		}
		bool IsStitchable(GraphChain<NodeType,EdgeType>* poChain,bool bIsPreStitching) const
		{
			poChain->ResetIterator();
			unsigned int i = 0;
			unsigned int j = 0;
			unsigned int iRemoteChainSize = poChain->GetSize();
			unsigned int iLocalChainSize = GetSize();
			GraphNode<NodeType,EdgeType>* poRemoteNode = NULL;
			GraphNode<NodeType,EdgeType>* poLocalNode = NULL;
			GraphNode<NodeType,EdgeType>* poRemoteFirstNode = poChain->GetFirst();
			GraphNode<NodeType,EdgeType>* poRemoteLastNode = poChain->GetLast();
			GraphNode<NodeType,EdgeType>* poLocalFirstNode = GetFirst();
			GraphNode<NodeType,EdgeType>* poLocalLastNode = GetLast();
			for(i = 0 ; i < iRemoteChainSize ; i++)
			{
				poRemoteNode = poChain->GetCurrentNode();
				if(poRemoteNode == NULL)
				{
					return false;
				}
				ResetIterator();
				for(j = 0 ; j < iLocalChainSize ; j++)
				{
					poLocalNode = GetCurrentNode();
					if(poRemoteNode == poLocalNode)
					{
						if((poRemoteNode == poRemoteFirstNode) && (poLocalNode == poLocalLastNode) && (!bIsPreStitching))
						{
							continue;
						}
						if((poRemoteNode == poRemoteLastNode) && (poLocalNode == poLocalFirstNode) && bIsPreStitching)
						{
							continue;
						}
						return false;
					}
					IncrementIterator();
				}
				poChain->IncrementIterator();
			}
		}
		GraphChain<NodeType,EdgeType>* Clone()
 		{
 			return new GraphChain<NodeType,EdgeType>(*this);
 		}
		bool IsEndNode(GraphNode<NodeType,EdgeType>* poNode) const
		{
			if(poNode == m_lpoChain.front())
			{
				return true;
			}
			else if(poNode == m_lpoChain.back())
			{
				return true;
			}
			return false;
		}
		GraphNode<NodeType,EdgeType>* GetOtherEnd(GraphNode<NodeType,EdgeType>* poNode) const
		{
			if(poNode == m_lpoChain.front())
			{
				return m_lpoChain.back();
			}
			else if(poNode == m_lpoChain.back())
			{
				return m_lpoChain.front();
			}
			return NULL;
		}
		void Reverse()
		{
			m_lpoChain.reverse();
		}
	private:
	
	protected:
		list< GraphNode<NodeType,EdgeType>* > m_lpoChain;
		typename list< GraphNode<NodeType,EdgeType>* >::iterator m_liIterator;
	};
}

#endif


