// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks.

#ifndef TREENODE_H_
#define TREENODE_H_

#include "list"

using namespace std;

template <typename Type> class TreeNode
{
	public:
		TreeNode()
		{
			Initialize();
		}
		TreeNode(const unsigned int iID)
		{
			Initialize();
			m_iID = iID;
		}
		TreeNode(const TreeNode& oNode)
		{
			*this = oNode;
		}
		~TreeNode()
		{
			Reset();
		}
		TreeNode& operator=(const TreeNode& oNode)
		{
			m_poParent = oNode.m_poParent;
			m_lpoChildren = oNode.m_lpoChildren;
			m_iLevel = oNode.m_iLevel;
			m_oData = oNode.m_oData;
			m_iID = oNode.m_iID;
			return *this;
		}
		void Reset()
		{
			Initialize();
		}
		void SetParent(TreeNode<Type>* poParent)
		{
			m_poParent = poParent;
		}
		TreeNode<Type>* GetParent()
		{
			return m_poParent;
		}
		TreeNode<Type>* AddChild(TreeNode<Type>* poChild)
		{
			if(poChild == NULL)
			{
				return NULL;
			}
			m_lpoChildren.push_back(poChild);
			poChild->m_poParent = this;
			poChild->m_iLevel = m_iLevel + 1;
			return poChild;
		}
		void RemoveChild(TreeNode<Type>* poChild)
		{
			typename list< TreeNode<Type>* >::iterator liChildren = m_lpoChildren.begin();
			while(liChildren != m_lpoChildren.end())
			{
				if((*liChildren) == poChild)
				{
					liChildren = m_lpoChildren.erase(liChildren);
				}
				else
				{
					liChildren++;
				}
			}
		}
		void DeleteChildren()
		{
			typename list< TreeNode<Type>* >::iterator liChildren;
			for(liChildren = m_lpoChildren.begin() ; liChildren != m_lpoChildren.end() ; liChildren++)
			{
				if((*liChildren) != NULL)
				{
					(*liChildren)->DeleteChildren();
					delete (*liChildren);
				}
			}
			m_lpoChildren.clear();
		}
		list< TreeNode<Type>* >* GetChildren()
		{
			return &m_lpoChildren;
		}
		void SetLevel(const unsigned int& iLevel)
		{
			m_iLevel = iLevel;
		}
		unsigned int GetLevel() const
		{
			return m_iLevel;
		}
		void SetData(Type oData)
		{
			m_oData = oData;
		}
		Type GetData() const
		{
			return m_oData;
		}
		Type* GetDataPointer()
		{
			return &m_oData;
		}
		unsigned int GetTreeSize() const
		{
			unsigned int iSize = 1;
			if(!m_lpoChildren.empty())
			{
				typename list< TreeNode<Type>* >::const_iterator liChildren;
				for(liChildren = m_lpoChildren.begin() ; liChildren != m_lpoChildren.end() ; liChildren++)
				{
					iSize = iSize + (*liChildren)->GetTreeSize();
				}
			}
			return iSize;
		}
		void SetID(const unsigned int& iID)
		{
			m_iID = iID;
		}
		unsigned int GetID() const
		{
			return m_iID;
		}
		bool IsLeaf() const
		{
			return (m_lpoChildren.empty());
		}
		unsigned int GetChildrenCount() const
		{
			return (unsigned int)m_lpoChildren.size();
		}

	private:

	protected:
		void Initialize()
		{
			m_poParent = NULL;
			m_lpoChildren.clear();
			m_iLevel = 0;
			m_iID = 0;
		}
		TreeNode<Type>* m_poParent;
		list< TreeNode<Type>* > m_lpoChildren;
		unsigned int m_iLevel;
		Type m_oData;
		unsigned int m_iID;
};

#endif



