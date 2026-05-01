//Ahmed Hussein
//ahussei4@jhu.edu

#ifndef BSPTREENODE_H_
#define BSPTREENODE_H_

#include "list"
#include "Plane.h"
#include "Tools.h"

using namespace std;
using namespace SupportSystem;

template <typename Type> class BSPTreeNode
{
	public:
		BSPTreeNode()
		{
			Initialize();
		}
		BSPTreeNode(const BSPTreeNode& oNode)
		{
			*this = oNode;
		}
		~BSPTreeNode()
		{
			Collapse();
		}
		BSPTreeNode& operator=(const BSPTreeNode& oNode)
		{
			m_poLeftChild = oNode.m_poLeftChild;
			m_poRightChild = oNode.m_poRightChild;
			m_loItems = oNode.m_loItems;
			m_dSeparationValue = oNode.m_dSeparationValue;
			m_iDimension = oNode.m_iDimension;
			return *this;
		}
		void Collapse()
		{
			if(m_poLeftChild != NULL)
			{
				m_poLeftChild->Collapse();
				delete m_poLeftChild;
			}
			if(m_poRightChild != NULL)
			{
				m_poRightChild->Collapse();
				delete m_poRightChild;
			}
			m_poLeftChild = NULL;
			m_poRightChild = NULL;
			m_loItems.clear();
			m_dSeparationValue = 0.0;
			m_iDimension = 3;
		}
		void Initialize()
		{
			m_poLeftChild = NULL;
			m_poRightChild = NULL;
			m_loItems.clear();
			m_dSeparationValue = 0.0;
			m_iDimension = 3;
		}
		void SetLeftChild(BSPTreeNode<Type>* poChild)
		{
			if(poChild == NULL)
			{
				return;
			}
			m_poLeftChild = poChild;
		}
		void SetRightChild(BSPTreeNode<Type>* poChild)
		{
			if(poChild == NULL)
			{
				return;
			}
			m_poRightChild = poChild;
		}
		BSPTreeNode<Type>* GetLeftChild() const
		{
			return m_poLeftChild;
		}
		BSPTreeNode<Type>* GetRightChild() const
		{
			return m_poRightChild;
		}
		void AddItem(Type oItem)
		{
			m_loItems.push_back(oItem);
		}
		list<Type>* GetItems()
		{
			return &m_loItems;
		}
		void SetSeparator(vector<Point>* pvoPoints,const unsigned int& iDimensionIndex)
		{
			vector<double> vdData;
			unsigned int iSize = (unsigned int)pvoPoints->size();
			vdData.resize(iSize);
			unsigned int i = 0;
			for(i = 0 ; i < iSize ; i++)
			{
				vdData[i] = pvoPoints->at(i).GetComponent(iDimensionIndex);
			}
			QuickSort(vdData);
			vector<Point> voPoints;
			voPoints.resize(iSize);
			if(iSize%2 == 0)
			{
				m_dSeparationValue = 0.5*(vdData[iSize/2] + vdData[iSize/2 + 1]);
			}
			else
			{
				m_dSeparationValue = vdData[iSize/2 + 1];
			}
			m_iDimension = iDimensionIndex;
		}
		int Classify(const Point& oPoint)
		{
			return Classify((Point*)&oPoint);
		}
		int Classify(Point* oPoint)
		{
			double dValue = oPoint->GetComponent(m_iDimension);
			if(dValue > m_dSeparationValue)
			{
				return 1;
			}
			else if(dValue < m_dSeparationValue)
			{
				return -1;
			}
			return 0;
		}
		int ClassifyBox(AxisAlignedBoundingBox* poBox)
		{
			double dXMin = poBox->GetXMin();
			double dXMax = poBox->GetXMax();
			double dYMin = poBox->GetYMin();
			double dYMax = poBox->GetYMax();
			double dZMin = poBox->GetZMin();
			double dZMax = poBox->GetZMax();
			int iSum = 0;
			iSum = iSum + Classify(Point(dXMin,dYMin,dZMin));
			iSum = iSum + Classify(Point(dXMax,dYMin,dZMin));
			iSum = iSum + Classify(Point(dXMax,dYMax,dZMin));
			iSum = iSum + Classify(Point(dXMin,dYMax,dZMin));
			iSum = iSum + Classify(Point(dXMin,dYMin,dZMax));
			iSum = iSum + Classify(Point(dXMax,dYMin,dZMax));
			iSum = iSum + Classify(Point(dXMax,dYMax,dZMax));
			iSum = iSum + Classify(Point(dXMin,dYMax,dZMax));
			
			if(iSum == 8)
			{
				return 1;
			}
			else if(iSum == -8)
			{
				return -1;
			}
			return 0;
		}
		bool IsLeaf()
		{
			if((m_poLeftChild == NULL) && (m_poRightChild == NULL))
			{
				return true;
			}
			return false;
		}
		double GetSeparator()
		{
			return m_dSeparationValue;
		}
		unsigned int GetDimension()
		{
			return m_iDimension;
		}
		
	private:
	
	protected:
		BSPTreeNode<Type>* m_poLeftChild;
		BSPTreeNode<Type>* m_poRightChild;
		list<Type> m_loItems;
		double m_dSeparationValue;
		unsigned int m_iDimension;
};

#endif



