// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks.

#ifndef FMMTREENODE_H_
#define FMMTREENODE_H_

#include "list"
#include "AxisAlignedBoundingBox.h"

using namespace std;
using namespace GeometrySystem;

namespace FMM
{
	class FMMTreeNode
	{
	public:
		FMMTreeNode();
		FMMTreeNode(const FMMTreeNode& oNode);
		~FMMTreeNode();
		FMMTreeNode& operator=(const FMMTreeNode& oNode);
		void Reset(const bool& bResetBox = true);
		void SetBox(const double& dXRange,const double& dYRange,const double& dZRange);
		void SetBox(const double& dXMin,const double& dXMax,const double& dYMin,const double& dYMax,const double& dZMin,const double& dZMax);
		void Build(const unsigned int& iDepth,const unsigned int& iStartLevel = 0);
		void SetIDs(unsigned int& iStartID);
		void SetParent(FMMTreeNode* poParent);
		const AxisAlignedBoundingBox* GetBox() const;
		void ClaimSources(list< Point* >* plpoSources);
		void ClaimTargets(list< Point* >* plpoTargets);
		
	private:

	protected:
		void Initialize();
		AxisAlignedBoundingBox m_oBox;
		unsigned int m_iLevel;
		unsigned int m_iID;
		FMMTreeNode* m_poParent;
		list< FMMTreeNode* > m_lpoChildren;
		list< Point* > m_lpoSources;
		list< Point* > m_lpoTargets;
	};
}


#endif




