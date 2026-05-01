// Paradis Processor Project
// Ahmed M. Hussein (mailto : am.hussin@gmail.com)
// June 2012

#ifndef DISLOCATIONCHAIN_H_
#define DISLOCATIONCHAIN_H_

#include "GraphChain.h"
#include "DislocationNode.h"
#include "DislocationSegment.h"

#define DislocationNetworkNode                                                 \
  GraphNode<DislocationNode, DislocationSegment> // ask
#define DislocationNetworkArm                                                  \
  GraphEdge<DislocationNode, DislocationSegment> // ask

using namespace GraphSystem;

namespace DislocationSystem {
enum ChainType {
  InactiveChain = -1,
  NonScrewChain = 0,
  BulkChain = 1,
  SurfaceChain = 2,
  IntersectionChain = 3,
  JunctionChain = 4
};

class DislocationChain {
public:
  DislocationChain();
  DislocationChain(const DislocationChain &oChain); // ask
  DislocationChain(GraphChain<DislocationNode, DislocationSegment> *poChain);
  ~DislocationChain();
  DislocationChain &operator=(const DislocationChain &oChain); // ask
  void Reset();
  void Collapse();
  void Set(GraphChain<DislocationNode, DislocationSegment> *poChain);
  DislocationChain *Split();
  DislocationChain *SplitEvenly();
  void Reverse();
  void SoftReverse();
  GraphChain<DislocationNode, DislocationSegment> *GetChain() const;
  Vector GetBurgersVector() const;
  Vector GetSlipPlaneNormal() const;
  Vector GetInitialLineDirection() const;
  ChainType GetType() const;
  unsigned int GetSize() const;
  void OverrideType(ChainType eType);
  DislocationNetworkNode *GetFirst() const;
  DislocationNetworkNode *GetLast() const;
  DislocationNetworkNode *GetOtherEnd(DislocationNetworkNode *poNode) const;
  bool IsEndNode(DislocationNetworkNode *poNode) const;
  bool IsFirstNode(DislocationNetworkNode *poNode) const;
  bool IsLastNode(DislocationNetworkNode *poNode) const;
  DislocationNetworkNode *GetCurrentNode() const;
  DislocationNetworkNode *GetPreviousNode();
  DislocationNetworkNode *GetNextNode();
  void IncrementIterator();
  void IncrementIterator(const unsigned int &iCount);
  void DecrementIterator();
  void DecrementIterator(const unsigned int &iCount);
  void ResetIterator();
  bool IsIteratorAtFirstNode() const;
  bool IsIteratorAtLastNode();
  bool IsLoop() const;
  void MoveIteratorToEnd();
  DislocationChain *Clone();
  double GetLength();
  DislocationChain *Replicate();
  DislocationChain *Dissociate(const Vector &oBurgersVector);
  void OverrideBurgersVector(const Vector &oBurgersVector);

private:
protected:
  void SetProperties();
  GraphChain<DislocationNode, DislocationSegment> *m_poChain;
  Vector m_oBurgersVector;
  Vector m_oSlipPlaneNormal;
  Vector m_oInitialLineDirection;
  ChainType m_eType;
};
} // namespace DislocationSystem

#endif
