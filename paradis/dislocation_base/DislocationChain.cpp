// Ahmed M. Hussein

#include "DislocationChain.h"

namespace DislocationSystem {
DislocationChain::DislocationChain() { Reset(); }
DislocationChain::DislocationChain(const DislocationChain &oChain) {
  Reset();
  *this = oChain;
}
DislocationChain::DislocationChain(
    GraphChain<DislocationNode, DislocationSegment> *poChain) {
  Reset();
  Set(poChain);
}
DislocationChain::~DislocationChain() { Reset(); }
DislocationChain &DislocationChain::operator=(const DislocationChain &oChain) {
  Reset();
  Set(oChain.m_poChain);
  return *this;
}
void DislocationChain::Reset() {
  m_poChain = NULL;
  m_oBurgersVector.Set(0.0, 0.0, 0.0);
  m_oSlipPlaneNormal.Set(0.0, 0.0, 0.0);
  m_oInitialLineDirection.Set(0.0, 0.0, 0.0);
  m_eType = InactiveChain;
}
void DislocationChain::Collapse() {
  if (m_poChain != NULL) {
    m_poChain->Clear();
    delete m_poChain;
  }
  Reset();
}
void DislocationChain::Set(
    GraphChain<DislocationNode, DislocationSegment> *poChain) {
  m_poChain = poChain;
  SetProperties();
}
DislocationChain *DislocationChain::Split() {
  GraphChain<DislocationNode, DislocationSegment> *poGraphNewChain =
      m_poChain->Split();
  if (poGraphNewChain == NULL) {
    return NULL;
  }
  DislocationChain *poNewDislocationChain = new DislocationChain();
  poNewDislocationChain->Set(poGraphNewChain);
  return poNewDislocationChain;
}
DislocationChain *DislocationChain::SplitEvenly() {
  GraphChain<DislocationNode, DislocationSegment> *poGraphNewChain =
      m_poChain->SplitEvenly();
  if (poGraphNewChain == NULL) {
    return NULL;
  }
  DislocationChain *poNewDislocationChain = new DislocationChain();
  poNewDislocationChain->Set(poGraphNewChain);
  return poNewDislocationChain;
}
void DislocationChain::Reverse() {
  SoftReverse();
  m_oBurgersVector.Reverse();
  m_oSlipPlaneNormal.Reverse();
  m_poChain->ResetIterator();
  DislocationNode *poNode1 = m_poChain->GetCurrentNode()->GetDataPointer();
  m_poChain->IncrementIterator();
  DislocationNode *poNode2 = m_poChain->GetCurrentNode()->GetDataPointer();
  m_poChain->ResetIterator();
  m_oInitialLineDirection.SetByPoints((*poNode1), (*poNode2));
  m_oInitialLineDirection.Normalize();
}
void DislocationChain::SoftReverse() { m_poChain->Reverse(); }
GraphChain<DislocationNode, DislocationSegment> *
DislocationChain::GetChain() const {
  return m_poChain;
}
Vector DislocationChain::GetBurgersVector() const { return m_oBurgersVector; }
Vector DislocationChain::GetSlipPlaneNormal() const {
  return m_oSlipPlaneNormal;
}
Vector DislocationChain::GetInitialLineDirection() const {
  return m_oInitialLineDirection;
}
ChainType DislocationChain::GetType() const { return m_eType; }
unsigned int DislocationChain::GetSize() const { return m_poChain->GetSize(); }
void DislocationChain::OverrideType(ChainType eType) { m_eType = eType; }
DislocationNetworkNode *DislocationChain::GetFirst() const {
  return m_poChain->GetFirst();
}
DislocationNetworkNode *DislocationChain::GetLast() const {
  return m_poChain->GetLast();
}
DislocationNetworkNode *
DislocationChain::GetOtherEnd(DislocationNetworkNode *poNode) const {
  return m_poChain->GetOtherEnd(poNode);
}
bool DislocationChain::IsEndNode(DislocationNetworkNode *poNode) const {
  return m_poChain->IsEndNode(poNode);
}
bool DislocationChain::IsFirstNode(DislocationNetworkNode *poNode) const {
  return (GetFirst() == poNode);
}
bool DislocationChain::IsLastNode(DislocationNetworkNode *poNode) const {
  return (GetLast() == poNode);
}
DislocationNetworkNode *DislocationChain::GetCurrentNode() const {
  return m_poChain->GetCurrentNode();
}
DislocationNetworkNode *DislocationChain::GetPreviousNode() {
  return m_poChain->GetPreviousNode();
}
DislocationNetworkNode *DislocationChain::GetNextNode() {
  return m_poChain->GetNextNode();
}
void DislocationChain::IncrementIterator() { m_poChain->IncrementIterator(); }
void DislocationChain::IncrementIterator(const unsigned int &iCount) {
  m_poChain->IncrementIterator(iCount);
}
void DislocationChain::DecrementIterator() { m_poChain->DecrementIterator(); }
void DislocationChain::DecrementIterator(const unsigned int &iCount) {
  m_poChain->DecrementIterator(iCount);
}
void DislocationChain::ResetIterator() { m_poChain->ResetIterator(); }
bool DislocationChain::IsIteratorAtFirstNode() const {
  return m_poChain->IsIteratorAtFirstNode();
}
bool DislocationChain::IsIteratorAtLastNode() {
  return m_poChain->IsIteratorAtLastNode();
}
bool DislocationChain::IsLoop() const { return m_poChain->IsLoop(); }
void DislocationChain::MoveIteratorToEnd() { m_poChain->MoveIteratorToEnd(); }
DislocationChain *DislocationChain::Clone() {
  return new DislocationChain(*this);
}
void DislocationChain::SetProperties() {
  m_poChain->ResetIterator();
  DislocationNetworkNode *poNode1 = m_poChain->GetCurrentNode();
  DislocationNetworkNode *poNode2 = m_poChain->GetNextNode();
  DislocationNetworkArm *poEdge = NULL;
  if (poNode2 != NULL) {
    if (poNode1->IsConnected(poNode2, poEdge)) {
      m_oSlipPlaneNormal = poEdge->GetDataPointer()->GetSlipPlaneNormal();
      m_oSlipPlaneNormal.Normalize();
      m_oBurgersVector = poEdge->GetDataPointer()->GetBurgersVector();
      m_oInitialLineDirection.SetByPoints((*poNode1->GetDataPointer()),
                                          (*poNode2->GetDataPointer()));
      m_oInitialLineDirection.Normalize();
    }
  }

  poNode1 = m_poChain->GetFirst();
  unsigned int iFirstNodeNeighboursCount = poNode1->GetAllNeighboursCount();
  if (iFirstNodeNeighboursCount > 4) {
    m_eType = InactiveChain;
  }
  poNode2 = m_poChain->GetLast();
  unsigned int iLastNodeNeighboursCount = poNode2->GetAllNeighboursCount();
  if (iLastNodeNeighboursCount > 4) {
    m_eType = InactiveChain;
  } else if (iFirstNodeNeighboursCount == 3 && iLastNodeNeighboursCount == 3) {
    m_eType = JunctionChain;
  } else if (iFirstNodeNeighboursCount > 2 || iLastNodeNeighboursCount > 2) {
    m_eType = IntersectionChain;
  } else if ((poNode1->GetDataPointer()->GetCategory() == SurfaceNode) &&
             (poNode2->GetDataPointer()->GetCategory() == SurfaceNode)) {
    m_eType = InactiveChain;
  } else if ((poNode1->GetDataPointer()->GetCategory() == SurfaceNode) ||
             (poNode2->GetDataPointer()->GetCategory() == SurfaceNode)) {
    m_eType = SurfaceChain;
  } else {
    m_eType = BulkChain;
  }
}
double DislocationChain::GetLength() {
  m_poChain->ResetIterator();
  unsigned int iSize = m_poChain->GetSize();
  unsigned int i = 0;
  DislocationNode *poCurrentNode = NULL;
  DislocationNode *poNextNode = NULL;
  double dLength = 0.0;
  for (i = 0; i < iSize - 1; i++) {
    poCurrentNode = m_poChain->GetCurrentNode()->GetDataPointer();
    poNextNode = m_poChain->GetNextNode()->GetDataPointer();
    if (poCurrentNode == NULL || poNextNode == NULL) {
      break;
    }
    dLength = dLength + poCurrentNode->Distance(*poNextNode);
    m_poChain->IncrementIterator();
  }
  return dLength;
}
DislocationChain *DislocationChain::Replicate() {
  GraphChain<DislocationNode, DislocationSegment> *poGraphChain =
      new GraphChain<DislocationNode, DislocationSegment>();
  unsigned int iSize = m_poChain->GetSize();
  unsigned int i = 0;
  m_poChain->ResetIterator();
  DislocationNode *poNode = NULL;
  for (i = 0; i < iSize; i++) {
    poGraphChain->AddNodeAtEnd(m_poChain->GetCurrentNode());
    m_poChain->IncrementIterator();
  }
  m_poChain->ResetIterator();
  DislocationChain *poChain = new DislocationChain();
  poChain->Set(poGraphChain);
  return poChain;
}
DislocationChain *DislocationChain::Dissociate(const Vector &oBurgersVector) {
  DislocationChain *poNewChain = Replicate();
  poNewChain->OverrideBurgersVector(m_oBurgersVector - oBurgersVector);
  OverrideBurgersVector(oBurgersVector);
  return poNewChain;
}
void DislocationChain::OverrideBurgersVector(const Vector &oBurgersVector) {
  m_oBurgersVector = oBurgersVector;
}
} // namespace DislocationSystem
