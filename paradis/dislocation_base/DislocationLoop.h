// Paradis Processor Project
// Ahmed M. Hussein (mailto : am.hussin@gmail.com)
// June 2012

#ifndef DISLOCATIONLOOP_H_
#define DISLOCATIONLOOP_H_

#include "DislocationChain.h"
#include "CircularLinkedList.h"

namespace DislocationSystem {
class DislocationLoop {
public:
  DislocationLoop();
  DislocationLoop(const DislocationLoop &oLoop);
  ~DislocationLoop();
  DislocationLoop &operator=(const DislocationLoop &oLoop);
  void Reset();
  bool GenerateLoop(DislocationChain *poStartingChain,
                    list<DislocationChain *> *plpoChains);
  bool IsPlanar();
  const list<DislocationNetworkNode *> &GetNodes();
  double GetSolidAngle(const Point &oPoint);
  Vector GetPointDisplacement(const Point &oPoint, const double &dPoissonRatio);
  void Print();
  DislocationNetworkNode *GetCurrentNode();
  DislocationNetworkNode *GetPreviousNode();
  DislocationNetworkNode *GetNextNode();
  void IncrementIterator();
  void DecrementIterator();
  void ResetIterator();
  bool IsAtBeginning() const;
  void CoarsenVirtualNodes();

private:
protected:
  void Initialize();
  void AppendChain(DislocationChain *poChain);
  void DropNode();
  CircularLinkedList<DislocationNetworkNode *> m_loNodesList;
  Vector m_oBurgersVector;
};
} // namespace DislocationSystem

#endif
