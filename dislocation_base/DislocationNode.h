// Paradis Processor Project
// Ahmed M. Hussein (mailto : am.hussin@gmail.com)
// June 2012

#ifndef DISLOCATIONNODE_H_
#define DISLOCATIONNODE_H_

#include "GenericNode.h"
#include "Vector.h"

using namespace EZ;

namespace DislocationSystem {
enum NodeCategory {
  UnconstrainedNode = 0,
  PinnedNode = 7,
  SurfaceNode = 8,
  VirtualNode = 10
};

class DislocationNode : public GenericNode {
public:
  DislocationNode();
  DislocationNode(const double &dX, const double &dY, const double &dZ);
  DislocationNode(const DislocationNode &oNode);
  DislocationNode(const GenericNode &oNode);
  DislocationNode(const Point &oNode);
  virtual ~DislocationNode();
  DislocationNode &operator=(const DislocationNode &oNode);
  DislocationNode &operator=(const GenericNode &oNode);
  DislocationNode &operator=(const Point &oPoint);
  virtual void Reset();
  void SetSurfaceNormal(const Vector &oNormal);
  Vector GetSurfaceNormal() const;

private:
protected:
  void Initialize();
  Vector m_oSurfaceNormal;
};
} // namespace DislocationSystem

#endif
