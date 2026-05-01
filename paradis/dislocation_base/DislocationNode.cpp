// Ahmed M. Hussein

#include "DislocationNode.h"

namespace DislocationSystem {
DislocationNode::DislocationNode() { Initialize(); }
DislocationNode::DislocationNode(const double &dX, const double &dY,
                                 const double &dZ) {
  Initialize();
  Set(dX, dY, dZ);
}
DislocationNode::DislocationNode(const DislocationNode &oNode) {
  *this = oNode;
}
DislocationNode::DislocationNode(const GenericNode &oNode) { *this = oNode; }
DislocationNode::DislocationNode(const Point &oNode) { *this = oNode; }
DislocationNode::~DislocationNode() { Reset(); }
DislocationNode &DislocationNode::operator=(const DislocationNode &oNode) {
  GenericNode::operator=(oNode);
  m_oSurfaceNormal = oNode.m_oSurfaceNormal;
  return *this;
}
DislocationNode &DislocationNode::operator=(const GenericNode &oNode) {
  Initialize();
  GenericNode::operator=(oNode);
  return *this;
}
DislocationNode &DislocationNode::operator=(const Point &oNode) {
  Initialize();
  GenericNode::operator=(oNode);
  return *this;
}
void DislocationNode::Reset() {
  GenericNode::Reset();
  m_oSurfaceNormal.Set(0.0, 0.0, 0.0);
}
void DislocationNode::SetSurfaceNormal(const Vector &oNormal) {
  m_oSurfaceNormal = oNormal;
}
Vector DislocationNode::GetSurfaceNormal() const { return m_oSurfaceNormal; }
void DislocationNode::Initialize() {
  GenericNode::Initialize();
  m_oSurfaceNormal.Set(0.0, 0.0, 0.0);
}
} // namespace DislocationSystem
