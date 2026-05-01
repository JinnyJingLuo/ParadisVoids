// Ahmed M. Hussein

#include "DislocationSegment.h"

namespace DislocationSystem {
DislocationSegment::DislocationSegment() { Reset(); }
DislocationSegment::DislocationSegment(const DislocationSegment &oSegment) {
  *this = oSegment;
}
DislocationSegment::~DislocationSegment() {}
DislocationSegment &DislocationSegment::
operator=(const DislocationSegment &oSegment) {
  m_oSlipPlaneNormal = oSegment.m_oSlipPlaneNormal;
  m_oBurgersVector = oSegment.m_oBurgersVector;
  return *this;
}
void DislocationSegment::Reset() {
  m_oSlipPlaneNormal.Set(0.0, 0.0, 0.0);
  m_oBurgersVector.Set(0.0, 0.0, 0.0);
}
void DislocationSegment::SetSlipPlaneNormal(const Vector &oNormal) {
  m_oSlipPlaneNormal = oNormal;
  m_oSlipPlaneNormal.Normalize();
}
void DislocationSegment::SetBurgersVector(const Vector &oBurgersVector) {
  m_oBurgersVector = oBurgersVector;
}
Vector DislocationSegment::GetSlipPlaneNormal() const {
  return m_oSlipPlaneNormal;
}
Vector DislocationSegment::GetBurgersVector() const { return m_oBurgersVector; }
} // namespace DislocationSystem
