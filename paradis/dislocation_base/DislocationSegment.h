// Paradis Processor Project
// Ahmed M. Hussein (mailto : am.hussin@gmail.com)
// June 2012

#ifndef DISLOCATIONSEGMENT_H_
#define DISLOCATIONSEGMENT_H_

#include "Vector.h"

using namespace EZ;

namespace DislocationSystem {
class DislocationSegment {
public:
  DislocationSegment();
  DislocationSegment(const DislocationSegment &oSegment);
  ~DislocationSegment();
  DislocationSegment &operator=(const DislocationSegment &oSegment);
  void Reset();
  void SetSlipPlaneNormal(const Vector &oNormal);
  void SetBurgersVector(const Vector &oBurgersVector);
  Vector GetSlipPlaneNormal() const;
  Vector GetBurgersVector() const;

private:
protected:
  Vector m_oSlipPlaneNormal;
  Vector m_oBurgersVector;
};
} // namespace DislocationSystem

#endif
