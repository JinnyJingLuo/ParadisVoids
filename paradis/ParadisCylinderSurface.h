#ifndef PARADISCYLINDERSURFACE_H_
#define PARADISCYLINDERSURFACE_H_

#include "ParadisSurface.h"
#include "Cylinder.h"

using namespace GeometrySystem;

class ParadisCylinderSurface : public ParadisSurface {
public:
  ParadisCylinderSurface();
  ParadisCylinderSurface(const ParadisCylinderSurface &oSurface);
  ~ParadisCylinderSurface();
  ParadisCylinderSurface &operator=(const ParadisCylinderSurface &oSurface);
  void Reset();
  virtual void Set(Home_t *poHome);

private:
protected:
  void Initialize();
  virtual bool IsPointInside(const Point &oPoint) const;
  virtual bool IsPointOnSurface(const Point &oPoint) const;
  virtual void GenerateTriangulations();

  Cylinder *m_poCylinder;
};

#endif
