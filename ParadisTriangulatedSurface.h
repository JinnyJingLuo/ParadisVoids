#ifndef PARADISTRIANGULATEDSURFACE_H_
#define PARADISTRIANGULATEDSURFACE_H_

#include "Home.h"
#include "ParadisSurface.h"

using namespace std;
using namespace EZ;
using namespace GeometrySystem;

class ParadisTriangulatdSurface : public ParadisSurface {
public:
  ParadisTriangulatdSurface();
  ParadisTriangulatdSurface(const ParadisTriangulatdSurface &oSurface);
  ~ParadisTriangulatdSurface();
  ParadisTriangulatdSurface &
  operator=(const ParadisTriangulatdSurface &oSurface);
  void Reset();
  virtual void Set(Home_t *poHome);

private:
protected:
  void Initialize();
  virtual bool IsPointInside(const Point &oPoint) const;
  virtual bool IsPointOnSurface(const Point &oPoint) const;
  virtual void GenerateTriangulations();
};

#endif
