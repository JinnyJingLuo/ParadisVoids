#ifndef PARADISBLOCKSURFACE_H_
#define PARADISBLOCKSURFACE_H_

#include "ParadisSurface.h"

class ParadisBlockSurface : public ParadisSurface {
public:
  ParadisBlockSurface();
  ParadisBlockSurface(const ParadisBlockSurface &oSurface);
  ~ParadisBlockSurface();
  ParadisBlockSurface &operator=(const ParadisBlockSurface &oSurface);
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
