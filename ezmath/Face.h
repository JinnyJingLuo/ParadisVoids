#ifndef FACE_H_
#define FACE_H_

#include "GeometricComponent.h"
#include "Point.h"
#include "string"

using namespace EZ;
using namespace std;

namespace GeometrySystem
{
	class Face : public GeometricComponent
	{
	public:
		Face();
		Face(const Face& oFace);
		virtual ~Face();
		Face& operator=(const Face& oFace);
		void Reset();
		virtual bool IsVisibleFromPoint(const Point& oPoint) const = 0;
	private:
	
	protected:
		void Initialize();
	};
}

#endif

