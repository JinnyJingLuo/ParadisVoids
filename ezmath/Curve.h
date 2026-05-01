#ifndef CURVE_H_
#define CURVE_H_

#include "Segment.h"
#include "CartesianOrthogonalCoordinateSystem.h"
#include "list"

using namespace std;

namespace GeometrySystem
{
	class Curve
	{
	public:
		Curve();
		Curve(const Curve& oCurve);
		~Curve();
		Curve& operator=(const Curve& oCurve);
		void Reset();
		void AddSegment(const Point& oStartPoint,const Point& oEndPoint);
		void FrontPushSegment(const Point& oStartPoint,const Point& oEndPoint);
		void WriteVTK(const string sFileName) const;
		void CyclicSort();
		void Localize(CartesianOrthogonalCoordinateSystem* poSystem);
		void Globalize(CartesianOrthogonalCoordinateSystem* poSystem);
		void Expand(const Vector& oPlaneNormal,const double& dDistance);
		list<Point> GetIntersectionPoints(const Segment& oSegment) const;
		list<Point> GetIntersectionPoints(const Curve& oCurve) const;
		void Split(const Point& oPoint);
		Curve* ExtractSubCurve(const Point& oStartPoint,const Point& oEndPoint);
		GenericNode* GetStartPoint() const;
		GenericNode* GetEndPoint() const;
		list<GenericNode*> GetPoints() const;
		void Flip();
		
	private:
	
	protected:
		void Initialize();
		list<Segment*> m_lpoSegments;
	};
}
#endif

