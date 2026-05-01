#ifndef SURFACELOOP_H_
#define SURFACELOOP_H_

#include "Vector.h"
#include "vector"
#include "list"
#include "CircularLinkedList.h"

using namespace std;
using namespace EZ;

class SurfaceLoop
{
public:
	SurfaceLoop();
	SurfaceLoop(const SurfaceLoop& oLoop);
	~SurfaceLoop();
	SurfaceLoop& operator=(const SurfaceLoop& oLoop);
	void Reset();
	static void SetVirtualNodeSpacing(const double& dSpacing);
	static void SetPoissonsRatio(const double& dRatio);
	static void SetGaussPointsCount(const unsigned int& iCount);
	void Set(const Point& oNode1,const Point& oNode2,const Vector& oBurgersVector,const Vector& oSlipNormal,const Vector& oSurfaceNormal1,const Vector& oSurfaceNormal2);
	Vector GetPointDisplacement(const Point& oPoint);
	string ToString();

private:

protected:
	void Initialize();
	double GetSolidAngle(const Point& oPoint);

	static unsigned int GaussPointsCount;
	static double VirtualNodeSpacing;
	static double PoissonsRatio;
	static vector<double> GaussPointsLocations;
	static vector<double> GaussPointsWeights;

	CircularLinkedList<Point> m_loNodesList;
	Vector m_oBurgersVector;
};


#endif





