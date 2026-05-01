#include "SurfaceLoop.h"
#include "MathServices.h"
#include "cmath"


unsigned int SurfaceLoop::GaussPointsCount = 0;
double SurfaceLoop::VirtualNodeSpacing = 0.0;
double SurfaceLoop::PoissonsRatio = 0.0;
vector<double> SurfaceLoop::GaussPointsLocations;
vector<double> SurfaceLoop::GaussPointsWeights;

SurfaceLoop::SurfaceLoop()
{
	Initialize();
}
SurfaceLoop::SurfaceLoop(const SurfaceLoop& oLoop)
{
	*this = oLoop;
}
SurfaceLoop::~SurfaceLoop()
{
	Reset();
}
SurfaceLoop& SurfaceLoop::operator=(const SurfaceLoop& oLoop)
{
	m_loNodesList = oLoop.m_loNodesList;
	m_oBurgersVector = oLoop.m_oBurgersVector;
	return *this;
}
void SurfaceLoop::Reset()
{
	m_loNodesList.Reset();
}
void SurfaceLoop::SetVirtualNodeSpacing(const double& dSpacing)
{
	VirtualNodeSpacing = dSpacing;
}
void SurfaceLoop::SetPoissonsRatio(const double& dRatio)
{
	PoissonsRatio = dRatio;
}
void SurfaceLoop::SetGaussPointsCount(const unsigned int& iCount)
{
	GaussPointsLocations.clear();
	GaussPointsWeights.clear();
	GaussPointsCount = iCount;
	MathServices::GenerateGaussPoints(GaussPointsLocations,GaussPointsWeights,GaussPointsCount);
}
void SurfaceLoop::Set(const Point& oNode1,const Point& oNode2,const Vector& oBurgersVector,const Vector& oSlipNormal,const Vector& oSurfaceNormal1,const Vector& oSurfaceNormal2)
{
	m_loNodesList.Reset();
	// the origin is always zero
	Point oOrigin(0.0,0.0,0.0);
	// the original nodes are the dislocation nodes
	Point oFirstOriginalNode = oNode1;
	Point oSecondOriginalNode = oNode2;

	// the virtual nodes must be projected in the surface normal direction but they have to stay on their slip plane
	m_oBurgersVector = oBurgersVector;

	Vector oPlaneNormal = Vector(oNode1)^Vector(oNode2);
	oPlaneNormal.Normalize();
	
	Vector oDisplacement = oSurfaceNormal1*VirtualNodeSpacing;
	double dDropHeight = oDisplacement*oPlaneNormal;
	oDisplacement = oDisplacement - oPlaneNormal*dDropHeight;
	Point oFirstVirtualNode = oFirstOriginalNode + oDisplacement;

	oDisplacement = oSurfaceNormal2*VirtualNodeSpacing;
	dDropHeight = oDisplacement*oPlaneNormal;
	oDisplacement = oDisplacement - oPlaneNormal*dDropHeight;
	Point oSecondVirtualNode = oSecondOriginalNode + oDisplacement;

	m_loNodesList.Append(oOrigin);
	m_loNodesList.Append(oFirstOriginalNode);
	m_loNodesList.Append(oFirstVirtualNode);
	m_loNodesList.Append(oSecondVirtualNode);
	m_loNodesList.Append(oSecondOriginalNode);
}
Vector SurfaceLoop::GetPointDisplacement(const Point& oPoint)
{
	double dSolidAngle = GetSolidAngle(oPoint);
	Point oThisNode;
	Point oNextNode;
	double dTolerance = 3.0E0;
	unsigned int i = 0;
	Vector oTangent;
	Point oFieldPoint;
	double dR = 0.0;
	double dRxx = 0.0;
	double dRxy = 0.0;
	double dRxz = 0.0;
	double dRyx = 0.0;
	double dRyy = 0.0;
	double dRyz = 0.0;
	double dRzx = 0.0;
	double dRzy = 0.0;
	double dRzz = 0.0;
	double dBx = m_oBurgersVector.GetX();
	double dBy = m_oBurgersVector.GetY();
	double dBz = m_oBurgersVector.GetZ();
	double dTx = 0.0;
	double dTy = 0.0;
	double dTz = 0.0;
	double dUx = 0.0;
	double dUy = 0.0;
	double dUz = 0.0;
	double dTemp = 0.0;
	double dJacobian = 0.0;
	double dMaterialParameter = 0.5/(1.0 - PoissonsRatio);
	Vector oTemp;
	m_loNodesList.ResetIterator();
	do
	{
		oThisNode = m_loNodesList.GetCurrentItem();
		oNextNode = m_loNodesList.GetNextItem();
		oTangent.SetByPoints(oThisNode,oNextNode);
		oTemp = oTangent;
		oTemp.Normalize();
		dTx = oTemp.GetX();
		dTy = oTemp.GetY();
		dTz = oTemp.GetZ();
		dJacobian = 0.5*oTangent.Length();
		for(i = 0 ; i < GaussPointsCount ; i++)
		{
			dTemp = 0.5*(GaussPointsLocations[i] + 1.0);
			oFieldPoint = oThisNode + oTangent*dTemp;
			dR = oPoint.Distance(oFieldPoint);
			if(dR < dTolerance)
			{
				continue;
			}
			dRxx = oPoint.GetRXX(oFieldPoint);
			dRxy = oPoint.GetRXY(oFieldPoint);
			dRxz = oPoint.GetRXZ(oFieldPoint);
			dRyx = oPoint.GetRYX(oFieldPoint);
			dRyy = oPoint.GetRYY(oFieldPoint);
			dRyz = oPoint.GetRYZ(oFieldPoint);
			dRzx = oPoint.GetRZX(oFieldPoint);
			dRzy = oPoint.GetRZY(oFieldPoint);
			dRzz = oPoint.GetRZZ(oFieldPoint);

			dTemp = dBz*dRyx*dTx - dBz*dRxx*dTy + dBy*dRxx*dTz - dBy*dRzx*dTx + dBx*dRzx*dTy - dBx*dRyx*dTz;
			dTemp = (dBz*dTy - dBy*dTz)/dR + dTemp*dMaterialParameter;
			dUx = dUx + dTemp*dJacobian*GaussPointsWeights[i];

			dTemp = dBz*dRyy*dTx - dBz*dRxy*dTy + dBy*dRxy*dTz - dBy*dRzy*dTx + dBx*dRzy*dTy - dBx*dRyy*dTz;
			dTemp = (dBx*dTz - dBz*dTx)/dR + dTemp*dMaterialParameter;
			dUy = dUy + dTemp*dJacobian*GaussPointsWeights[i];

			dTemp = dBz*dRyz*dTx - dBz*dRxz*dTy + dBy*dRxz*dTz - dBy*dRzz*dTx + dBx*dRzz*dTy - dBx*dRyz*dTz;
			dTemp = (dBy*dTx - dBx*dTy)/dR + dTemp*dMaterialParameter;
			dUz = dUz + dTemp*dJacobian*GaussPointsWeights[i];
		}
		m_loNodesList.IncrementIterator();
	}
	while(!m_loNodesList.IsAtBeginning());

	Vector oDisplacement(dUx - dBx*dSolidAngle,dUy - dBy*dSolidAngle,dUz - dBz*dSolidAngle);
	oDisplacement = oDisplacement*(1.0/4.0/PI);
	return oDisplacement;
}
void SurfaceLoop::Initialize()
{
	m_loNodesList.Reset();
	m_oBurgersVector.Set(0.0,0.0,0.0);
}
double SurfaceLoop::GetSolidAngle(const Point& oPoint)
{
	// the projection is assumed to be RIGHT bounded by the loop, that is, the area projection
	// lies to the LEFT of the bounding curve. a traveller moving along the bounding curve would
	// see the projection at his LEFT at all times
	CircularLinkedList<Point> loProjections;
	Vector oTempVector;

	unsigned int i = 0;
	unsigned int iSize = m_loNodesList.GetSize();
	m_loNodesList.ResetIterator();
	for(i = 0 ; i < iSize ; i++)
	{
		oTempVector.SetByPoints(oPoint,m_loNodesList.GetCurrentItem());
		oTempVector.Normalize();
		loProjections.Append(oTempVector);
		m_loNodesList.IncrementIterator();
	}

	Point* poPreviousNode = NULL;
	Point* poThisNode = NULL;
	Point* poNextNode = NULL;
	Vector oPreviousVector;
	Vector oCurrentVector;
	Vector oNextVector;
	// see if the loop is right oriented
	bool bRightOriented = true;
	loProjections.ResetIterator();
	poPreviousNode = loProjections.GetPreviousItemPointer();
	poThisNode = loProjections.GetCurrentItemPointer();
	poNextNode = loProjections.GetNextItemPointer();
	double dTolerance = 1.0E-18;

	// constrain the first node to be convex and determine the loop orientation based
	// on it, if it turns out to be concave according to the right loop assumption, reverse
	// the loop and set the flag to be left oriented
	oPreviousVector.SetByPoints(*poPreviousNode,*poThisNode);
	oNextVector.SetByPoints(*poThisNode,*poNextNode);
	oTempVector = oPreviousVector^oNextVector;
	oTempVector.Normalize();
	oCurrentVector = *poThisNode;
	if(oTempVector*oCurrentVector > 0.0)
	{
		// concave first node
		bRightOriented = false;
		loProjections.Reverse();
	}

	bool bKeepIterating = false;
	double dArea = 0.0;
	double dArcLength1 = 0.0;
	double dArcLength2 = 0.0;
	double dArcLength3 = 0.0;
	double dSemiPerimeter = 0.0;
	double dAngularExcess = 0.0;

	// now do the actual triangle cutting, start by taking out the concave angles
	loProjections.ResetIterator();
	do
	{
		bKeepIterating = false;
		poPreviousNode = loProjections.GetPreviousItemPointer();
		poThisNode = loProjections.GetCurrentItemPointer();
		poNextNode = loProjections.GetNextItemPointer();

		oPreviousVector.SetByPoints(*poPreviousNode,*poThisNode);
		oNextVector.SetByPoints(*poThisNode,*poNextNode);
		oTempVector = oPreviousVector^oNextVector;
		oTempVector.Normalize();
		oCurrentVector = *poThisNode;

		if(oTempVector*oCurrentVector > 0.0)
		{
			oPreviousVector = *poPreviousNode;
			oCurrentVector = *poThisNode;
			oNextVector = *poNextNode;
			oPreviousVector = oPreviousVector^oCurrentVector;
			oNextVector = oCurrentVector^oNextVector;
			oPreviousVector.Normalize();
			oNextVector.Normalize();
			if(!oPreviousVector.IsSame(oNextVector) && !oPreviousVector.IsOpposite(oNextVector))
			{
				oPreviousVector = *poPreviousNode;
				oCurrentVector = *poThisNode;
				oNextVector = *poNextNode;
				dArcLength1 = acos(oPreviousVector*oCurrentVector);
				dArcLength2 = acos(oCurrentVector*oNextVector);
				dArcLength3 = acos(oNextVector*oPreviousVector);
				dSemiPerimeter = 0.5*(dArcLength1 + dArcLength2 + dArcLength3);
				dAngularExcess = tan(0.5*dSemiPerimeter)*tan(0.5*(dSemiPerimeter - dArcLength1));
				dAngularExcess = dAngularExcess*tan(0.5*(dSemiPerimeter - dArcLength2))*tan(0.5*(dSemiPerimeter - dArcLength3));
				if(dAngularExcess > dTolerance)
				{
					dAngularExcess = 4.0*atan(sqrt(dAngularExcess));
					dArea = dArea - dAngularExcess;
				}
			}
			// now, cut that triangle out by simply dropping the current point from the list
			loProjections.DropItem();
			bKeepIterating = true;
		}
		else
		{
			loProjections.IncrementIterator();
		}
	}
	while(!loProjections.IsAtBeginning() || bKeepIterating);

	// now we have a convex polygon, work on it
	loProjections.ResetIterator();
	while(loProjections.GetSize() >= 3)
	{
		poPreviousNode = loProjections.GetPreviousItemPointer();
		poThisNode = loProjections.GetCurrentItemPointer();
		poNextNode = loProjections.GetNextItemPointer();

		// this node is guaranteed to be convex, no need to check on anything
		// compute the area of that triangle, take into account that the triangle lies on
		// a spherical surface
		// compute the first arc length, which is simply the angle between the lines connecting
		// the center of the sphere (the origin in this case) to the end pojnts of the arc
		oPreviousVector = *poPreviousNode;
		oCurrentVector = *poThisNode;
		oNextVector = *poNextNode;
		oPreviousVector = oPreviousVector^oCurrentVector;
		oNextVector = oCurrentVector^oNextVector;
		oPreviousVector.Normalize();
		oNextVector.Normalize();
		if(!oPreviousVector.IsSame(oNextVector) && !oPreviousVector.IsOpposite(oNextVector))
		{
			oPreviousVector = *poPreviousNode;
			oCurrentVector = *poThisNode;
			oNextVector = *poNextNode;
			dArcLength1 = acos(oPreviousVector*oCurrentVector);
			dArcLength2 = acos(oCurrentVector*oNextVector);
			dArcLength3 = acos(oNextVector*oPreviousVector);
			dSemiPerimeter = 0.5*(dArcLength1 + dArcLength2 + dArcLength3);
			dAngularExcess = tan(0.5*dSemiPerimeter)*tan(0.5*(dSemiPerimeter - dArcLength1));
			dAngularExcess = dAngularExcess*tan(0.5*(dSemiPerimeter - dArcLength2))*tan(0.5*(dSemiPerimeter - dArcLength3));
			if(dAngularExcess > dTolerance)
			{
				dAngularExcess = 4.0*atan(sqrt(dAngularExcess));
				dArea = dArea + dAngularExcess;
			}
		}
		// now, cut that triangle out by simply dropping the current point from the list
		loProjections.DropItem();
	}

	// fix left oriented loops
	if(!bRightOriented)
	{
		dArea = -dArea;
	}
	return dArea;
}
string SurfaceLoop::ToString()
{
	m_loNodesList.ResetIterator();
	char cString[256];
	string sString = "";
	unsigned int i = 0;
	unsigned int iSize = m_loNodesList.GetSize();
	Point oThisNode;
	for(i = 0 ; i < iSize ; i++)
	{
		oThisNode = m_loNodesList.GetCurrentItem();
		sprintf(cString,"(%lf,%lf,%lf)",oThisNode.GetX(),oThisNode.GetY(),oThisNode.GetZ());
		sString = sString + cString;
		m_loNodesList.IncrementIterator();
	}
	return sString;
}



