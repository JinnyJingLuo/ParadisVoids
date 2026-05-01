#include "stdio.h"
#include "PrecipitateStructure.h"
#include "Plane.h"
#include "string"
#include "math.h"

using namespace std;
using namespace GeometrySystem;

PrecipitateStructure GenerateStructure(int argc,char** argv)
{
	PrecipitateStructure oStructure;
	oStructure.Set(string(argv[1]));
	return oStructure;
}

void RemoveNonIntersecting(PrecipitateStructure* poStructure,int argc,char** argv)
{
	double dA = atof(argv[2]);
	double dB = atof(argv[3]);
	double dC = atof(argv[4]);
	double dD = atof(argv[5]);
	Vector oNormal(dA,dB,dC);
	oNormal.Normalize();
	// get a point on the plane
	double dX = 0.0;
	double dY = 0.0;
	double dZ = 0.0;
	if(fabs(oNormal.GetX()) > fabs(oNormal.GetY()))
	{
		if(fabs(oNormal.GetX()) > fabs(oNormal.GetZ()))
		{
			// x is the largest
			dX = dD/dA;
		}
		else
		{
			// z is the largest
			dZ = dD/dC;
		}
	}
	else
	{
		if(fabs(oNormal.GetY()) > fabs(oNormal.GetZ()))
		{
			// y is the largest
			dY = dD/dB;
		}
		else
		{
			// z is the largest
			dZ = dD/dC;
		}
	}
	Point oPlanePoint(dX,dY,dZ);
	Plane oPlane(oNormal,oPlanePoint);
	printf("before plane intersection : %d\n",poStructure->GetPrecipitatesCount());
	poStructure->RemoveNonIntersectingPrecipitates(oPlane);
	printf("after plane intersection : %d\n",poStructure->GetPrecipitatesCount());
	poStructure->WritePrecipitateStructure("red_precip.txt");
}

void RemoveExternal(PrecipitateStructure* poStructure)
{
	double dXDim = 10000.0;
	double dYDim = 10000.0;
	double dZDim = 10000.0;
	double ChannelWidth = 200.0;
	
	Vector oXAxis(1.0,0.0,-1.0);
	Vector oZAxis(1.0,1.0,1.0);
	Vector oYAxis = oZAxis^oXAxis;
	oXAxis.Normalize();
	oYAxis.Normalize();
	oZAxis.Normalize();
	
	Point oPXPlanePoint = oXAxis*(0.5*dXDim - ChannelWidth);
	Point oNXPlanePoint = oXAxis*(-0.5*dXDim + ChannelWidth);
	Point oPYPlanePoint = oYAxis*(0.5*dYDim - ChannelWidth);
	Point oNYPlanePoint = oYAxis*(-0.5*dYDim + ChannelWidth);
	Point oPZPlanePoint = oZAxis*(0.5*dZDim - ChannelWidth);
	Point oNZPlanePoint = oZAxis*(-0.5*dZDim + ChannelWidth);
	
	Plane oPXPlane(oXAxis,oPXPlanePoint);
	Plane oNXPlane(oXAxis*(-1.0),oNXPlanePoint);
	Plane oPYPlane(oYAxis,oPYPlanePoint);
	Plane oNYPlane(oYAxis*(-1.0),oNYPlanePoint);
	Plane oPZPlane(oZAxis,oPZPlanePoint);
	Plane oNZPlane(oZAxis*(-1.0),oNZPlanePoint);
	
	unsigned int iPerturbationPointsCount = 10;
	double dPerturnationMean = 30.0;
	
	poStructure->PlaneCut(oPXPlane,true,iPerturbationPointsCount,dPerturnationMean);
	printf("after px cut %d\n",poStructure->GetPrecipitatesCount());
	poStructure->WritePrecipitateStructure("final_px_precip.txt");
	
	poStructure->PlaneCut(oNXPlane,true,iPerturbationPointsCount,dPerturnationMean);
	printf("after nx cut %d\n",poStructure->GetPrecipitatesCount());
	poStructure->WritePrecipitateStructure("final_nx_precip.txt");
	
	poStructure->PlaneCut(oPYPlane,true,iPerturbationPointsCount,dPerturnationMean);
	printf("after py cut %d\n",poStructure->GetPrecipitatesCount());
	poStructure->WritePrecipitateStructure("final_py_precip.txt");
	
	poStructure->PlaneCut(oNYPlane,true,iPerturbationPointsCount,dPerturnationMean);
	printf("after ny cut %d\n",poStructure->GetPrecipitatesCount());
	poStructure->WritePrecipitateStructure("final_ny_precip.txt");
	
	poStructure->PlaneCut(oPZPlane,true,iPerturbationPointsCount,dPerturnationMean);
	printf("after pz cut %d\n",poStructure->GetPrecipitatesCount());
	poStructure->WritePrecipitateStructure("final_pz_precip.txt");
	
	poStructure->PlaneCut(oNZPlane,true,iPerturbationPointsCount,dPerturnationMean);
	printf("after nz cut %d\n",poStructure->GetPrecipitatesCount());
	poStructure->WritePrecipitateStructure("final_nz_precip.txt");
	
	poStructure->WritePrecipitateStructure("final_precip.txt");
}

void Trim()
{
	printf("trimmer running\n");
	AxisAlignedBoundingBox oBox;
	oBox.SetXMin(-1000.0);
	oBox.SetXMax(1000.0);
	oBox.SetYMin(-1000.0);
	oBox.SetYMax(1000.0);
	oBox.SetZMin(-1000.0);
	oBox.SetZMax(1000.0);
	
	unsigned int iSize = 100;
	unsigned int i = 0;
	list< Point > loPoints;

	for(i = 0 ; i < iSize ; i++)
	{
		loPoints.push_back(oBox.GenerateRandomPoint(1.0));
	}
	Polyhedron oOriginal;
	oOriginal.CreateAsHullFromPoints(&loPoints);
	Vector oNormal(0.0,0.0,1.0);
	Point oPoint(0.0,0.0,0.0);
	Polyhedron oCut = oOriginal.PlaneCut(Plane(oNormal,oPoint),true,10,50.0);
	oOriginal.WriteParaview("original.vtk");
	oCut.WriteParaview("cut.vtk");
}

int main(int argc,char** argv)
{	
	if(argc < 6)
	{
		printf("error: too few input arguments\n");
		printf("usage: trim precip_input_file a b c d (such that the plane equation is ax + by + cz = d)\n");
		return 1;
	}
	PrecipitateStructure oStructure = GenerateStructure(argc,argv);
	RemoveNonIntersecting(&oStructure,argc,argv);
	RemoveExternal(&oStructure);
	return 0;
}

