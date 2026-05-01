#include "iostream"
#include "Tools.h"
#include "Vector.h"
#include "Matrix.h"
#include "cmath"

using namespace SupportSystem;
using namespace EZ;

void CompareDisplacements(const Vector& oDisplacement1,const Vector& oDisplacement2)
{
	Vector oDifference = oDisplacement1 - oDisplacement2;
	double dLength1 = oDisplacement1.Length();
	double dLength2 = oDisplacement2.Length();
	if(dLength1 < 1.0E-10)
	{
		if(dLength2 < 1.0E-10)
		{
			return;
		}
	}
	double dReferenceLength = min(dLength1,dLength2);
	double dTolerance = 1.0E-6;
	if(fabs(oDifference.GetX()/dReferenceLength) > dTolerance)
	{
		//printf("displacements differ : (%lf,%lf,%lf) vs (%lf,%lf,%lf)\n",oDisplacement1.GetX(),oDisplacement1.GetY(),oDisplacement1.GetZ(),oDisplacement2.GetX(),oDisplacement2.GetY(),oDisplacement2.GetZ());
		printf("displacements differ : %lf\n",oDifference.Length()/dReferenceLength);
		return;
	}
	if(fabs(oDifference.GetY()/dReferenceLength) > dTolerance)
	{
		//printf("displacements differ : (%lf,%lf,%lf) vs (%lf,%lf,%lf)\n",oDisplacement1.GetX(),oDisplacement1.GetY(),oDisplacement1.GetZ(),oDisplacement2.GetX(),oDisplacement2.GetY(),oDisplacement2.GetZ());
		printf("displacements differ : %lf\n",oDifference.Length()/dReferenceLength);
		return;
	}
	if(fabs(oDifference.GetZ()/dReferenceLength) > dTolerance)
	{
		//printf("displacements differ : (%lf,%lf,%lf) vs (%lf,%lf,%lf)\n",oDisplacement1.GetX(),oDisplacement1.GetY(),oDisplacement1.GetZ(),oDisplacement2.GetX(),oDisplacement2.GetY(),oDisplacement2.GetZ());
		printf("displacements differ : %lf\n",oDifference.Length()/dReferenceLength);
		return;
	}
}
void CompareForces(const Vector& oForce1,const Vector& oForce2)
{
	Vector oDifference = oForce1 - oForce2;
	double dLength1 = oForce1.Length();
	double dLength2 = oForce2.Length();
	if(dLength1 < 1.0E-6)
	{
		if(dLength2 < 1.0E-6)
		{
			return;
		}
	}
	double dReferenceLength = min(dLength1,dLength2);
	double dTolerance = 1.0E-6;
	if(fabs(oDifference.GetX()/dReferenceLength) > dTolerance)
	{
		//printf("forces differ : (%lf,%lf,%lf) vs (%lf,%lf,%lf)\n",oForce1.GetX(),oForce1.GetY(),oForce1.GetZ(),oForce2.GetX(),oForce2.GetY(),oForce2.GetZ());
		printf("forces differ : %lf\n",oDifference.Length()/dReferenceLength);
		return;
	}
	if(fabs(oDifference.GetY()/dReferenceLength) > dTolerance)
	{
		//printf("forces differ : (%lf,%lf,%lf) vs (%lf,%lf,%lf)\n",oForce1.GetX(),oForce1.GetY(),oForce1.GetZ(),oForce2.GetX(),oForce2.GetY(),oForce2.GetZ());
		printf("forces differ : %lf\n",oDifference.Length()/dReferenceLength);
		return;
	}
	if(fabs(oDifference.GetZ()/dReferenceLength) > dTolerance)
	{
		//printf("forces differ : (%lf,%lf,%lf) vs (%lf,%lf,%lf)\n",oForce1.GetX(),oForce1.GetY(),oForce1.GetZ(),oForce2.GetX(),oForce2.GetY(),oForce2.GetZ());
		printf("forces differ : %lf\n",oDifference.Length()/dReferenceLength);
		return;
	}
}
void CompareStresses(const Matrix& oStress1,const Matrix& oStress2)
{
	Matrix oDifference = oStress1 - oStress2;
	double dReference1 = oStress1.GetAbsoluteMaximum();
	double dReference2 = oStress2.GetAbsoluteMaximum();
	if(dReference1 < 1.0E-6)
	{
		if(dReference2 < 1.0E-6)
		{
			return;
		}
	}
	double dReference = min(dReference1,dReference2);
	double dTolerance = 1.0E-6;
	unsigned int i = 0;
	unsigned int j = 0;
	for(i = 1 ; i <= 3 ; i++)
	{
		for(j = 1 ; j <= 3 ; j++)
		{
			if(fabs(oDifference.Get(i,j)/dReference) > dTolerance)
			{
				//printf("stresses differ : (%d,%d) : %lf vs %lf : ref : %lf\n",i,j,oStress1.Get(i,j),oStress2.Get(i,j),dReference);
				printf("stresses differ : %lf\n",oDifference.GetNorm()/dReference);
				return;
			}
		}
	}
}


int main(int argc,char** argv)
{
	if(argc < 3)
	{
		printf("error: missing arguments\n");
		printf("usage: femdiff file1 file2\n");
		return 1;
	}
	FILE* fpFile1 = fopen(argv[1],"r");
	FILE* fpFile2 = fopen(argv[2],"r");
	if(fpFile1 == NULL)
	{
		printf("error: cannot find file %s\n",argv[1]);
		return 1;
	}
	if(fpFile2 == NULL)
	{
		printf("error: cannot find file %s\n",argv[2]);
		return 1;
	}
	unsigned int iNodesCount1 = 0;
	unsigned int iNodesCount2 = 0;
	string sRead = "";
	
	sRead = GetRealString(512,fpFile1);
	sscanf(sRead.c_str(),"%d\n",&iNodesCount1);
	sRead = GetRealString(512,fpFile2);
	sscanf(sRead.c_str(),"%d\n",&iNodesCount2);
	
	if(iNodesCount1 != iNodesCount2)
	{
		printf("different node counts : %d - %d\n",iNodesCount1,iNodesCount2);
		fclose(fpFile1);
		fclose(fpFile2);
		return 0;
	}
	
	unsigned int i = 0;
	double dTemp1 = 0.0;
	double dTemp2 = 0.0;
	double dTemp3 = 0.0;
	double dTemp4 = 0.0;
	double dTemp5 = 0.0;
	double dTemp6 = 0.0;
	Vector oDisplacement1;
	Vector oDisplacement2;
	Vector oForce1;
	Vector oForce2;
	Vector oVelocity1;
	Vector oVelocity2;
	Vector Acceleration1;
	Vector Acceleration2;
	Matrix oStress1(3,3);
	Matrix oStress2(3,3);
	for(i = 0 ; i < iNodesCount1 ; i++)
	{
		sRead = GetRealString(512,fpFile1);
		sscanf(sRead.c_str(),"%*f\t\t%*f\t\t%*f\t\t%lf\t\t%lf\t\t%lf\t\t%lf\t\t%lf\t\t%lf\n",&dTemp1,&dTemp2,&dTemp3,&dTemp4,&dTemp5,&dTemp6);
		oDisplacement1.Set(dTemp1,dTemp2,dTemp3);
		oForce1.Set(dTemp4,dTemp5,dTemp6);
		sRead = GetRealString(512,fpFile1);
		sscanf(sRead.c_str(),"%lf\t\t%lf\t\t%lf\t\t%lf\t\t%lf\t\t%lf\n",&dTemp1,&dTemp2,&dTemp3,&dTemp4,&dTemp5,&dTemp6);
		oStress1.Set(1,1,dTemp1);
		oStress1.Set(2,2,dTemp2);
		oStress1.Set(3,3,dTemp3);
		oStress1.Set(1,2,dTemp4);
		oStress1.Set(2,3,dTemp5);
		oStress1.Set(3,1,dTemp6);
		oStress1.Set(2,1,dTemp4);
		oStress1.Set(3,2,dTemp5);
		oStress1.Set(1,3,dTemp6);
		oStress1.Filter();
		sRead = GetRealString(512,fpFile1);
		sscanf(sRead.c_str(),"%lf\t\t%lf\t\t%lf\t\t%lf\t\t%lf\t\t%lf\n",&dTemp1,&dTemp2,&dTemp3,&dTemp4,&dTemp5,&dTemp6);
		oVelocity1.Set(dTemp1,dTemp2,dTemp3);
		Acceleration1.Set(dTemp4,dTemp5,dTemp6);
		
		sRead = GetRealString(512,fpFile2);
		sscanf(sRead.c_str(),"%*f\t\t%*f\t\t%*f\t\t%lf\t\t%lf\t\t%lf\t\t%lf\t\t%lf\t\t%lf\n",&dTemp1,&dTemp2,&dTemp3,&dTemp4,&dTemp5,&dTemp6);
		oDisplacement2.Set(dTemp1,dTemp2,dTemp3);
		oForce2.Set(dTemp4,dTemp5,dTemp6);
		sRead = GetRealString(512,fpFile2);
		sscanf(sRead.c_str(),"%lf\t\t%lf\t\t%lf\t\t%lf\t\t%lf\t\t%lf\n",&dTemp1,&dTemp2,&dTemp3,&dTemp4,&dTemp5,&dTemp6);
		oStress2.Set(1,1,dTemp1);
		oStress2.Set(2,2,dTemp2);
		oStress2.Set(3,3,dTemp3);
		oStress2.Set(1,2,dTemp4);
		oStress2.Set(2,3,dTemp5);
		oStress2.Set(3,1,dTemp6);
		oStress2.Set(2,1,dTemp4);
		oStress2.Set(3,2,dTemp5);
		oStress2.Set(1,3,dTemp6);
		oStress2.Filter();
		sRead = GetRealString(512,fpFile2);
		sscanf(sRead.c_str(),"%lf\t\t%lf\t\t%lf\t\t%lf\t\t%lf\t\t%lf\n",&dTemp1,&dTemp2,&dTemp3,&dTemp4,&dTemp5,&dTemp6);
		oVelocity2.Set(dTemp1,dTemp2,dTemp3);
		Acceleration2.Set(dTemp4,dTemp5,dTemp6);
		
		CompareDisplacements(oDisplacement1,oDisplacement2);
		CompareDisplacements(oVelocity1,oVelocity2);
		CompareDisplacements(Acceleration1,Acceleration2);
		CompareForces(oForce1,oForce2);
		CompareStresses(oStress1,oStress2);
	}
	
	fclose(fpFile1);
	fclose(fpFile2);

	return 0;
}

