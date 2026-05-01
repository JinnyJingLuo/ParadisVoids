#include "iostream"
#include "Tools.h"
#include "float.h"
#include "cmath"

using namespace std;
using namespace SupportSystem;


int main(int argc,char** argv)
{
	if(argc < 2)
	{
		printf("error: output file missing\n");
		return 1;
	}
	FILE* fpFile = fopen(argv[1],"r");
	if(fpFile == NULL)
	{
		printf("error: cannot input output file\n");
		return 1;
	}
	unsigned int i = 0;
	string sRead;
	unsigned int iStartLine = 1010964 - 12;
	for(i = 0 ; i < iStartLine ; i++)
	{
		sRead = GetRealString(500,fpFile);
	}
	// done skipping
	// loop till the end line
	unsigned int iEndLine = 1013269 - 12;
	char cTest;
	double dTemp = 0.0;
	double dMax = -DBL_MAX;
	unsigned int iMax = 0;
	dMax = 0.0;
	for(i = iStartLine ; i < iEndLine ; i++)
	{
		sRead = GetRealString(1024,fpFile);
		cTest = sRead[0];
		
		if(sRead.find("dDeltaPlasticStrain") != sRead.npos)
		{
			sscanf(sRead.c_str(),"%*s%lf\n",&dTemp);
			dMax = dMax + dTemp;
			printf("current %25.20f\n",dTemp);
			if(fabs(dTemp) > dMax)
			{
				dMax = fabs(dTemp);
				iMax = i;
			}
		}
	}
	
			printf("max is %25.20f at %d\n",dMax,iMax);
	fclose(fpFile);
	return 0;
}


