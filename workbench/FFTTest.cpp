#include "iostream"
#include "Randomizer.h"
#include "cmath"
#include "Point.h"
#include "SurfaceProfile.h"

using namespace EZ;

int main(int argc,char** argv)
{
	SurfaceProfile oSurface;
	oSurface.Read(argv[1]);
	FourierTransformable* poAutoCorrelation = oSurface.GetAutoCorrelation();
	poAutoCorrelation->WriteData("auto_corr.txt");
	oSurface.GetAmplitudeCharacteristics();
	delete poAutoCorrelation;
	
// 	oSurface.WriteTransform("transform.txt");
// 	oSurface.Shift();
// 	oSurface.WriteTransform("shift.txt");
	// unsigned int iSpectrumSize = 0;
// 	double* pdSpectrum = oSurface.GetPowerSpectrum(iSpectrumSize,true);
// 	unsigned int i = 0;
// 	// write the file
// 	FILE* fpFile = fopen(argv[2],"w");
// 	for(i = 0 ; i < iSpectrumSize ; i++)
// 	{
// 		fprintf(fpFile,"%d\t\t%lf\n",i,pdSpectrum[i]);
// 	}
// 	fclose(fpFile);
// 	delete [] pdSpectrum;
}

