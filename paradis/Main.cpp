/***************************************************************************
 *
 *  Function    : Main
 *  Description : main routine for ParaDiS simulation
 *
 **************************************************************************/
#include <stdio.h>
#include <time.h>
#include "Home.h"
#include "Init.h"
#include "ParadisBlockSurface.h"
#include "ParadisCylinderSurface.h"
#include "ParadisExternalLoadServer.h"
#include "ParadisCrossSlipServer.h"
#include "LoopFilter.h"
#include "DSMPI.h"
#include "Tools.h"
#include "ParadisPrecipitateServer.h"
#include "queue"

using namespace SupportSystem;
using namespace std;

/*
//struct Pixel
//{
//	int m_iX;
//	int m_iY;
//	int m_iColor;
//	int m_iType;
//};
//
//Pixel** ExtractVolume(Pixel** ppoData,const unsigned int iStartPixelX,const
unsigned int iStartPixelY)
//{
//	return NULL;
//}
//
//int main(int argc,char** argv)
//{
//	FILE* fpFile = fopen("precip42.txt","r");
//	string sRead;
//	sRead = GetRealString(512,fpFile);
//	unsigned int iXResolution = 0;
//	unsigned int iYResolution = 0;
//	sscanf(sRead.c_str(),"%d X %d\n",&iXResolution,&iYResolution);
//	Pixel** ppoData = new Pixel*[iXResolution];
//	unsigned int i = 0;
//	unsigned int j = 0;
//	for(i = 0 ; i < iXResolution ; i++)
//	{
//		ppoData[i] = new Pixel[iYResolution];
//	}
//
//	unsigned int iIndex = 0;
//	int iRead = 0;
//	while(!feof(fpFile))
//	{
//		iRead = fgetc(fpFile);
//		if((iRead == 48) || (iRead == 49))		// read 0 or 1
//		{
//			i = iIndex/iXResolution;
//			j = iIndex%iYResolution;
//			ppoData[i][j].m_iX = i;
//			ppoData[i][j].m_iY = j;
//			ppoData[i][j].m_iColor = iRead - 48;
//			ppoData[i][j].m_iType = 0;
//			iIndex++;
//		}
//	}
//	fclose(fpFile);
//
//	// a pixel can only be of type 0 : unprocessed, or 1 : processed and
internal, 2 :  processed and external
//	// or 3 : in queue for processing
//	queue< Pixel > qoPixels;
//	int iInitialPixelX = 180;
//	int iIntiialPixelY = 280;
//	qoPixels.push(ppoData[iInitialPixelX][iIntiialPixelY]);
//	int iTargetColor = 0;
//	int iMaxI = 0;
//	int iMinI = iXResolution;
//	int iMaxJ = 0;
//	int iMinJ = iYResolution;
//	Pixel oPixel;
//	while(!qoPixels.empty())
//	{
//		// read the first pixel coordinates
//		oPixel = qoPixels.front();
//		// remove it from the queue
//		qoPixels.pop();
//		i = oPixel.m_iX;
//		j = oPixel.m_iY;
//		// see if it is of the target color
//		if(oPixel.m_iColor == iTargetColor)
//		{
//			// mark this node as processed and internal
//			oPixel.m_iType = 1;
//			if(i < iMinI)
//			{
//				iMinI = i;
//			}
//			if(i > iMaxI)
//			{
//				iMaxI = i;
//			}
//			if(j < iMinJ)
//			{
//				iMinJ = j;
//			}
//			if(j > iMaxJ)
//			{
//				iMaxJ = j;
//			}
//			// add its unprocessed neighbours in the 4 principal
directions
//			if((i > 0) && (ppoData[i - 1][j].m_iType == 0))
//			{
//				qoPixels.push(ppoData[i - 1][j]);
//				ppoData[i - 1][j].m_iType = 3;
//			}
//			if((i < iXResolution - 1) && (ppoData[i + 1][j].m_iType
== 0))
//			{
//				qoPixels.push(ppoData[i + 1][j]);
//				ppoData[i + 1][j].m_iType = 3;
//			}
//			if((j > 0) && (ppoData[i][j - 1].m_iType == 0))
//			{
//				qoPixels.push(ppoData[i][j - 1]);
//				ppoData[i][j - 1].m_iType = 3;
//			}
//			if((j < iYResolution - 1) && (ppoData[i][j + 1].m_iType
== 0))
//			{
//				qoPixels.push(ppoData[i][j + 1]);
//				ppoData[i][j + 1].m_iType = 3;
//			}
//		}
//		else
//		{
//			// mark it as processed and external
//			oPixel.m_iType = 2;
//		}
//	}
//
//	unsigned int iXRange = iMaxI - iMinI + 1;
//	unsigned int iYRange = iMaxJ - iMinJ + 1;
//	// write the output
//	FILE* fpOutput = fopen("precipitate_01.txt","w");
//	fprintf(fpOutput,"%d X %d\n",iXRange,iYRange);
//	for(i = iMinI ; i <= iMaxI ; i++)
//	{
//		for(j = iMinJ ; j <= iMaxJ ; j++)
//		{
//			if(piTypes[i][j] == 1)
//			{
//				fprintf(fpOutput,"0 ");
//			}
//			else
//			{
//				fprintf(fpOutput,"1 ");
//			}
//			// new line every 40 characters
//			if(((j - iMinJ + 1)%40) == 0)
//			{
//				fprintf(fpOutput,"\n");
//			}
//		}
//		fprintf(fpOutput,"\n");
//	}
//	fclose(fpOutput);
//
//	for(i = 0 ; i < iXResolution ; i++)
//	{
//		if(piData[i] != NULL)
//		{
//			delete [] piData[i];
//		}
//		if(piTypes[i] != NULL)
//		{
//			delete [] piTypes[i];
//		}
//	}
//	delete [] piData;
//	delete [] piTypes;
//}
*/

/*
void GeneratePoints()
{
        FILE* fpFile = fopen("input.txt","w");
        double dX = 0.0;
        double dY = 0.0;
        double dZ = 0.0;
        unsigned int iSize = 1000;
        unsigned int i = 0;
        fprintf(fpFile,"%d\n",iSize);
        for(i = 0 ; i < iSize ; i++)
        {
                dX = Randomizer::RandomNormal(0,2000);
                dY = Randomizer::RandomNormal(0,2000);
                dZ = Randomizer::RandomNormal(0,2000);
                fprintf(fpFile,"%e\t\t%e\t\t%e\n",dX,dY,dZ);
        }
        fclose(fpFile);
}

void ReadInput(const string& sFileName,list< Point >* ploPoints)
{
        // read the input nodes, notice that not all of the nodes read here will
end up in the final convex hull, a
        // cleaning step will be carried out later to get rid of these nodes
        ploPoints->clear();
        FILE* fpFile = fopen(sFileName.c_str(),"r");
        string sRead;
        sRead = GetRealString(512,fpFile);
        unsigned int iVerticesCount = 0;
        sscanf(sRead.c_str(),"%d\n",&iVerticesCount);
        unsigned int i = 0;
        double dX = 0.0;
        double dY = 0.0;
        double dZ = 0.0;
        for(i = 0 ; i < iVerticesCount ; i++)
        {
                sRead = GetRealString(512,fpFile);
                sscanf(sRead.c_str(),"%lf\t%lf\t%lf\n",&dX,&dY,&dZ);
                ploPoints->push_back(Point(dX,dY,dZ));
        }
        fclose(fpFile);
}

Matrix DoScrew(const Point& oPoint)
{
        double dMu = 8.0e10;
        double dNu = 0.3;
        double dX = oPoint.GetX();
        double dY = oPoint.GetY();
        double dZ = oPoint.GetZ();
        double dR = (dX*dX + dY*dY);
        double dF1 = dMu/2.0/PI;
        Matrix oStress(3,3);
        oStress.Set(1,1,0.0);
        oStress.Set(1,2,0.0);
        oStress.Set(1,3,-dF1*dY/dR);
        oStress.Set(2,1,0.0);
        oStress.Set(2,2,0.0);
        oStress.Set(2,3,dF1*dX/dR);
        oStress.Set(3,1,oStress.Get(1,3));
        oStress.Set(3,2,oStress.Get(2,3));
        oStress.Set(3,3,0.0);
        return oStress;
}

Matrix DoEdge(const Point& oPoint)
{
        double dMu = 8.0e10;
        double dNu = 0.3;
        double dX = oPoint.GetX();
        double dY = oPoint.GetY();
        double dZ = oPoint.GetZ();
        double dR2 = (dX*dX + dY*dY)*(dX*dX + dY*dY);
        double dF1 = dMu/2.0/PI/(1.0 - dNu);
        double dF2 = dMu*dNu/PI/(1.0 - dNu);
        Matrix oStress(3,3);
        oStress.Set(1,1,-dF1*dY*(3.0*dX*dX + dY*dY)/dR2);
        oStress.Set(1,2,dF1*dX*(dX*dX - dY*dY)/dR2);
        oStress.Set(1,3,0.0);
        oStress.Set(2,1,oStress.Get(1,2));
        oStress.Set(2,2,dF1*dY*(dX*dX - dY*dY)/dR2);
        oStress.Set(2,3,0.0);
        oStress.Set(3,1,0.0);
        oStress.Set(3,2,0.0);
        oStress.Set(3,3,-dF2*dY/sqrt(dR2));
        return oStress;
}

void StressTest()
{
        Point oQ(0.0,0.0,-100000.0);
        Point oP(0.0,0.0,100000.0);
        Vector oBurgers(0.0,0.0,1.0);
        oBurgers = oBurgers;
        double dMu = 8.0e10;
        double dNu = 0.3;

        double dR = 0.01;
        Point oY;
        double dTheta = 0.0;
        unsigned int iSize = 360000;
        unsigned int i = 0;
        double dThetaIncrement = 2.0*PI/iSize;
        Matrix oStress1;
        Matrix oStress2;
        double dMax = 0.0;
        double dError = 0.0;
        for(i = 0 ; i < iSize ; i++)
        {
                dTheta = i*dThetaIncrement;
                oY.Set(dR*cos(dTheta),dR*sin(dTheta),0.0);
                oStress1 =
ParadisExternalLoadServer::ComputeSegmentStressAtPoint(oQ,oP,oBurgers,oY,dMu,dNu);
                oStress2 = DoScrew(oY);
                dError = fabs(100.0*((oStress1 -
oStress2).GetNorm()/oStress2.GetNorm())); printf("error percentage :
%e\n",dError); if(dError > dMax)
                {
                        dMax = dError;
                }
        }
        printf("max error : %e\n",dMax);
}

void SingularStressTest()
{
        Point oQ(-5.0,0.0,0.0);
        Point oP(5.0,0.0,0.0);
        Vector oBurgers(0.0,3.0,1.0);
        double dMu = 8.0e10;
        double dNu = 0.3;

        Point oY(19.0,5.5,04);
        Matrix oStress =
ParadisExternalLoadServer::ComputeSegmentStressAtPoint(oQ,oP,oBurgers,oY,dMu,dNu);

        oStress.Print();
}
*/

void Do() {
  Point oP1(-4900.0, -5050.0, 0.0);
  Point oP2(4900.0, -5050.0, 0.0);
  Point oP3(-4900.0, -5500.0, 0.0);
  Point oP4(4900.0, -5500.0, 0.0);

  CartesianOrthogonalCoordinateSystem oSystem;
  oSystem.SetXZ(Vector(1.0, 0.0, -1.0), Vector(1.0, 1.0, 1.0));
  printf("%s\n", oSystem.GetInGlobalCoordinates(oP1).ToString().c_str());
  printf("%s\n", oSystem.GetInGlobalCoordinates(oP2).ToString().c_str());
  printf("%s\n", oSystem.GetInGlobalCoordinates(oP3).ToString().c_str());
  printf("%s\n", oSystem.GetInGlobalCoordinates(oP4).ToString().c_str());
}

int main(int argc, char **argv) {
  if (argc != 2) {
    printf("error: incorrect number of arguments\n");
    return 1;
  }
  string sControlFileName = argv[1];

  // initialize the system and all the servers
  Home_t *home = (Home_t *)calloc(1, sizeof(Home_t));
  DSMPI::Initialize(argc, argv);
  home->myDomain = DSMPI::GetProcessID();
  home->numDomains = DSMPI::GetProcessesCount();
  if (home->myDomain == 0) {
    printf("initializing paradis\n");
    fflush(NULL);
  }
  DSMPI::Barrier();
  Initialize(home, sControlFileName, true);
  home->cycle = home->param->cycleStart;
  unsigned int iCycleEnd = home->param->cycleStart + home->param->maxstep;
  DSMPI::Barrier();

  // intitialize collision server
  if (home->myDomain == 0) {
    printf("initializing collision server\n");
    fflush(NULL);
  }
  home->poCollisionServer = ParadisCollisionServer::CreateInstance();
  home->poCollisionServer->Set(home);
  DSMPI::Barrier();

  // intitialize precipitates
  if (home->myDomain == 0) {
    printf("initializing precipitates\n");
    fflush(NULL);
  }
  home->poPrecipitateServer = ParadisPrecipitateServer::CreateInstance();
  home->poPrecipitateServer->Set(home);

  // set the surface server
  if (home->myDomain == 0) {
    printf("initializing surface\n");
    fflush(NULL);
  }
  if (home->param->GeometryType == 1) {
    home->poSurface = new ParadisBlockSurface;
  } else if (home->param->GeometryType == 4) {
    home->poSurface = new ParadisCylinderSurface;
  } else {
    printf("error: unknown surface geometry : %d\n", home->param->GeometryType);
    return 2;
  }
  home->poSurface->Set(home);

  // set the FEM server
  if (home->myDomain == 0) {
    printf("initializing FEM\n");
    fflush(NULL);
  }
  home->poExternalLoadServer = ParadisExternalLoadServer::CreateInstance();
  if (!home->poExternalLoadServer->Set(home)) {
    printf("error: couldn't set paradis external load server\n");
    return 1;
  }
  double dFEMTimeStep = home->poExternalLoadServer->GetTimeStep();
  double dLastTimeStepQuotient =
      home->poExternalLoadServer->GetCurrentTime() / dFEMTimeStep;
  double dTimeStepQuotient = 0.0;

  // remove initial surface arms if required
  if (home->myDomain == 0) {
    printf("removing initial surface arms\n");
    fflush(NULL);
  }
  DSMPI::Barrier();
  if (home->param->BoundaryType != PERIODIC_BOUNDARY) {
    // junjie: Some nodes are close to the rigid boundary, when restart, those
    // nodes
    // can potentially be considered a node out of the surface and placed to its
    // old position. But it's your first step, you don't have old positions. We
    // better not do surface check when we restart.
    if (home->cycle <= 1) {
      home->poSurface->CheckNodes(home);
      home->poSurface->StoreSurfaceArms(home);
    }
  }
  DSMPI::Barrier();

  // the main paradis loop
  unsigned int iNodesCount = GetNodeCount(home);
  if (home->myDomain == 0) {
    printf("entering paradis main loop, initial node count is %d\n",
           iNodesCount);
    fflush(NULL);
  }
  DSMPI::Barrier();
  while (home->cycle < iCycleEnd) {
    // apply FEM if required
    if (home->param->EnableFEM != 0) {
      dTimeStepQuotient = floor(home->param->timeNow / dFEMTimeStep);
      if (dTimeStepQuotient > dLastTimeStepQuotient) {
        dLastTimeStepQuotient = dTimeStepQuotient;
        home->poExternalLoadServer->RunFEM(home->param->timeNow, home);
      }
    }
    DSMPI::Barrier();
    ParadisStep(home);
    DSMPI::Barrier();
  }

  ParadisCrossSlipServer::Free();
  delete home->poPrecipitateServer;
  delete home->poSurface;
  delete home->poExternalLoadServer;
  ParadisFinish(home);
  return 0;
}
