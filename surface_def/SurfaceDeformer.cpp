#include "iostream"
#include "SurfaceDeformer.h"
#include "Tools.h"
#include "CartesianOrthogonalCoordinateSystem.h"
#include "SurfaceSegment.h"
#include "SurfaceProfile.h"
#include "mpi.h"
#include "math.h"

using namespace GeometrySystem;
using namespace SupportSystem;


TimePoint::TimePoint()
{
	Initialize();
}
TimePoint::TimePoint(const unsigned int& iTimeStep,const string& sFileName)
{
	Initialize();
	SetTimeStep(iTimeStep);
	OpenFile(sFileName);
}
TimePoint::TimePoint(const TimePoint& oPoint)
{
	*this = oPoint;
}
TimePoint::~TimePoint()
{
	Reset();
}
TimePoint& TimePoint::operator=(const TimePoint& oPoint)
{
	m_iTimeStep = oPoint.m_iTimeStep;
	fpFilterFile = oPoint.fpFilterFile;
	return *this;
}
void TimePoint::Reset()
{
	CloseFile();
	Initialize();
}
void TimePoint::OpenFile(const string& sFileName)
{
	fpFilterFile = fopen(sFileName.c_str(),"w");
}
void TimePoint::CloseFile()
{
	if(fpFilterFile != NULL)			fclose(fpFilterFile);
	fpFilterFile = NULL;
}
void TimePoint::FlushFile()
{
	if(fpFilterFile != NULL)			fflush(fpFilterFile);
}
void TimePoint::Write(const string& sString)
{
	fprintf(fpFilterFile,sString.c_str());
}
bool TimePoint::IsBefore(const unsigned int& iTimeStep)
{
	return (iTimeStep < m_iTimeStep);
}
void TimePoint::Initialize()
{
	m_iTimeStep = 0;
	fpFilterFile = NULL;
}
void TimePoint::SetTimeStep(const unsigned int& iTimeStep)
{
	m_iTimeStep = iTimeStep;
}

SurfaceDeformer* SurfaceDeformer::m_poInstance = NULL;
SurfaceDeformer* SurfaceDeformer::GetInstance()
{
	if(m_poInstance == NULL)
	{
		m_poInstance = new SurfaceDeformer;
	}
	return m_poInstance;
}
SurfaceDeformer::SurfaceDeformer()
{
	Initialize();
}
SurfaceDeformer::~SurfaceDeformer()
{
	Reset();
}
void SurfaceDeformer::Reset()
{
	ClearSurfaceSegments();
	ClearTimePoints();
	Initialize();
}
void SurfaceDeformer::SetParallelizationParameters(const int& iProcessesCount,const int& iProcessID)
{
	m_iProcessesCount = iProcessesCount;
	m_iProcessID = iProcessID;
}
void SurfaceDeformer::Read(const string& sInputFile)
{
	FILE* fpInput = fopen(sInputFile.c_str(),"r");
	
	m_sSurfaceSegmentsFile = GetRealString(512,fpInput);
	m_sOutputFile = GetRealString(512,fpInput);
	
	string sRead = "";
	sRead = GetRealString(512,fpInput);
	sscanf(sRead.c_str(),"%lf\n",&m_dPoissonsRatio);
	sRead = GetRealString(512,fpInput);
	sscanf(sRead.c_str(),"%lf\n",&m_dVirtualNodeSpacing);
	sRead = GetRealString(512,fpInput);
	sscanf(sRead.c_str(),"%d\n",&m_iGaussPointsCount);
	
	double dX = 0.0;
	double dY = 0.0;
	double dZ = 0.0;
	sRead = GetRealString(512,fpInput);
	sscanf(sRead.c_str(),"%lf\t%lf\t%lf\n",&dX,&dY,&dZ);
	m_oPlanePoint.Set(dX,dY,dZ);
	sRead = GetRealString(512,fpInput);
	sscanf(sRead.c_str(),"%lf\t%lf\t%lf\n",&dX,&dY,&dZ);
	m_oPlaneNormal.Set(dX,dY,dZ);
	sRead = GetRealString(512,fpInput);
	sscanf(sRead.c_str(),"%lf\t%lf\t%lf\n",&dX,&dY,&dZ);
	m_oPlaneX.Set(dX,dY,dZ);
	
	sRead = GetRealString(512,fpInput);
	sscanf(sRead.c_str(),"%lf\n",&m_dLength);
	sRead = GetRealString(512,fpInput);
	sscanf(sRead.c_str(),"%lf\n",&m_dWidth);
	
	sRead = GetRealString(512,fpInput);
	sscanf(sRead.c_str(),"%d\n",&m_iLengthResolution);
	sRead = GetRealString(512,fpInput);
	sscanf(sRead.c_str(),"%d\n",&m_iWidthResolution);
	
	fclose(fpInput);
	
	// set the surface segment class parameters
	SurfaceSegment::SetVirtualNodeSpacing(m_dVirtualNodeSpacing);
	SurfaceSegment::SetPoissonsRatio(m_dPoissonsRatio);
	SurfaceSegment::SetGaussPointsCount(m_iGaussPointsCount);

	SurfaceLoop::SetVirtualNodeSpacing(m_dVirtualNodeSpacing);
	SurfaceLoop::SetPoissonsRatio(m_dPoissonsRatio);
	SurfaceLoop::SetGaussPointsCount(m_iGaussPointsCount);
}
void SurfaceDeformer::Print() const
{
	printf("input : %s\n",m_sSurfaceSegmentsFile.c_str());
	printf("output : %s\n",m_sOutputFile.c_str());
	printf("plane point : %lf\t%lf\t%lf\n",m_oPlanePoint.GetX(),m_oPlanePoint.GetY(),m_oPlanePoint.GetZ());
	printf("plane normal : %lf\t%lf\t%lf\n",m_oPlaneNormal.GetX(),m_oPlaneNormal.GetY(),m_oPlaneNormal.GetZ());
	printf("plane X : %lf\t%lf\t%lf\n",m_oPlaneX.GetX(),m_oPlaneX.GetY(),m_oPlaneX.GetZ());
	printf("length : %lf\n",m_dLength);
	printf("width : %lf\n",m_dWidth);
	printf("Poisson's ratio : %lf\n",m_dPoissonsRatio);
	printf("virtual node spacing : %lf\n",m_dVirtualNodeSpacing);
	printf("length resolution : %d\n",m_iLengthResolution);
	printf("width resolution : %d\n",m_iWidthResolution);
	printf("Gauss points count : %d\n",m_iGaussPointsCount);
}
void SurfaceDeformer::GenerateSurfacePoints()
{
	m_vvoSurfacePoints.clear();
	m_vvoDisplacements.clear();
	unsigned int i = 0;
	m_vvoSurfacePoints.resize(m_iLengthResolution);
	m_vvoDisplacements.resize(m_iLengthResolution);
	for(i = 0 ; i < m_iLengthResolution ; i++)
	{
		m_vvoSurfacePoints[i].resize(m_iWidthResolution);
		m_vvoDisplacements[i].resize(m_iWidthResolution);
	}
	
	unsigned int j = 0;
	double dX = 0.0;
	double dY = 0.0;
	double dZ = 0.0;
	double dXMin = -0.5*m_dWidth;
	double dYMin = -0.5*m_dLength;
	double dXIncrement = m_dWidth/(m_iWidthResolution - 1);
	double dYIncrement = m_dLength/(m_iLengthResolution - 1);
	CartesianOrthogonalCoordinateSystem oSystem;
	oSystem.SetOrigin(m_oPlanePoint);
	oSystem.SetXZ(m_oPlaneX,m_oPlaneNormal);
	for(i = 0 ; i < m_iLengthResolution ; i++)
	{
		dY = dYMin + i*dYIncrement;
		for(j = 0 ; j < m_iWidthResolution ; j++)
		{
			dX = dXMin + j*dXIncrement;
			m_vvoSurfacePoints[i][j] = oSystem.GetInGlobalCoordinates(Point(dX,dY,dZ));
		}
	}
}
void SurfaceDeformer::PrintSurfacePoints() const
{
	unsigned int i = 0;
	unsigned int j = 0;
	Point oPoint;
	for(i = 0 ; i < m_iLengthResolution ; i++)
	{
		for(j = 0 ; j < m_iWidthResolution ; j++)
		{
			oPoint = m_vvoSurfacePoints[i][j];
			printf("%lf,%lf,%lf\n",oPoint.GetX(),oPoint.GetY(),oPoint.GetZ());
		}
	}
}
void SurfaceDeformer::ReadSurfaceSegments()
{
	ClearSurfaceSegments();
	FILE* fpFile = fopen(m_sSurfaceSegmentsFile.c_str(),"r");
	string sRead;
	double dX1 = 0.0;
	double dY1 = 0.0;
	double dZ1 = 0.0;
	double dX2 = 0.0;
	double dY2 = 0.0;
	double dZ2 = 0.0;
	double dNX1 = 0.0;
	double dNY1 = 0.0;
	double dNZ1 = 0.0;
	double dNX2 = 0.0;
	double dNY2 = 0.0;
	double dNZ2 = 0.0;
	Point oNode1;
	Point oNode2;
	Vector oBurgersVector;
	Vector oSlipPlaneNormal;
	Vector oSurfaceNormal1;
	Vector oSurfaceNormal2;
	SurfaceLoop* poSurfaceLoop = NULL;
	double dTolerance = 1.0E-2;

	while(!feof(fpFile))
	{
		sRead = GetRealString(1024,fpFile);
		if(sRead.empty())
		{
			break;
		}
		sscanf(sRead.c_str(),"%*d,%*d : (%*d,%*d) -> (%*d,%*d) : (%lf,%lf,%lf) / (%lf,%lf,%lf)\n",&dX1,&dY1,&dZ1,&dX2,&dY2,&dZ2);
		oSlipPlaneNormal.Set(dX1,dY1,dZ1);
		oBurgersVector.Set(dX2,dY2,dZ2);
		sRead = GetRealString(1024,fpFile);
		sscanf(sRead.c_str(),"(%lf,%lf,%lf) -> (%lf,%lf,%lf) : (%lf,%lf,%lf) -> (%lf,%lf,%lf)\n",&dX1,&dY1,&dZ1,&dX2,&dY2,&dZ2,&dNX1,&dNY1,&dNZ1,&dNX2,&dNY2,&dNZ2);
		oNode1.Set(dX1,dY1,dZ1);
		oNode2.Set(dX2,dY2,dZ2);
		oSurfaceNormal1.Set(dNX1,dNY1,dNZ1);
		oSurfaceNormal2.Set(dNX2,dNY2,dNZ2);
		if(oNode1.Distance(oNode2) < dTolerance)
		{
			continue;
		}
		oSurfaceNormal1.Normalize();
		oSurfaceNormal2.Normalize();

		poSurfaceLoop = new SurfaceLoop;
		poSurfaceLoop->Set(oNode1,oNode2,oBurgersVector,oSlipPlaneNormal,oSurfaceNormal1,oSurfaceNormal2);
		m_lpoLoops.push_back(poSurfaceLoop);
	}
	fclose(fpFile);
	printf("read %d segs\n",(unsigned int)m_lpoLoops.size());
}
void SurfaceDeformer::ReadSurfaceSegmentsForReduction(FILE* fpFile)
{
	ClearSurfaceSegments();
	string sRead;
	double dX1 = 0.0;
	double dY1 = 0.0;
	double dZ1 = 0.0;
	double dX2 = 0.0;
	double dY2 = 0.0;
	double dZ2 = 0.0;
	double dNX1 = 0.0;
	double dNY1 = 0.0;
	double dNZ1 = 0.0;
	double dNX2 = 0.0;
	double dNY2 = 0.0;
	double dNZ2 = 0.0;
	DislocationNode oNode1;
	DislocationNode oNode2;
	DislocationSegment oSegment;
	SurfaceSegment* poSurfaceSegment = NULL;
	unsigned int iDomain1 = 0;
	unsigned int iDomain2 = 0;
	unsigned int iIndex1 = 0;
	unsigned int iIndex2 = 0;
	
	while(!feof(fpFile))
	{
		sRead = GetRealString(1024,fpFile);
		if(sRead.empty())
		{
			break;
		}
		sscanf(sRead.c_str(),"%*d,%*d : (%d,%d) -> (%d,%d) : (%lf,%lf,%lf) / (%lf,%lf,%lf)\n",&iDomain1,&iIndex1,&iDomain2,&iIndex2,&dX1,&dY1,&dZ1,&dX2,&dY2,&dZ2);
		// skip real segments
		if((iDomain1 != iDomain2) || (iIndex1 != iIndex2))
		{
			sRead = GetRealString(1024,fpFile);
			continue;
		}
		oSegment.SetSlipPlaneNormal(Vector(dX1,dY1,dZ1));
		oSegment.SetBurgersVector(Vector(dX2,dY2,dZ2));
		sRead = GetRealString(1024,fpFile);
		sscanf(sRead.c_str(),"(%lf,%lf,%lf) -> (%lf,%lf,%lf) : (%lf,%lf,%lf) -> (%lf,%lf,%lf)\n",&dX1,&dY1,&dZ1,&dX2,&dY2,&dZ2,&dNX1,&dNY1,&dNZ1,&dNX2,&dNY2,&dNZ2);
		oNode1.Set(dX1,dY1,dZ1);
		oNode2.Set(dX2,dY2,dZ2);
		oNode1.SetSurfaceNormal(Vector(dNX1,dNY1,dNZ1));
		oNode2.SetSurfaceNormal(Vector(dNX2,dNY2,dNZ2));
		poSurfaceSegment = new SurfaceSegment;
		poSurfaceSegment->Set(oNode1,oNode2,oSegment);
		m_lpoSegments.push_back(poSurfaceSegment);
	}
}
void SurfaceDeformer::ReadSurfaceSegmentsForReduction(const string& sFileName)
{
	FILE* fpFile = fopen(sFileName.c_str(),"r");
	ReadSurfaceSegmentsForReduction(fpFile);
	fclose(fpFile);
}
bool SurfaceDeformer::ReduceSurfaceMotion()
{
	unsigned int iInitialSegmentsCount = (unsigned int)m_lpoSegments.size();
	printf(" : %d ",iInitialSegmentsCount);
	list<SurfaceSegment*>::iterator liOuterSegments = m_lpoSegments.begin();
	list<SurfaceSegment*>::iterator liInnerSegments;
	double dTolerance = 1.0E-6;
	liOuterSegments = m_lpoSegments.begin();
	while(liOuterSegments != m_lpoSegments.end())
	{
		if((*liOuterSegments)->GetLength() < dTolerance)
		{
			delete (*liOuterSegments);
			liOuterSegments = m_lpoSegments.erase(liOuterSegments);
		}
		else
		{
			liInnerSegments = liOuterSegments;
			liInnerSegments++;
			while(liInnerSegments != m_lpoSegments.end())
			{
				if((*liOuterSegments)->IsPostMergeable((*liInnerSegments)))
				{
					(*liOuterSegments)->PostMerge((*liInnerSegments));
					delete (*liInnerSegments);
					liInnerSegments = m_lpoSegments.erase(liInnerSegments);
				}
				else
				{
					liInnerSegments++;
				}
			}
			liOuterSegments++;
		}
	}
	unsigned int iFinalSegmentsCount = (unsigned int)m_lpoSegments.size();
	printf("---> %d\n",iFinalSegmentsCount);
	if(iInitialSegmentsCount == iFinalSegmentsCount)
	{
		return true;
	}
	return false;
}
unsigned int SurfaceDeformer::GenerateMotionChunks(const string& sFileName,const unsigned int& iChunkSize)
{
	FILE* fpFile = fopen(sFileName.c_str(),"r");
	if(fpFile == NULL)
	{
		return 0;
	}
	string sRead;
	unsigned int iDomain1 = 0;
	unsigned int iDomain2 = 0;
	unsigned int iIndex1 = 0;
	unsigned int iIndex2 = 0;

	unsigned int iChunkOrder = 0;
	string sRealSegmentsFileName = "realsegs.seg";
	string sCurrentChunkFileName = GetChunkFileName(iChunkOrder);

	FILE* fpRealSegmentsFile = fopen(sRealSegmentsFileName.c_str(),"w");
	FILE* fpCurrentChunkFile = fopen(sCurrentChunkFileName.c_str(),"w");
	unsigned int iSegmentsCount = 0;

	printf("in chunk file %d",iChunkOrder);
	while(!feof(fpFile))
	{
		sRead = GetRealString(1024,fpFile);
		if(sRead.empty())
		{
			printf("done, breaking the loop\n");
			fflush(NULL);
			break;
		}
		try
		{
			sscanf(sRead.c_str(),"%*d,%*d : (%d,%d) -> (%d,%d) : (%l*f,%l*f,%l*f) / (%l*f,%l*f,%l*f)\n",&iDomain1,&iIndex1,&iDomain2,&iIndex2);
		}
		catch(...)
		{
			// couldn't read first line, loop till we hit a valid first line
			while(true)
			{
				try
				{
					sRead = GetRealString(1024,fpFile);
					sscanf(sRead.c_str(),"%*d,%*d : (%d,%d) -> (%d,%d) : (%l*f,%l*f,%l*f) / (%l*f,%l*f,%l*f)\n",&iDomain1,&iIndex1,&iDomain2,&iIndex2);
				}
				catch(...)
				{
					continue;
				}
				break;
			}
			// then continue executing from this point on
		}
		if((iDomain1 == iDomain2) && (iIndex1 == iIndex2))
		{
			// the segment is actually a node motion, put it in the chunk file
			// first of all, make sure that the current chunk file has space for this segment
			if(iSegmentsCount >= iChunkSize)
			{
				// close the current chunk file
				fclose(fpCurrentChunkFile);
				// get a new chunk file name
				iChunkOrder = iChunkOrder + 1;
				sCurrentChunkFileName = GetChunkFileName(iChunkOrder);
				// open a new chunk file
				fpCurrentChunkFile = fopen(sCurrentChunkFileName.c_str(),"w");
				// reset the segments count
				iSegmentsCount = 0;
				printf(" ... done\n");
				printf("in chunk file %d",iChunkOrder);
			}
			// put the first line in the chunk file
			fprintf(fpCurrentChunkFile,"%s\n",sRead.c_str());
			sRead = GetRealString(1024,fpFile);
			// put the second line in the chunk file
			fprintf(fpCurrentChunkFile,"%s\n",sRead.c_str());
			iSegmentsCount = iSegmentsCount + 1;
		}
		else
		{
			// write real segments to the segments file
			// put the first line in the real segments file
			fprintf(fpRealSegmentsFile,"%s\n",sRead.c_str());
			sRead = GetRealString(1024,fpFile);
			// put the second line in the real segments file
			fprintf(fpRealSegmentsFile,"%s\n",sRead.c_str());
		}
	}
	fclose(fpFile);
	fclose(fpRealSegmentsFile);
	fclose(fpCurrentChunkFile);
	printf(" ... done\n");
	fflush(NULL);
	// return the number of the generated chunks
	return (iChunkOrder + 1);
}
bool SurfaceDeformer::ProcessChunks()
{
	int iProcessID = 0;
	int iChunkOrderStep = 1;
 	MPI_Comm_rank(MPI_COMM_WORLD,&iProcessID);
 	MPI_Comm_size(MPI_COMM_WORLD,&iChunkOrderStep);
 	int iChunkOrder = iProcessID;
	FILE* fpChunkFile = NULL;
	string sChunkFileName = "";
	char sOutputFileName[512];
	sprintf(sOutputFileName,"reduced_segs_p_%d.txt",iProcessID);
	FILE* fpSurfaceMotionFile = fopen(sOutputFileName,"w");
	bool bChunkDone = true;
	while(true)
	{
		sChunkFileName = GetChunkFileName(iChunkOrder);
		fpChunkFile = fopen(sChunkFileName.c_str(),"r");
		if(fpChunkFile == NULL)
		{
			// set the chunk order to the last file processed by this process, if the process has no
			// files to process, the chunk order will be negative
			iChunkOrder = iChunkOrder - iChunkOrderStep;
			break;
		}
		ReadSurfaceSegmentsForReduction(fpChunkFile);
		printf("in chunk number %d ",iChunkOrder);
		bChunkDone = ReduceSurfaceMotion();
		WriteSurfaceSegments(fpSurfaceMotionFile);
		fclose(fpChunkFile);
		// remove the file
		remove(sChunkFileName.c_str());
		iChunkOrder = iChunkOrder + iChunkOrderStep;
	}
	fclose(fpSurfaceMotionFile);
	// if this process didn't process any chunks, mark it as done
	if(iChunkOrder < 0)
	{
		bChunkDone = true;
	}
	MPI_Barrier(MPI_COMM_WORLD);
//	// get the maximum processed chunk order, this is the number of chunks in the system
// 	int iMaxChunkOrder = 0;
// 	MPI_Allreduce(&iChunkOrder,&iMaxChunkOrder,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
	// cat files
	if(iProcessID == 0)
	{
		system("cat reduced_segs_p_*.txt > reduced_segs.txt");
		remove("reduced_segs_p_*.txt");
	}
	MPI_Barrier(MPI_COMM_WORLD);
// 	if((iMaxChunkOrder % iChunkOrderStep) < iProcessID)
// 	{
// 		return false;
// 	}
// 	return true;
	// return the status of the last processed chunk, if the system is processing more than one
	// chunk at this time, then this status will never be used, otherwise, if it is processing 
	// only one status, all the processes will return true except for process 0 which will return true 
	// or false depending on the input and output segments counts. This will be processed in the calling
	// method and the decision to continue the recursive merging will be made accordingly
	return bChunkDone;
}
void SurfaceDeformer::CalculateDisplacements()
{
	unsigned int i = 0;
	unsigned int j = 0;
	list<SurfaceLoop*>::iterator liLoops;
	Vector oDisplacement;
	// determine the rows that will be calculated by this process
	unsigned int iNumberOfRows = m_iLengthResolution/m_iProcessesCount;
	unsigned int iFirstRow = 0;
	if(m_iLengthResolution%m_iProcessesCount > m_iProcessID)
	{
		iNumberOfRows = iNumberOfRows + 1;
		iFirstRow = m_iProcessID*iNumberOfRows;
	}
	else
	{
		iFirstRow = m_iProcessID*iNumberOfRows + m_iLengthResolution%m_iProcessesCount;
	}
	unsigned int iLastRow = iFirstRow + iNumberOfRows;
	for(i = iFirstRow ; i < iLastRow ; i++)
	{
		printf("proc %d : i : %d of %d\n",m_iProcessID,i,m_iLengthResolution);
		for(j = 0 ; j < m_iWidthResolution ; j++)
		{	
			oDisplacement.Set(0.0,0.0,0.0);
			for(liLoops = m_lpoLoops.begin() ; liLoops != m_lpoLoops.end() ; liLoops++)
			{
				oDisplacement = oDisplacement + (*liLoops)->GetPointDisplacement(m_vvoSurfacePoints[i][j]);
			}
			m_vvoDisplacements[i][j] = oDisplacement;
		}
	}
}
void SurfaceDeformer::WriteOutput() const
{
	FILE* fpFile = NULL;
	unsigned int i = 0;
	unsigned int j = 0;
	unsigned int k = 0;
	Point oPoint;
	Vector oDisplacement;
	// determine the rows that will be calculated by this process
	unsigned int iNumberOfRows = m_iLengthResolution/m_iProcessesCount;
	unsigned int iFirstRow = 0;
	if(m_iLengthResolution%m_iProcessesCount > m_iProcessID)
	{
		iNumberOfRows = iNumberOfRows + 1;
		iFirstRow = m_iProcessID*iNumberOfRows;
	}
	else
	{
		iFirstRow = m_iProcessID*iNumberOfRows + m_iLengthResolution%m_iProcessesCount;
	}
	unsigned int iLastRow = iFirstRow + iNumberOfRows;

	for(k = 0 ; k < m_iProcessesCount ; k++)
	{
		if(k == m_iProcessID)
		{
			if(k == 0)
			{
				fpFile = fopen(m_sOutputFile.c_str(),"w");
			}
			else
			{
				fpFile = fopen(m_sOutputFile.c_str(),"a");
			}
			for(i = iFirstRow ; i < iLastRow ; i++)
			{
				for(j = 0 ; j < m_iWidthResolution ; j++)
				{	
					oPoint = m_vvoSurfacePoints[i][j];
					oDisplacement = m_vvoDisplacements[i][j];
					fprintf(fpFile,"%e,%e,%e : %e,%e,%e\n",oPoint.GetX(),oPoint.GetY(),oPoint.GetZ(),oDisplacement.GetX(),oDisplacement.GetY(),oDisplacement.GetZ());
				}
			}
			fclose(fpFile);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}
void SurfaceDeformer::Initialize()
{
	m_sSurfaceSegmentsFile = "";
	m_sOutputFile = "";
	m_oPlanePoint.Set(0.0,0.0,0.0);
	m_oPlaneNormal.Set(0.0,0.0,0.0);
	m_oPlaneX.Set(0.0,0.0,0.0);
	m_dLength = 0.0;
	m_dWidth = 0.0;
	m_dVirtualNodeSpacing = 0.0;
	m_dPoissonsRatio = 0.0;
	m_iLengthResolution = 0;
	m_iWidthResolution = 0;
	m_iGaussPointsCount = 0;
	m_iProcessesCount = 1;
	m_iProcessID = 0;
	m_vvoSurfacePoints.clear();
	m_vvoDisplacements.clear();
	m_lpoSegments.clear();
	m_lpoLoops.clear();
	m_lpoTimePoints.clear();
}
void SurfaceDeformer::ClearSurfaceSegments()
{
	list<SurfaceSegment*>::iterator liSegments;
	for(liSegments = m_lpoSegments.begin() ; liSegments != m_lpoSegments.end() ; liSegments++)
	{
		if((*liSegments) != NULL)
		{
			delete (*liSegments);
		}
	}
	m_lpoSegments.clear();

	list<SurfaceLoop*>::iterator liLoops;
	for(liLoops = m_lpoLoops.begin() ; liLoops != m_lpoLoops.end() ; liLoops++)
	{
		if((*liLoops) != NULL)
		{
			delete (*liLoops);
		}
	}
	m_lpoLoops.clear();
}
void SurfaceDeformer::WriteSurfaceSegments(FILE* fpFile) const
{
	list<SurfaceSegment*>::const_iterator liSegments;
	for(liSegments = m_lpoSegments.begin() ; liSegments != m_lpoSegments.end() ; liSegments++)
	{
		(*liSegments)->Write(fpFile);
	}
}
string SurfaceDeformer::GetChunkFileName(const unsigned int& iChunkOrder) const
{
	char cWrite[256];
	sprintf(cWrite,"motion_%d.mtn",iChunkOrder);
	string sName = cWrite;
	return sName;
}
string SurfaceDeformer::AddTimePoint(const unsigned int& iTimeStep,const string& sBaseFileName)
{
	char cFileName[256];
	sprintf(cFileName,"%s_%d",sBaseFileName.c_str(),iTimeStep);
	string sFileName = string(cFileName);
	TimePoint* poPoint = new TimePoint(iTimeStep,sFileName);
	m_lpoTimePoints.push_back(poPoint);
	return sFileName;
}
void SurfaceDeformer::ReduceSegmentGroup(string& sInputFileName,const unsigned int& iOutputOrder)
{
	unsigned int iChunkSizeIncrement = 3000;
	unsigned int iInitialChunkSize = 3000;
	unsigned int iChunkSize = iInitialChunkSize - iChunkSizeIncrement;
	bool bIsDone = false;
	int iProcessID = 0;
	MPI_Comm_rank(MPI_COMM_WORLD,&iProcessID);
	int iIsDone = 0;
	unsigned int iChunksCount = 0;
	while(true)
	{
		iChunkSize = iChunkSize + iChunkSizeIncrement;
		if(iProcessID == 0)
		{
			if(iChunkSize > iInitialChunkSize)
			{
				// this is not the first iteration
				// rename the last output file and use it as the new input
				rename("reduced_segs.txt","reduced_segs.txt_old");
				sInputFileName = "reduced_segs.txt_old";
			}
			printf("processing file : %s\n",sInputFileName.c_str());
			fflush(NULL);
			iChunksCount = GenerateMotionChunks(sInputFileName,iChunkSize);
		}
		MPI_Bcast(&iChunksCount,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		bIsDone = ProcessChunks();
		MPI_Barrier(MPI_COMM_WORLD);
		if(iProcessID == 0)
		{
			if(iChunkSize == iInitialChunkSize)
			{
				rename("realsegs.seg","realsegs.seg_00");
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		// if there is only one chunk (no matter what the number of processes is) and this chunk is done
		// (meaning, the input and output segments counts are the same) break out of the loop
		// to determine this :
		// 1. if there are more than one file to process, don't bother checking because we are not
		// done merging all the files yet
		// 2. if there is only one chunk, all the processes will have nothing to process except for process 0
		// 3. if process zero returns true for the input and output segments count, then we are done
		// 4. otherwise, continue
		
		if(iChunksCount <= 1)
		{
			// there is only one file to process, see if it was fully processed or not, since we only
			// care about process 0 here, we just need to broadcast this value from process 0 to all
			// the other processors
			iIsDone = 1;
			if(!bIsDone)
			{
				iIsDone = 0;
			}
			MPI_Bcast(&iIsDone,1,MPI_INT,0,MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			if(iIsDone == 1)
			{
				break;
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(iProcessID == 0)
	{
		// now since we are done with the reductions, remove the old reduced segs file
		remove("reduced_segs.txt_old");
		// rename the final results files
		char cOutputFileName[512];
		sprintf(cOutputFileName,"process_%d.mtn",iOutputOrder);
		rename("reduced_segs.txt",cOutputFileName);
		sprintf(cOutputFileName,"process_%d.seg",iOutputOrder);
		rename("realsegs.seg_00",cOutputFileName);
		remove("realsegs.seg");
	}
	MPI_Barrier(MPI_COMM_WORLD);
}
void SurfaceDeformer::CollectSurfaceMotion(const string& sInputFileName)
{
	FILE* fpInput = fopen(sInputFileName.c_str(),"r");
	string sRead;
	unsigned int iTimeStep = 0;
	list< TimePoint* >::iterator liPoint;
	while(!feof(fpInput))
	{
		sRead = GetRealString(1024,fpInput);
		if(sRead.empty())
		{
			break;
		}
		sscanf(sRead.c_str(),"%d,%*d : (%*d,%*d) -> (%*d,%*d) : (%*f,%*f,%*f) / (%*f,%*f,%*f)\n",&iTimeStep);
		for(liPoint = m_lpoTimePoints.begin() ; liPoint != m_lpoTimePoints.end() ; liPoint++)
		{
			if((*liPoint)->IsBefore(iTimeStep))
			{
				// write the two lines in the time point file if the time step is before the time point
				sRead = sRead + "\n";
				(*liPoint)->Write(sRead);
				sRead = GetRealString(1024,fpInput);
				sRead = sRead + "\n";
				(*liPoint)->Write(sRead);
				break;
			}
		}
	}
	// when done with this input file, flush all the file buffers
	for(liPoint = m_lpoTimePoints.begin() ; liPoint != m_lpoTimePoints.end() ; liPoint++)
	{
		(*liPoint)->FlushFile();
	}
	fclose(fpInput);
}
unsigned int SurfaceDeformer::GetProcessID() const
{
	return m_iProcessID;
}
void SurfaceDeformer::ClearTimePoints()
{
	list< TimePoint* >::iterator liPoint;
	for(liPoint = m_lpoTimePoints.begin() ; liPoint != m_lpoTimePoints.end() ; liPoint++)
	{
		(*liPoint)->CloseFile();
		delete (*liPoint);
	}
	m_lpoTimePoints.clear();
}

void Filter(SurfaceDeformer* poDeformer,unsigned int iGroupOrder,const unsigned int& iMaxGroupOrder,const char* cBaseFileName)
{
	char cFileName[512];
	string sFileName;
	if(poDeformer->GetProcessID() == 0)
	{
		while(iGroupOrder < iMaxGroupOrder)
		{
			sprintf(cFileName,"%s_%d.txt",cBaseFileName,iGroupOrder);
			sFileName = cFileName;
			printf("processing file %s\n",cFileName);
			poDeformer->CollectSurfaceMotion(sFileName);
			iGroupOrder = iGroupOrder + 1;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

void Reduce(SurfaceDeformer* poDeformer,unsigned int iGroupOrder,const unsigned int& iMaxGroupOrder,const char* cBaseFileName)
{
	char cFileName[512];
	string sFileName;
	while(iGroupOrder < iMaxGroupOrder)
	{
		sprintf(cFileName,"%s_%d.txt",cBaseFileName,iGroupOrder);
		sFileName = cFileName;
		poDeformer->ReduceSegmentGroup(sFileName,iGroupOrder);
		MPI_Barrier(MPI_COMM_WORLD);
		iGroupOrder = iGroupOrder + 1;
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

void CalculateDisplacement(SurfaceDeformer* poDeformer)
{
	poDeformer->GenerateSurfacePoints();
	poDeformer->ReadSurfaceSegments();
	poDeformer->CalculateDisplacements();
	poDeformer->WriteOutput();
}


int main(int argc,char** argv)
{	
	if(argc < 2)
	{
		printf("error: few arguments\n");
		printf("usage: surfdef -[sf | sr | sd | sac | ssc | sspd | gf | gac | gsc | gspd] ...\n");
		return 1;
	}
	
	int iProcessID = 0;
	int iProcessesCount = 1;
 	MPI_Init(&argc,&argv);
 	MPI_Comm_rank(MPI_COMM_WORLD,&iProcessID);
 	MPI_Comm_size(MPI_COMM_WORLD,&iProcessesCount);
	
	SurfaceDeformer* poDeformer = SurfaceDeformer::GetInstance();
	poDeformer->SetParallelizationParameters(iProcessesCount,iProcessID);
	string sRunType = string(argv[1]);
	if(sRunType.compare("-sf") == 0)
	{
		if(argc != 6)
		{
			printf("error: wrong number of arguments\n");
			printf("usage: surfdef -sf max_group_order first_group_order base_file_name max_time_step\n");
			return 2;
		}
		poDeformer->AddTimePoint(atoi(argv[5]),"filtered_surfseg");
		Filter(poDeformer,atoi(argv[3]),atoi(argv[2]),argv[4]);
	}
	else if(sRunType.compare("-sr") == 0)
	{
		if(argc != 5)
		{
			printf("error: wrong number of arguments\n");
			printf("usage: surfdef -sr max_group_order first_group_order base_file_name\n");
			return 3;
		}
		Reduce(poDeformer,atoi(argv[3]),atoi(argv[2]),argv[4]);
	}
	else if(sRunType.compare("-sd") == 0)
	{
		if(argc != 3)
		{
			printf("error: wrong number of arguments\n");
			printf("usage: surfdef -sd input_file\n");
			return 4;
		}
		poDeformer->Read(argv[2]);
		CalculateDisplacement(poDeformer);
	}
	else if(sRunType.compare("-sac") == 0)
	{
		if(argc != 3)
		{
			printf("error: wrong number of arguments\n");
			printf("usage: surfdef -sac input_file\n");
			return 5;
		}
		SurfaceProfile oProfile;
		oProfile.Read(argv[2]);
		oProfile.GetAmplitudeCharacteristics();
	}
	else if(sRunType.compare("-ssc") == 0)
	{
		if(argc != 3)
		{
			printf("error: wrong number of arguments\n");
			printf("usage: surfdef -ssc input_file\n");
			return 5;
		}
		SurfaceProfile oProfile;
		oProfile.Read(argv[2]);
		oProfile.GetSpatialCharacteristics();
	}
	else if(sRunType.compare("-sspd") == 0)
	{
		if(argc != 3)
		{
			printf("error: wrong number of arguments\n");
			printf("usage: surfdef -sspd input_file\n");
			return 5;
		}
		SurfaceProfile oProfile;
		oProfile.Read(argv[2]);
		oProfile.Shift();
		unsigned int iSize = 0;
		double* pdSPD = oProfile.GetPowerSpectrum(iSize);
		unsigned int i = 0;
		printf("%d\n",iSize);
		for(i = 0 ; i < iSize ; i++)
		{
			printf("%d\t\t%25.20f\n",i,pdSPD[i]);
		}
		delete [] pdSPD;
	}
	else if(sRunType.compare("-shd") == 0)
	{
		if((argc != 4) && (argc != 5))
		{
			printf("error: wrong number of arguments\n");
			printf("usage: surfdef -shd input_file xsize [ysize]\n");
			return 5;
		}
		SurfaceProfile oProfile;
		oProfile.Read(argv[2]);
		double dXHalfSize = 0.5*atof(argv[3]);
		double dYHalfSize = dXHalfSize;
		if(argc == 5)
		{
			dYHalfSize = 0.5*atof(argv[4]);
		}
		oProfile.SetXBounds(-dXHalfSize,dXHalfSize);
		oProfile.SetYBounds(-dYHalfSize,dYHalfSize);
		printf("%lf\n",oProfile.GetHausdorffDimension());
	}
	else if(sRunType.compare("-gf") == 0)
	{
		if(argc != 5)
		{
			printf("error: wrong number of arguments\n");
			printf("usage: surfdef -gf cycles_file max_group_order base_file_path\n");
			return 6;
		}
		FILE* fpCyclesFile = fopen(argv[2],"r");
		unsigned int iMaxGroupOrder = atoi(argv[3]);
		string sBaseFilePath = string(argv[4]);
		// read through the cycles file
		string sRead;
		int iTimePoint = 0;
		char cCommandString[1024];
		
		string sFileName;
		list<string> lsFileNames;
		while(!feof(fpCyclesFile))
		{
			sRead = GetRealString(1024,fpCyclesFile);
			if(sRead.empty())					break;
			sscanf(sRead.c_str(),"%*s : %d : %*f\n",&iTimePoint);
			// add the time point and save the file name for later moving
			sFileName = poDeformer->AddTimePoint(iTimePoint,"filtered_surfseg");
			lsFileNames.push_back(sFileName);
		}
		fclose(fpCyclesFile);
		// call the filtering function
		Filter(poDeformer,0,iMaxGroupOrder,sBaseFilePath.c_str());
		poDeformer->ClearTimePoints();
		// when done, create directories and move the filtered files to them
		unsigned int iCount = 0;
		list<string>::iterator liFiles;
		for(liFiles = lsFileNames.begin() ; liFiles != lsFileNames.end() ; liFiles++)
		{
			iCount++;
			sprintf(cCommandString,"mkdir reduced_%d",iCount);
			system(cCommandString);
			sprintf(cCommandString,"mv %s reduced_%d//filtered_0.txt",(*liFiles).c_str(),iCount);
			system(cCommandString);
		}
		lsFileNames.clear();
	}
	else if(sRunType.compare("-gac") == 0)
	{
		if(argc != 3)
		{
			printf("error: wrong number of arguments\n");
			printf("usage: surfdef -gac count\n");
			return 7;
		}
		SurfaceProfile oProfile;
		unsigned int i = 0;
		unsigned int j = 0;
		unsigned int iCount = atoi(argv[2]);
		char cTemp[256];
		string surfaces[6] = {"nx","ny","nz","px","py","pz"};
		string sOutputFile = "";
		for(i = iProcessID ; i < 6 ; i = i + iProcessesCount)
		{
			sprintf(cTemp,"%s_amp.dat",surfaces[i].c_str());
			sOutputFile = string(cTemp);
			for(j = 1 ; j <= iCount ; j++)
			{
				sprintf(cTemp,"%s_%d.top",surfaces[i].c_str(),j);
				printf("processing %s\n",cTemp);
				oProfile.Read(cTemp);
				oProfile.GetAmplitudeCharacteristics(sOutputFile);
			}
		}
	}
	else if(sRunType.compare("-gsc") == 0)
	{
		if(argc != 3)
		{
			printf("error: wrong number of arguments\n");
			printf("usage: surfdef -gsc count\n");
			return 8;
		}
		SurfaceProfile oProfile;
		unsigned int i = 0;
		unsigned int j = 0;
		unsigned int iCount = atoi(argv[2]);
		char cTemp[256];
		string surfaces[6] = {"nx","ny","nz","px","py","pz"};
		string sOutputFile = "";
		for(i = iProcessID ; i < 6 ; i = i + iProcessesCount)
		{
			sprintf(cTemp,"%s_spa.dat",surfaces[i].c_str());
			sOutputFile = string(cTemp);
			for(j = 1 ; j <= iCount ; j++)
			{
				sprintf(cTemp,"%s_%d.top",surfaces[i].c_str(),j);
				printf("processing %s\n",cTemp);
				oProfile.Read(cTemp);
				oProfile.GetSpatialCharacteristics(sOutputFile);
			}
		}
	}
	else if(sRunType.compare("-gspd") == 0)
	{
	
	}
	else
	{
		printf("error: unknown run type %s\n",sRunType.c_str());
		printf("usage: surfdef -[sf | sr | sd | gf | gr | gd] ...\n");
		return 10;
	}
	
	delete poDeformer;
 	MPI_Finalize();
	return 0;
}

