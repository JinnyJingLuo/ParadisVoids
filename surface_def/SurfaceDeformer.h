#ifndef SURFACEDEFORMER_H_
#define SURFACEDEFORMER_H_

#include "string"
#include "vector"
#include "list"
#include "Vector.h"
#include "SurfaceSegment.h"
#include "SurfaceLoop.h"

using namespace std;
using namespace EZ;

class TimePoint
{
public:
	TimePoint();
	TimePoint(const unsigned int& iTimeStep,const string& sFileName);
	TimePoint(const TimePoint& oPoint);
	~TimePoint();
	TimePoint& operator=(const TimePoint& oPoint);
	void Reset();
	void OpenFile(const string& sFileName);
	void CloseFile();
	void FlushFile();
	void Write(const string& sString);
	bool IsBefore(const unsigned int& iTimeStep);
	void SetTimeStep(const unsigned int& iTimeStep);
	
private:

protected:
	void Initialize();
	unsigned int m_iTimeStep;
	FILE* fpFilterFile;
};

class SurfaceDeformer
{
public:
	static SurfaceDeformer* GetInstance();
	~SurfaceDeformer();
	void Reset();
	void SetParallelizationParameters(const int& iProcessesCount,const int& iProcessID);
	void Read(const string& sInputFile);
	void Print() const;
	void GenerateSurfacePoints();
	void PrintSurfacePoints() const;
	void ReadSurfaceSegments();
	void ReadSurfaceSegmentsForReduction(const string& sFileName);
	void ReadSurfaceSegmentsForReduction(FILE* fpFile);
	unsigned int GenerateMotionChunks(const string& sFileName,const unsigned int& iChunkSize = 10000);
	bool ReduceSurfaceMotion();
	void CalculateDisplacements();
	void WriteOutput() const;
	bool ProcessChunks();
	void ReduceSegmentGroup(string& sInputFileName,const unsigned int& iOutputOrder);
	string AddTimePoint(const unsigned int& iTimeStep,const string& sBaseFileName);
	void CollectSurfaceMotion(const string& sInputFileName);
	unsigned int GetProcessID() const;
	void ClearTimePoints();
	
private:

protected:
	static SurfaceDeformer* m_poInstance;
	SurfaceDeformer();
	void Initialize();
	void ClearSurfaceSegments();
	void WriteSurfaceSegments(FILE* fpFile) const;
	string GetChunkFileName(const unsigned int& iChunkOrder) const;
	
	string m_sSurfaceSegmentsFile;
	string m_sOutputFile;
	Point m_oPlanePoint;
	Vector m_oPlaneNormal;
	Vector m_oPlaneX;
	double m_dLength;
	double m_dWidth;
	double m_dVirtualNodeSpacing;
	double m_dPoissonsRatio;
	unsigned int m_iLengthResolution;
	unsigned int m_iWidthResolution;
	unsigned int m_iGaussPointsCount;
	int m_iProcessesCount;
	int m_iProcessID;
	list< TimePoint* > m_lpoTimePoints;

	vector< vector<Point> > m_vvoSurfacePoints;
	vector< vector<Point> > m_vvoDisplacements;
	list<SurfaceSegment*> m_lpoSegments;
	list<SurfaceLoop*> m_lpoLoops;
};


#endif


