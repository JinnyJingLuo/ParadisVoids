// Paradis Processor Project
// Ahmed M. Hussein (mailto : am.hussin@gmail.com)
// June 2012

#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include "TriPatch.h"

namespace GeometrySystem
{
	class Boundary
	{
	public:
		Boundary();
		~Boundary();
		Boundary(const Boundary& oBoundary);
		Boundary& operator=(const Boundary& oBoundary);
		void Reset();
		void Set(const string& sFileName);
		void Set(list<GenericNode*>* plpoPoints,list<TriPatch*>* plpoTriangles);
		list<TriPatch*>* GetTriangles();
		list<GenericNode*>* GetTriangulationPoints();
		void Refine(const unsigned int& iRefinementsCount = 1);
		void WriteTriangles(FILE* fpFile);
		void WriteTriangles(const string& sFileName);
		void WriteParaviewTriangles(const string& sFilePath,const unsigned int& iPointsCategory,const unsigned int& iTrianglesCategory);
		list<Point*> GenerateSurfacePoints(const unsigned int& iResolution);
	private:
		
	protected:
	 	void Initialize();
		void RefineTriangulation();
		list<TriPatch*> m_lpoTriangles;
 		list<GenericNode*> m_lpoTriangulationPoints;
	};
}

#endif


