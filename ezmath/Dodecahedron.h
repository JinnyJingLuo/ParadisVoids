// Paradis Processor Project
// Ahmed M. Hussein (mailto : am.hussin@gmail.com)
// June 2012

#ifndef DODECAHEDRON_H_
#define DODECAHEDRON_H_

#include "stdio.h"
#include "list"
#include "TriPatch.h"
#include "Geometry.h"

using namespace std;

namespace GeometrySystem
{
	class Dodecahedron : public Geometry
	{
	public:
		Dodecahedron();
		Dodecahedron(const Dodecahedron& oDodecahedron);
		~Dodecahedron();
		Dodecahedron& operator=(const Dodecahedron& oDodecahedron);
		void Reset();
		void Set(const double& dDiameter,const bool& bHalfGrain = false);
		double GetDiameter() const;
		list<GenericNode*>* GetTriangulationPoints();
		list<TriPatch*>* GetTriangles();
		bool IsPointInside(const Point& oPoint,const double& dTolerance = 1.0E-6) const;
		void WriteParaviewFile(const string& sFileName) const;
		double GetVolume() const;
		void WriteTriangulation(const string& sFileName) const;
		void SetHalfGrain(const bool& bHalfGrain = false);
		void SetOrigin(const Point& oOrigin);
		void ReleaseTriangulations();
		virtual Geometry* Clone();
		
	private:
	
	protected:
		void Initialize();
		void GenerateTriangulations();
		double m_dDiameter;
		list<GenericNode*> m_lpoTriangulationPoints;
		list<TriPatch*> m_lpoTriangles;
		bool m_bIsHalfGrain;
	};
}


#endif



