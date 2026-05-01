// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef PATCH_H_
#define PATCH_H_

#include "vector"
using namespace std;

namespace GeometrySystem
{
	class Patch
	{
	public:
		~Patch();
		virtual vector<Patch*> GenerateTriangulation() const = 0;
	protected:
	
	private:
	
	};
}

#endif

