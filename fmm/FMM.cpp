#include "iostream"
#include "FMMTreeNode.h"
#include "Randomizer.h"

using namespace FMM;
using namespace EZ;


int main(int argc,char** argv)
{
	printf("starting FMM calculation\n");
	FMMTreeNode oRoot;
	double dXRange = 1000.0;
	double dYRange = 1000.0;
	double dZRange = 1000.0;
	unsigned int iDepth = 3;
	oRoot.SetBox(dXRange,dYRange,dZRange);
	oRoot.Build(iDepth);
	
	// generate the sources and targets
	unsigned int iPointsCount = 10000;
	unsigned int i = 0;
	list< Point > loPoints;
	list< Point* > lpoSources;
	list< Point* > lpoTargets;
	double dX = 0.0;
	double dY = 0.0;
	double dZ = 0.0;
	for(i = 0 ; i < iPointsCount ; i++)
	{
		dX = Randomizer::Random(-0.5*dXRange,0.5*dXRange);
		dY = Randomizer::Random(-0.5*dYRange,0.5*dYRange);
		dZ = Randomizer::Random(-0.5*dZRange,0.5*dZRange);
		loPoints.push_back(Point(dX,dY,dZ));
		lpoSources.push_back(&loPoints.back());
		lpoTargets.push_back(&loPoints.back());
	}
	
	// distribute the sources and targets over the tree nodes
	oRoot.ClaimSources(&lpoSources);
	oRoot.ClaimTargets(&lpoTargets);
	return 0;
}


