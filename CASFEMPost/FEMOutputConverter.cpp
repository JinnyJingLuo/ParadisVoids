#include "iostream"
#include "Tools.h"
#include "list"
#include "vector"


using namespace std;
using namespace SupportSystem;

int main(int argc,char** argv)
{
	if(argc != 4)
	{
		printf("incorrect number of arguments\n");
		exit(1);
	}

	string sRead;
	string sInputFileName = argv[1];
	string sOutputFileName = argv[2];
	double dDisplacementMagnificationFactor = (double)atof(argv[3]);
	
	FILE* fpInput = fopen(sInputFileName.c_str(),"r");
	FILE* fpOutput = fopen(sOutputFileName.c_str(),"w");

	// write file header
	fprintf(fpOutput,"# vtk DataFile Version 3.0\n");
	fprintf(fpOutput,"output from file %s\n",sInputFileName.c_str());
	fprintf(fpOutput,"ASCII\n");
	fprintf(fpOutput,"DATASET UNSTRUCTURED_GRID\n");

	// read and write the nodes
	unsigned int iNodesCount = 0;
	sRead = GetRealString(500,fpInput);
	sscanf(sRead.c_str(),"%d\n",&iNodesCount);
	fprintf(fpOutput,"POINTS %d float\n",iNodesCount);
	unsigned int i = 0;
	// it is far more efficient (memory wise) to use lists instead of arrays here
	double dX = 0.0;
	double dY = 0.0;
	double dZ = 0.0;
	double dU = 0.0;
	double dQ = 0.0;
	double dQX = 0.0;
	double dQY = 0.0;
	double dQZ = 0.0;

	list<double> ldU;
	list<double> ldQ;
	list<double> ldQX;
	list<double> ldQY;
	list<double> ldQZ;

	//unsigned int iNodesCategory = 1;
	
	for(i = 0 ; i < iNodesCount ; i++)
	{
		sRead = GetRealString(500,fpInput);
		sscanf(sRead.c_str(),"%lf\t\t%lf\t\t%lf\t\t%lf\t\t%lf\n",&dX,&dY,&dZ,&dU,&dQ);
		sRead = GetRealString(500,fpInput);
		sscanf(sRead.c_str(),"%lf\t\t%lf\t\t%lf\n",&dQX,&dQY,&dQZ);
		
		ldU.push_back(dU);
		ldQ.push_back(dQ);
		ldQX.push_back(dQX);
		ldQY.push_back(dQY);
		ldQZ.push_back(dQZ);
	
		fprintf(fpOutput,"%e\t\t%e\t\t%e\n",dX,dY,dZ);
	}
	
	// read and write the elements
	unsigned int iElementsCount = 0;
	sRead = GetRealString(500,fpInput);
	sscanf(sRead.c_str(),"%d\n",&iElementsCount);
	unsigned int j = 0;
	
	// write the elements connectivity
	unsigned int iNodesPerElement = 20;
	fprintf(fpOutput,"CELLS %d %d\n",iElementsCount,(iNodesPerElement + 1)*iElementsCount);
	vector<unsigned int> viNodeIndices;
	viNodeIndices.resize(iNodesPerElement);
	unsigned int k = 0;

	// go over the elements data, read and write them
	for(i = 0 ; i < iElementsCount ; i++)
	{
		sRead = GetRealString(500,fpInput);
		sscanf(sRead.c_str(),"%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\n",&viNodeIndices[0],&viNodeIndices[1],&viNodeIndices[2],&viNodeIndices[3],&viNodeIndices[4],&viNodeIndices[5],&viNodeIndices[6],&viNodeIndices[7],&viNodeIndices[8],&viNodeIndices[9],&viNodeIndices[10],&viNodeIndices[11],&viNodeIndices[12],&viNodeIndices[13],&viNodeIndices[14],&viNodeIndices[15],&viNodeIndices[16],&viNodeIndices[17],&viNodeIndices[18],&viNodeIndices[19]);
		fprintf(fpOutput,"20\t\t");			
		for(k = 0 ; k < iNodesPerElement ; k++)
		{
			fprintf(fpOutput,"%d\t\t",viNodeIndices[k] - 1);
		}			
		fprintf(fpOutput,"\n");
	}
	
	// write the elements types
	fprintf(fpOutput,"CELL_TYPES %d\n",iElementsCount);
	for(j = 0 ; j < iElementsCount ; j++)
	{
		fprintf(fpOutput,"25\n");
	}

	// now write the nodal displacements and forces

	fprintf(fpOutput,"POINT_DATA %d\n",iNodesCount);
	fprintf(fpOutput,"SCALARS node_id float 1\n");
	fprintf(fpOutput,"LOOKUP_TABLE default\n");
	for(i = 0 ; i < iNodesCount ; i++)
	{
		fprintf(fpOutput,"%f\n",(float)i);
	}

	fprintf(fpOutput,"SCALARS potential float 1\n");
	fprintf(fpOutput,"LOOKUP_TABLE default\n");
	for(i = 0 ; i < iNodesCount ; i++)
	{
		dU = ldU.front();
		ldU.pop_front();
		fprintf(fpOutput,"%f\n",dU);
	}
	ldU.clear();

	fprintf(fpOutput,"SCALARS volume_flux float 1\n");
	fprintf(fpOutput,"LOOKUP_TABLE default\n");
	for(i = 0 ; i < iNodesCount ; i++)
	{
		dQ = ldQ.front();
		ldQ.pop_front();
		fprintf(fpOutput,"%f\n",dQ);
	}
	ldQ.clear();
	
	fprintf(fpOutput,"VECTORS flux_vector float\n");
	for(i = 0 ; i < iNodesCount ; i++)
	{
		dQX = ldQX.front();
		dQY = ldQY.front();
		dQZ = ldQZ.front();
		ldQX.pop_front();
		ldQY.pop_front();
		ldQZ.pop_front();
		fprintf(fpOutput,"%e\t\t%e\t\t%e\n",dQX,dQY,dQZ);
	}
	ldQX.clear();
	ldQY.clear();
	ldQZ.clear();

	fclose(fpInput);
	fclose(fpOutput);

	return 0;
}




