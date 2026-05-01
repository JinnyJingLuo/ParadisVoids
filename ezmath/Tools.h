// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef TOOLS_H_
#define TOOLS_H_

#include "string"
#include "Vector.h"
#include "Matrix.h"
#include "Randomizer.h"

using namespace std;
using namespace EZ;

namespace SupportSystem
{
	string GetRealString(unsigned int iMaxSize,FILE* fpFile);
	void PrintOnScreen(const string& sStatement);
	template<typename Type> int Partition(vector<Type>& vtArray,unsigned const int& iStartIndex,unsigned const int& iEndIndex,const unsigned int& iPivotIndex)
	{
		// get the pivot value (comparison value)
		Type tPivotValue = vtArray[iPivotIndex];
		// place it at the end of the array
		Type tTempValue = tPivotValue;
		vtArray[iPivotIndex] = vtArray[iEndIndex];
		vtArray[iEndIndex] = tTempValue;

		unsigned int iNewPivotIndex = iStartIndex;
		unsigned int i = 0;
		for(i = iStartIndex ; i < iEndIndex ; i++)
		{
			if(vtArray[i] < tPivotValue)
			{
				// swap elements if required
				tTempValue = vtArray[iNewPivotIndex];
				vtArray[iNewPivotIndex] = vtArray[i];
				vtArray[i] = tTempValue;
				iNewPivotIndex = iNewPivotIndex + 1;
			}
		}
		// put the pivot in its new location
		tTempValue = vtArray[iNewPivotIndex];
		vtArray[iNewPivotIndex] = vtArray[iEndIndex];
		vtArray[iEndIndex] = tTempValue;
		// return the new pivot index
		return iNewPivotIndex;
	}
	template<typename Type> int Partition(vector<Type>& vtArray,vector<unsigned int>& viSortedIndices,unsigned const int& iStartIndex,unsigned const int& iEndIndex,const unsigned int& iPivotIndex)
	{
		// get the pivot value (comparison value)
		Type tPivotValue = vtArray[iPivotIndex];
		
		unsigned int iPivotValue = viSortedIndices[iPivotIndex];
		
		// place it at the end of the array
		Type tTempValue = tPivotValue;
		vtArray[iPivotIndex] = vtArray[iEndIndex];
		vtArray[iEndIndex] = tTempValue;
		
		unsigned int iTempValue = iPivotValue;
		viSortedIndices[iPivotIndex] = viSortedIndices[iEndIndex];
		viSortedIndices[iEndIndex] = iTempValue;
		
		unsigned int iNewPivotIndex = iStartIndex;
		unsigned int i = 0;
		for(i = iStartIndex ; i < iEndIndex ; i++)
		{
			if(vtArray[i] < tPivotValue)
			{
				// swap elements if required
				tTempValue = vtArray[iNewPivotIndex];
				vtArray[iNewPivotIndex] = vtArray[i];
				vtArray[i] = tTempValue;
				
				iTempValue = viSortedIndices[iNewPivotIndex];
				viSortedIndices[iNewPivotIndex] = viSortedIndices[i];
				viSortedIndices[i] = iTempValue;
				
				iNewPivotIndex = iNewPivotIndex + 1;
			}
		}
		// put the pivot in its new location
		tTempValue = vtArray[iNewPivotIndex];
		vtArray[iNewPivotIndex] = vtArray[iEndIndex];
		vtArray[iEndIndex] = tTempValue;
		
		iTempValue = viSortedIndices[iNewPivotIndex];
		viSortedIndices[iNewPivotIndex] = viSortedIndices[iEndIndex];
		viSortedIndices[iEndIndex] = iTempValue;
		// return the new pivot index
		return iNewPivotIndex;
	}
	template<typename Type> void QuickSort(vector<Type>& vtArray,unsigned const int& iStartIndex,unsigned const int& iEndIndex)
	{
		// check the problem validity
		unsigned int iSize = (unsigned int)vtArray.size();
		if(iSize <= 1)
		{
			return;
		}
		
		if(iStartIndex >= iEndIndex)
		{
			return;
		}
		
		if(iEndIndex >= iSize)
		{
			return;
		}
		
		// pick a random pivot
		unsigned int iPivotIndex = (unsigned int)(Randomizer::RandomInteger(iStartIndex,iEndIndex));
 		// partition the array into 2 arrays based on the chosen pivot
 		unsigned int iNewPivotIndex = Partition(vtArray,iStartIndex,iEndIndex,iPivotIndex);
 		// sort each array recursively
 		QuickSort(vtArray,iStartIndex,iNewPivotIndex - 1);
 		QuickSort(vtArray,iNewPivotIndex + 1,iEndIndex);
	}
	template<typename Type> void QuickSort(vector<Type>& vtArray,vector<unsigned int>& viSortedIndices,unsigned const int& iStartIndex,unsigned const int& iEndIndex)
	{
		// check the problem validity
		unsigned int iSize = (unsigned int)vtArray.size();
		if(iSize <= 1)
		{
			return;
		}
		
		if(iStartIndex >= iEndIndex)
		{
			return;
		}
		
		if(iEndIndex >= iSize)
		{
			return;
		}
		
		// pick a random pivot
		unsigned int iPivotIndex = (unsigned int)(Randomizer::RandomInteger(iStartIndex,iEndIndex));
 		// partition the array into 2 arrays based on the chosen pivot
 		unsigned int iNewPivotIndex = Partition(vtArray,viSortedIndices,iStartIndex,iEndIndex,iPivotIndex);
 		// sort each array recursively
 		QuickSort(vtArray,viSortedIndices,iStartIndex,iNewPivotIndex - 1);
 		QuickSort(vtArray,viSortedIndices,iNewPivotIndex + 1,iEndIndex);
	}
	template<typename Type> void QuickSort(vector<Type>& vtArray)
	{
		unsigned int iEndIndex = (unsigned int)vtArray.size() - 1;
		QuickSort(vtArray,0,iEndIndex);
	}
	template<typename Type> void QuickSort(vector<Type>& vtArray,vector<unsigned int>& viSortedIndices)
	{
		unsigned int iEndIndex = (unsigned int)vtArray.size() - 1;
		viSortedIndices.clear();
		viSortedIndices.resize(iEndIndex + 1);
		unsigned int i = 0;
		for(i = 0 ; i <= iEndIndex ; i++)
		{
			viSortedIndices[i] = i;
		}
		QuickSort(vtArray,viSortedIndices,0,iEndIndex);
	}
}

#endif

