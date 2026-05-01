// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#include "Vector.h"
#include "Tools.h"
#include "math.h"

using namespace EZ;

namespace SupportSystem
{
	string GetRealString(unsigned int iMaxSize,FILE* fpFile)
	{
		bool bIsReal = false;
		char cReadLine[1000];
		string sString = "";

		while(!bIsReal)
		{
			cReadLine[0] = 0;
			sString = "";
			if(feof(fpFile))
			{
				return sString;
			}
			fgets(cReadLine,iMaxSize,fpFile);
			sString = string(cReadLine);
			if(sString.size() < 1)
			{
				continue;
			}
			// if it is a new line or a blank space, skip it
			if(sString.size() == 1)
			{
				if(sString[0] == 0 || sString[0] == 10 || sString[0] == 13 || sString[0] == 32)
				{
					continue;
				}
			}
			sString = sString.substr(sString.find_first_not_of(' '),sString.size() - 1);
			if(sString[0] != 0 && sString[0] != '*' && sString[0] != 10 && sString[0] != 13 && sString[0] != '#' && sString[0] != '!')
			{
				bIsReal = true;
			}
		}
		return sString;
	}
	void PrintOnScreen(const string& sStatement)
	{
		printf("%s\n",sStatement.c_str());
		fflush(NULL);
	}
}

