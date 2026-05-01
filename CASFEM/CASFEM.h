#ifndef CASFEM_H_
#define CASFEM_H_

#include "FEMSolver.h"
 
using namespace FEMSystem;
 
class CASFEM
{
public:
 	static CASFEM* CreateInstance();
 	~CASFEM();
 	void Reset();
 	bool Initialize(const string& sFileName);
 	void Run();
 	void WriteInLog(const string& sLogString);
 	FEMSolver* GetFEMSolver() const;
 	MainDataStructure* GetMainDataStructure() const;
 	
 	
private:
 
protected:
	static CASFEM* CASFEMInstance;
	CASFEM();
 	FILE* m_fpLog;
 	MainDataStructure* m_poData;
 	FEMSolver* m_poFEMSolver;
};
 

#endif


