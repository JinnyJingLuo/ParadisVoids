/****************************************************************************
 *
 *      Restart.h  Contains various prototypes for functions related
 *                 to reading/writing restart files.
 *
 ***************************************************************************/
#ifndef _Restart_h
#define _Restart_h

/*
 *      Prototypes for functions involved in reading the restart files
 */
void AssignNodesToDomains(Home_t *home, InData_t *inData, int nodeCount,
                          int ***nodeLists, int **listCounts);
void FreeNodeLists(Home_t *home, int ***nodeLists, int **listCounts);
void ReadControlFile(Home_t *home, const char *ctrlFileName);
void ReadNodeDataFile(Home_t *home, InData_t *inData, const char *dataFile);
void ReadStressDataFile(Home_t *home, InData_t *inData, const char *stressFile);

/*
 *      Prototypes for functions involved in writing the restart files
 */
void SetLatestRestart(char *fileName);
void WriteRestart(Home_t *home, char *baseFileName, int ioGroup,
                  int firstInGroup, int writePrologue, int writeEpilogue);

#endif
