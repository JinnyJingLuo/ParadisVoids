#ifndef DSMPI_H_
#define DSMPI_H_

#include "mpi.h"

class DSMPI {
public:
  static void Initialize(int argc, char **argv);
  static void Finalize();
  static int GetProcessID();
  static int GetProcessesCount();
  static void Barrier();
  static void Broadcast(void *buffer, int count, MPI_Datatype datatype,
                        int root);
  static void Reduce(void *sendbuf, void *recvbuf, int count,
                     MPI_Datatype datatype, MPI_Op op, int root);
  static void AllReduce(void *sendbuf, void *recvbuf, int count,
                        MPI_Datatype datatype, MPI_Op op);

private:
protected:
  static int GetDataSize(MPI_Datatype datatype);
  static int m_iProcessesCount;
  static int m_iProcessID;
};

#endif
