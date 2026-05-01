#include "DSMPI.h"

int DSMPI::m_iProcessesCount = 1;
int DSMPI::m_iProcessID = 0;
void DSMPI::Initialize(int argc, char **argv) {
#ifdef PARALLEL
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &m_iProcessID);
  MPI_Comm_size(MPI_COMM_WORLD, &m_iProcessesCount);
#endif
}
void DSMPI::Finalize() {
#ifdef PARALLEL
  MPI_Finalize();
#endif
}
int DSMPI::GetProcessID() { return m_iProcessID; }
int DSMPI::GetProcessesCount() { return m_iProcessesCount; }
void DSMPI::Barrier() {
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}
void DSMPI::Broadcast(void *buffer, int count, MPI_Datatype datatype,
                      int root) {
#ifdef PARALLEL
  MPI_Bcast(buffer, count, datatype, root, MPI_COMM_WORLD);
#endif
}
void DSMPI::Reduce(void *sendbuf, void *recvbuf, int count,
                   MPI_Datatype datatype, MPI_Op op, int root) {
#ifdef PARALLEL
  MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, MPI_COMM_WORLD);
#else
  // copy the copy array into the receive array
  memcpy(recvbuf, sendbuf, count * GetDataSize(datatype));
#endif
}
void DSMPI::AllReduce(void *sendbuf, void *recvbuf, int count,
                      MPI_Datatype datatype, MPI_Op op) {
#ifdef PARALLEL
  MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, MPI_COMM_WORLD);
#else
  // copy the copy array into the receive array
  memcpy(recvbuf, sendbuf, count * GetDataSize(datatype));
#endif
}
int DSMPI::GetDataSize(MPI_Datatype datatype) {
  if (datatype == MPI_CHAR) {
    return sizeof(char);
  }
  if (datatype == MPI_INT) {
    return sizeof(int);
  }
  if (datatype == MPI_DOUBLE) {
    return sizeof(int);
  }
  return 1;
}
