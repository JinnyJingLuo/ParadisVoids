#include "APBEventCommunicator.h"
#include "ParadisPrecipitateServer.h"
#include "mpi.h"
#include "math.h"

APBEvent::~APBEvent() { Reset(); }
APBEvent &APBEvent::operator=(const APBEvent &oEvent) {
  m_eType = oEvent.m_eType;
  m_iSlipPlaneID = oEvent.m_iSlipPlaneID;
  m_dPlaneSpacing = oEvent.m_dPlaneSpacing;
  return *this;
}
void APBEvent::Reset() { Initialize(); }
void APBEvent::SetAPB(APB *poAPB) {
  m_iSlipPlaneID = poAPB->GetSlipPlaneID();
  m_dPlaneSpacing = poAPB->GetPlaneSpacing();
}
void APBEvent::SetType(const APBEventType &eType) { m_eType = eType; }
APBEventType APBEvent::GetType() const { return m_eType; }
void APBEvent::Unpack(list<double> *pldData) {
  m_iSlipPlaneID = (unsigned int)floor(pldData->front() + 0.5);
  pldData->pop_front();
  m_dPlaneSpacing = pldData->front();
  pldData->pop_front();
}
void APBEvent::Initialize() {
  m_eType = NULLAPBEventType;
  m_iSlipPlaneID = 0;
  m_dPlaneSpacing = 0.0;
}
APBEvent *APBEvent::CreateEvent(const APBEventType &eType) {
  if (eType == APBCreation) {
    return new APBCreationEvent;
  } else if (eType == APBCellShearing) {
    return new APBCellShearingEvent;
  }
  return NULL;
}

APBCreationEvent::APBCreationEvent() { Initialize(); }
APBCreationEvent::~APBCreationEvent() { Reset(); }
APBCreationEvent::APBCreationEvent(const APBCreationEvent &oEvent) {
  *this = oEvent;
}
APBCreationEvent &APBCreationEvent::operator=(const APBCreationEvent &oEvent) {
  APBEvent::operator=(oEvent);
  return *this;
}
void APBCreationEvent::Reset() { APBEvent::Reset(); }
void APBCreationEvent::SetAPB(APB *poAPB) { APBEvent::SetAPB(poAPB); }
void APBCreationEvent::Pack(vector<list<double>> *pvldData,
                            const unsigned int &iDomainID) {
  // APB creation events should be communicated to all the domains
  unsigned int iDomainsCount = (unsigned int)pvldData->size();
  unsigned int i = 0;

  for (i = 0; i < iDomainsCount; i++) {
    // don't self send the creation data
    if (i == iDomainID)
      continue;
    // push the event type followed by the event data
    pvldData->at(i).push_back((double)APBCreation);
    pvldData->at(i).push_back((double)m_iSlipPlaneID);
    pvldData->at(i).push_back(m_dPlaneSpacing);
  }
}
void APBCreationEvent::Unpack(list<double> *pldData) {
  APBEvent::Unpack(pldData);
}
void APBCreationEvent::Process(ParadisPrecipitateServer *poServer) {
  APB *poAPB = poServer->GetAPBByPlane(m_iSlipPlaneID, m_dPlaneSpacing);
  if (poAPB != NULL) {
    // the APB is already there, don't add it
    return;
  }
  poServer->AddAPBByPlane(m_iSlipPlaneID, m_dPlaneSpacing);
}
void APBCreationEvent::Initialize() {
  APBEvent::Initialize();
  m_eType = APBCreation;
}

APBCellShearingEvent::APBCellShearingEvent() { Initialize(); }
APBCellShearingEvent::~APBCellShearingEvent() { Reset(); }
APBCellShearingEvent::APBCellShearingEvent(const APBCellShearingEvent &oEvent) {
  *this = oEvent;
}
APBCellShearingEvent &APBCellShearingEvent::
operator=(const APBCellShearingEvent &oEvent) {
  APBEvent::operator=(oEvent);
  m_iCellX = oEvent.m_iCellX;
  m_iCellY = oEvent.m_iCellY;
  m_iOwningDomain = oEvent.m_iOwningDomain;
  m_iXShear = oEvent.m_iXShear;
  m_iYShear = oEvent.m_iYShear;
  m_iZShear = oEvent.m_iZShear;
  return *this;
}
void APBCellShearingEvent::Reset() {
  APBEvent::Reset();
  Initialize();
}
void APBCellShearingEvent::SetCellPoint(APBPoint *poPoint) {
  m_iCellX = poPoint->m_iX;
  m_iCellY = poPoint->m_iY;
  m_iOwningDomain = poPoint->GetOwningDomain();
  m_iXShear = poPoint->GetXShear();
  m_iYShear = poPoint->GetYShear();
  m_iZShear = poPoint->GetZShear();
}
void APBCellShearingEvent::Pack(vector<list<double>> *pvldData,
                                const unsigned int &iDomainID) {
  // APB cell flipping events should only be communicated to the owning domain
  // 	if(m_iOwningDomain == iDomainID)				return;
  // 	// push the event type followed by the event data
  // 	pvldData->at(m_iOwningDomain).push_back((double)APBCellFlipping);
  // 	pvldData->at(m_iOwningDomain).push_back((double)m_iSlipSystemID);
  // 	pvldData->at(m_iOwningDomain).push_back(m_dPlaneSpacing);
  // 	pvldData->at(m_iOwningDomain).push_back((double)m_iCellX);
  // 	pvldData->at(m_iOwningDomain).push_back((double)m_iCellY);

  unsigned int iDomainsCount = (unsigned int)pvldData->size();
  unsigned int i = 0;

  for (i = 0; i < iDomainsCount; i++) {
    // don't self send the creation data
    if (i == iDomainID)
      continue;
    // push the event type followed by the event data
    pvldData->at(i).push_back((double)APBCellShearing);
    pvldData->at(i).push_back((double)m_iSlipPlaneID);
    pvldData->at(i).push_back(m_dPlaneSpacing);
    pvldData->at(i).push_back((double)m_iCellX);
    pvldData->at(i).push_back((double)m_iCellY);
    pvldData->at(i).push_back((double)m_iXShear);
    pvldData->at(i).push_back((double)m_iYShear);
    pvldData->at(i).push_back((double)m_iZShear);
  }
}
void APBCellShearingEvent::Unpack(list<double> *pldData) {
  APBEvent::Unpack(pldData);
  m_iCellX = (unsigned int)floor(pldData->front() + 0.5);
  pldData->pop_front();
  m_iCellY = (unsigned int)floor(pldData->front() + 0.5);
  pldData->pop_front();
  m_iXShear = (unsigned int)floor(pldData->front() + 0.5);
  pldData->pop_front();
  m_iYShear = (unsigned int)floor(pldData->front() + 0.5);
  pldData->pop_front();
  m_iZShear = (unsigned int)floor(pldData->front() + 0.5);
  pldData->pop_front();
}
void APBCellShearingEvent::Process(ParadisPrecipitateServer *poServer) {
  APB *poAPB = poServer->GetAPBByPlane(m_iSlipPlaneID, m_dPlaneSpacing);
  if (poAPB == NULL) {
    // something is wrong, this processor doesn't have the APB that owns the
    // flipped cell
    printf("error: processing cell shearing event failed - APB not found\n");
    fflush(NULL);
    return;
  }
  poAPB->ShearCell(m_iCellX, m_iCellY, m_iXShear, m_iYShear, m_iZShear);
}
void APBCellShearingEvent::Initialize() {
  APBEvent::Initialize();
  m_eType = APBCellShearing;
  m_iCellX = 0;
  m_iCellY = 0;
  m_iOwningDomain = 0;
  m_iXShear = 0;
  m_iYShear = 0;
  m_iZShear = 0;
}

APBEventCommunicator *APBEventCommunicator::m_poCommunicator = NULL;
APBEventCommunicator *APBEventCommunicator::GetCommunicator() {
  if (m_poCommunicator == NULL) {
    m_poCommunicator = new APBEventCommunicator;
  }
  return m_poCommunicator;
}
APBEventCommunicator::APBEventCommunicator() { Initialize(); }
APBEventCommunicator::~APBEventCommunicator() { Reset(); }
void APBEventCommunicator::Reset() { ClearEvents(); }
void APBEventCommunicator::SetDomainID(const unsigned int &iDomainID) {
  m_iDomainID = iDomainID;
}
void APBEventCommunicator::SetDomainsCount(const unsigned int &iDomainsCount) {
  m_iDomainsCount = iDomainsCount;
}
unsigned int APBEventCommunicator::GetDomainID() const { return m_iDomainID; }
void APBEventCommunicator::AddEvent(APBEvent *poEvent) {
  if (poEvent == NULL)
    return;
  m_lpoEvents.push_back(poEvent);
}
void APBEventCommunicator::Communicate(ParadisPrecipitateServer *poServer) {
  vector<list<double>> vldData;
  vldData.resize(m_iDomainsCount);
  // pack the data
  Pack(&vldData);

  // loop over all the packed data vectors and determine the data size to be
  // sent to each process, self sending/receiving is being automatically avoided
  unsigned int *piSentDataSizes = new unsigned int[m_iDomainsCount];
  unsigned int *piReceivedDataSizes = new unsigned int[m_iDomainsCount];
  MPI_Request *poSendRequests = new MPI_Request[m_iDomainsCount];
  MPI_Request *poReceiveRequests = new MPI_Request[m_iDomainsCount];

  // communicate the data sizes
  // issue receives
  unsigned int i = 0;
  for (i = 0; i < m_iDomainsCount; i++) {
    piReceivedDataSizes[i] = 0;
    MPI_Irecv(&piReceivedDataSizes[i], 1, MPI_INT, i,
              APBDataLengthCommunication, MPI_COMM_WORLD,
              &poReceiveRequests[i]);
  }
  // send
  for (i = 0; i < m_iDomainsCount; i++) {
    piSentDataSizes[i] = (unsigned int)vldData[i].size();
    MPI_Isend(&piSentDataSizes[i], 1, MPI_INT, i, APBDataLengthCommunication,
              MPI_COMM_WORLD, &poSendRequests[i]);
  }
  // wait for the communications to finish, ignore the statuses
  MPI_Waitall(m_iDomainsCount, poSendRequests, MPI_STATUSES_IGNORE);
  MPI_Waitall(m_iDomainsCount, poReceiveRequests, MPI_STATUSES_IGNORE);

  // now every processor knows the amount of data it will receive, communicate
  // the data
  double **ppdSentData = new double *[m_iDomainsCount];
  double **ppdReceivedData = new double *[m_iDomainsCount];
  // issue receives
  for (i = 0; i < m_iDomainsCount; i++) {
    if (piReceivedDataSizes[i] == 0)
      continue;
    ppdReceivedData[i] = new double[piReceivedDataSizes[i]];
    MPI_Irecv(ppdReceivedData[i], piReceivedDataSizes[i], MPI_DOUBLE, i,
              APBDataCommunication, MPI_COMM_WORLD, &poReceiveRequests[i]);
  }
  // send
  unsigned int j = 0;
  list<double>::iterator liData;
  for (i = 0; i < m_iDomainsCount; i++) {
    if (piSentDataSizes[i] == 0)
      continue;
    ppdSentData[i] = new double[piSentDataSizes[i]];
    j = 0;
    for (liData = vldData[i].begin(); liData != vldData[i].end(); liData++) {
      ppdSentData[i][j] = (*liData);
      j++;
    }
    MPI_Isend(ppdSentData[i], piSentDataSizes[i], MPI_DOUBLE, i,
              APBDataCommunication, MPI_COMM_WORLD, &poSendRequests[i]);
    // clear the data list
    vldData[i].clear();
  }
  // wait for the communications to finish, ignore the statuses,
  // copy to data list once received and free send and receive arrays
  for (i = 0; i < m_iDomainsCount; i++) {
    if (piSentDataSizes[i] == 0)
      continue;
    MPI_Wait(&poSendRequests[i], MPI_STATUSES_IGNORE);
    delete[] ppdSentData[i];
  }
  for (i = 0; i < m_iDomainsCount; i++) {
    if (piReceivedDataSizes[i] == 0)
      continue;
    MPI_Wait(&poReceiveRequests[i], MPI_STATUSES_IGNORE);
    for (j = 0; j < piReceivedDataSizes[i]; j++) {
      vldData[i].push_back(ppdReceivedData[i][j]);
    }
    delete[] ppdReceivedData[i];
  }
  delete[] ppdSentData;
  delete[] ppdReceivedData;
  delete[] piSentDataSizes;
  delete[] piReceivedDataSizes;
  delete[] poSendRequests;
  delete[] poReceiveRequests;

  Unpack(&vldData);
  vldData.clear();
  ProcessEvents(poServer);
}
void APBEventCommunicator::Initialize() {
  m_lpoEvents.clear();
  m_iDomainID = 0;
  m_iDomainsCount = 1;
}
void APBEventCommunicator::Pack(vector<list<double>> *pvldData) {
  // pack all the events into the data list vector
  // the function loops over all of the events and calls the packing functions
  // for the events, once done, the list of events gets cleared to make space
  // for the events unpacking later
  list<APBEvent *>::iterator liEvents;
  for (liEvents = m_lpoEvents.begin(); liEvents != m_lpoEvents.end();
       liEvents++) {
    (*liEvents)->Pack(pvldData, m_iDomainID);
  }
  ClearEvents();
}
void APBEventCommunicator::Unpack(vector<list<double>> *pvldData) {
  // we start with an empty list of events and unpack all the received data into
  // it
  unsigned int i = 0;
  APBEvent *poEvent = NULL;
  APBEventType eEventType = NULLAPBEventType;
  list<double> *pldData = NULL;

  for (i = 0; i < m_iDomainsCount; i++) {
    pldData = &pvldData->at(i);
    while (!pldData->empty()) {
      eEventType = (APBEventType)((int)floor(pldData->front() + 0.5));
      pldData->pop_front();
      poEvent = APBEvent::CreateEvent(eEventType);
      poEvent->Unpack(pldData);
      m_lpoEvents.push_back(poEvent);
    }
  }
}
void APBEventCommunicator::ProcessEvents(ParadisPrecipitateServer *poServer) {
  list<APBEvent *>::iterator liEvents;
  for (liEvents = m_lpoEvents.begin(); liEvents != m_lpoEvents.end();
       liEvents++) {
    (*liEvents)->Process(poServer);
  }
  ClearEvents();
}
void APBEventCommunicator::ClearEvents() {
  list<APBEvent *>::iterator liEvents;
  for (liEvents = m_lpoEvents.begin(); liEvents != m_lpoEvents.end();
       liEvents++) {
    if ((*liEvents) != NULL)
      delete (*liEvents);
  }
  m_lpoEvents.clear();
}
