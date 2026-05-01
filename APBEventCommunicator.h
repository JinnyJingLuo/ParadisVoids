#ifndef APBEVENTCOMMUNICATOR_H_
#define APBEVENTCOMMUNICATOR_H_

#include "vector"
#include "list"

using namespace std;

enum APBEventType {
  NULLAPBEventType = 0,
  APBCreation = 1,
  APBCellShearing = 2
};

enum APBCommunicationMessages {
  NULLAPBCommunicationMessage = 0,
  APBDataLengthCommunication = 1,
  APBDataCommunication = 2
};

class APB;
class APBPoint;
class ParadisPrecipitateServer;

class APBEvent {
public:
  virtual ~APBEvent();
  APBEvent &operator=(const APBEvent &oEvent);
  virtual void Reset();
  virtual void SetAPB(APB *poAPB);
  void SetType(const APBEventType &eType);
  APBEventType GetType() const;
  virtual void Pack(vector<list<double>> *pvldData,
                    const unsigned int &iDomainID) = 0;
  virtual void Unpack(list<double> *pldData);
  virtual void Process(ParadisPrecipitateServer *poServer) = 0;
  static APBEvent *CreateEvent(const APBEventType &eType);

private:
protected:
  virtual void Initialize();
  unsigned int m_iSlipPlaneID;
  double m_dPlaneSpacing;
  APBEventType m_eType;
};

class APBCreationEvent : public APBEvent {
public:
  APBCreationEvent();
  ~APBCreationEvent();
  APBCreationEvent(const APBCreationEvent &oEvent);
  APBCreationEvent &operator=(const APBCreationEvent &oEvent);
  void Reset();
  void SetAPB(APB *poAPB);
  void Pack(vector<list<double>> *pvldData, const unsigned int &iDomainID);
  void Unpack(list<double> *pldData);
  void Process(ParadisPrecipitateServer *poServer);

private:
protected:
  void Initialize();
};

class APBCellShearingEvent : public APBEvent {
public:
  APBCellShearingEvent();
  ~APBCellShearingEvent();
  APBCellShearingEvent(const APBCellShearingEvent &oEvent);
  APBCellShearingEvent &operator=(const APBCellShearingEvent &oEvent);
  void Reset();
  void SetCellPoint(APBPoint *poPoint);
  void Pack(vector<list<double>> *pvldData, const unsigned int &iDomainID);
  void Unpack(list<double> *pldData);
  void Process(ParadisPrecipitateServer *poServer);

private:
protected:
  void Initialize();
  unsigned int m_iCellX;
  unsigned int m_iCellY;
  unsigned int m_iOwningDomain;
  int m_iXShear;
  int m_iYShear;
  int m_iZShear;
};

class APBEventCommunicator {
public:
  static APBEventCommunicator *GetCommunicator();
  ~APBEventCommunicator();
  void Reset();
  void SetDomainID(const unsigned int &iDomainID);
  void SetDomainsCount(const unsigned int &iDomainsCount);
  unsigned int GetDomainID() const;
  void AddEvent(APBEvent *poEvent);
  void Communicate(ParadisPrecipitateServer *poServer);

private:
protected:
  static APBEventCommunicator *m_poCommunicator;
  APBEventCommunicator();
  void Initialize();
  void Pack(vector<list<double>> *pvldData);
  void Unpack(vector<list<double>> *pvldData);
  void ProcessEvents(ParadisPrecipitateServer *poServer);
  void ClearEvents();
  list<APBEvent *> m_lpoEvents;
  unsigned int m_iDomainID;
  unsigned int m_iDomainsCount;
};

#endif
