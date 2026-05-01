// Paradis Processor Project
// Ahmed M. Hussein (mailto : am.hussin@gmail.com)
// June 2012

#include "DislocationLoop.h"
#include "cmath"
#include "MathServices.h"

using namespace EZ;

namespace DislocationSystem {
DislocationLoop::DislocationLoop() { Initialize(); }
DislocationLoop::DislocationLoop(const DislocationLoop &oLoop) {
  *this = oLoop;
}
DislocationLoop::~DislocationLoop() { Reset(); }
DislocationLoop &DislocationLoop::operator=(const DislocationLoop &oLoop) {
  Reset();
  m_loNodesList = oLoop.m_loNodesList;
  m_loNodesList.ResetIterator();
  m_oBurgersVector = oLoop.m_oBurgersVector;
  return *this;
}
void DislocationLoop::Reset() {
  m_loNodesList.Reset();
  m_oBurgersVector.Set(0.0, 0.0, 0.0);
}
bool DislocationLoop::GenerateLoop(DislocationChain *poStartingChain,
                                   list<DislocationChain *> *plpoChains) {
  Reset();
  if (poStartingChain == NULL) {
    return false;
  }
  list<DislocationChain *> lpoChains;
  lpoChains.clear();
  lpoChains.push_back((poStartingChain));
  list<DislocationChain *>::iterator liChains;
  if (poStartingChain->IsLoop()) {
    // remove it from chains list
    liChains = plpoChains->begin();
    while (liChains != plpoChains->end()) {
      if (poStartingChain == (*liChains)) {
        liChains = plpoChains->erase(liChains);
      } else {
        liChains++;
      }
    }
    AppendChain(lpoChains.front());
    m_oBurgersVector = lpoChains.front()->GetBurgersVector();
    return true;
  }
  DislocationChain *poCurrentChain = NULL;
  DislocationNetworkNode *poFirstLoopNode = poStartingChain->GetFirst();
  DislocationNetworkNode *poLastLoopNode = NULL;
  unsigned int iEndNeighboursCount = 0;
  unsigned int i = 0;
  list<DislocationChain *> lpoTempChains;
  Vector oBurgersVector;
  Vector oTempVector;

  list<DislocationChain *>::iterator liTempChains;
  DislocationNetworkNode *poTempNode = NULL;
  bool bNextChainFound = false;
  DislocationChain *poTempChain = NULL;
  while (true) {
    // get current chain
    poCurrentChain = lpoChains.back();
    // remove it from chains list
    liChains = plpoChains->begin();
    while (liChains != plpoChains->end()) {
      if (poCurrentChain == (*liChains)) {
        liChains = plpoChains->erase(liChains);
      } else {
        liChains++;
      }
    }
    // see if this chain touches the loop start node
    if (poCurrentChain->GetLast() == poFirstLoopNode) {
      break;
    }
    // get last node and initialize touching chains search
    poLastLoopNode = poCurrentChain->GetLast();
    lpoTempChains.clear();
    iEndNeighboursCount = poLastLoopNode->GetAllNeighboursCount() - 1;
    i = 0;
    oBurgersVector = poCurrentChain->GetBurgersVector();
    for (liChains = plpoChains->begin(); liChains != plpoChains->end();
         liChains++) {
      if ((*liChains)->IsEndNode(poLastLoopNode)) {
        if ((*liChains) == poCurrentChain) {
          continue;
        }
        i = i + 1;
        lpoTempChains.push_back((*liChains));
        if (i == iEndNeighboursCount) {
          break;
        }
      }
    }
    bNextChainFound = false;
    for (liChains = lpoTempChains.begin(); liChains != lpoTempChains.end();
         liChains++) {
      oTempVector = (*liChains)->GetBurgersVector();
      if (oBurgersVector.IsSame(oTempVector)) {
        if ((*liChains)->IsLastNode(poLastLoopNode)) {
          (*liChains)->SoftReverse();
        }
        bNextChainFound = true;
      } else if (oBurgersVector.IsOpposite(oTempVector)) {
        (*liChains)->Reverse();
        if ((*liChains)->IsLastNode(poLastLoopNode)) {
          (*liChains)->SoftReverse();
        }
        bNextChainFound = true;
      }
      if (bNextChainFound) {
        lpoChains.push_back((*liChains));
        break;
      }
    }
    if (!bNextChainFound) {
      // chain with the required burgers vector was not found, now, this is
      // complicated first of all, look for a chian which closes the loop
      poTempChain = NULL;
      for (liChains = lpoTempChains.begin(); liChains != lpoTempChains.end();
           liChains++) {
        if ((*liChains)->IsEndNode(poFirstLoopNode)) {
          poTempChain = (*liChains);
          bNextChainFound = true;
          break;
        }
      }
      // if not found, look for a chain whose end node was NOT visited before in
      // that loop, this prevents the self intersecting loops from being formed
      if (!bNextChainFound) {
        for (liChains = lpoTempChains.begin(); liChains != lpoTempChains.end();
             liChains++) {
          if (bNextChainFound) {
            break;
          }
          // get the node which is NOT the current chain's end node
          poTempNode = (*liChains)->GetOtherEnd(poLastLoopNode);
          // now loop over all of the chains that are already in that loop and
          // get a chain whose ends were not visited before
          for (liTempChains = lpoChains.begin();
               liTempChains != lpoChains.end(); liTempChains++) {
            if (!(*liTempChains)->IsEndNode(poTempNode)) {
              poTempChain = (*liChains);
              bNextChainFound = true;
              break;
            }
          }
        }
      }
      // now, the next chain should be available, in case it is not, there is
      // something really wrong going on
      if (!bNextChainFound) {
        printf("chain not found, let's crash !!!!\n");
        return false;
      }
      if (poTempChain->IsLastNode(poLastLoopNode)) {
        poTempChain->Reverse();
      }
      lpoChains.push_back(poTempChain);
      poTempChain = poTempChain->Dissociate(oBurgersVector);
      plpoChains->push_back(poTempChain);
    }
  }
  // set burgers vector
  m_oBurgersVector = oBurgersVector;
  // now, add the nodes of the chains forming this loop
  for (liChains = lpoChains.begin(); liChains != lpoChains.end(); liChains++) {
    AppendChain((*liChains));
  }
  lpoChains.clear();
  return true;
}
bool DislocationLoop::IsPlanar() {
  DislocationNetworkNode *poCurrentNode = NULL;
  DislocationNetworkNode *poNextNode = NULL;
  DislocationNetworkArm *poConnectingArm = NULL;
  Vector oInitialNormal;
  Vector oTempNormal;
  ResetIterator();
  poCurrentNode = GetCurrentNode();
  poNextNode = GetNextNode();
  if (poCurrentNode->IsConnected(poNextNode, poConnectingArm)) {
    oInitialNormal = poConnectingArm->GetDataPointer()->GetSlipPlaneNormal();
  } else {
    return false;
  }

  ResetIterator();
  do {
    poCurrentNode = GetCurrentNode();
    poNextNode = GetNextNode();
    if (poCurrentNode->IsConnected(poNextNode, poConnectingArm)) {
      oTempNormal = poConnectingArm->GetDataPointer()->GetSlipPlaneNormal();
      if ((!oTempNormal.IsSame(oInitialNormal))) {
        if (!(oTempNormal.IsOpposite(oInitialNormal))) {
          return false;
        }
      }
    } else {
      return false;
    }
    IncrementIterator();
  } while (!IsAtBeginning());
  return true;
}
const list<DislocationNetworkNode *> &DislocationLoop::GetNodes() {
  return m_loNodesList.GetList();
}
double DislocationLoop::GetSolidAngle(const Point &oPoint) {
  // the projection is assumed to be RIGHT bounded by the loop, that is, the
  // area projection lies to the LEFT of the bounding curve. a traveller moving
  // along the bounding curve would see the projection at his LEFT at all times
  list<DislocationNetworkNode *> lpoNodes = GetNodes();
  CircularLinkedList<Point> loProjections;
  list<DislocationNetworkNode *>::iterator liNodes;
  Vector oTempVector;

  for (liNodes = lpoNodes.begin(); liNodes != lpoNodes.end(); liNodes++) {
    oTempVector.SetByPoints(oPoint, *((*liNodes)->GetDataPointer()));
    oTempVector.Normalize();
    loProjections.Append(oTempVector);
  }

  Point *poPreviousNode = NULL;
  Point *poThisNode = NULL;
  Point *poNextNode = NULL;
  Vector oPreviousVector;
  Vector oCurrentVector;
  Vector oNextVector;
  // see if the loop is right oriented
  bool bRightOriented = true;
  loProjections.ResetIterator();
  poPreviousNode = loProjections.GetPreviousItemPointer();
  poThisNode = loProjections.GetCurrentItemPointer();
  poNextNode = loProjections.GetNextItemPointer();
  double dTolerance = 1.0E-18;

  // constrain the first node to be convex and determine the loop orientation
  // based on it, if it turns out to be concave according to the right loop
  // assumption, reverse the loop and set the flag to be left oriented
  oPreviousVector.SetByPoints(*poPreviousNode, *poThisNode);
  oNextVector.SetByPoints(*poThisNode, *poNextNode);
  oTempVector = oPreviousVector ^ oNextVector;
  oTempVector.Normalize();
  oCurrentVector = *poThisNode;
  if (oTempVector * oCurrentVector > 0.0) {
    // concave first node
    bRightOriented = false;
    loProjections.Reverse();
  }

  bool bKeepIterating = false;
  double dArea = 0.0;
  double dArcLength1 = 0.0;
  double dArcLength2 = 0.0;
  double dArcLength3 = 0.0;
  double dSemiPerimeter = 0.0;
  double dAngularExcess = 0.0;

  // now do the actual triangle cutting, start by taking out the concave angles
  loProjections.ResetIterator();
  do {
    bKeepIterating = false;
    poPreviousNode = loProjections.GetPreviousItemPointer();
    poThisNode = loProjections.GetCurrentItemPointer();
    poNextNode = loProjections.GetNextItemPointer();

    oPreviousVector.SetByPoints(*poPreviousNode, *poThisNode);
    oNextVector.SetByPoints(*poThisNode, *poNextNode);
    oTempVector = oPreviousVector ^ oNextVector;
    oTempVector.Normalize();
    oCurrentVector = *poThisNode;

    if (oTempVector * oCurrentVector > 0.0) {
      oPreviousVector = *poPreviousNode;
      oCurrentVector = *poThisNode;
      oNextVector = *poNextNode;
      oPreviousVector = oPreviousVector ^ oCurrentVector;
      oNextVector = oCurrentVector ^ oNextVector;
      oPreviousVector.Normalize();
      oNextVector.Normalize();
      if (!oPreviousVector.IsSame(oNextVector) &&
          !oPreviousVector.IsOpposite(oNextVector)) {
        oPreviousVector = *poPreviousNode;
        oCurrentVector = *poThisNode;
        oNextVector = *poNextNode;
        dArcLength1 = acos(oPreviousVector * oCurrentVector);
        dArcLength2 = acos(oCurrentVector * oNextVector);
        dArcLength3 = acos(oNextVector * oPreviousVector);
        dSemiPerimeter = 0.5 * (dArcLength1 + dArcLength2 + dArcLength3);
        dAngularExcess = tan(0.5 * dSemiPerimeter) *
                         tan(0.5 * (dSemiPerimeter - dArcLength1));
        dAngularExcess = dAngularExcess *
                         tan(0.5 * (dSemiPerimeter - dArcLength2)) *
                         tan(0.5 * (dSemiPerimeter - dArcLength3));
        if (dAngularExcess > dTolerance) {
          dAngularExcess = 4.0 * atan(sqrt(dAngularExcess));
          dArea = dArea - dAngularExcess;
        }
      }
      // now, cut that triangle out by simply dropping the current point from
      // the list
      loProjections.DropItem();
      bKeepIterating = true;
    } else {
      loProjections.IncrementIterator();
    }
  } while (!loProjections.IsAtBeginning() || bKeepIterating);

  // now we have a convex polygon, work on it
  loProjections.ResetIterator();
  while (loProjections.GetSize() >= 3) {
    poPreviousNode = loProjections.GetPreviousItemPointer();
    poThisNode = loProjections.GetCurrentItemPointer();
    poNextNode = loProjections.GetNextItemPointer();

    // this node is guaranteed to be convex, no need to check on anything
    // compute the area of that triangle, take into account that the triangle
    // lies on a spherical surface compute the first arc length, which is simply
    // the angle between the lines connecting the center of the sphere (the
    // origin in this case) to the end pojnts of the arc
    oPreviousVector = *poPreviousNode;
    oCurrentVector = *poThisNode;
    oNextVector = *poNextNode;
    oPreviousVector = oPreviousVector ^ oCurrentVector;
    oNextVector = oCurrentVector ^ oNextVector;
    oPreviousVector.Normalize();
    oNextVector.Normalize();
    if (!oPreviousVector.IsSame(oNextVector) &&
        !oPreviousVector.IsOpposite(oNextVector)) {
      oPreviousVector = *poPreviousNode;
      oCurrentVector = *poThisNode;
      oNextVector = *poNextNode;
      dArcLength1 = acos(oPreviousVector * oCurrentVector);
      dArcLength2 = acos(oCurrentVector * oNextVector);
      dArcLength3 = acos(oNextVector * oPreviousVector);
      dSemiPerimeter = 0.5 * (dArcLength1 + dArcLength2 + dArcLength3);
      dAngularExcess =
          tan(0.5 * dSemiPerimeter) * tan(0.5 * (dSemiPerimeter - dArcLength1));
      dAngularExcess = dAngularExcess *
                       tan(0.5 * (dSemiPerimeter - dArcLength2)) *
                       tan(0.5 * (dSemiPerimeter - dArcLength3));
      if (dAngularExcess > dTolerance) {
        dAngularExcess = 4.0 * atan(sqrt(dAngularExcess));
        dArea = dArea + dAngularExcess;
      }
    }
    // now, cut that triangle out by simply dropping the current point from the
    // list
    loProjections.DropItem();
  }

  // fix left oriented loops
  if (!bRightOriented) {
    dArea = -dArea;
  }
  return dArea;
}
Vector DislocationLoop::GetPointDisplacement(const Point &oPoint,
                                             const double &dPoissonRatio) {
  double dSolidAngle = GetSolidAngle(oPoint);
  Point *poThisNode = NULL;
  Point *poNextNode = NULL;
  // the tolerance is too high to prevent displacement singularities
  double dTolerance = 1.0E-2;
  unsigned int i = 0;
  unsigned int iGaussPointsCount = 150;
  vector<double> vdGaussPointsLocations;
  vector<double> vdGaussPointsWeights;
  Vector oTangent;
  MathServices::GenerateGaussPoints(vdGaussPointsLocations,
                                    vdGaussPointsWeights, iGaussPointsCount);
  Point oFieldPoint;
  double dR = 0.0;
  double dRxx = 0.0;
  double dRxy = 0.0;
  double dRxz = 0.0;
  double dRyx = 0.0;
  double dRyy = 0.0;
  double dRyz = 0.0;
  double dRzx = 0.0;
  double dRzy = 0.0;
  double dRzz = 0.0;
  double dBx = m_oBurgersVector.GetX();
  double dBy = m_oBurgersVector.GetY();
  double dBz = m_oBurgersVector.GetZ();
  double dTx = 0.0;
  double dTy = 0.0;
  double dTz = 0.0;
  double dUx = 0.0;
  double dUy = 0.0;
  double dUz = 0.0;
  double dTemp = 0.0;
  double dJacobian = 0.0;
  double dMaterialParameter = 0.5 / (1.0 - dPoissonRatio);
  Vector oTemp;
  ResetIterator();
  do {
    poThisNode = GetCurrentNode()->GetDataPointer();
    poNextNode = GetNextNode()->GetDataPointer();
    oTangent.SetByPoints(*poThisNode, *poNextNode);
    oTemp = oTangent;
    oTemp.Normalize();
    dTx = oTemp.GetX();
    dTy = oTemp.GetY();
    dTz = oTemp.GetZ();
    dJacobian = 0.5 * oTangent.Length();
    for (i = 0; i < iGaussPointsCount; i++) {
      dTemp = 0.5 * (vdGaussPointsLocations[i] + 1.0);
      oFieldPoint = *poThisNode + oTangent * dTemp;
      dR = oPoint.Distance(oFieldPoint);
      if (dR < dTolerance) {
        continue;
      }
      dRxx = oPoint.GetRXX(oFieldPoint);
      dRxy = oPoint.GetRXY(oFieldPoint);
      dRxz = oPoint.GetRXZ(oFieldPoint);
      dRyx = oPoint.GetRYX(oFieldPoint);
      dRyy = oPoint.GetRYY(oFieldPoint);
      dRyz = oPoint.GetRYZ(oFieldPoint);
      dRzx = oPoint.GetRZX(oFieldPoint);
      dRzy = oPoint.GetRZY(oFieldPoint);
      dRzz = oPoint.GetRZZ(oFieldPoint);

      dTemp = dBz * dRyx * dTx - dBz * dRxx * dTy + dBy * dRxx * dTz -
              dBy * dRzx * dTx + dBx * dRzx * dTy - dBx * dRyx * dTz;
      dTemp = (dBz * dTy - dBy * dTz) / dR + dTemp * dMaterialParameter;
      dUx = dUx + dTemp * dJacobian * vdGaussPointsWeights[i];

      dTemp = dBz * dRyy * dTx - dBz * dRxy * dTy + dBy * dRxy * dTz -
              dBy * dRzy * dTx + dBx * dRzy * dTy - dBx * dRyy * dTz;
      dTemp = (dBx * dTz - dBz * dTx) / dR + dTemp * dMaterialParameter;
      dUy = dUy + dTemp * dJacobian * vdGaussPointsWeights[i];

      dTemp = dBz * dRyz * dTx - dBz * dRxz * dTy + dBy * dRxz * dTz -
              dBy * dRzz * dTx + dBx * dRzz * dTy - dBx * dRyz * dTz;
      dTemp = (dBy * dTx - dBx * dTy) / dR + dTemp * dMaterialParameter;
      dUz = dUz + dTemp * dJacobian * vdGaussPointsWeights[i];
    }
    IncrementIterator();
  } while (!IsAtBeginning());

  Vector oDisplacement(dUx - dBx * dSolidAngle, dUy - dBy * dSolidAngle,
                       dUz - dBz * dSolidAngle);
  oDisplacement = oDisplacement * (1.0 / 4.0 / PI);
  return oDisplacement;
}
void DislocationLoop::Print() {
  printf("Burgers Vector : (%lf,%lf,%lf) - nodes : ", m_oBurgersVector.GetX(),
         m_oBurgersVector.GetY(), m_oBurgersVector.GetZ());
  ResetIterator();
  do {
    printf("%d\t", GetCurrentNode()->GetID());
    IncrementIterator();
  } while (!IsAtBeginning());
  printf("\n");
}
DislocationNetworkNode *DislocationLoop::GetCurrentNode() {
  return m_loNodesList.GetCurrentItem();
}
DislocationNetworkNode *DislocationLoop::GetPreviousNode() {
  return m_loNodesList.GetPreviousItem();
}
DislocationNetworkNode *DislocationLoop::GetNextNode() {
  return m_loNodesList.GetNextItem();
}
void DislocationLoop::IncrementIterator() { m_loNodesList.IncrementIterator(); }
void DislocationLoop::DecrementIterator() { m_loNodesList.DecrementIterator(); }
void DislocationLoop::ResetIterator() { m_loNodesList.ResetIterator(); }
bool DislocationLoop::IsAtBeginning() const {
  return m_loNodesList.IsAtBeginning();
}
void DislocationLoop::CoarsenVirtualNodes() {
  DislocationNetworkNode *poPreviousNode = NULL;
  DislocationNetworkNode *poThisNode = NULL;
  DislocationNetworkNode *poNextNode = NULL;
  Vector oPreviousSurfaceNormal;
  Vector oThisSurfaceNormal;
  Vector oNextSurfaceNormal;
  bool bNodeDropped = false;
  ResetIterator();
  do {
    poPreviousNode = GetPreviousNode();
    poThisNode = GetCurrentNode();
    poNextNode = GetNextNode();
    bNodeDropped = false;
    if (poPreviousNode->GetDataPointer()->GetCategory() == VirtualNode) {
      oPreviousSurfaceNormal =
          poPreviousNode->GetDataPointer()->GetSurfaceNormal();
      if (poThisNode->GetDataPointer()->GetCategory() == VirtualNode) {
        oThisSurfaceNormal = poThisNode->GetDataPointer()->GetSurfaceNormal();
        if (poNextNode->GetDataPointer()->GetCategory() == VirtualNode) {
          oNextSurfaceNormal = poNextNode->GetDataPointer()->GetSurfaceNormal();
          if (oPreviousSurfaceNormal.IsSame(oThisSurfaceNormal)) {
            if (oThisSurfaceNormal.IsSame(oNextSurfaceNormal)) {
              DropNode();
              bNodeDropped = true;
            }
          }
        }
      }
    }
    if (!bNodeDropped) {
      IncrementIterator();
    }
  } while (!IsAtBeginning() || bNodeDropped);
}
void DislocationLoop::Initialize() {
  m_loNodesList.Reset();
  m_oBurgersVector.Set(0.0, 0.0, 0.0);
}
void DislocationLoop::AppendChain(DislocationChain *poChain) {
  unsigned int i = 0;
  unsigned int iSize = poChain->GetSize();
  poChain->ResetIterator();
  for (i = 0; i < iSize - 1; i++) {
    m_loNodesList.Append(poChain->GetCurrentNode());
    poChain->IncrementIterator();
  }
}
void DislocationLoop::DropNode() { m_loNodesList.DropItem(); }
} // namespace DislocationSystem
