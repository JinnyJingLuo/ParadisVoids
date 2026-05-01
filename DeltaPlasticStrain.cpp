/***************************************************************************
 *
 *      Module:     DeltaPlasticStrain.
 *
 *      Description:  This contains a simple generic dispatch function that
 *                    will invoke the version of DeltaPlasticStrain*()
 *                    appropriate to the type of material being simulated
 *
 ***************************************************************************/
#include "Home.h"
#include "Vector.h"
#include "DSMPI.h"
#include "CartesianOrthogonalCoordinateSystem.h"

using namespace EZ;

void FixBurgersVector(Vector &oBurgersVector) {
  double dTolerance = 1.0E-3;
  bool bFlip = false;
  if (fabs(oBurgersVector.GetX()) < dTolerance) {
    if (fabs(oBurgersVector.GetY()) < dTolerance) {
      if (oBurgersVector.GetZ() <= dTolerance) {
        bFlip = true;
      }
    } else if (oBurgersVector.GetY() < 0.0) {
      bFlip = true;
    }
  } else if (oBurgersVector.GetX() < 0.0) {
    bFlip = true;
  }
  if (bFlip) {
    oBurgersVector.Reverse();
  }
}
double GetTriSweptArea(const Point &oP1, const Point &oP2, const Point &oP3,
                       Vector oNormal) {
  CartesianOrthogonalCoordinateSystem oSystem;
  oSystem.SetOrigin(oP1);
  oSystem.SetXZ(Vector(oP1, oP2), oNormal);
  Point oPoint1 = oSystem.GetInLocalCoordinates(oP1);
  Point oPoint2 = oSystem.GetInLocalCoordinates(oP2);
  Point oPoint3 = oSystem.GetInLocalCoordinates(oP3);
  // see if the area is positive or negative
  double dArea = 0.5 * oPoint3.GetX() * oPoint3.GetY();
  return dArea;
}
double GetSweptArea(Param_t *poParam, const Vector &oNormal, Node_t *poNode,
                    Node_t *poNeighbour) {
  Point oNode1Old(poNode->oldx, poNode->oldy, poNode->oldz);
  Point oNode1New =
      GetNearestImage(poParam, poNode->oldx, poNode->oldy, poNode->oldz,
                      poNode->x, poNode->y, poNode->z);
  Point oNode2Old =
      GetNearestImage(poParam, poNode->oldx, poNode->oldy, poNode->oldz,
                      poNeighbour->oldx, poNeighbour->oldy, poNeighbour->oldz);
  Point oNode2New =
      GetNearestImage(poParam, poNode->oldx, poNode->oldy, poNode->oldz,
                      poNeighbour->x, poNeighbour->y, poNeighbour->z);

  CartesianOrthogonalCoordinateSystem oSystem;
  oSystem.SetOrigin(oNode1Old);
  oSystem.SetXZ(Vector(oNode1Old, oNode2Old), Vector(oNormal));
  Point oP1 = oSystem.GetInLocalCoordinates(oNode1Old);
  Point oP2 = oSystem.GetInLocalCoordinates(oNode2Old);
  Point oP3 = oSystem.GetInLocalCoordinates(oNode1New);
  Point oP4 = oSystem.GetInLocalCoordinates(oNode2New);

  double dOriginalLength = oP2.GetX();
  Vector oNewSegment = Vector(oP3, oP4);
  double dNewLength = oNewSegment.Length();

  double dL = dOriginalLength;
  double dAlpha = dNewLength - dOriginalLength;
  double dGamma = 0.5 * (oP3.GetX() + oP4.GetX() - dOriginalLength);
  double dDelta = 0.5 * (oP3.GetY() + oP4.GetY());
  double dBeta = acos(oNewSegment.GetX() / dNewLength);
  double dTheta = 0.5 * PI;
  double dAngleSum = dTheta + dBeta;

  double dTolerance = 1.0E-6;
  double dArea = 0.0;
  if (fabs(dBeta) > dTolerance) {
    double dI1 = dL * dGamma / dBeta * (sin(dAngleSum) - sin(dTheta));
    double dI2 = dAlpha * dGamma / dBeta / dBeta *
                 (dBeta * sin(dAngleSum) + cos(dAngleSum) - cos(dTheta));
    double dI3 = -dL * dDelta / dBeta * (cos(dAngleSum) - cos(dTheta));
    double dI4 = dAlpha * dDelta / dBeta / dBeta *
                 (sin(dAngleSum) - dBeta * cos(dAngleSum) - sin(dTheta));
    dArea = dI1 + dI2 + dI3 + dI4;
  } else {
    dArea = (dGamma * cos(dTheta) + dDelta * sin(dTheta)) * (dL + dAlpha / 2);
  }
  return dArea;
}
void DeltaPlasticStrain(Home_t *home) {
  unsigned int i = 0;
  unsigned int j = 0;
  double dDeltaPlasticStrain[3][3];
  double dDeltaPlasticSpin[3][3];
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      dDeltaPlasticStrain[i][j] = 0.0;
      dDeltaPlasticSpin[i][j] = 0.0;
    }
  }
  home->param->disloDensity = 0.0;

  // loop over all the segments
  Node_t *poNode = NULL;
  Node_t *poNeighbour = NULL;
  Vector oBurgersVector;
  Vector oTemp;
  Vector oNormal;
  double dArea = 0.0;
  double dDiadicProduct[3][3];
  unsigned int k = 0;
  unsigned int l = 0;
  double dTolerance = 1.0E-6;
  for (i = 0; i < home->newNodeKeyPtr; i++) {
    poNode = home->nodeKeys[i];
    if (poNode == NULL) {
      continue;
    }
    for (j = 0; j < poNode->numNbrs; j++) {
      poNeighbour = GetNeighborNode(home, poNode, j);
      if (poNeighbour == NULL) {
        printf("WARNING: Neighbor not found at %s line %d\n", __FILE__,
               __LINE__);
        continue;
      }
      if ((poNode->constraint >=SURFACE_NODE) &&
          (poNeighbour->constraint >=SURFACE_NODE)) {
        continue;
      }
      // avoid adding the segment plastic strain twice
      if (fabs(poNode->burgX[j]) < dTolerance) {
        if (fabs(poNode->burgY[j]) < dTolerance) {
          if (poNode->burgZ[j] <= dTolerance) {
            continue;
          }
        } else if (poNode->burgY[j] < 0.0) {
          continue;
        }
      } else if (poNode->burgX[j] < 0.0) {
        continue;
      }

      // make sure that the segment is planar
      oNormal.Set(poNode->nx[j], poNode->ny[j], poNode->nz[j]);
      oTemp.Set(poNeighbour->x - poNode->x, poNeighbour->y - poNode->y,
                poNeighbour->z - poNode->z);
      oNormal.Normalize();
      oTemp.Normalize();
      if (fabs(oNormal * oTemp) > dTolerance) {
        //				printf("@ %d : segment (%d,%d) -->
        //(%d,%d)
        // is off plane by %e :
        //%e,%e,%e\n",home->myDomain,poNode->myTag.domainID,poNode->myTag.index,poNeighbour->myTag.domainID,poNeighbour->myTag.index,fabs(oNormal*oTemp),poNode->nx[j],poNode->ny[j],poNode->nz[j]);
        //				printf("@ %d : current positions :
        //%e,%e,%e
        //-->
        //%e,%e,%e\n",home->myDomain,poNode->x,poNode->y,poNode->z,poNeighbour->x,poNeighbour->y,poNeighbour->z);
        //				printf("@ %d : velocities : %e,%e,%e -->
        //%e,%e,%e\n",home->myDomain,poNode->vX,poNode->vY,poNode->vZ,poNeighbour->vX,poNeighbour->vY,poNeighbour->vZ);
        //				printf("@ %d : old positions : %e,%e,%e
        //-->
        //%e,%e,%e\n",home->myDomain,poNode->oldx,poNode->oldy,poNode->oldz,poNeighbour->oldx,poNeighbour->oldy,poNeighbour->oldz);
        continue;
      }
      oBurgersVector.Set(poNode->burgX[j], poNode->burgY[j], poNode->burgZ[j]);
      dArea = GetSweptArea(home->param, oNormal, poNode, poNeighbour);

      for (k = 0; k < 3; k++) {
        for (l = 0; l < 3; l++) {
          dDiadicProduct[k][l] = dArea * oNormal.GetComponentZeroBased(k) *
                                 oBurgersVector.GetComponentZeroBased(l);
        }
      }

      for (k = 0; k < 3; k++) {
        for (l = 0; l < 3; l++) {
          dDeltaPlasticStrain[k][l] +=
              (dDiadicProduct[k][l] + dDiadicProduct[l][k]) * 0.5 /
              home->param->simVol;
          dDeltaPlasticSpin[k][l] +=
              (dDiadicProduct[k][l] - dDiadicProduct[l][k]) * 0.5 /
              home->param->simVol;
        }
      }

      // calculate the density contribution of the segment
      oTemp = GetNearestImage(home->param, poNode, poNeighbour);
      oTemp = oTemp - Vector(poNode->x, poNode->y, poNode->z);
      home->param->disloDensity += oTemp.Length() * home->param->burgVolFactor;
    }
  }

  // sum the density, plastic strain and spin from all processors
  home->param->delpStrain[0] = dDeltaPlasticStrain[0][0];
  home->param->delpStrain[1] = dDeltaPlasticStrain[1][1];
  home->param->delpStrain[2] = dDeltaPlasticStrain[2][2];
  home->param->delpStrain[3] = dDeltaPlasticStrain[1][2];
  home->param->delpStrain[4] = dDeltaPlasticStrain[0][2];
  home->param->delpStrain[5] = dDeltaPlasticStrain[0][1];

  home->param->delpSpin[0] = dDeltaPlasticSpin[0][0];
  home->param->delpSpin[1] = dDeltaPlasticSpin[1][1];
  home->param->delpSpin[2] = dDeltaPlasticSpin[2][2];
  home->param->delpSpin[3] = dDeltaPlasticSpin[1][2];
  home->param->delpSpin[4] = dDeltaPlasticSpin[0][2];
  home->param->delpSpin[5] = dDeltaPlasticSpin[0][1];

  double dTempPlasticStrain[6];
  double dTempPlasticSpin[6];
  double dTempDensity = 0.0;
  DSMPI::AllReduce(home->param->delpStrain, dTempPlasticStrain, 6, MPI_DOUBLE,
                   MPI_SUM);
  DSMPI::AllReduce(home->param->delpSpin, dTempPlasticSpin, 6, MPI_DOUBLE,
                   MPI_SUM);
  DSMPI::AllReduce(&home->param->disloDensity, &dTempDensity, 1, MPI_DOUBLE,
                   MPI_SUM);
  for (i = 0; i < 6; i++) {
    home->param->delpStrain[i] = dTempPlasticStrain[i];
    home->param->delpSpin[i] = dTempPlasticSpin[i];
  }
  home->param->disloDensity = dTempDensity;

  // update the total plastic strain and spin
  for (int i = 0; i < 6; i++) {
    home->param->totpStn[i] +=
        home->param->delpStrain[i] + home->param->mergedelpStrain[i];
    home->param->totpSpn[i] +=
        home->param->delpSpin[i] + home->param->mergedelpSpin[i];
  }
}
void StoreMergePlasticSweep(Home_t *home, const Point &oP1, const Point &oP2,
                            const Point &oP3, const Vector &oBurgersVector,
                            const Vector &oNormal) {
  Vector oWorkingBurgersVector = oBurgersVector;
  // FixBurgersVector(oWorkingBurgersVector);
  double dArea = GetTriSweptArea(oP1, oP2, oP3, oNormal);

  unsigned int i = 0;
  unsigned int j = 0;
  double dDeltaPlasticStrain[3][3];
  double dDeltaPlasticSpin[3][3];
  double dFactor = 0.5 * dArea / home->param->simVol;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      dDeltaPlasticStrain[i][j] =
          dFactor * oNormal.GetComponentZeroBased(i) *
              oWorkingBurgersVector.GetComponentZeroBased(j) +
          oNormal.GetComponentZeroBased(j) *
              oWorkingBurgersVector.GetComponentZeroBased(i);
      dDeltaPlasticSpin[i][j] =
          dFactor * oNormal.GetComponentZeroBased(i) *
              oWorkingBurgersVector.GetComponentZeroBased(j) -
          oNormal.GetComponentZeroBased(j) *
              oWorkingBurgersVector.GetComponentZeroBased(i);
    }
  }

  // store the plastic strain and spin
  home->param->mergedelpStrain[0] += dDeltaPlasticStrain[0][0];
  home->param->mergedelpStrain[1] += dDeltaPlasticStrain[1][1];
  home->param->mergedelpStrain[2] += dDeltaPlasticStrain[2][2];
  home->param->mergedelpStrain[3] += dDeltaPlasticStrain[1][2];
  home->param->mergedelpStrain[4] += dDeltaPlasticStrain[0][2];
  home->param->mergedelpStrain[5] += dDeltaPlasticStrain[0][1];

  home->param->mergedelpSpin[0] += dDeltaPlasticSpin[0][0];
  home->param->mergedelpSpin[1] += dDeltaPlasticSpin[1][1];
  home->param->mergedelpSpin[2] += dDeltaPlasticSpin[2][2];
  home->param->mergedelpSpin[3] += dDeltaPlasticSpin[1][2];
  home->param->mergedelpSpin[4] += dDeltaPlasticSpin[0][2];
  home->param->mergedelpSpin[5] += dDeltaPlasticSpin[0][1];
}
