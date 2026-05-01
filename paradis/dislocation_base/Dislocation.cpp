// Paradis Processor Project
// Ahmed M. Hussein (mailto : am.hussin@gmail.com)
// June 2012

#include "Dislocation.h"
#include "MathServices.h"
#include "stdio.h"

namespace DislocationSystem {
Vector Dislocation::GetPointDisplacement(
    GraphEdge<DislocationNode, DislocationSegment> *poDislocationSegment,
    Point *poPoint, bool bIsDebug) {
  //  		Vector oDisplacement(0.0,0.0,0.0);
  //  		unsigned int i = 0;
  //  		unsigned int iGaussPointsCount = 56;
  //  		vector<double> vdGaussPointsPositions;
  //  		vector<double> vdGaussPointsWeights;
  //  		MathServices::GenerateGaussPoints(vdGaussPointsPositions,vdGaussPointsWeights,iGaussPointsCount);
  //  		double dEta = 0.0;
  //  		double dWeight = 0.0;
  //  		double dFactor = 0.0;
  //  		Point oGaussPoint;
  //  		Vector oTangent;
  //  		double dR = 0.0;
  //  		double dTolerance = 1.0E-10;
  //  		double dMinR = 2.0*dTolerance;
  //  		Vector oRVector;
  //  		double dTemp1 = 0.0;
  //  		double dTemp2 = 0.0;
  //  		double dTemp3 = 0.0;
  //  		double dTemp4 = 0.0;
  //  		double dTemp5 = 0.0;
  //
  //   		Point oFirstPoint = poDislocationSegment->GetFirst()->GetData();
  //   		Point oLastPoint = poDislocationSegment->GetLast()->GetData();
  //   		oTangent.SetByPoints(oFirstPoint,oLastPoint);
  //  		double dTangentMagnitude = oTangent.Length();
  //  		oTangent.Normalize();
  //
  //  		Vector oBurgersVector =
  //  poDislocationSegment->GetDataPointer()->GetBurgersVector();
  //  Vector
  //  oBurgersDirection = oBurgersVector.GetDirection(); 		double
  //  dBurgersMagnitude = oBurgersVector.Length();
  //  		//double dNu = poSlipSystem->GetPoissonRatio();
  //
  //  		// TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP
  //  *****
  //  TEMP ***** TEMP *****
  //  		// TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP
  //  *****
  //  TEMP ***** TEMP *****
  //  		// TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP
  //  *****
  //  TEMP ***** TEMP *****
  //  		// TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP
  //  *****
  //  TEMP ***** TEMP *****
  //  		// TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP
  //  *****
  //  TEMP ***** TEMP *****
  //  		// TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP
  //  *****
  //  TEMP ***** TEMP ***** 		double dNu = 0.3;
  //  		// TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP
  //  *****
  //  TEMP ***** TEMP *****
  //  		// TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP
  //  *****
  //  TEMP ***** TEMP *****
  //  		// TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP
  //  *****
  //  TEMP ***** TEMP *****
  //  		// TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP
  //  *****
  //  TEMP ***** TEMP *****
  //  		// TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP
  //  *****
  //  TEMP ***** TEMP *****
  //  		// TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP ***** TEMP
  //  *****
  //  TEMP ***** TEMP *****
  //
  //  		Vector oArbitraryVector(1.0,1.0,1.0);
  //  		oArbitraryVector.Normalize();
  //
  // 		for(i = 0; i < iGaussPointsCount ; i++)
  // 		{
  // 			dEta = 0.5*(vdGaussPointsPositions[i] + 1.0);
  // 			dWeight = vdGaussPointsWeights[i];
  // 			oGaussPoint = oFirstPoint*(1.0 - dEta) +
  // oLastPoint*dEta; 			oRVector = Vector(*poPoint,oGaussPoint);
  // dFactor
  // = 1.0*dWeight; 			dR = oRVector.Length(); 			if(dR <
  // dTolerance)
  // 			{
  // 				continue;
  // 			}
  // 			if(dR < dMinR)			// minimum allowable distance,
  // probably
  // equal to the minimum dislocation radius
  // 			{
  // 				dR = dMinR;
  // 			}
  // 			oRVector.Normalize();
  // 			dTemp1 = (oRVector^oTangent)*oBurgersDirection;
  // 			if(dTemp1 < 1.0E-12)
  // 			{
  // 				continue;
  // 			}
  // 			dTemp2 = (oArbitraryVector^oRVector)*oTangent;
  // 			dTemp3 =
  // dFactor*dTangentMagnitude*dBurgersMagnitude*dTemp1/(8.0*PI*dR)/(1.0 - dNu);
  // 			dTemp4 = 2.0*(1.0 - dNu)*dTemp2/dTemp1/(1.0 +
  // oArbitraryVector*oRVector); 			dTemp5 = (1.0 - 2.0*dNu)/dTemp1;
  // oDisplacement
  // = oDisplacement + (oBurgersDirection*dTemp4 +
  // (oTangent^oBurgersDirection)*dTemp5 + oRVector)*dTemp3;
  // 		}
  //   		return oDisplacement;
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  Vector oDisplacement(0.0, 0.0, 0.0);
  unsigned int i = 0;
  unsigned int iGaussPointsCount = 20;
  vector<double> vdGaussPointsPositions;
  vector<double> vdGaussPointsWeights;
  MathServices::GenerateGaussPoints(vdGaussPointsPositions,
                                    vdGaussPointsWeights, iGaussPointsCount);
  double dEta = 0.0;
  double dWeight = 0.0;
  Point oGaussPoint;
  Vector oTangent;
  double dR = 0.0;
  double dTolerance = 1.0E-6;
  // double dMinR = 2.0*dTolerance;
  Vector oRVector;
  double dX = 0.0;
  double dY = 0.0;
  double dZ = 0.0;
  Vector oBurgersVector =
      poDislocationSegment->GetDataPointer()->GetBurgersVector();
  double dBX = oBurgersVector.GetX();
  double dBY = oBurgersVector.GetY();
  double dBZ = oBurgersVector.GetZ();
  double dTemp1 = 0.0;
  double dTemp2 = 0.0;
  double dTemp3 = 0.0;
  double dTemp4 = 0.0;
  double dUX = 0.0;
  double dUY = 0.0;
  double dUZ = 0.0;
  double dNu = 0.3;
  double dK = 1 / 2.0 / (1.0 - dNu);
  Point oFirstPoint = poDislocationSegment->GetFirst()->GetData();
  Point oLastPoint = poDislocationSegment->GetLast()->GetData();
  oTangent.SetByPoints(oFirstPoint, oLastPoint);
  // solution according to the expression in Peach paper
  for (i = 0; i < iGaussPointsCount; i++) {
    dEta = 0.5 * (vdGaussPointsPositions[i] + 1.0);
    dWeight = vdGaussPointsWeights[i];
    oGaussPoint = oFirstPoint * (1.0 - dEta) + oLastPoint * dEta;
    oRVector = Vector(oGaussPoint, *poPoint);

    dX = oRVector.GetX();
    dY = oRVector.GetY();
    dZ = oRVector.GetZ();
    dR = oRVector.Length();
    if (dR < dTolerance) {
      continue;
    }
    // displacement x component

    if (bIsDebug) {
      printf("starter @ %d : %25.20f\t\t%25.20f\t\t%25.20f\t\t%25.20f\n", i + 1,
             dX, dY, dZ, dR);
    }

    dTemp1 = dY * dZ / 3.0 / dR *
             (1 / (dX * dX + dZ * dZ) - 1 / (dX * dX + dY * dY)) * dBX;
    dTemp2 = dX * dZ / 3.0 / dR *
             (1 / (dX * dX + dY * dY) - 1 / (dZ * dZ + dY * dY)) * dBX;
    dTemp3 = dY * dX / 3.0 / dR *
             (1 / (dY * dY + dZ * dZ) - 1 / (dX * dX + dZ * dZ)) * dBX;
    // check for NaN
    if (dTemp1 != dTemp1) {
      dTemp1 = 0.0;
    }
    if (dTemp2 != dTemp2) {
      dTemp2 = 0.0;
    }
    if (dTemp3 != dTemp3) {
      dTemp3 = 0.0;
    }
    if (bIsDebug) {
      printf("1st X temps @ %d : %25.20f\t\t%25.20f\t\t%25.20f\n", i + 1,
             dTemp1, dTemp2, dTemp3);
    }

    dTemp1 = dTemp1 + dX * dK / dR / dR / dR * (-dBY * dZ + dBZ * dY);
    dTemp2 = dTemp2 + dX * dK / dR / dR / dR * (-dBZ * dX + dBX * dZ) -
             (1.0 - dK) / dR * dBZ;
    dTemp3 = dTemp3 + dX * dK / dR / dR / dR * (-dBX * dY + dBY * dX) +
             (1.0 - dK) / dR * dBY;

    if (bIsDebug) {
      printf("2nd X temps @ %d : %25.20f\t\t%25.20f\t\t%25.20f\n", i + 1,
             dTemp1, dTemp2, dTemp3);
    }

    dTemp1 = dTemp1 * oTangent.GetX();
    dTemp2 = dTemp2 * oTangent.GetY();
    dTemp3 = dTemp3 * oTangent.GetZ();

    if (bIsDebug) {
      printf("3rd X temps @ %d : %25.20f\t\t%25.20f\t\t%25.20f\n", i + 1,
             dTemp1, dTemp2, dTemp3);
      // printf("%25.20f\t\t%25.20f\t\t%25.20f\n",dTemp1,dTemp2,dTemp3);
    }

    dTemp4 = dTemp1 + dTemp2 + dTemp3;
    dTemp4 = dTemp4 * dWeight;
    dUX = dUX + dTemp4;

    // displacement y component
    dTemp1 = dZ * dX / 3.0 / dR *
             (1 / (dY * dY + dX * dX) - 1 / (dY * dY + dZ * dZ)) * dBY;
    dTemp2 = dY * dX / 3.0 / dR *
             (1 / (dY * dY + dZ * dZ) - 1 / (dX * dX + dZ * dZ)) * dBY;
    dTemp3 = dZ * dY / 3.0 / dR *
             (1 / (dZ * dZ + dX * dX) - 1 / (dY * dY + dX * dX)) * dBY;
    // check for NaN
    if (dTemp1 != dTemp1) {
      dTemp1 = 0.0;
    }
    if (dTemp2 != dTemp2) {
      dTemp2 = 0.0;
    }
    if (dTemp3 != dTemp3) {
      dTemp3 = 0.0;
    }
    dTemp1 = dTemp1 + dY * dK / dR / dR / dR * (-dBZ * dX + dBX * dZ);
    dTemp2 = dTemp2 + dY * dK / dR / dR / dR * (-dBX * dY + dBY * dX) -
             (1.0 - dK) / dR * dBX;
    dTemp3 = dTemp3 + dY * dK / dR / dR / dR * (-dBY * dZ + dBZ * dY) +
             (1.0 - dK) / dR * dBZ;

    dTemp1 = dTemp1 * oTangent.GetY();
    dTemp2 = dTemp2 * oTangent.GetZ();
    dTemp3 = dTemp3 * oTangent.GetX();
    dTemp4 = dTemp1 + dTemp2 + dTemp3;
    dTemp4 = dTemp4 * dWeight;
    dUY = dUY + dTemp4;

    // displacement z component
    dTemp1 = dX * dY / 3.0 / dR *
             (1 / (dZ * dZ + dY * dY) - 1 / (dZ * dZ + dX * dX)) * dBZ;
    dTemp2 = dZ * dY / 3.0 / dR *
             (1 / (dZ * dZ + dX * dX) - 1 / (dY * dY + dX * dX)) * dBZ;
    dTemp3 = dX * dZ / 3.0 / dR *
             (1 / (dX * dX + dY * dY) - 1 / (dZ * dZ + dY * dY)) * dBZ;
    // check for NaN
    if (dTemp1 != dTemp1) {
      dTemp1 = 0.0;
    }
    if (dTemp2 != dTemp2) {
      dTemp2 = 0.0;
    }
    if (dTemp3 != dTemp3) {
      dTemp3 = 0.0;
    }
    dTemp1 = dTemp1 + dZ * dK / dR / dR / dR * (-dBX * dY + dBY * dX);
    dTemp2 = dTemp2 + dZ * dK / dR / dR / dR * (-dBY * dZ + dBZ * dY) -
             (1.0 - dK) / dR * dBY;
    dTemp3 = dTemp3 + dZ * dK / dR / dR / dR * (-dBZ * dX + dBX * dZ) +
             (1.0 - dK) / dR * dBX;

    dTemp1 = dTemp1 * oTangent.GetZ();
    dTemp2 = dTemp2 * oTangent.GetX();
    dTemp3 = dTemp3 * oTangent.GetY();
    dTemp4 = dTemp1 + dTemp2 + dTemp3;
    dTemp4 = dTemp4 * dWeight;
    dUZ = dUZ + dTemp4;
  }
  double dLength = oTangent.Length();
  oDisplacement.Set(dUX / 4.0 / PI / dLength, dUY / 4.0 / PI / dLength,
                    dUZ / 4.0 / PI / dLength);
  return oDisplacement;
}
Matrix Dislocation::GetSegmentPointStress(const Point &oSegmentStart,
                                          const Point &oSegmentEnd,
                                          const Vector &oBurgersVector,
                                          const double &dShearModulus,
                                          const double &dPoissonRatio,
                                          const Point &oFieldPoint) {
  Matrix oStress(3, 3);

  unsigned int iGaussPointsCount = 20;
  vector<double> vdGaussPointsPositions;
  vector<double> vdGaussPointsWeights;
  MathServices::GenerateGaussPoints(vdGaussPointsPositions,
                                    vdGaussPointsWeights, iGaussPointsCount);
  double dEta = 0.0;
  double dWeight = 0.0;
  Point oGaussPoint;
  Vector oRVector;
  Vector oSegmentVector(oSegmentStart, oSegmentEnd);
  double dR = 0.0;
  double dJacobians[4] = {0.0, oSegmentVector.GetX(), oSegmentVector.GetY(),
                          oSegmentVector.GetZ()};
  Matrix oTempStress(3, 3);
  unsigned int p = 0;
  unsigned int i = 0;
  unsigned int j = 0;
  unsigned int m = 0;
  unsigned int n = 0;
  unsigned int k = 0;
  double dTerm1 = 0.0;
  double dTerm2 = 0.0;
  double dTerm3 = 0.0;
  double dTerm4 = 0.0;
  double dTemp = 0.0;
  double dTolerance = 1E-10;
  double dFactor = 1.0 / (1.0 - dPoissonRatio);
  double dStressComponent = 0.0;

  for (p = 0; p < iGaussPointsCount; p++) {
    dEta = 0.5 * (vdGaussPointsPositions[p] + 1.0);
    dWeight = vdGaussPointsWeights[p];
    oGaussPoint = oSegmentStart * (1.0 - dEta) + oSegmentEnd * dEta;
    oRVector = Vector(oGaussPoint, oFieldPoint);
    dR = oRVector.Length();
    if (dR < dTolerance) {
      dR = dTolerance;
    }

    for (i = 1; i <= 3; i++) {
      for (j = 1; j <= 3; j++) {
        dStressComponent = 0.0;
        for (n = 1; n <= 3; n++) {
          dTerm1 = 0.0;
          dTerm2 = 0.0;
          dTerm3 = 0.0;
          dTerm4 = 0.0;
          for (m = 1; m <= 3; m++) {
            dTerm1 = dTerm1 -
                     MathServices::PermutationSymbol(j, m, n) *
                         oRVector.GetComponent(m) / dR / dR / dR *
                         dJacobians[i];
            dTerm2 = dTerm2 -
                     MathServices::PermutationSymbol(i, m, n) *
                         oRVector.GetComponent(m) / dR / dR / dR *
                         dJacobians[j];
            for (k = 1; k <= 3; k++) {
              dTerm3 = dTerm3 +
                       2.0 * MathServices::PermutationSymbol(k, m, n) *
                           oRVector.GetComponent(m) / dR / dR / dR *
                           dJacobians[k] * MathServices::KroneckerDelta(i, j);
              dTemp = 3.0 * oRVector.GetComponent(i) *
                      oRVector.GetComponent(j) * oRVector.GetComponent(k) / dR /
                      dR;
              dTemp =
                  dTemp -
                  MathServices::KroneckerDelta(i, j) * oRVector.GetComponent(m);
              dTemp =
                  dTemp -
                  MathServices::KroneckerDelta(j, m) * oRVector.GetComponent(i);
              dTemp =
                  dTemp -
                  MathServices::KroneckerDelta(i, m) * oRVector.GetComponent(j);
              dTemp = MathServices::PermutationSymbol(k, m, n) * dTerm4 / dR /
                      dR / dR * dJacobians[k];
              dTerm4 = dTerm4 + dTemp;
            }
          }
          dStressComponent = dStressComponent +
                             (dTerm1 + dTerm2 + dFactor * (dTerm4 - dTerm3)) *
                                 oBurgersVector.GetComponent(n);
        }
        oTempStress.Set(i, j, dStressComponent);
      }
    }
    oStress = oStress + oTempStress * dWeight;
  }
  oStress = oStress * (dShearModulus / 4.0 / PI);
  return oStress;
}
Matrix Dislocation::GetSegmentPointStrain(const Point &oSegmentStart,
                                          const Point &oSegmentEnd,
                                          const Vector &oBurgersVector,
                                          const double &dPoissonRatio,
                                          const Point &oFieldPoint) {
  Matrix oStrain(3, 3);

  unsigned int iGaussPointsCount = 20;
  vector<double> vdGaussPointsPositions;
  vector<double> vdGaussPointsWeights;
  MathServices::GenerateGaussPoints(vdGaussPointsPositions,
                                    vdGaussPointsWeights, iGaussPointsCount);
  double dEta = 0.0;
  double dWeight = 0.0;
  Point oGaussPoint;
  Vector oRVector;
  Vector oSegmentVector(oSegmentStart, oSegmentEnd);
  double dR = 0.0;
  double dJacobians[4] = {0.0, oSegmentVector.GetX(), oSegmentVector.GetY(),
                          oSegmentVector.GetZ()};
  Matrix oTempStrain(3, 3);
  unsigned int p = 0;
  unsigned int i = 0;
  unsigned int j = 0;
  unsigned int m = 0;
  unsigned int l = 0;
  unsigned int k = 0;
  double dTerm1 = 0.0;
  double dTerm2 = 0.0;
  double dTerm3 = 0.0;
  double dTemp = 0.0;
  double dTolerance = 1E-10;
  double dFactor = 2.0 / (1.0 - dPoissonRatio);
  double dStrainComponent = 0.0;

  for (p = 0; p < iGaussPointsCount; p++) {
    dEta = 0.5 * (vdGaussPointsPositions[p] + 1.0);
    dWeight = vdGaussPointsWeights[p];
    oGaussPoint = oSegmentStart * (1.0 - dEta) + oSegmentEnd * dEta;
    oRVector = Vector(oGaussPoint, oFieldPoint);
    dR = oRVector.Length();
    if (dR < dTolerance) {
      dR = dTolerance;
    }

    for (i = 1; i <= 3; i++) {
      for (j = 1; j <= 3; j++) {
        dStrainComponent = 0.0;
        for (k = 1; k <= 3; k++) {
          dTerm1 = 0.0;
          dTerm2 = 0.0;
          dTerm3 = 0.0;
          for (l = 1; l <= 3; l++) {
            dTerm1 = dTerm1 -
                     2.0 * MathServices::PermutationSymbol(i, k, l) *
                         oRVector.GetComponent(j) / dR / dR / dR *
                         dJacobians[i];
            dTerm2 = dTerm2 -
                     2.0 * MathServices::PermutationSymbol(j, k, l) *
                         oRVector.GetComponent(i) / dR / dR / dR *
                         dJacobians[j];
            for (m = 1; m <= 3; m++) {
              dTemp = 3.0 * oRVector.GetComponent(i) *
                      oRVector.GetComponent(j) * oRVector.GetComponent(m) / dR /
                      dR;
              dTemp =
                  dTemp -
                  MathServices::KroneckerDelta(i, j) * oRVector.GetComponent(m);
              dTemp =
                  dTemp -
                  MathServices::KroneckerDelta(j, m) * oRVector.GetComponent(i);
              dTemp =
                  dTemp -
                  MathServices::KroneckerDelta(i, m) * oRVector.GetComponent(j);
              dTemp = MathServices::PermutationSymbol(k, m, l) * dTemp / dR /
                      dR / dR * dJacobians[k];
              dTerm3 = dTerm3 + dFactor * dTemp;
            }
          }
          dStrainComponent =
              dStrainComponent +
              (dTerm1 + dTerm2 + dTerm3) * oBurgersVector.GetComponent(l);
        }
        oTempStrain.Set(i, j, dStrainComponent);
      }
    }
    oStrain = oStrain + oTempStrain * dWeight;
  }
  oStrain = oStrain * (1.0 / 8.0 / PI);
  return oStrain;
}
} // namespace DislocationSystem
