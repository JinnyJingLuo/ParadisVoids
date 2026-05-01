# CEMEL-Paradis

In house rewrite of paradis by Ahmed Hussien, and used, maintained and modified
by the El-Awady group (CEMEL) at JHU.

## Notable differences:

Ahmed's code, except that the collisions(Hinge-Joint collisions) between a
sessile dislocation and a mobile one is enabled

## Change History
This code started from Ahmed's version.

2020/04/30: Merged Junjie & Yejun's changes

Flags to activate the functions:

Junjie:
  1. Peierls Stress:
    1. EnablePeierls
    2. PeierlsScrew
    3. PeierlsEdge
    4. PeierlsOthers
  2. Twin Plane Cross Slip:
    1. EnableTwinPlaneCrossSlip
    2. A, B, C, X_1, Y_1, Z_1

Yejun:
  1. ThermoStress
    1. SType: External stress field only, no Temperature gradient
    2. Ttype: Temperature-affected mobility law and cross slip events
