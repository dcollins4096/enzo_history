*****************************************************************************
*                                                                           *
* Copyright 2004 Greg Bryan                                                 *
* Copyright 2004 Laboratory for Computational Astrophysics                  *
* Copyright 2004 Board of Trustees of the University of Illinois            *
* Copyright 2004 Regents of the University of California                    *
*                                                                           *
* This software is released under the terms of the "Enzo Public License"    *
* in the accompanying LICENSE file.                                         *
*                                                                           *
*****************************************************************************

C
      PARAMETER(NUME=12,MAXLN=220)
      PARAMETER(NIT=200,TEMMIN=3.0,DELT=0.03)
      PARAMETER(NID=300,DENMIN=-12.0,DELD=0.05,NID2=240)
      PARAMETER(NIB=400,FREQDEL=0.02,FREQMIN=1.0)
C
      COMMON/ATM/NJ(12),ABUNJ(12)
     .          ,WJ(12,MAXLN),E3J(12,MAXLN),FJ(12,MAXLN)
     .          ,EJ(12,30),EAJ(12,30)
     .          ,S2J(12,30),LLJ(12,30),S3J(12,30)
     .          ,S4J(12,30),S5J(12,30)
     .          ,AN(NID),ABIN(NIB),IPHOT
C
C
