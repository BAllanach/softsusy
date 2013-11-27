C-----------------------------------------------------------------------
C     THIS SUBROUTINE MUST BE CALLED AT THE BEGINNING OF A RUN - ONCE
C     PER PROGRAM. INITIALISES HERWIG.
C-----------------------------------------------------------------------
      SUBROUTINE INITHWG
C---COMMON BLOCKS ARE INCLUDED AS FILE HERWIG65.INC
      INCLUDE 'HERWIG65.INC'
      INTEGER N
      EXTERNAL HWUDAT
C---MAX NUMBER OF EVENTS THIS RUN
      MAXEV=0
C---BEAM PARTICLES
      PART1='P'
      PART2='P'
C---BEAM MOMENTA
      PBEAM1=7000.
      PBEAM2=PBEAM1
C---PROCESS
      IPROC=3000
C---INITIALISE OTHER COMMON BLOCKS
      CALL HWIGIN
C---USER CAN RESET PARAMETERS AT
C   THIS POINT, OTHERWISE DEFAULT
C   VALUES IN HWIGIN WILL BE USED.
      syspin=.false.
      PRVTX=.FALSE.
      MAXER=MAXEV/100
      PTMIN=50.

      RETURN
      END


C-----------------------------------------------------------------------
C     THIS SUBROUTINE MUST BE CALLED AT THE BEGINNING OF A RUN - ONCE
C     PER PROGRAM. INITIALISES HERWIG. FOR E+ E-
C-----------------------------------------------------------------------
      SUBROUTINE INITHWG2
C---COMMON BLOCKS ARE INCLUDED AS FILE HERWIG65.INC
      INCLUDE 'HERWIG65.INC'
      INTEGER N
      EXTERNAL HWUDAT
C---MAX NUMBER OF EVENTS THIS RUN
      MAXEV=0
C---BEAM PARTICLES
      PART1='E-'
      PART2='E+'
C---BEAM MOMENTA
      PBEAM1=1000.
      PBEAM2=PBEAM1
C---PROCESS
      IPROC=700
C C---INITIALISE OTHER COMMON BLOCKS
      CALL HWIGIN
C     C---USER CAN RESET PARAMETERS AT
C     C   THIS POINT, OTHERWISE DEFAULT
C     C   VALUES IN HWIGIN WILL BE USED.
      syspin=.false.
      PRVTX=.FALSE.
      MAXER=MAXEV/100
      PTMIN=PBEAM1/70.D0

      RETURN
      END

C-----------------------------------------------------------------------
C     READS IN HERWIG INPUT FILE+INITIALISES ITSELF
C-----------------------------------------------------------------------
      SUBROUTINE READINP
      INCLUDE 'HERWIG65.INC'

C---READ IN SUSY INPUT FILE, IN THIS CASE LHC SUGRA POINT 2
      OPEN(UNIT=LRSUSY,FORM='FORMATTED',STATUS='OLD',ERR=999,
     $     FILE='hwg.in')
      CALL HWISSP
      CLOSE(UNIT=LRSUSY)
C---COMPUTE PARAMETER-DEPENDENT CONSTANTS
      CALL HWUINC
C---CALL HWUSTA TO MAKE ANY PARTICLE STABLE
      CALL HWUSTA('PI0     ')

      RETURN
 999  WRITE (6,*)
      WRITE (6,*) 'SUSY input file did not open correctly.'
      WRITE (6,*) 'Please check that it is in the right place.'
      WRITE (6,*) 'Examples can be obtained from the ISAWIG web page.'
      WRITE (6,*)

      END

C-----------------------------------------------------------------------
C     PURELY FOR DEBUGGING PURPOSES
C-----------------------------------------------------------------------
      SUBROUTINE READINP2
      INCLUDE 'HERWIG65.INC'

C---READ IN SUSY INPUT FILE, IN THIS CASE LHC SUGRA POINT 2
      OPEN(UNIT=LRSUSY,FORM='FORMATTED',STATUS='OLD',ERR=999,
     $     FILE='hwg.in2')
      CALL HWISSP
      CLOSE(UNIT=LRSUSY)
C---COMPUTE PARAMETER-DEPENDENT CONSTANTS
      CALL HWUINC
C---CALL HWUSTA TO MAKE ANY PARTICLE STABLE
      CALL HWUSTA('PI0     ')

      RETURN

 999  WRITE (6,*)
      WRITE (6,*) 'SUSY input file did not open correctly.'
      WRITE (6,*) 'Please check that it is in the right place.'
      WRITE (6,*) 'Examples can be obtained from the ISAWIG web page.'
      WRITE (6,*)

      END


C-----------------------------------------------------------------------
C     CALCULATES CROSS-SECTION FOR A GIVEN PROCESS CODE IPP. READS IN A
C     HERWIG INPUT FILE hwg.in FIRST. INITIALISES ITSELF.
C-----------------------------------------------------------------------
      SUBROUTINE CALCSIG(IPP,sigmaGlu,sigmaUpL,sigmaDnL,sigmaChL
     $     ,sigmaStL)
C---COMMON BLOCKS ARE INCLUDED AS FILE HERWIG65.INC
      INCLUDE 'HERWIG65.INC'
      include "particleid.inc"
      INTEGER IPP,N,gluinoNum, upLNum, dnLNum, stLNum, chLNum
      DOUBLE PRECISION SIGMA, sigmaGlu,sigmaUpL,sigmaDnL,sigmaChL
     $     ,sigmaStL

      gluinoNum = 0
      upLNum = 0
      dnLNum = 0
      chLNum = 0
      stLNum = 0
c     First work out the number of squarks and gluinos produced
      IPROC=3010 ! sparton production
C---INITIALISE RUN
      CALL BENIN
C---INITIALISE ELEMENTARY PROCESS
      CALL HWEINI

      MAXEV=100
      DO 100 N=1,MAXEV
C---  INITIALISE EVENT
         CALL HWUINE
C---  GENERATE HARD SUBPROCESS
         CALL HWEPRO
C---  FINISH EVENT
         CALL HWUFNE
      if (idhep(7).eq.gluino) gluinoNum = gluinoNum + 1
      if (idhep(8).eq.gluino) gluinoNum = gluinoNum + 1
      if (idhep(7).eq.abs(down_squark_L)) dnLNum = dnLNum + 1
      if (idhep(7).eq.abs(up_squark_L)) upLNum = upLNum + 1
      if (idhep(7).eq.abs(charm_squark_L)) chLNum = chLNum + 1
      if (idhep(7).eq.abs(strange_squark_L)) stLNum = stLNum + 1
      if (idhep(8).eq.abs(down_squark_L)) dnLNum = dnLNum + 1
      if (idhep(8).eq.abs(up_squark_L)) upLNum = upLNum + 1
      if (idhep(8).eq.abs(charm_squark_L)) chLNum = chLNum + 1
      if (idhep(8).eq.abs(strange_squark_L)) stLNum = stLNum + 1
 100  CONTInUE
C---  TERMINATE ELEMENTARY PROCESS
      CALL HWEFIN
C---  USER'S TERMINAL CALCULATIONS
      CALL HWAEND

      SIGMA=1000.*AVWGT ! total sparton cross-section
      
      sigmaDnL = dble(dnLNum) / dble(MAXEV) * sigma
      sigmaUpL = dble(upLNum) / dble(MAXEV) * sigma
      sigmaChL = dble(chLNum) / dble(MAXEV) * sigma
      sigmaStL = dble(stLNum) / dble(MAXEV) * sigma
      sigmaGlu = dble(gluinoNum) / dble(MAXEV) * sigma

      print *,"Number of gluinos=",sigmaGlu
c      print *,"Number of 1-2 gen squarks=",sigmaUpL,sigmaDnL,sigmaChL
c     $     ,sigmaStL

      END


C-----------------------------------------------------------------------
      SUBROUTINE UPINIT
C-----------------------------------------------------------------------
C     DUMMY UPINIT ROUTINE DELETE AND REPLACE IF USING LES HOUCHES
C     INTERFACE
C-----------------------------------------------------------------------
      WRITE (6,10)
 10   FORMAT(/10X,'UPINIT CALLED BUT NOT LINKED')
c      STOP
      END
C-----------------------------------------------------------------------
CDECK  ID>,  UPEVNT.
*CMZ :-        -16/07/02  10.30.48  by  Peter Richardson
*-- Author :    Peter Richardson
C-----------------------------------------------------------------------
      SUBROUTINE UPEVNT
C-----------------------------------------------------------------------
C     DUMMY UPEVNT ROUTINE DELETE AND REPLACE IF USING LES HOUCHES
C     INTERFACE
C-----------------------------------------------------------------------
      WRITE (6,10)
 10   FORMAT(/10X,'UPEVNT CALLED BUT NOT LINKED')
      STOP
      END


      SUBROUTINE HWAEND
C     USER'S ROUTINE FOR TERMINAL CALCULATIONS, HISTOGRAM OUTPUT, ETC
      END


C     INITIALISES IN BETWEEN DIFFERENT PROCESS CALCULATION
      SUBROUTINE BENIN
C-----------------------------------------------------------------------
C     SETS INPUT PARAMETERS
C----------------------------------------------------------------------
      INCLUDE 'HERWIG65.INC'
      IPRINT=0
C---UNIT FOR WRITING EVENT DATA IN HWANAL (IF ZERO THEN NOT WRITTEN)
      LWEVT=0
C---ASSUMED MAXIMUM WEIGHT (ZERO TO RECOMPUTE)
      WGTMAX=0.
C---CURRENT NO OF EVENTS
      NEVHEP=0
C---CURRENT NO OF ENTRIES IN /HEPEVT/
      NHEP=0

C---CHECK WHETHER SUSY DATA INPUTTED
      SUSYIN = .FALSE.

C---SPIN CORRELATIONS IN TOP/TAU/SUSY DECAYS
      SYSPIN = .FALSE.
  999 END


