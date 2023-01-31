*      
*     SOME AUXILIARY FUNCTIONS CALLED BY THE DSZ, BDSZ & DDS ROUTINES
*
*     Last updates:
*
*                 13/07/2004: 1) unused variables removed
*                 24/11/2003: 1) phi expansion around lambda=0 removed
*
*                 06/08/2003: 1) phi expansion around x/z=0,y/z=0
*
*                 13/05/2003: 1) new routine for the complex dilogarithm
*                             2) phi expansion around lambda=0
*                             3) more Passarino-Veltman functions
*
***********************************************************************
*

      real*8 function myAA(m,q)      
      real*8 m,q

      if(m.ne.0d0) then
         myAA = m*(1d0-Log(m/q))
      else
         myAA = 0d0
      endif

      return
      end

*
***********************************************************************
*


      real*8 function myB0(q,m1,m2,mu2) 

c     from Degrassi and Sirlin, Phys. Rev. D46 (1992) 3104.
      
      real*8 q,m1,m2,Omega,mu2

      if(q.eq.0d0) then

         if(m1.eq.0d0.and.m2.ne.0d0) then
            myB0 = 1d0-Log(m2/mu2)
         elseif(m1.ne.0d0.and.m2.eq.0d0) then
            myB0 = 1d0-Log(m1/mu2)
         elseif(m1.eq.m2) then
            myB0 = -Log(m1/mu2)
         else
            myB0 = 1d0 - Log(m2/mu2) + m1/(m1-m2)*Log(m2/m1)
         endif
         
      else

         if(m1.eq.0d0.and.m2.ne.0d0) then
            
            if(m2.ne.q) then
               myB0 = -(Log(m2/mu2)-2-(m2/q-1d0)*Log(dabs(1d0-q/m2))) 
            else 
               myB0 = -(Log(m2/mu2) - 2)
            endif
            
         elseif(m2.eq.0d0.and.m1.ne.0d0) then
            
            if(m1.ne.q) then
               myB0 = -(Log(m1/mu2)-2-(m1/q-1d0)*Log(dabs(1d0-q/m1))) 
            else
               myB0 = -(Log(m1/mu2) - 2)
            endif
            
         elseif(m2.eq.0d0.and.m1.eq.0d0) then
            
            myB0 = -(Log(q/mu2) - 2) ! cut the imaginary part (I Pi)
            
         else
            
            myB0 = -( dlog(q/mu2)-2.d0 + 
     1           1.d0/2.d0*( 1.d0 + (m1/q-m2/q))*dlog(m1/q) +
     2           1.d0/2.d0*( 1.d0 - (m1/q-m2/q))*dlog(m2/q) +
     3           2.d0*Omega(m1/q,m2/q))
            
         endif
         
      endif

      return
      end
      
c     function Omega(a,b) contained in myB0
      real*8 function Omega(a,b)
      real*8 a,b,cbig
      Cbig = (a+b)/2.d0 - (a-b)**2.d0/4.d0 -1.d0/4.d0
      if(Cbig.gt.0.d0) then
         Omega = dsqrt(Cbig)*
     1        (datan((1.d0 + a - b)/(2.d0*dsqrt(Cbig))) +
     2        datan((1.d0 - a + b)/(2.d0*dsqrt(Cbig))) )
      elseif(Cbig.lt.0d0) then
         Cbig = - Cbig
         Omega = 1.d0/2.d0*dsqrt(Cbig)*
     1        dlog((a/2.d0 +b/2.d0 -1.d0/2.d0 -dsqrt(Cbig))/
     2        (a/2.d0 + b/2.d0 -1.d0/2.d0 + dsqrt(Cbig)))
      else
         Omega = 0         
      endif

      return
      end

*
**********************************************************************
*
      
      real*8 function myB1(p,m1,m2,q)

      implicit none

      real*8 p,m1,m2,q
      real*8 myAA,myB0
      
      if(p.eq.0d0) then
         myB1 = (1d0-Log(m2/q)+m1**2/(m1-m2)**2*Log(m2/m1)
     $        +(m1+m2)/(m1-m2)/2d0)/2d0
      else
         myB1 = (myAA(m2,q)-myAA(m1,q)+(p+m1-m2)*myB0(p,m1,m2,q))/2d0/p
      endif
      
      return
      end

*
**********************************************************************
*
      
      real*8 function myG(p,m1,m2,q)

      implicit none

      real*8 p,m1,m2,q
      real*8 myAA,myB0

      myG = (p-m1-m2)*myB0(p,m1,m2,q) - myAA(m1,q) - myAA(m2,q)

      return
      end
      
*
**********************************************************************
*

      function phi(x,y,z)

c     from Davydychev and Tausk, Nucl. Phys. B397 (1993) 23

      implicit none
      real*8 x,y,z,phi,pphi,myphi
      
      if(x.le.z.and.y.le.z) then
         pphi = myphi(x,y,z)
      elseif(z.le.x.and.y.le.x) then
         pphi = z/x*myphi(z,y,x)
      elseif(z.le.y.and.x.le.y) then
         pphi = z/y*myphi(z,x,y)
      endif

      phi = pphi
      
      end
      
      function myphi(x,y,z)
      
      implicit none

      real*8 x,y,z,myphi
      real*8 u,v
      real*8 Pi,Li2
      complex*16 clam,cxp,cxm,CLI2,ccphi

      Pi = 3.14159265358979d0

c     auxiliary variables

      u = x/z
      v = y/z
      
      if(u.le.1d-8) then
         
         if(v.ne.1d0) then
            myphi = (log(u)*log(v)+2d0*Li2(1d0-v))/(1d0-v)
         else
            myphi = 2d0-log(u)
         endif

      elseif(v.le.1d-8) then

         if(u.ne.1d0) then
            myphi = (log(v)*log(u)+2d0*Li2(1d0-u))/(1d0-u)
         else
            myphi = 2d0-log(v)
         endif

      else
         
         if((1d0-u-v)**2.ge.4d0*u*v) then         
            clam = dcmplx(dsqrt((1d0-u-v)**2 - 4d0*u*v),0d0)
         else
            clam = dcmplx(0d0,dsqrt(4d0*u*v - (1d0-u-v)**2))
         endif
         
         cxp = (1d0+(u-v)-clam)/2d0
         cxm = (1d0-(u-v)-clam)/2d0
         
c     phi function from eq. (A4)
            
         ccphi = (2d0*log(cxp)*log(cxm) - log(u)*log(v) - 
     &        2d0*(CLI2(cxp) + CLI2(cxm)) + Pi**2/3d0)/clam
         myphi = dreal(ccphi)
                     
      endif
      
      return
      end


*
***********************************************************************
*

      function delt(x,y,z)

      implicit none

      real*8 x,y,z,delt

      delt = x**2 + y**2 + z**2 - 2d0*(x*y + x*z + y*z)

      return
      end

*
***********************************************************************
*

      function Li2(x)

      implicit none

      complex*16 CLI2,z
      real*8 x,Li2

      z = DCMPLX(x,0d0)
      Li2 = dreal(CLI2(z))

      return
      end

*
***********************************************************************
*

      COMPLEX*16 FUNCTION CLI2(Z)

c     just call the CSPEN routine
      
      COMPLEX*16 Z,CSPEN

      CLI2 = CSPEN(Z)

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SPENCE-FUNKTION KOMPLEX, FREI NACH HOLLIK                     C
C---------------------------------------------------------------------C
C       20.07.83    LAST CHANGED 10.05.89        ANSGAR DENNER        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        FUNCTION CSPEN(ZZ)

        COMPLEX*16 CSPEN,W,SUM,ZZ,Z,U
        REAL*8 RZ,AZ,A1
        REAL*8 B(9)
C     BEACHTE:                 B(N)=B2N
C     B(1)=1./6.
C     B(2)=-1./30.
C     B(3)=1./42.
C     B(4)=-1./30.
C     B(5)=5./66.
C     B(6)=-691./2730.
C     B(7)=7./6.
C     B(8)=-3617./510.
C     B(9)=43867./798.
C     PI=3.1415926535897932384
C     PI*PI/6.=1.6449..., PI*PI/3=3.28986...

      B(1)=0.1666666666666666666666666667d0
      B(2)=-0.0333333333333333333333333333d0
      B(3)=0.0238095238095238095238095238d0
      B(4)=-0.0333333333333333333333333333d0
      B(5)=0.0757575757575757575757575758d0
      B(6)=-0.2531135531135531135531135531d0
      B(7)=1.1666666666666666666666666667d0
      B(8)=-7.09215686274509804d0
      B(9)=54.97117794486215539d0

c      write(*,*) 'z:',z
      Z =ZZ*DCMPLX(1D0)
      RZ=DREAL(Z)
      AZ=CDABS(Z)
      A1=CDABS(1D0-Z)
c      write(*,*)'z, rz, az, a1:',z,rz,az,a1
C     IF((SNGL(RZ) .EQ. 0.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN
C ---> CHANGED  10.5.89
      IF(AZ .LT. 1D-20) THEN
        CSPEN=-CDLOG(1D0-Z)
c        write(*,*) 'cspen:', cspen
        RETURN
      END IF
      IF((SNGL(RZ) .EQ. 1.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN
        CSPEN=1.64493406684822643D0
c        write(*,*) 'cspen:', cspen
        RETURN
      END IF
      IF(RZ.GT.5D-1) GOTO 20
      IF(AZ.GT.1D0) GOTO 10
      W=-CDLOG(1D0-Z)
      SUM=W-0.25D0*W*W
      U=W
      IF(CDABS(U).LT.1D-10) GOTO 2
c      write(*,*) 'u:',u
c      write(*,*) 'sum:',sum
      DO 1 K=1,9
      U=U*W*W/DFLOAT(2*K*(2*K+1))
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 2
      SUM=SUM+U*B(K)
 1    CONTINUE
 2    CSPEN=SUM
c        write(*,*) 'cspen:', cspen
      RETURN
10    W=-CDLOG(1D0-1D0/Z)
      SUM=W-0.25D0*W*W
      U=W
      IF(CDABS(U).LT.1D-10) GOTO 12

      DO 11 K=1,9
      U=U*W*W/DFLOAT(2*K*(2*K+1))
      IF(CDABS(B(K)*U/SUM).LT.1D-20) GOTO 12
      SUM=SUM+U*B(K)
11    CONTINUE
12    CSPEN=-SUM-1.64493406684822643D0-.5D0*CDLOG(-Z)**2
c        write(*,*) 'cspen:', cspen
      RETURN
20    IF(A1.GT.1D0) GOTO 30
      W=-CDLOG(Z)
      SUM=W-0.25D0*W*W
      U=W
      IF(CDABS(U).LT.1D-10) GOTO 22
      DO 21 K=1,9
      U=U*W*W/DFLOAT(2*K*(2*K+1))
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 22
      SUM=SUM+U*B(K)
21    CONTINUE
22    CSPEN=-SUM+1.64493406684822643D0-CDLOG(Z)*CDLOG(1D0-Z)
c        write(*,*) 'cspen:', cspen
      RETURN
30    W=CDLOG(1D0-1D0/Z)
      SUM=W-0.25D0*W*W
      U=W
      IF(CDABS(U).LT.1D-10) GOTO 32
      DO 31 K=1,9
      U=U*W*W/DFLOAT(2*K*(2*K+1))
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 32
      SUM=SUM+U*B(K)
31    CONTINUE
32    CSPEN=SUM+3.28986813369645287D0
     *               +.5D0*CDLOG(Z-1D0)**2-CDLOG(Z)*CDLOG(1D0-Z)
c        write(*,*) 'cspen:', cspen
      END