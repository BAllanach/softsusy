
      subroutine dszhiggs(t,mg,T1,T2,st,ct,q,mu,tanb,v2,gs,
     $     OS,S11,S22,S12) bind(C, name="dszhiggs_")
      use iso_c_binding

c     Two-loop O(a_t a_s) corrections to the CP-even Higgs mass matrix. 
c     Routine written by P. Slavich (e-mail: slavich@mppmu.mpg.de).
c     Based on G. Degrassi, P. Slavich and F. Zwirner, 
c     Nucl. Phys. B611 (2001) 403 [hep-ph/0105096].
c
c     Last update:  24/02/2004: mg is given as input instead of mg^2;
c                               value of pi corrected (10th digit);
c                               unused variables cleaned up.
c                   22/10/2002: gs is given as input.
c     
c
c     I/O PARAMETERS:
c     t = m_top^2, mg = m_gluino, T1 = m_stop1^2, T2 = m_stop2^2,
c     st = sin(theta_stop), ct = cos(theta_stop), q = Q^2 (ren. scale),
c     mu = Higgs mixing parameter, tanb = tan(beta), v2 = v^2, 
c     gs = strong coupling constant,
c     OS = renormalization scheme for 1-loop (0 = DRbar, 1 = On-Shell),
c     Sij = 2-loop corrections to the CP-even Higgs mass matrix elements.

      implicit none

      integer OS
      real*8 ht,gs,k,mt,pi,v2
      real*8 t,mg,T1,T2,st,ct,q,A,X,mu,tanb,sb,s2t,c2t
      real*8 F1,F2,F3,sF2,sF3
      real*8 DF1,DF2,DF3,DsF2,DsF3
      real*8 F2_s,sF2_A,sF3_A
      real*8 S11,S22,S12,osdr

c$$$      pi = 3.14159265897d0
      pi = 3.1415926535898d0

      mt = dsqrt(t)

      s2t = 2d0*ct*st
      c2t = ct**2 - st**2

      X = (T1-T2)*s2t/2d0/mt    ! eq. (19) of DSZ
      A = X - mu/tanb           ! notice the sign convention for mu

      sb = dsin(datan(tanb))
      ht = dsqrt(2d0/v2)*mt/sb
      
      k = 4d0*gs**2/(16d0*Pi**2)**2 ! gs^2/(16 Pi^2)^2 CF Nc
      
      call strfuncs(t,mg,T1,T2,s2t,c2t,q,F1,F2,F3)
      call strsfuncs(mg,T1,T2,q,A,sF2,sF3)
      call strdfuncs(t,mg,T1,T2,s2t,c2t,q,A,X,
     $     DF1,DF2,DF3,DsF2,DsF3)

      osdr = 1d0*OS

      if(s2t.ne.0d0.and.A.ne.0d0) then
         
         S11 = .5d0 * ht**2 * mu**2 * s2t**2 * (F3 + osdr*DF3) ! eq. (25)
         
         S12 = .5d0 * ht**2 * mu * A  * s2t**2 * (F3 + sF3 +   ! eq. (26)
     $        osdr*(DF3 + DsF3)) + 
     $        ht**2 * mt * mu * s2t * (F2 + osdr*DF2)
         
         S22 = .5d0 * ht**2 * A**2 * s2t**2 * (F3 + 2d0*sF3 +  ! eq. (27)
     $        osdr*(DF3 + 2d0*DsF3)) + 
     $        2d0 * ht**2 * mt * A * s2t * (F2 + sF2 + 
     $        osdr*(DF2 + DsF2)) + 
     $        2d0 * ht**2 * mt**2 * (F1 + osdr*DF1)

c     some of the functions have poles in s2t=0 or in A=0. 
c     when necessary we consider the residues:
         
      elseif(s2t.eq.0.and.A.eq.0) then
         
         S11 = 0d0
         S12 = 0d0
         S22 = 2 * ht**2 * mt**2 * (F1 + osdr*DF1)

      elseif(s2t.eq.0.and.A.ne.0d0) then 

         call strresfuncs(t,mg,T1,T2,q,F2_s,sF2_A,sF3_A)     

         S11 = 0d0
         S12 = ht**2 * mt * mu * (F2_s + osdr*DF2)
         S22 = 2d0 * ht**2 * mt**2 * (F1 + osdr*DF1) +
     $        2d0 * ht**2 * mt * A * (F2_s + osdr*DF2)

      elseif(s2t.ne.0d0.and.A.eq.0) then

         call strresfuncs(t,mg,T1,T2,q,F2_s,sF2_A,sF3_A)     

         S11 = .5d0 * ht**2 * mu**2 * s2t**2 * (F3 + osdr*DF3)
         S12 = .5d0 * ht**2 * mu * s2t**2 * (sF3_A + osdr*DsF3) + 
     $        ht**2 * mt * mu * s2t * (F2 + osdr*DF2)
         S22 = 2d0 * ht**2 * mt**2 * (F1 + osdr*DF1) +
     $        2d0 * ht**2 * mt * s2t * (sF2_A + osdr*DsF2)
 
      endif

      S11 = k*S11
      S12 = k*S12
      S22 = k*S22

      return
      end

*
***********************************************************************
*
      
      subroutine strfuncs(t,mg,T1,T2,s2t,c2t,q,F1,F2,F3)

      implicit none
      real*8 t,mg,T1,T2,s2t,c2t,q,F1,F2,F3
      real*8 strF1ab,strF1c,strF2ab,strF2c,strF3ab,strF3c

      F1 = strF1ab(t,T1,T2,s2t,c2t,q) 
     $     + strF1c(t,mg,T1,s2t,q)
     $     + strF1c(t,mg,T2,-s2t,q)
      
      F2 = strF2ab(T1,T2,s2t,c2t,q) 
     $     + strF2c(t,mg,T1,T2,s2t,q)
     $     - strF2c(t,mg,T2,T1,-s2t,q)

      F3 = strF3ab(T1,T2,s2t,c2t,q) 
     $     + strF3c(t,mg,T1,T2,s2t,q)
     $     + strF3c(t,mg,T2,T1,-s2t,q)

      return
      end

*
*********************************************************************
*
            
      function strF1ab(t,T1,T2,s2t,c2t,q)

      implicit none
      real*8 t,T1,T2,s2t,c2t,q
      real*8 strF1ab

      strF1ab =                    ! eq. (32) 
     $     -6*(1-Log(t/q))+5*Log(T1*T2/t**2)+Log(T1*T2/t**2)**2
     $     +8*Log(t/q)**2-4*Log(T1/q)**2-4*Log(T2/q)**2
     $     -c2t**2*(2-Log(T1/q)-Log(T2/q)-Log(T1/T2)**2)
     $     -s2t**2*(T1/T2*(1-Log(T1/q))+T2/T1*(1-Log(T2/q)))
      
      return
      end
      
*
*********************************************************************
*

      function strF1c(t,mg,T1,s2t,q)

      implicit none
      real*8 t,g,mt,mg,T1,s2t,q
      real*8 strF1c,phi,del

      mt = sqrt(t)
      g = mg**2

      del = g**2 + t**2 + T1**2 - 2*(g*t + g*T1 + t*T1)
      
      strF1c =                     ! eq. (A1)
     $     +4*(t+g-mg*mt*s2t)/T1*(1-Log(g/q))
     $     +4*Log(t/g) - 2*Log(T1/g)
     $     +2d0/del*(4*g**2*Log(T1/g)
     $     +(g**2-T1**2+t*(10*g+3*t+2*t*g/T1-2*t**2/T1))*Log(t/g))
     $     +2*mg/mt*s2t*(Log(T1/q)**2+2*Log(t/q)*Log(T1/q))
     $     +4*mg/mt*s2t/del*(g*(T1-t-g)
     $     *Log(T1/g)+t*(T1-3*g-2*t-(t*g-t**2)/T1)*Log(t/g))
     $     +(4*g*(t+g-T1-2*mg*mt*s2t)/del
     $     -4*mg/mt*s2t)*phi(t,T1,g)

      return
      end

*
*********************************************************************
*

      function strF2ab(T1,T2,s2t,c2t,q)
      
      implicit none
      real*8 T1,T2,s2t,c2t,q
      real*8 strF2ab
      
      strF2ab =                    ! eq. (33)
     $     5*Log(T1/T2)-3*(Log(T1/q)**2-Log(T2/q)**2)
     $     +c2t**2*(5*Log(T1/T2)
     $     -(T1+T2)/(T1-T2)*Log(T1/T2)**2
     $     -2/(T1-T2)*(T1*Log(T1/q)-T2*Log(T2/q))*Log(T1/T2))
     $     +s2t**2*(T1/T2*(1-Log(T1/q))-T2/T1*(1-Log(T2/q)))
      
      return
      end

*
*********************************************************************
*

      function strF2c(t,mg,T1,T2,s2t,q)

      implicit none
      real*8 t,g,mg,mt,T1,T2,s2t,q
      real*8 strF2c,phi,del

      mt = sqrt(t)
      g = mg**2

      del = g**2 + t**2 + T1**2 - 2*(g*t + g*T1 + t*T1)
      
      strF2c =                     ! eq. (A2)
     $     4*(t+g)/T1-4*mg/mt*s2t/(T1-T2)*(3*T1-t*T2/T1)
     $     +2*mg/mt*s2t/(T1-T2)*(
     $     (4*t+5*T1+T2)*Log(T1/q)-2*t*T2/T1*Log(g/q))
     $     -4*(g+t)/T1*Log(g/q) - 2*Log(T1/g) 
     $     +2/del*(2*g*(g+t-T1)*Log(T1/g)
     $     +2*t*(3*g+2*t-T1+(g*t-t**2)/T1)*Log(t/g))
     $     -4*mg*mt*s2t/del/T1*(2*g*T1*Log(T1/g)-
     $     ((t-T1)**2-g*(t+T1))*Log(t/g))
     $     -8*mg*mt/s2t/(T1-T2)*(Log(T1/q)-Log(t/q)*Log(T1/q))
     $     -mg/mt*s2t/(T1-T2)*((T1+T2)*Log(T1/q)**2
     $     +(10*t-2*g+T1+T2)*Log(t/q)*Log(T1/q)
     $     +(2*g-2*t+T1+T2)*Log(T1/q)*Log(g/q))
     $     +(8*g*t/del-8*mg*mt/s2t/(T1-T2)
     $     +2*s2t/mg/mt/(T1-T2)*(4*g*t-del)
     $     +s2t/mg/mt/del*(T1-g-t)**3)*phi(t,T1,g)

      return
      end

*
*********************************************************************
*
      
      function strF3ab(T1,T2,s2t,c2t,q)

      implicit none
      real*8 T1,T2,s2t,c2t,q
      real*8 strF3ab

      strF3ab =                    ! eq. (34)
     $     (3+9*c2t**2)*(2-(T1+T2)/(T1-T2)*Log(T1/T2))
     $     +4-(3+13*c2t**2)/(T1-T2)*(T1*Log(T1/q)-T2*Log(T2/q))
     $     +3*(T1+T2)/(T1-T2)*(Log(T1/q)**2-Log(T2/q)**2)
     $     -c2t**2*(4-((T1+T2)/(T1-T2))**2*Log(T1/T2)**2
     $     -6*(T1+T2)/(T1-T2)**2
     $     *(T1*Log(T1/q)-T2*Log(T2/q))*Log(T1/T2))
     $     -s2t**2*(T1/T2+T2/T1 + 2*Log(T1*T2/q**2)
     $     -T1**2/T2/(T1-T2)*Log(T1/q)
     $     +T2**2/T1/(T1-T2)*Log(T2/q))

      return
      end
      
*
*********************************************************************
*

      function strF3c(t,mg,T1,T2,s2t,q)

      implicit none
      real*8 t,g,mt,mg,T1,T2,s2t,q
      real*8 strF3c,phi,del

      mt = sqrt(t)
      g = mg**2

      del = g**2 + t**2 + T1**2 - 2*(g*t + g*T1 + t*T1)
      
      strF3c =                     ! eq. (A3)
     $     -4*T2/T1/(T1-T2)*(g+t)
     $     +4*mg*mt*s2t/(T1-T2)**2*(21*T1-T2**2/T1)
     $     +4/(T1-T2)*(g*T2/T1*Log(g/q)-2*(t+g)*Log(T1/q))
     $     -24*mg*mt*s2t/(T1-T2)**2*(3*T1+T2)*Log(T1/q)
     $     +4*t/T1/del*(2*g*T1*Log(T1/q)-g*(g-t+T1)*Log(g/q)+
     $     (g*(T+T1)-(t-T1)**2)*Log(t/q))
     $     -4*mg*mt*s2t/T1/del*(t*(g-t+T1)*Log(t/q)
     $     -g*(g-t-T1)*Log(g/q)+T1*(g+t-T1)*Log(T1/q))
     $     +2*(2*g+2*t-T1-T2)/(T1-T2)*Log(g*t/q**2)*Log(T1/q)
     $     +12*mg*mt*s2t/(T1-T2)**2*(2*(g-t)*Log(g/t)*Log(T1/q)
     $     +(T1+T2)*Log(t*g/q**2)*Log(T1/q))
     $     +8*mg*mt/s2t/(T1-T2)**2*
     $     (-8*T1+2*(3*T1+T2)*Log(T1/q)-2*(g-t)*Log(g/t)*Log(T1/q)
     $     -(T1+T2)*Log(t*g/q**2)*Log(T1/q))
     $     -((8/s2t-12*s2t)*mt/mg/(T1-T2)**2
     $     *(2*del+(g+t-T1)*(T1-T2))
     $     +(4*del+8*g*t)/g/(T1-T2)+2*(g+t-T1)/g
     $     -4*t*(g+t-T1- 2*mg*mt*s2t)/del)*phi(t,T1,g)

      return
      end

*
*********************************************************************
*

      subroutine strsfuncs(mg,T1,T2,q,A,sF2,sF3)
      
c     shift to the Fi functions due to the renormalization of A
      
      implicit none
      real*8 mg,T1,T2,q,A,sF2,sF3
      
      sF2  = mg/A *       ! eq. (35)
     $     2d0*(dlog(T2/q)**2 - dlog(T1/q)**2)
      
      sF3  = mg/A *       ! eq. (36)
     $     (8d0 - 2d0*(T1+T2)/(T1-T2)*(dlog(T2/q)**2 - dlog(T1/q)**2)
     $     + 8d0/(T1-T2)*(T2*dlog(T2/q) - T1*dlog(T1/q)))
      
      return
      end

*
*********************************************************************
*

      subroutine strresfuncs(t,mg,T1,T2,q,F2_s,sF2_A,sF3_A)
      
c     residues of some singular functions for s2t=0 and for A=0
      
      implicit none      
      real*8 t,g,mt,mg,T1,T2,q,sF2_A,sF3_A,F2_s,phi
      
      mt = sqrt(t)
      g = mg**2

      F2_s =  -8*mg*mt/(T1-T2)*(
     $     (Log(T1/q)-Log(t/q)*Log(T1/q)+phi(t,T1,g))-
     $     (Log(T2/q)-Log(t/q)*Log(T2/q)+phi(t,T2,g)))
      
      sF2_A  = mg*
     $     2d0*(dlog(T2/q)**2 - dlog(T1/q)**2)
      
      sF3_A  = mg*
     $     (8d0 - 2d0*(T1+T2)/(T1-T2)*(dlog(T2/q)**2 - dlog(T1/q)**2)
     $     + 8d0/(T1-T2)*(T2*dlog(T2/q) - T1*dlog(T1/q)))

      return
      end
      
*     
***********************************************************************
*

      subroutine strdfuncs(t,mg,T1,T2,s2t,c2t,q,At,X,
     $     DF1,DF2,DF3,DsF2,DsF3)
      
c     shift of the parameters from DRbar to On-Shell scheme
 
      implicit none      
      real*8 t,g,mt,mg,T1,T2,s2t,c2t,q,At,X,myB0
      real*8 DF1,DF2,DF3,DsF2,DsF3
      real*8 msdr
      real*8 F1o,F2o,F3o,dm1,dm2,dmt,dAt,dth,ds2t

      msdr = -5d0
      mt = Sqrt(t)
      g = mg**2

      F1o = Log(T1/q) + Log(T2/q) - 2d0*Log(t/q) ! eq. (31)
      F2o = Log(T1/q) - Log(T2/q) 
      F3o = 2d0 - (T1+T2)/(T1-T2)*(Log(T1/q) - Log(T2/q))

      dmt =                     ! eq. (B2)
     $     mt*(3*Log(t/q) + msdr + .5d0*(2*g/t*(Log(g/q)-1)
     $      -T1/t*(Log(T1/q)-1) - T2/t*(Log(T2/q)-1)
     $      +(g+t-T1 - 2*s2t*mg*mt)/t*myB0(t,g,T1,q)
     $      +(g+t-T2 + 2*s2t*mg*mt)/t*myB0(t,g,T2,q)))

      dm1 =                     ! eq. (B3)
     $     T1*(3*Log(T1/q) - 7 - c2t**2*(Log(T1/q)-1)
     $      -s2t**2*T2/T1*(Log(T2/q)-1) + 2*(
     $      g/T1*(Log(g/q)-1) + t/T1*(Log(t/q)-1)
     $      +(T1-g-t + 2*s2t*mg*mt)/T1*myB0(T1,t,g,q)))
      
      dm2 =                     ! eq. (B4)
     $     T2*(3*Log(T2/q) - 7 - c2t**2*(Log(T2/q)-1)
     $      -s2t**2*T1/T2*(Log(T1/q)-1) + 2*(
     $      g/T2*(Log(g/q)-1) + t/T2*(Log(t/q)-1)
     $      +(T2-g-t - 2*s2t*mg*mt)/T2*myB0(T2,t,g,q)))

c     On-Shell theta-stop: asymmetric definition used in FeynHiggs
      dth = (4d0*mg*mt*c2t*myB0(T1,t,g,q) +
     $     c2t*s2t*(T2*(1d0-Log(T2/q))-T1*(1d0-Log(T1/q))))/(T1-T2)      

c$$$c     On-Shell theta-stop: eq. (B6)-(B7) of DSZ 
c$$$      dth = (4d0*mg*mt*c2t*(myB0(T1,t,g,q)+myB0(T2,t,g,q)) +
c$$$     $     2d0*c2t*s2t*(T2*(1d0-Log(T2/q))-T1*(1d0-Log(T1/q))))/
c$$$     $     2d0/(T1-T2)      

      ds2t = 2d0*c2t*dth

      dAt = ((dm1-dm2)/(T1-T2) + ds2t/s2t - dmt/mt)*X ! eq. (B8)

      DF1 = dm1/T1 + dm2/T2 - 4d0*dmt/mt + 4d0*dmt/mt*F1o ! eq. (37)
      DF2 = dm1/T1 - dm2/T2 + (3d0*dmt/mt + ds2t/s2t)*F2o ! eq. (38)
      DF3 = (2d0*T1*T2/(T1-T2)**2*Log(T1/T2) - (T1+T2)/(T1-T2)) 
     $     *(dm1/T1-dm2/T2) + (2d0*dmt/mt + 2d0*ds2t/s2t)*F3o ! eq. (39)

      DsF2 = dAt/At * F2o       ! eq. (40)
      DsF3 = dAt/At * F3o       ! eq. (41)

c     residues of some singular functions for s2t=0 and for A=0

      if(s2t.eq.0d0) then
         DF2  = ds2t*F2o
         DsF2 = ds2t*X/At * F2o
      endif

      if(At.eq.0d0) then
         DsF2 = dAt * F2o
         DsF3 = dAt * F3o
      endif

      return
      end












