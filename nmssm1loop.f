
      subroutine gettadS(g,gp,ll,kk,ht,hb,htau,v1,v2,xx,
     $     Ak,Al,At,Ab,Atau,Q,tadS)
      
      implicit none

      double precision g,gp,ll,kk,ht,hb,htau,v1,v2,xx,
     $     Ak,Al,At,Ab,Atau,Q,tadS(3)
      
      double precision mstop(2),msbot(2),mstau(2),msup(2),msdown(2),
     $     msel(2),msnutau,msnue,Rt(2,2),Rb(2,2),Rtau(2,2),
     $     NN(5,5),UU(2,2),VV(2,2),mch(2),mne(5),
     $     mhh(3),maa(3),mhc(2),RS(3,3),RP(3,3),RC(2,2),myA0,pi,
     $     mw2,mz2,mt,mb,mtau,gb2,sq2

      double precision lsstt(3,3,2,2),lssbb(3,3,2,2),lsstata(3,3,2,2),
     $     lssntnt(3,3),lssuu(3,3,2,2),lssdd(3,3,2,2),lssee(3,3,2,2),
     $     lssnn(3,3),lstt(3,2,2),lsbb(3,2,2),lstata(3,2,2),lsntnt(3),
     $     lsuu(3,2,2),lsdd(3,2,2),lsee(3,2,2),lsnn(3),
     $     lsshh(3,3,3,3),lssaa(3,3,3,3),lshh(3,3,3),lsaa(3,3,3),
     $     lsscc(3,3,2,2),lscc(3,2,2),lschch(3,2,2),lsnene(3,5,5)

      integer i,k

      common/TREE_MASSES/mhh,maa,mhc,RS,RP,RC,mch,UU,VV,mne,NN,
     $     mstop,msbot,mstau,Rt,Rb,Rtau,msup,msdown,msel,msnutau,msnue
      
      pi = 4d0*atan(1.d0)
      sq2 = sqrt(2d0)

c     the gauge and fermion contributions

      gb2 = (g**2+gp**2)/2d0

      mw2 = g**2/2d0*(v1**2+v2**2)
      mz2 = gb2*(v1**2+v2**2)

      mt = ht*v2
      mb = hb*v1
      mtau = htau*v1

      tadS(1) = 0d0
     $     -6*sq2*hb*mb*myA0(mb**2,Q**2)
     $     -2*sq2*htau*mtau*myA0(mtau**2,Q**2)
     $     +3d0*v1/sq2*g**2*myA0(mw2,Q**2)+3d0*v1/sq2*gb2*myA0(mz2,Q**2)

      tadS(2) = 0d0
     $     -6*sq2*ht*mt*myA0(mt**2,Q**2)
     $     +3d0*v2/sq2*g**2*myA0(mw2,Q**2)+3d0*v2/sq2*gb2*myA0(mz2,Q**2)

      tadS(3) = 0d0

c     the sfermion contributions

      call coupl_s_sf(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     $     Rt,Rb,Rtau,lsstt,lssbb,lsstata,lssntnt,lssuu,lssdd,
     $     lssee,lssnn,lstt,lsbb,lstata,lsntnt,lsuu,lsdd,lsee,lsnn)

      do i=1,3

         do k = 1,2
            
            tadS(i) = tadS(i)
     $           + 3*lstt(i,k,k)*myA0(mstop(k)**2,Q**2)               
     $           + 3*lsbb(i,k,k)*myA0(msbot(k)**2,Q**2)               
     $           + lstata(i,k,k)*myA0(mstau(k)**2,Q**2)               
     $           + 6*lsuu(i,k,k)*myA0(msup(k)**2,Q**2)
     $           + 6*lsdd(i,k,k)*myA0(msdown(k)**2,Q**2)
     $           + 2*lsee(i,k,k)*myA0(msel(k)**2,Q**2)

         enddo

         tadS(i) = tadS(i)      ! add sneutrinos (no sum)
     $        + 2*lsnn(i)*myA0(msnue**2,Q**2)
     $        + lsntnt(i)*myA0(msnutau**2,Q**2)               

      enddo

c     the Higgs contributions

      call coupl_s_hh(g,gp,ll,kk,v1,v2,xx,Al,Ak,RS,RP,
     $     lsshh,lssaa,lshh,lsaa,lsscc,lscc)

      do i=1,3

         do k = 1,3             ! neutral Higgses
            
            tadS(i) = tadS(i)
     $           + lshh(i,k,k)*myA0(mhh(k)**2,Q**2)               
     $           + lsaa(i,k,k)*myA0(maa(k)**2,Q**2)               
            
         enddo

         do k = 1,2             ! charged Higgses
            
            tadS(i) = tadS(i)
     $           + lscc(i,k,k)*myA0(mhc(k)**2,Q**2)

         enddo

      enddo

c     the chargino and neutralino contributions

      call coupl_s_ino(g,gp,ll,kk,NN,UU,VV,lsnene,lschch)

      do i=1,3
         
         do k = 1,5             ! neutralinos
            tadS(i) = tadS(i)
     $           - 4*lsnene(i,k,k)*mne(k)*myA0(mne(k)**2,Q**2)
         enddo

         do k = 1,2             ! charginos
            tadS(i) = tadS(i)
     $           - 4*lschch(i,k,k)*mch(k)*myA0(mch(k)**2,Q**2)
         enddo

      enddo

      do i=1,3
         tadS(i) = tadS(i)/16d0/pi**2
      enddo

      return
      end
      
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     

      subroutine getPiSS(g,gp,ll,kk,ht,hb,htau,v1,v2,xx,
     $     Ak,Al,At,Ab,Atau,p,Q,piSS)
      
      implicit none

      double precision g,gp,ll,kk,ht,hb,htau,v1,v2,xx,
     $     Ak,Al,At,Ab,Atau,p,Q,piSS(3,3)
      
      double precision mstop(2),msbot(2),mstau(2),msup(2),msdown(2),
     $     msel(2),msnutau,msnue,Rt(2,2),Rb(2,2),Rtau(2,2),
     $     NN(5,5),UU(2,2),VV(2,2),mch(2),mne(5),
     $     mhh(3),maa(3),mhc(2),RS(3,3),RP(3,3),RC(2,2),
     $     myB0NM,myA0,myF,myG,pi,cb,sb,mw2,mz2,mt2,mb2,mtau2,gb2

      double precision lsstt(3,3,2,2),lssbb(3,3,2,2),lsstata(3,3,2,2),
     $     lssntnt(3,3),lssuu(3,3,2,2),lssdd(3,3,2,2),lssee(3,3,2,2),
     $     lssnn(3,3),lstt(3,2,2),lsbb(3,2,2),lstata(3,2,2),lsntnt(3),
     $     lsuu(3,2,2),lsdd(3,2,2),lsee(3,2,2),lsnn(3),
     $     lsshh(3,3,3,3),lssaa(3,3,3,3),lshh(3,3,3),lsaa(3,3,3),
     $     lsscc(3,3,2,2),lscc(3,2,2),lschch(3,2,2),lsnene(3,5,5)

      integer i,j,k,l
      
      common/TREE_MASSES/mhh,maa,mhc,RS,RP,RC,mch,UU,VV,mne,NN,
     $     mstop,msbot,mstau,Rt,Rb,Rtau,msup,msdown,msel,msnutau,msnue
      
      pi = 4d0*atan(1.d0)

      do i=1,3                  ! initialize
         do j = 1,3
            PiSS(i,j) = 0d0
         enddo
      enddo

c     the gauge and fermion contributions

      cb = v1/sqrt(v1**2+v2**2)
      sb = v2/sqrt(v1**2+v2**2)
      gb2 = (g**2+gp**2)/2d0

      mw2 = g**2/2d0*(v1**2+v2**2)
      mz2 = gb2*(v1**2+v2**2)

      mt2 = ht**2*v2**2
      mb2 = hb**2*v1**2
      mtau2 = htau**2*v1**2

      PiSS(1,1) = 
     $     3*hb**2*((p**2-4*mb2)*myB0NM(p**2,mb2,mb2,Q**2)
     $     -2*myA0(mb2,Q**2))
     $     +htau**2*((p**2-4*mtau2)*myB0NM(p**2,mtau2,mtau2,Q**2)
     $     -2*myA0(mtau2,Q**2))
     $     +7d0/2d0*cb**2*(g**2*mw2*myB0NM(p**2,mw2,mw2,Q**2)
     $     +gb2*mz2*myB0NM(p**2,mz2,mz2,Q**2))
     $     +2*g**2*myA0(mw2,Q**2)+2*gb2*myA0(mz2,Q**2)         

      PiSS(1,2) =
     $     7d0/2d0*cb*sb*(g**2*mw2*myB0NM(p**2,mw2,mw2,Q**2)
     $     +gb2*mz2*myB0NM(p**2,mz2,mz2,Q**2))
      
      PiSS(2,1) = PiSS(1,2)

      PiSS(2,2) = 3*ht**2*((p**2-4*mt2)*myB0NM(p**2,mt2,mt2,Q**2)
     $     -2*myA0(mt2,Q**2))
     $     +7d0/2d0*sb**2*(g**2*mw2*myB0NM(p**2,mw2,mw2,Q**2)
     $     +gb2*mz2*myB0NM(p**2,mz2,mz2,Q**2))
     $     +2*g**2*myA0(mw2,Q**2)+2*gb2*myA0(mz2,Q**2)         

c     pseudoscalar-gauge contribution

      do i = 1,2
         do j = 1,2

            do k = 1,3
               PiSS(i,j) = PiSS(i,j) + (-1)**(i+j)*
     $              gb2/2d0*RP(k,i)*RP(k,j)*myF(p**2,maa(k)**2,mz2,Q**2)
            enddo

            do k = 1,2
               PiSS(i,j) = PiSS(i,j) + (-1)**(i+j)*
     $             g**2/2d0*RC(k,i)*RC(k,j)*myF(p**2,mhc(k)**2,mw2,Q**2)
            enddo

         enddo
      enddo

c     the sfermion contributions

      call coupl_s_sf(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     $     Rt,Rb,Rtau,lsstt,lssbb,lsstata,lssntnt,lssuu,lssdd,
     $     lssee,lssnn,lstt,lsbb,lstata,lsntnt,lsuu,lsdd,lsee,lsnn)
      
      
      do i=1,3
         do j = 1,3
            do k = 1,2

               PiSS(i,j) = PiSS(i,j)
     $              + 6*lsstt(i,j,k,k)*myA0(mstop(k)**2,Q**2)               
     $              + 6*lssbb(i,j,k,k)*myA0(msbot(k)**2,Q**2)               
     $              + 2*lsstata(i,j,k,k)*myA0(mstau(k)**2,Q**2)               
     $              + 12*lssuu(i,j,k,k)*myA0(msup(k)**2,Q**2)
     $              + 12*lssdd(i,j,k,k)*myA0(msdown(k)**2,Q**2)
     $              + 4*lssee(i,j,k,k)*myA0(msel(k)**2,Q**2)
               
               do l = 1,2

                  PiSS(i,j) = PiSS(i,j)
     $                 + 3*lstt(i,k,l)*lstt(j,l,k)
     $                 *myB0NM(p**2,mstop(k)**2,mstop(l)**2,Q**2)               
     $                 + 3*lsbb(i,k,l)*lsbb(j,l,k)
     $                 *myB0NM(p**2,msbot(k)**2,msbot(l)**2,Q**2)               
     $                 + lstata(i,k,l)*lstata(j,l,k)
     $                 *myB0NM(p**2,mstau(k)**2,mstau(l)**2,Q**2)               
     $                 + 6*lsuu(i,k,l)*lsuu(j,l,k)
     $                 *myB0NM(p**2,msup(k)**2,msup(l)**2,Q**2)
     $                 + 6*lsdd(i,k,l)*lsdd(j,l,k)
     $                 *myB0NM(p**2,msdown(k)**2,msdown(l)**2,Q**2)
     $                 + 2*lsee(i,k,l)*lsee(j,l,k)
     $                 *myB0NM(p**2,msel(k)**2,msel(l)**2,Q**2)
               enddo
            enddo

            PiSS(i,j) = PiSS(i,j) ! add sneutrinos (no sum)
     $           + 4*lssnn(i,j)*myA0(msnue**2,Q**2)
     $           + 2*lsnn(i)*lsnn(j)
     $           *myB0NM(p**2,msnue**2,msnue**2,Q**2)
     $           + 2*lssntnt(i,j)*myA0(msnutau**2,Q**2)               
     $           + lsntnt(i)*lsntnt(j)
     $           *myB0NM(p**2,msnutau**2,msnutau**2,Q**2)            

         enddo
      enddo
      
c     the Higgs contributions

      call coupl_s_hh(g,gp,ll,kk,v1,v2,xx,Al,Ak,RS,RP,
     $     lsshh,lssaa,lshh,lsaa,lsscc,lscc)

      do i=1,3
         do j = 1,3

            do k = 1,3          ! neutral Higgses

               PiSS(i,j) = PiSS(i,j)
     $              + 2*lsshh(i,j,k,k)*myA0(mhh(k)**2,Q**2)               
     $              + 2*lssaa(i,j,k,k)*myA0(maa(k)**2,Q**2)               

               do l = 1,3

                  PiSS(i,j) = PiSS(i,j)
     $                 + 2*lshh(i,k,l)*lshh(j,k,l)
     $                 *myB0NM(p**2,mhh(k)**2,mhh(l)**2,Q**2)               
     $                 + 2*lsaa(i,k,l)*lsaa(j,k,l)
     $                 *myB0NM(p**2,maa(k)**2,maa(l)**2,Q**2)               

               enddo
            enddo

            do k = 1,2          ! charged Higgses
               
               PiSS(i,j) = PiSS(i,j)
     $              + 2*lsscc(i,j,k,k)*myA0(mhc(k)**2,Q**2)

               do l = 1,2
                  
                  PiSS(i,j) = PiSS(i,j)
     $                 + lscc(i,k,l)*lscc(j,l,k)
     $                 *myB0NM(p**2,mhc(k)**2,mhc(l)**2,Q**2)

               enddo
               
            enddo
         enddo
      enddo

c     the chargino and neutralino contributions

      call coupl_s_ino(g,gp,ll,kk,NN,UU,VV,lsnene,lschch)
      
      do i=1,3
         do j = 1,3

            do k = 1,5          ! neutralinos
               do l = 1,5
                  
                  PiSS(i,j) = PiSS(i,j)
     $                 + 4*lsnene(i,k,l)*lsnene(j,k,l)*(
     $                 myG(p**2,mne(k)**2,mne(l)**2,Q**2)               
     $                 -2*mne(k)*mne(l)
     $                 *myB0NM(p**2,mne(k)**2,mne(l)**2,Q**2))               
               enddo
            enddo

            do k = 1,2          ! charginos
               do l = 1,2
                  
                  PiSS(i,j) = PiSS(i,j)
     $                 + 2*(lschch(i,k,l)*lschch(j,k,l)*
     $                 myG(p**2,mch(k)**2,mch(l)**2,Q**2)               
     $                 -2*lschch(i,k,l)*lschch(j,l,k)*mch(k)*mch(l)
     $                 *myB0NM(p**2,mch(k)**2,mch(l)**2,Q**2))               
               enddo
            enddo

         enddo
      enddo

      do i=1,3
         do j = 1,3
            PiSS(i,j) = PiSS(i,j)/16d0/pi**2
         enddo
      enddo

      return
      end
      
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     

c$$$      subroutine getpipp_(g,gp,ll,kk,ht,hb,htau,v1,v2,xx,
c$$$     $     Ak,Al,At,Ab,Atau,p,Q,piPP)
 
c Temp change to input
   
      subroutine getPiPP(g,gp,ll,kk,ht,hb,htau,v1,v2,xx,
     $     Ak,Al,At,Ab,Atau,p,Q,piPP)
  
      implicit none

      double precision g,gp,ll,kk,ht,hb,htau,v1,v2,xx,
     $     Ak,Al,At,Ab,Atau,p,Q,piPP(3,3)
      
c New variables for temp test
      double precision  piPP11,piPP12,piPP13,piPP22,piPP23,piPP33

      double precision mstop(2),msbot(2),mstau(2),msup(2),msdown(2),
     $     msel(2),msnutau,msnue,Rt(2,2),Rb(2,2),Rtau(2,2),
     $     NN(5,5),UU(2,2),VV(2,2),mch(2),mne(5),
     $     mhh(3),maa(3),mhc(2),RS(3,3),RP(3,3),RC(2,2),
     $     myB0NM,myA0,myF,myG,pi,cb,sb,mw2,mz2,mt2,mb2,mtau2,gb2,ghost

      double precision lpptt(3,3,2,2),lppbb(3,3,2,2),lpptata(3,3,2,2),
     $     lppntnt(3,3),lppuu(3,3,2,2),lppdd(3,3,2,2),lppee(3,3,2,2),
     $     lppnn(3,3),lptt(3,2,2),lpbb(3,2,2),lptata(3,2,2),lpntnt(3),
     $     lpuu(3,2,2),lpdd(3,2,2),lpee(3,2,2),lpnn(3),
     $     lpphh(3,3,3,3),lppaa(3,3,3,3),lpah(3,3,3),
     $     lppcc(3,3,2,2),lpcc(3,2,2),lpchch(3,2,2),lpnene(3,5,5)

      integer i,j,k,l
      
      common/TREE_MASSES/mhh,maa,mhc,RS,RP,RC,mch,UU,VV,mne,NN,
     $     mstop,msbot,mstau,Rt,Rb,Rtau,msup,msdown,msel,msnutau,msnue
      
      pi = 4d0*atan(1.d0)

      do i=1,3                  ! initialize
         do j = 1,3
            PiPP(i,j) = 0d0
         enddo
      enddo

c     the gauge and fermion contributions

      cb = v1/sqrt(v1**2+v2**2)
      sb = v2/sqrt(v1**2+v2**2)
      gb2 = (g**2+gp**2)/2d0

      mw2 = g**2/2d0*(v1**2+v2**2)
      mz2 = gb2*(v1**2+v2**2)

      mt2 = ht**2*v2**2
      mb2 = hb**2*v1**2
      mtau2 = htau**2*v1**2

      ghost = g**2/2d0*mw2*myB0NM(p**2,mw2,mw2,Q**2)

      PiPP(1,1) = 
     $     3*hb**2*(p**2*myB0NM(p**2,mb2,mb2,Q**2)-2*myA0(mb2,Q**2))
     $     +htau**2*(p**2*myB0NM(p**2,mtau2,mtau2,Q**2)
     $     -2*myA0(mtau2,Q**2))
     $     +2*g**2*myA0(mw2,Q**2)+2*gb2*myA0(mz2,Q**2)
     $     +cb**2*ghost

      PiPP(1,2) = 
     $     -cb*sb*ghost
      
      PiPP(2,1) = PiPP(1,2)

      PiPP(2,2) = 
     $     3*ht**2*(p**2*myB0NM(p**2,mt2,mt2,Q**2)-2*myA0(mt2,Q**2))
     $     +2*g**2*myA0(mw2,Q**2)+2*gb2*myA0(mz2,Q**2)         
     $     +sb**2*ghost

c     scalar-Z contribution

      do i = 1,2
         do j = 1,2

            do k = 1,3
               PiPP(i,j) = PiPP(i,j) + (-1)**(i+j)*
     $              gb2/2d0*RS(k,i)*RS(k,j)*myF(p**2,mhh(k)**2,mz2,Q**2)
            enddo

            do k = 1,2
               PiPP(i,j) = PiPP(i,j) + 
     $             g**2/2d0*RC(k,i)*RC(k,j)*myF(p**2,mhc(k)**2,mw2,Q**2)
            enddo

         enddo
      enddo

c     the sfermion contributions

      call coupl_p_sf(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     $     Rt,Rb,Rtau,lpptt,lppbb,lpptata,lppntnt,lppuu,lppdd,
     $     lppee,lppnn,lptt,lpbb,lptata,lpntnt,lpuu,lpdd,lpee,lpnn)
      
      
      do i=1,3
         do j = 1,3
            do k = 1,2

               PiPP(i,j) = PiPP(i,j)
     $              + 6*lpptt(i,j,k,k)*myA0(mstop(k)**2,Q**2)               
     $              + 6*lppbb(i,j,k,k)*myA0(msbot(k)**2,Q**2)               
     $              + 2*lpptata(i,j,k,k)*myA0(mstau(k)**2,Q**2)               
     $              + 12*lppuu(i,j,k,k)*myA0(msup(k)**2,Q**2)
     $              + 12*lppdd(i,j,k,k)*myA0(msdown(k)**2,Q**2)
     $              + 4*lppee(i,j,k,k)*myA0(msel(k)**2,Q**2)
               
               do l = 1,2

                  PiPP(i,j) = PiPP(i,j)
     $                 - 3*lptt(i,k,l)*lptt(j,l,k)
     $                 *myB0NM(p**2,mstop(k)**2,mstop(l)**2,Q**2)               
     $                 - 3*lpbb(i,k,l)*lpbb(j,l,k)
     $                 *myB0NM(p**2,msbot(k)**2,msbot(l)**2,Q**2)               
     $                 - lptata(i,k,l)*lptata(j,l,k)
     $                 *myB0NM(p**2,mstau(k)**2,mstau(l)**2,Q**2)               
     $                 - 6*lpuu(i,k,l)*lpuu(j,l,k)
     $                 *myB0NM(p**2,msup(k)**2,msup(l)**2,Q**2)
     $                 - 6*lpdd(i,k,l)*lpdd(j,l,k)
     $                 *myB0NM(p**2,msdown(k)**2,msdown(l)**2,Q**2)
     $                 - 2*lpee(i,k,l)*lpee(j,l,k)
     $                 *myB0NM(p**2,msel(k)**2,msel(l)**2,Q**2)
               enddo
            enddo

            PiPP(i,j) = PiPP(i,j) ! add sneutrinos (no sum)
     $           + 4*lppnn(i,j)*myA0(msnue**2,Q**2)
     $           - 2*lpnn(i)*lpnn(j)
     $           *myB0NM(p**2,msnue**2,msnue**2,Q**2)
     $           + 2*lppntnt(i,j)*myA0(msnutau**2,Q**2)               
     $           - lpntnt(i)*lpntnt(j)
     $           *myB0NM(p**2,msnutau**2,msnutau**2,Q**2)            

         enddo
      enddo
      
c     the Higgs contributions

      call coupl_p_hh(g,gp,ll,kk,v1,v2,xx,Al,Ak,RS,RP,
     $     lpphh,lppaa,lpah,lppcc,lpcc)

      do i=1,3
         do j = 1,3

            do k = 1,3         ! neutral Higgses

               PiPP(i,j) = PiPP(i,j)
     $              + 2*lpphh(i,j,k,k)*myA0(mhh(k)**2,Q**2)               
     $              + 2*lppaa(i,j,k,k)*myA0(maa(k)**2,Q**2)               

               do l = 1,3

                  PiPP(i,j) = PiPP(i,j)
     $                 + lpah(i,k,l)*lpah(j,k,l)
     $                 *myB0NM(p**2,maa(k)**2,mhh(l)**2,Q**2)               

               enddo
            enddo

            do k = 1,2          ! charged Higgses
               
               PiPP(i,j) = PiPP(i,j)
     $              + 2*lppcc(i,j,k,k)*myA0(mhc(k)**2,Q**2)

               do l = 1,2
                  
                  PiPP(i,j) = PiPP(i,j)
     $                 - lpcc(i,k,l)*lpcc(j,l,k)
     $                 *myB0NM(p**2,mhc(k)**2,mhc(l)**2,Q**2)

               enddo
               
            enddo
         enddo
      enddo

c     the chargino and neutralino contributions

      call coupl_p_ino(g,gp,ll,kk,NN,UU,VV,lpnene,lpchch)
      
      do i=1,3
         do j = 1,3

            do k = 1,5          ! neutralinos
               do l = 1,5
                  
                  PiPP(i,j) = PiPP(i,j)
     $                 + 4*lpnene(i,k,l)*lpnene(j,k,l)*(
     $                 myG(p**2,mne(k)**2,mne(l)**2,Q**2)               
     $                 +2*mne(k)*mne(l)
     $                 *myB0NM(p**2,mne(k)**2,mne(l)**2,Q**2))               
               enddo
            enddo

            do k = 1,2          ! charginos
               do l = 1,2
                  
                  PiPP(i,j) = PiPP(i,j)
     $                 + 2*(lpchch(i,k,l)*lpchch(j,k,l)*
     $                 myG(p**2,mch(k)**2,mch(l)**2,Q**2)               
     $                 + 2*lpchch(i,k,l)*lpchch(j,l,k)*mch(k)*mch(l)
     $                 *myB0NM(p**2,mch(k)**2,mch(l)**2,Q**2))               
               enddo
            enddo

         enddo
      enddo

      do i=1,3
         do j = 1,3
            PiPP(i,j) = PiPP(i,j)/16d0/pi**2
         enddo
      enddo

      return
      end
      
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     

      subroutine getPiZZ(g,gp,ht,hb,htau,v1,v2,p,Q,piZZ)
      
      implicit none

      double precision g,gp,ht,hb,htau,v1,v2,p,Q,piZZ
      
      double precision mstop(2),msbot(2),mstau(2),msup(2),msdown(2),
     $     msel(2),msnutau,msnue,Rt(2,2),Rb(2,2),Rtau(2,2),del(2,2),
     $     NN(5,5),UU(2,2),VV(2,2),mch(2),mne(5),
     $     mhh(3),maa(3),mhc(2),RS(3,3),RP(3,3),RC(2,2),
     $     myB0NM,myB22T,myH,pi,cb,sb,mw2,mz2,mt2,mb2,mtau2,gb2,
     $     sw2,cw2,guL,guR,gdL,gdR,geL,geR,gnu,lznene(5,5),
     $     azchch(2,2),bzchch(2,2)

      integer i,j
      
      common/TREE_MASSES/mhh,maa,mhc,RS,RP,RC,mch,UU,VV,mne,NN,
     $     mstop,msbot,mstau,Rt,Rb,Rtau,msup,msdown,msel,msnutau,msnue
      
      pi = 4d0*atan(1.d0)

      PiZZ = 0d0

      cb = v1/sqrt(v1**2+v2**2)
      sb = v2/sqrt(v1**2+v2**2)
      gb2 = (g**2+gp**2)/2d0
      sw2 = gp**2/2d0/gb2
      cw2 = 1-sw2

      mw2 = g**2/2d0*(v1**2+v2**2)
      mz2 = gb2*(v1**2+v2**2)

c     the higgs and gauge contributions
   
      do i = 1,3                ! neutral scalars

         PiZZ = PiZZ            ! scalar-Z
     $        + 2*gb2*mz2*(cb*RS(i,1)+sb*RS(i,2))**2
     $        *myB0NM(p**2,mhh(i)**2,mz2,Q**2)
            
         do j = 1,3

            PiZZ = PiZZ
     $           - 2*gb2*(RS(i,1)*RP(j,1)-RS(i,2)*RP(j,2))**2
     $           *myB22T(p**2,mhh(i)**2,maa(j)**2,Q**2)

         enddo
      enddo
       
     
      do i = 1,2                ! charged scalars

         PiZZ = PiZZ            
     $        - 2*gb2*(cw2-sw2)**2*myB22T(p**2,mhc(i)**2,mhc(i)**2,Q**2)

      enddo


      PiZZ = PiZZ               ! pure gauge
     $     -4*gb2*cw2**2*(2*p**2+mw2-mz2*sw2**2/cw2)
     $     *myB0NM(p**2,mw2,mw2,Q**2)
     $     -16*gb2*cw2**2*myB22T(p**2,mw2,mw2,Q**2)
    

c     the fermion contributions

      guL =  .5d0 - 2d0/3d0*sw2
      gdL = -.5d0 + 1d0/3d0*sw2
      gnu =  .5d0
      geL = -.5d0 + sw2
      guR =  2d0/3d0*sw2
      gdR = -1d0/3d0*sw2
      geR = -sw2

      mt2 = ht**2*v2**2
      mb2 = hb**2*v1**2
      mtau2 = htau**2*v1**2
      
      PiZZ = PiZZ + 2*gb2*(
     $     (6*(guL**2+guR**2)+6*(gdL**2+gdR**2)+2*(geL**2+geR**2)
     $     +3*gnu**2)*myH(p**2,0d0,0d0,Q**2)
     $     +3*((guL**2+guR**2)*myH(p**2,mt2,mt2,Q**2)
     $     -4*guL*guR*mt2*myB0NM(p**2,mt2,mt2,Q**2))
     $     +3*((gdL**2+gdR**2)*myH(p**2,mb2,mb2,Q**2)
     $     -4*gdL*gdR*mb2*myB0NM(p**2,mb2,mb2,Q**2))
     $     +((geL**2+geR**2)*myH(p**2,mtau2,mtau2,Q**2)
     $     -4*geL*geR*mtau2*myB0NM(p**2,mtau2,mtau2,Q**2)))

c     the sfermion contributions

      del(1,1) = 1d0
      del(1,2) = 0d0
      del(2,1) = 0d0
      del(2,2) = 1d0

      do i=1,2
         do j = 1,2

            PiZZ = PiZZ - 8*gb2*(
     $           6*(guL*del(i,1)*del(j,1)-guR*del(i,2)*del(j,2))**2
     $           *myB22T(p**2,msup(i)**2,msup(j)**2,Q**2)
     $           +6*(gdL*del(i,1)*del(j,1)-gdR*del(i,2)*del(j,2))**2
     $           *myB22T(p**2,msdown(i)**2,msdown(j)**2,Q**2)
     $           +2*(geL*del(i,1)*del(j,1)-geR*del(i,2)*del(j,2))**2
     $           *myB22T(p**2,msel(i)**2,msel(j)**2,Q**2)
     $           +3*(guL*Rt(i,1)*Rt(j,1)-guR*Rt(i,2)*Rt(j,2))**2
     $           *myB22T(p**2,mstop(i)**2,mstop(j)**2,Q**2)
     $           +3*(gdL*Rb(i,1)*Rb(j,1)-gdR*Rb(i,2)*Rb(j,2))**2
     $           *myB22T(p**2,msbot(i)**2,msbot(j)**2,Q**2)
     $           +(geL*Rtau(i,1)*Rtau(j,1)-geR*Rtau(i,2)*Rtau(j,2))**2
     $           *myB22T(p**2,mstau(i)**2,mstau(j)**2,Q**2))
            
         enddo
      enddo

      PiZZ = PiZZ - 8*gb2*gnu**2*(
     $     2*myB22T(p**2,msnue**2,msnue**2,Q**2)
     $     +myB22T(p**2,msnutau**2,msnutau**2,Q**2))

c     the chargino and neutralino contributions

      call coupl_Z_ino(g,gp,NN,UU,VV,lznene,azchch,bzchch)

      do i = 1,5                ! neutralinos
         do j = 1,5
            
            PiZZ = PiZZ
     $           +4*lznene(i,j)**2*(myH(p**2,mne(i)**2,mne(j)**2,Q**2)
     $           -2*mne(i)*mne(j)*myB0NM(p**2,mne(i)**2,mne(j)**2,Q**2))

         enddo
      enddo
      
      do i = 1,2                ! charginos
         do j = 1,2
            
            PiZZ = PiZZ
     $           +(azchch(i,j)**2+bzchch(i,j)**2)
     $           *myH(p**2,mch(i)**2,mch(j)**2,Q**2)
     $           +4*azchch(i,j)*bzchch(i,j)*
     $           mch(i)*mch(j)*myB0NM(p**2,mch(i)**2,mch(j)**2,Q**2)
            
         enddo
      enddo
      PiZZ = PiZZ/16d0/pi**2

      return
      end
      
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     

      subroutine getPiWW(g,gp,ht,hb,htau,v1,v2,p,Q,piWW)
      
      implicit none

      double precision g,gp,ht,hb,htau,v1,v2,p,Q,piWW
      
      double precision mstop(2),msbot(2),mstau(2),msup(2),msdown(2),
     $     msel(2),msnutau,msnue,Rt(2,2),Rb(2,2),Rtau(2,2),del(2,2),
     $     NN(5,5),UU(2,2),VV(2,2),mch(2),mne(5),
     $     mhh(3),maa(3),mhc(2),RS(3,3),RP(3,3),RC(2,2),
     $     myB0NM,myB22T,myH,pi,cb,sb,mw2,mz2,mt2,mb2,mtau2,gb2,
     $     sw2,cw2,awnech(5,2),bwnech(5,2)

      integer i,j
      
      common/TREE_MASSES/mhh,maa,mhc,RS,RP,RC,mch,UU,VV,mne,NN,
     $     mstop,msbot,mstau,Rt,Rb,Rtau,msup,msdown,msel,msnutau,msnue
      
      pi = 4d0*atan(1.d0)

      PiWW = 0d0

      cb = v1/sqrt(v1**2+v2**2)
      sb = v2/sqrt(v1**2+v2**2)
      gb2 = (g**2+gp**2)/2d0
      sw2 = gp**2/2d0/gb2
      cw2 = 1-sw2

      mw2 = g**2/2d0*(v1**2+v2**2)
      mz2 = gb2*(v1**2+v2**2)

c     the higgs and gauge contributions

      do i = 1,3                

         PiWW = PiWW            ! scalar-W
     $        + g**2*mw2*(cb*RS(i,1)+sb*RS(i,2))**2
     $        *myB0NM(p**2,mhh(i)**2,mw2,Q**2)
         
         do j = 1,2

            PiWW = PiWW         ! scalar-charged
     $           - g**2*(RS(i,1)*RC(j,1)-RS(i,2)*RC(j,2))**2
     $           *myB22T(p**2,mhh(i)**2,mhc(j)**2,Q**2)

            PiWW = PiWW         ! pseudo-charged
     $           - g**2*(RP(i,1)*RC(j,1)+RP(i,2)*RC(j,2))**2
     $           *myB22T(p**2,maa(i)**2,mhc(j)**2,Q**2)

         enddo
      enddo

      PiWW = PiWW               ! pure gauge
     $     -8*g**2*cw2*myB22T(p**2,mw2,mz2,Q**2)
     $     -g**2*((4*p**2+mw2+mz2)*cw2-mz2*sw2**2)
     $     *myB0NM(p**2,mz2,mw2,Q**2)
     $     -sw2*g**2*(8*myB22T(p**2,mw2,0d0,Q**2)
     $     +4*p**2*myB0NM(p**2,mw2,0d0,Q**2))

c     the fermion contributions

      mt2 = ht**2*v2**2
      mb2 = hb**2*v1**2
      mtau2 = htau**2*v1**2
      
      PiWW = PiWW + g**2/2d0*(
     $     8*myH(p**2,0d0,0d0,Q**2)
     $     +3*myH(p**2,mt2,mb2,Q**2)
     $     +myH(p**2,mtau2,0d0,Q**2))

c     the sfermion contributions

      del(1,1) = 1d0
      del(1,2) = 0d0
      del(2,1) = 0d0
      del(2,2) = 1d0

      do i=1,2
         do j = 1,2

            PiWW = PiWW - 2*g**2*(
     $           6*(del(i,1)*del(j,1))**2
     $           *myB22T(p**2,msup(i)**2,msdown(j)**2,Q**2)
     $           +3*(Rt(i,1)*Rb(j,1))**2
     $           *myB22T(p**2,mstop(i)**2,msbot(j)**2,Q**2))
            
         enddo
         
         PiWW = PiWW - 2*g**2*(
     $        2*del(i,1)**2*myB22T(p**2,msel(i)**2,msnue**2,Q**2)
     $        +Rtau(i,1)**2*myB22T(p**2,mstau(i)**2,msnutau**2,Q**2))

      enddo

c     the chargino/neutralino contribution

      call coupl_W_ino(g,NN,UU,VV,awnech,bwnech)

      do i = 1,5
         do j = 1,2
            
            PiWW = PiWW
     $           +(awnech(i,j)**2+bwnech(i,j)**2)
     $           *myH(p**2,mne(i)**2,mch(j)**2,Q**2)
     $           +4*awnech(i,j)*bwnech(i,j)*
     $           mne(i)*mch(j)*myB0NM(p**2,mne(i)**2,mch(j)**2,Q**2)
            
         enddo
      enddo

      PiWW = PiWW/16d0/pi**2

      return
      end
      
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     

      subroutine treemasses(g,gp,ll,kk,ht,hb,htau,v1,v2,xx,M1,M2,
     $     Ak,Al,At,Ab,Atau,mQ3,mtr,mbr,mQ,mur,mdr,mL3,mtaur,mL,mer,
     $     Q,errmass)

      double precision g,gp,ll,kk,ht,hb,htau,v1,v2,xx,M1,M2,
     $     Ak,Al,At,Ab,Atau,mQ3,mtr,mbr,mQ,mur,mdr,mL3,mtaur,mL,mer,Q

      logical errmass
      
      double precision mstop(2),msbot(2),mstau(2),msup(2),msdown(2),
     $     msel(2),msnutau,msnue,Rt(2,2),Rb(2,2),Rtau(2,2),
     $     NN(5,5),UU(2,2),VV(2,2),mch(2),mne(5),
     $     mhh(3),maa(3),mhc(2),RS(3,3),RP(3,3),RC(2,2)

      logical errsfer,errhiggs

      common/TREE_MASSES/mhh,maa,mhc,RS,RP,RC,mch,UU,VV,mne,NN,
     $     mstop,msbot,mstau,Rt,Rb,Rtau,msup,msdown,msel,msnutau,msnue

      common/FOREFFPOT/mt,T1,T2,st,ct,Q2
      double precision mt,T1,T2,st,ct,Q2

c     compute all the tree-level masses

      call tree_sfermions(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     $     mQ3,mtr,mbr,mQ,mur,mdr,mL3,mtaur,mL,mer,mstop,msbot,
     $     mstau,Rt,Rb,Rtau,msup,msdown,msel,msnutau,msnue,errsfer)

      mt = ht*v2
      T1 = mstop(1)**2          ! for the Higgs mass calculation
      T2 = mstop(2)**2
      st = Rt(1,2)              ! always true???
      ct = Rt(1,1)
      Q2 = Q**2
    

      call tree_higgses(g,gp,ll,kk,v1,v2,xx,Ak,Al,mhh,maa,mhc,RS,RP,RC,
     $     errhiggs)
      
      
      call tree_charginos(g,ll,v1,v2,M2,xx,mch,UU,VV)

      call tree_neutralinos(g,gp,ll,kk,v1,v2,xx,M1,M2,mne,NN)

      errmass = errhiggs.or.errsfer

c$$$      write(*,*) 'scalar',mhh
c$$$      write(*,*) RS(1,1),RS(1,2),RS(1,3)
c$$$      write(*,*) RS(2,1),RS(2,2),RS(2,3)
c$$$      write(*,*) RS(3,1),RS(3,2),RS(3,3)
c$$$
c$$$      write(*,*) 'pseudo',maa
c$$$      write(*,*) RP(1,1),RP(1,2),RP(1,3)
c$$$      write(*,*) RP(2,1),RP(2,2),RP(2,3)
c$$$      write(*,*) RP(3,1),RP(3,2),RP(3,3)
c$$$
c$$$      write(*,*) 'charged',mhc
c$$$      
c$$$      write(*,*) 'sup',msup,mstop
c$$$      write(*,*) Rt(1,1),Rt(1,2)
c$$$      write(*,*) Rt(2,1),Rt(2,2)
c$$$
c$$$      write(*,*) 'sdown',msdown,msbot
c$$$      write(*,*) Rb(1,1),Rb(1,2)
c$$$      write(*,*) Rb(2,1),Rb(2,2)
c$$$
c$$$      write(*,*) 'slep',msel,mstau
c$$$      write(*,*) Rtau(1,1),Rtau(1,2)
c$$$      write(*,*) Rtau(2,1),Rtau(2,2)
c$$$
c$$$      write(*,*) 'sneut',msnue,msnutau
c$$$
c$$$      write(*,*) 'charg',mch
c$$$      write(*,*) UU(1,1),UU(1,2),VV(1,1),VV(1,2)
c$$$      write(*,*) UU(2,1),UU(2,2),VV(2,1),VV(2,2)
c$$$
c$$$      write(*,*) 'neut',mne
c$$$      write(*,*) NN(1,1),NN(1,2),NN(1,3),NN(1,4),NN(1,5)
c$$$      write(*,*) NN(2,1),NN(2,2),NN(2,3),NN(2,4),NN(2,5)
c$$$      write(*,*) NN(3,1),NN(3,2),NN(3,3),NN(3,4),NN(3,5)
c$$$      write(*,*) NN(4,1),NN(4,2),NN(4,3),NN(4,4),NN(4,5)
c$$$      write(*,*) NN(5,1),NN(5,2),NN(5,3),NN(5,4),NN(5,5)

      return
      end
      
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     

      subroutine tree_charginos(g,ll,v1,v2,M2,xx,xmc,u,v)

      implicit none

      double precision g,ll,v1,v2,M2,xx,xmc(2),u(2,2),v(2,2)
      double precision mc(2,2),xtx(2,2),xxt(2,2),ut(2,2),vt(2,2),
     $     bu(2),bv(2),umc(2),vmc(2),mtemp

      integer i,j,k,nrot
      
      mc(1,1) = M2
      mc(1,2) = g*v2
      mc(2,1) = g*v1
      mc(2,2) = ll*xx

      do i=1,2
         do j=1,2
            xxt(i,j) = 0d0
            xtx(i,j) = 0d0
         enddo
      enddo
      
      do i=1,2
         do j=1,2
            do k=1,2
               xxt(i,j) = xxt(i,j) + mc(i,k)*mc(j,k)
               xtx(i,j) = xtx(i,j) + mc(k,i)*mc(k,j)
            enddo
         enddo
      enddo

      call jacobi(xxt,2,2,umc,ut,nrot)
      call jacobi(xtx,2,2,vmc,vt,nrot)

      if(abs(vmc(1)-umc(1)).gt.1d-6) then ! swap eigenstates 
         do j=1,2               
            bv(j)=vt(j,1)
            vt(j,1)=vt(j,2)
            vt(j,2)=bv(j)
         enddo
      endif
      
      do i=1,2
         do j=1,2
            u(i,j) = ut(j,i)
            v(i,j) = vt(j,i)
         enddo
      enddo
      
      xmc(1) = 0d0
      xmc(2) = 0d0

      do i=1,2
         do j=1,2
            xmc(1) = xmc(1) + u(1,i)*mc(i,j)*v(1,j)
            xmc(2) = xmc(2) + u(2,i)*mc(i,j)*v(2,j)
         enddo
      enddo

c     order the eigenstates
      
      if(abs(xmc(1)).gt.abs(xmc(2))) then 
         mtemp=xmc(1)
         xmc(1)=xmc(2)
         xmc(2)=mtemp
         do j=1,2
            bu(j)=u(1,j)
            u(1,j)=u(2,j)
            u(2,j)=bu(j)
            bv(j)=v(1,j)
            v(1,j)=v(2,j)
            v(2,j)=bv(j)
         enddo
      endif

      return 
      end

c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     

      subroutine tree_neutralinos(g,gp,ll,kk,v1,v2,xx,M1,M2,xmn,Z)

      implicit none
      
      double precision g,gp,ll,kk,v1,v2,xx,M1,M2,xmn(5),Z(5,5)

      double precision xn(5,5),mn(5),ymn(5),xx0,xx1,sq2,zx(5,5)
      
      integer i,i1,j,idummy,iord(5),irem(5),nrot
      

      sq2 = sqrt(2d0)

      xn(1,1)=M1
      xn(1,2)=0d0
      xn(1,3)=-gp*v1/sq2
      xn(1,4)=gp*v2/sq2
      xn(1,5)=0d0
      xn(2,1)=0d0
      xn(2,2)=M2
      xn(2,3)=g*v1/sq2
      xn(2,4)=-g*v2/sq2
      xn(2,5)=0d0
      xn(3,1)=xn(1,3)
      xn(3,2)=xn(2,3)
      xn(3,3)=0d0
      xn(3,4)=-ll*xx
      xn(3,5)=-ll*v2
      xn(4,1)=xn(1,4)
      xn(4,2)=xn(2,4)
      xn(4,3)=xn(3,4)
      xn(4,4)=0d0
      xn(4,5)=-ll*v1
      xn(5,1)=xn(1,5)
      xn(5,2)=xn(2,5)
      xn(5,3)=xn(3,5)
      xn(5,4)=xn(4,5)
      xn(5,5)=2*kk*xx
      
      call jacobi(xn,5,5,ymn,zx,nrot)
      
c     ordering the disorder 
      
      do i=1,5
         mn(i) = abs(ymn(i))
      enddo
      
      xx0 = dmin1(mn(1),mn(2),mn(3),mn(4),mn(5))
      xx1 = dmax1(mn(1),mn(2),mn(3),mn(4),mn(5))
      idummy = 1
      do i = 1,5
         if(mn(i).eq.xx0)then
            iord(1) = i
         elseif(mn(i).eq.xx1)then
            iord(5) = i
         else
            irem(idummy) = i
            idummy = idummy+1
         endif
      enddo

      xx0 = dmin1(mn(irem(1)),mn(irem(2)),mn(irem(3)))
      xx1 = dmax1(mn(irem(1)),mn(irem(2)),mn(irem(3)))

      do i = 1,3
         if(mn(irem(i)).eq.xx0)then
            iord(2) = irem(i)
         elseif(mn(irem(i)).eq.xx1)then
            iord(4) = irem(i)
         else
            iord(3) = irem(i)
         endif
      enddo
c     
      do j=1,5
         i=iord(j)
         xmn(j)=ymn(i)
         do i1=1,5
            z(j,i1)=zx(i1,i)    ! note that ZX ~ Z^T
         enddo
      enddo
      
      return
      end
     
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     

      subroutine tree_higgses(g,gp,ll,kk,v1,v2,xx,Ak,Al,
     $     mss,maa,mhc,RS,RP,RC,errhiggs)

      implicit none

      double precision g,gp,ll,kk,v1,v2,xx,Ak,Al,
     $     mss(3),maa(3),mhc(2),RS(3,3),RP(3,3),RC(2,2)

      logical errhiggs

      double precision gb2,As,MS(3,3),MP(3,3),ZS(3,3),ZP(3,3),
     $     ms2(3),msd(3),ma2(3),mad(3)

      double precision cb,sb,gold

      integer i,j,k,nrot

      common/FOREFFPOT/mt,T1,T2,st,ct,Q2
      double precision mt,T1,T2,st,ct,Q2,DMS(3,3),DMP(3,3)

      gb2 = (g**2+gp**2)/2d0
      As = Al+kk*xx
            
      MS(1,1) = gb2*v1**2 + ll*xx*v2/v1*As
      MS(1,2) = (2*ll**2-gb2)*v1*v2-ll*xx*As
      MS(1,3) = 2*ll**2*v1*xx-ll*v2*(As+kk*xx)
      MS(2,2) = gb2*v2**2 + ll*xx*v1/v2*As
      MS(2,3) = 2*ll**2*v2*xx-ll*v1*(As+kk*xx)
      MS(3,3) = ll*Al*v1*v2/xx+kk*xx*(Ak+4*kk*xx)
      MS(2,1) = MS(1,2)
      MS(3,1) = MS(1,3)
      MS(3,2) = MS(2,3)

  
      gold = gb2*(v1**2+v2**2)  ! gauge-fixing mass
      sb = v2/sqrt(v1**2+v2**2)
      cb = v1/sqrt(v1**2+v2**2)
c$$$      gold = 0.0

     

      MP(1,1) = ll*xx*v2/v1*As + cb**2*gold
      MP(1,2) = ll*xx*As - cb*sb*gold
      MP(1,3) = ll*v2*(As-3*kk*xx)
      MP(2,2) = ll*xx*v1/v2*As + sb**2*gold
      MP(2,3) = ll*v1*(As-3*kk*xx)
      MP(3,3) = 4*ll*kk*v1*v2+ll*Al*v1*v2/xx-3*kk*Ak*xx
      MP(2,1) = MP(1,2)
      MP(3,1) = MP(1,3)
      MP(3,2) = MP(2,3)

      call jacobi(MS,3,3,ms2,ZS,nrot)
      call jacobi(MP,3,3,ma2,ZP,nrot)
      
c     take the square roots and check that it's all right

      errhiggs = .false.

      do i = 1,3         
         if(ms2(i).ge.0d0) then 
            msd(i) = sqrt(ms2(i))
         else
            errhiggs = .true.
            write(*,*) 'tachionic scalar',i,' at tree level'
         endif
      enddo

      do i = 1,3         
         if(ma2(i).ge.-1d-12) then ! allow for tiny nonzero Goldstone mass 
            mad(i) = sqrt(abs(ma2(i)))
         else
            errhiggs = .true.
            write(*,*) 'tachionic pseudoscalar',i,' at tree level'
         endif
      enddo

c     if there was a tachionic mass, try again with corrected matrices

c$$$      if(errhiggs) then
c$$$         
c$$$         call effpot(1,mt,0d0,T1,T2,st,ct,Q2,v2/v1,
c$$$     $        sqrt(v1**2+v2**2),ll,xx,0d0,DMS,DMP) 
c$$$
c$$$         do i=1,3
c$$$            do j=1,3               
c$$$               MS(i,j) = MS(i,j) + DMS(i,j)
c$$$               MP(i,j) = MP(i,j) + DMP(i,j)
c$$$            enddo
c$$$         enddo
c$$$
c$$$         call jacobi(MS,3,3,ms2,ZS,nrot)
c$$$         call jacobi(MP,3,3,ma2,ZP,nrot)
c$$$         
c$$$         errhiggs = .false.
c$$$
c$$$         do i = 1,3         
c$$$            if(ms2(i).ge.0d0) then 
c$$$               msd(i) = sqrt(ms2(i))
c$$$            else
c$$$               errhiggs = .true.
c$$$               write(*,*) 'tachionic scalar',i,' at one loop'
c$$$               msd(i) = -sqrt(abs(ms2(i)))
c$$$            endif
c$$$         enddo
c$$$         
c$$$         do i = 1,3         
c$$$            if(ma2(i).ge.-1d-12) then ! allow for tiny nonzero Goldstone mass 
c$$$               mad(i) = sqrt(abs(ma2(i)))
c$$$            else
c$$$               errhiggs = .true.
c$$$               write(*,*) 'tachionic pseudoscalar',i,' at one loop'
c$$$               mad(i) = -sqrt(abs(ma2(i)))
c$$$            endif
c$$$         enddo
c$$$
c$$$      endif

c     order the disorder

      mss(1) = dmin1(msd(1),msd(2),msd(3))
      mss(3) = dmax1(msd(1),msd(2),msd(3))

      do i=1,3         
         if(msd(i).gt.mss(1).and.msd(i).lt.mss(3)) then
            mss(2) = msd(i)
         endif         
      enddo

      do i = 1,3
         do j = 1,3
            if(mss(i).eq.msd(j)) then
               do k = 1,3                  
                  RS(i,k) = ZS(k,j)
               enddo
            endif
         enddo
      enddo

      maa(1) = dmin1(mad(1),mad(2),mad(3))
      maa(3) = dmax1(mad(1),mad(2),mad(3))

      do i=1,3         
         if(mad(i).gt.maa(1).and.mad(i).lt.maa(3)) then
            maa(2) = mad(i)
         endif         
      enddo

      do i = 1,3
         do j = 1,3
            if(maa(i).eq.mad(j)) then
               do k = 1,3                  
                  RP(i,k) = ZP(k,j)
               enddo
            endif
         enddo
      enddo

c$$$c     add gauge-fixing mass for G0 (improve later)
c$$$
c$$$      maa(1) = sqrt(gb2*(v1**2+v2**2))
 
c     the charged Higgses

      mhc(1) = sqrt(g**2*(v1**2+v2**2)/2d0)
      
      mhc(2) = sqrt((ll*xx*As-ll**2*v1*v2)*(v1**2+v2**2)/v1/v2
     $     + g**2*(v1**2+v2**2)/2d0)

      RC(1,1) = -v1/sqrt(v1**2+v2**2)
      RC(1,2) = v2/sqrt(v1**2+v2**2)
      RC(2,1) = RC(1,2)
      RC(2,2) = v1/sqrt(v1**2+v2**2)

      return
      end
     
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     

      subroutine tree_sfermions(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     $     mQ3,mtr,mbr,mQ,mur,mdr,mL3,mtaur,mL,mer,mstop,msbot,
     $     mstau,Rt,Rb,Rtau,msup,msdown,msel,msnutau,msnue,errsfer)

      implicit none

      double precision g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     $     mQ3,mtr,mbr,mQ,mur,mdr,mL3,mtaur,mL,mer,mstop(2),msbot(2),
     $     mstau(2),Rt(2,2),Rb(2,2),Rtau(2,2),msup(2),msdown(2),msel(2),
     $     msnutau,msnue

      logical errsfer,errt,errb,errta

      double precision msuL,msuR,msdL,msdR,mslL,mslR

c     first the ones that mix

      call diagsfe(1,ht,g,gp,ll,v1,v2,xx,At,mQ3,mtr,mstop,Rt,errt)
      call diagsfe(2,hb,g,gp,ll,v1,v2,xx,Ab,mQ3,mbr,msbot,Rb,errb)
      call diagsfe(3,htau,g,gp,ll,v1,v2,xx,Atau,mL3,mtaur,mstau,Rtau,
     $     errta)
      
      errsfer = errt.or.errb.or.errta

c     then the others

      msuL = sqrt(mQ**2 + ( g**2/2d0 - gp**2/6d0)*(v1**2-v2**2)/2d0)
      msdL = sqrt(mQ**2 + (-g**2/2d0 - gp**2/6d0)*(v1**2-v2**2)/2d0)
      mslL = sqrt(mL**2 + (-g**2/2d0 + gp**2/2d0)*(v1**2-v2**2)/2d0)

      msuR = sqrt(mur**2 + ( 2*gp**2/3d0)*(v1**2-v2**2)/2d0)
      msdR = sqrt(mdr**2 + (  -gp**2/3d0)*(v1**2-v2**2)/2d0)
      mslR = sqrt(mer**2 + (  -gp**2)*(v1**2-v2**2)/2d0)

      msnutau = sqrt(mL3**2 + (g**2/2d0 + gp**2/2d0)*(v1**2-v2**2)/2d0)
      msnue   = sqrt(mL**2 + (g**2/2d0 + gp**2/2d0)*(v1**2-v2**2)/2d0)

      msup(1) = msuL
      msup(2) = msuR
      msdown(1) = msdL
      msdown(2) = msdR
      msel(1) = mslL
      msel(2) = mslR

      return
      end

      
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     

      subroutine diagsfe(n,hf,g,gp,ll,v1,v2,xx,Af,mL,mR,mass,R,error)

      implicit none

      integer n
      double precision hf,g,gp,ll,v1,v2,xx,Af,mL,mR,mass(2),R(2,2)
      logical error

      double precision dL,dR,X,mf,MS(2,2),RT(2,2),mass2(2)
      integer nrot,i,j

      if(n.eq.1) then
         dL = ( g**2/2d0 - gp**2/6d0)*(v1**2-v2**2)/2d0
         dR = ( 2*gp**2/3d0)*(v1**2-v2**2)/2d0
         mf = hf*v2
         X = Af-ll*xx*v1/v2
      elseif(n.eq.2) then
         dL = (-g**2/2d0 - gp**2/6d0)*(v1**2-v2**2)/2d0
         dR = (  -gp**2/3d0)*(v1**2-v2**2)/2d0
         mf = hf*v1
         X = Af-ll*xx*v2/v1
      elseif(n.eq.3) then
         dL = (-g**2/2d0 + gp**2/2d0)*(v1**2-v2**2)/2d0
         dR = (  -gp**2)*(v1**2-v2**2)/2d0
         mf = hf*v1
         X = Af-ll*xx*v2/v1
      endif

      MS(1,1) = mL**2 + mf**2 + dL
      MS(1,2) = mf*X
      MS(2,1) = MS(1,2)
      MS(2,2) = mR**2 + mf**2 + dR

      call jacobi(MS,2,2,mass2,Rt,nrot)

      error = .false.

      if(mass2(1).ge.0d0) then
         mass(1) = sqrt(mass2(1))
      else
         error = .true.
         mass(1) = -sqrt(abs(mass2(1)))
      endif

      if(mass2(2).ge.0d0) then
         mass(2) = sqrt(mass2(2))
      else
         error = .true.
         mass(2) = -sqrt(abs(mass2(2)))
      endif

      do i=1,2
         do j=1,2
            R(i,j) = Rt(j,i)
         enddo
      enddo

      return
      end


      
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     

      subroutine coupl_s_sf(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     $     Rt,Rb,Rtau,lsstt,lssbb,lsstata,lssntnt,lssuu,lssdd,
     $     lssee,lssnn,lstt,lsbb,lstata,lsntnt,lsuu,lsdd,lsee,lsnn)

      implicit none

      double precision g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     $     Rt(2,2),Rb(2,2),Rtau(2,2),lsstt(3,3,2,2),lssbb(3,3,2,2),
     $     lsstata(3,3,2,2),lssntnt(3,3),lssuu(3,3,2,2),lssdd(3,3,2,2),
     $     lssee(3,3,2,2),lssnn(3,3),lstt(3,2,2),lsbb(3,2,2),
     $     lstata(3,2,2),lsntnt(3),lsuu(3,2,2),lsdd(3,2,2),lsee(3,2,2),
     $     lsnn(3),lsstt_int(3,3,2,2),lssbb_int(3,3,2,2),
     $     lsstata_int(3,3,2,2),lstt_int(3,2,2),lsbb_int(3,2,2),
     $     lstata_int(3,2,2)

      double precision sq2,gb2,sw2,guL,guR,gdL,gdR,geL,geR,gnu

      integer i,j,k,l,a,b

      sq2 = sqrt(2d0)
      gb2 = (g**2+gp**2)/2d0
      sw2 = gp**2/2d0/gb2

      guL =  .5d0 - 2d0/3d0*sw2
      gdL = -.5d0 + 1d0/3d0*sw2
      gnu =  .5d0
      geL = -.5d0 + sw2
      guR =  2d0/3d0*sw2
      gdR = -1d0/3d0*sw2
      geR = -sw2

c     FIRST TWO GENERATIONS

c     quartic, up squarks

      do i=1,3
         do j=1,3
            do k=1,2
               do l=1,2
                  lssuu(i,j,k,l) = 0d0
               enddo
            enddo
         enddo
      enddo

      lssuu(1,1,1,1) = gb2*guL/2d0
      lssuu(1,1,2,2) = gb2*guR/2d0
      lssuu(2,2,1,1) = -gb2*guL/2d0
      lssuu(2,2,2,2) = -gb2*guR/2d0
      
c     trilinear, up squarks

      do i=1,3
         do k=1,2
            do l=1,2
               lsuu(i,k,l) = 0d0
            enddo
         enddo
      enddo

      lsuu(1,1,1) = sq2*gb2*guL*v1
      lsuu(1,2,2) = sq2*gb2*guR*v1
      lsuu(2,1,1) = -sq2*gb2*guL*v2
      lsuu(2,2,2) = -sq2*gb2*guR*v2
      
c     quartic, down squarks

      do i=1,3
         do j=1,3
            do k=1,2
               do l=1,2
                  lssdd(i,j,k,l) = 0d0
               enddo
            enddo
         enddo
      enddo

      lssdd(1,1,1,1) = gb2*gdL/2d0
      lssdd(1,1,2,2) = gb2*gdR/2d0
      lssdd(2,2,1,1) = -gb2*gdL/2d0
      lssdd(2,2,2,2) = -gb2*gdR/2d0
      
c     trilinear, down squarks

      do i=1,3
         do k=1,2
            do l=1,2
               lsdd(i,k,l) = 0d0
            enddo
         enddo
      enddo

      lsdd(1,1,1) = sq2*gb2*gdL*v1
      lsdd(1,2,2) = sq2*gb2*gdR*v1
      lsdd(2,1,1) = -sq2*gb2*gdL*v2
      lsdd(2,2,2) = -sq2*gb2*gdR*v2

c     quartic, charged sleptons

      do i=1,3
         do j=1,3
            do k=1,2
               do l=1,2
                  lssee(i,j,k,l) = 0d0
               enddo
            enddo
         enddo
      enddo

      lssee(1,1,1,1) = gb2*geL/2d0
      lssee(1,1,2,2) = gb2*geR/2d0
      lssee(2,2,1,1) = -gb2*geL/2d0
      lssee(2,2,2,2) = -gb2*geR/2d0

c     trilinear, charged sleptons

      do i=1,3
         do k=1,2
            do l=1,2
               lsee(i,k,l) = 0d0
            enddo
         enddo
      enddo

      lsee(1,1,1) = sq2*gb2*geL*v1
      lsee(1,2,2) = sq2*gb2*geR*v1
      lsee(2,1,1) = -sq2*gb2*geL*v2
      lsee(2,2,2) = -sq2*gb2*geR*v2

c     quartic, sneutrinos

      do i=1,3
         do j=1,3
            lssnn(i,j) = 0d0
         enddo
      enddo

      lssnn(1,1) = gb2*gnu/2d0
      lssnn(2,2) = -gb2*gnu/2d0
      
c     trilinear, sneutrinos

      lsnn(1) = sq2*gb2*gnu*v1
      lsnn(2) = -sq2*gb2*gnu*v2
      lsnn(3) = 0d0

c     THIRD GENERATION

c     quartic, top squarks

      do i=1,3
         do j=1,3
            do k=1,2
               do l=1,2
                  lsstt_int(i,j,k,l) = lssuu(i,j,k,l)
               enddo
            enddo
         enddo
      enddo

      lsstt_int(2,2,1,1) = lsstt_int(2,2,1,1) + ht**2/2d0
      lsstt_int(2,2,2,2) = lsstt_int(2,2,2,2) + ht**2/2d0

      lsstt_int(1,3,1,2) = lsstt_int(1,3,1,2) - ht*ll/4d0
      lsstt_int(1,3,2,1) = lsstt_int(1,3,1,2)
      lsstt_int(3,1,1,2) = lsstt_int(1,3,1,2)
      lsstt_int(3,1,2,1) = lsstt_int(1,3,1,2)

c     trilinear, top squarks

      do i=1,3
         do k=1,2
            do l=1,2
               lstt_int(i,k,l) = lsuu(i,k,l)
            enddo
         enddo
      enddo

      lstt_int(2,1,1) = lstt_int(2,1,1) + sq2*ht**2*v2
      lstt_int(2,2,2) = lstt_int(2,2,2) + sq2*ht**2*v2

      lstt_int(1,1,2) = lstt_int(1,1,2) - ht*ll*xx/sq2
      lstt_int(2,1,2) = lstt_int(2,1,2) + ht*At/sq2
      lstt_int(3,1,2) = lstt_int(3,1,2) - ht*ll*v1/sq2

      lstt_int(1,2,1) = lstt_int(1,1,2) 
      lstt_int(2,2,1) = lstt_int(2,1,2) 
      lstt_int(3,2,1) = lstt_int(3,1,2) 

c     quartic, bottom squarks

      do i=1,3
         do j=1,3
            do k=1,2
               do l=1,2
                  lssbb_int(i,j,k,l) = lssdd(i,j,k,l)
               enddo
            enddo
         enddo
      enddo

      lssbb_int(1,1,1,1) = lssbb_int(1,1,1,1) + hb**2/2d0
      lssbb_int(1,1,2,2) = lssbb_int(1,1,2,2) + hb**2/2d0

      lssbb_int(2,3,1,2) = lssbb_int(2,3,1,2) - hb*ll/4d0
      lssbb_int(2,3,2,1) = lssbb_int(2,3,1,2)
      lssbb_int(3,2,1,2) = lssbb_int(2,3,1,2)
      lssbb_int(3,2,2,1) = lssbb_int(2,3,1,2)

c     trilinear, bottom squarks

      do i=1,3
         do k=1,2
            do l=1,2
               lsbb_int(i,k,l) = lsdd(i,k,l)
            enddo
         enddo
      enddo

      lsbb_int(1,1,1) = lsbb_int(1,1,1) + sq2*hb**2*v1
      lsbb_int(1,2,2) = lsbb_int(1,2,2) + sq2*hb**2*v1

      lsbb_int(1,1,2) = lsbb_int(1,1,2) + hb*Ab/sq2
      lsbb_int(2,1,2) = lsbb_int(2,1,2) - hb*ll*xx/sq2
      lsbb_int(3,1,2) = lsbb_int(3,1,2) - hb*ll*v2/sq2

      lsbb_int(1,2,1) = lsbb_int(1,1,2) 
      lsbb_int(2,2,1) = lsbb_int(2,1,2) 
      lsbb_int(3,2,1) = lsbb_int(3,1,2) 

c     quartic, staus

      do i=1,3
         do j=1,3
            do k=1,2
               do l=1,2
                   lsstata_int(i,j,k,l) = lssee(i,j,k,l)
               enddo
            enddo
         enddo
      enddo

      lsstata_int(1,1,1,1) = lsstata_int(1,1,1,1) + htau**2/2d0
      lsstata_int(1,1,2,2) = lsstata_int(1,1,2,2) + htau**2/2d0

      lsstata_int(2,3,1,2) = lsstata_int(2,3,1,2) - htau*ll/4d0
      lsstata_int(2,3,2,1) = lsstata_int(2,3,1,2)
      lsstata_int(3,2,1,2) = lsstata_int(2,3,1,2)
      lsstata_int(3,2,2,1) = lsstata_int(2,3,1,2)
      
c     trilinear, staus

      do i=1,3
         do k=1,2
            do l=1,2
               lstata_int(i,k,l) = lsee(i,k,l) 
            enddo
         enddo
      enddo

      lstata_int(1,1,1) = lstata_int(1,1,1) + sq2*htau**2*v1
      lstata_int(1,2,2) = lstata_int(1,2,2) + sq2*htau**2*v1

      lstata_int(1,1,2) = lstata_int(1,1,2) + htau*Atau/sq2
      lstata_int(2,1,2) = lstata_int(2,1,2) - htau*ll*xx/sq2
      lstata_int(3,1,2) = lstata_int(3,1,2) - htau*ll*v2/sq2

      lstata_int(1,2,1) = lstata_int(1,1,2) 
      lstata_int(2,2,1) = lstata_int(2,1,2) 
      lstata_int(3,2,1) = lstata_int(3,1,2) 

c     sneutrinos (same as first two generations)

      do i=1,3
         lsntnt(i) = lsnn(i)
         do j=1,3
            lssntnt(i,j) = lssnn(i,j)
         enddo
      enddo

c     now rotate the third-generation couplings

c     quartic
      do i = 1,3
         do j = 1,3
            do k=1,2
               do l=1,2
                  lsstt(i,j,k,l) = 0d0
                  lssbb(i,j,k,l) = 0d0
                  lsstata(i,j,k,l) = 0d0
                  do a=1,2
                     do b=1,2
                        lsstt(i,j,k,l) = lsstt(i,j,k,l)
     $                       +Rt(k,a)*Rt(l,b)*lsstt_int(i,j,a,b)
                        lssbb(i,j,k,l) = lssbb(i,j,k,l)
     $                       +Rb(k,a)*Rb(l,b)*lssbb_int(i,j,a,b)
                        lsstata(i,j,k,l) = lsstata(i,j,k,l)
     $                       +Rtau(k,a)*Rtau(l,b)*lsstata_int(i,j,a,b)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
           
c     trilinear
      do i = 1,3         
         do k=1,2
            do l=1,2
               lstt(i,k,l) = 0d0
               lsbb(i,k,l) = 0d0
               lstata(i,k,l) = 0d0
               do a=1,2
                  do b=1,2
                     lstt(i,k,l) = lstt(i,k,l)
     $                    +Rt(k,a)*Rt(l,b)*lstt_int(i,a,b)
                     lsbb(i,k,l) = lsbb(i,k,l)
     $                    +Rb(k,a)*Rb(l,b)*lsbb_int(i,a,b)
                     lstata(i,k,l) = lstata(i,k,l)
     $                    +Rtau(k,a)*Rtau(l,b)*lstata_int(i,a,b)
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
      end

c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     

      subroutine coupl_s_hh(g,gp,ll,kk,v1,v2,xx,Al,Ak,RS,RP,
     $     lsshh,lssaa,lshh,lsaa,lsscc,lscc)

      implicit none

      double precision g,gp,ll,kk,v1,v2,xx,Al,Ak,RS(3,3),RP(3,3),
     $     lsshh(3,3,3,3),lssaa(3,3,3,3),lshh(3,3,3),lsaa(3,3,3),
     $     lsscc(3,3,2,2),lscc(3,2,2)

      double precision lssss(3,3,3,3),lsspp(3,3,3,3),
     $     lsss(3,3,3),lspp(3,3,3),gb2,sq2,c2b,s2b

      integer i,j,k,l,a,b
      
      gb2 = (g**2+gp**2)/2d0
      sq2 = sqrt(2d0)
      
      do i=1,3                  ! initialize
         do j=1,3
            do k=1,3
               lsss(i,j,k)=0d0
               lspp(i,j,k)=0d0
               lshh(i,j,k)=0d0
               lsaa(i,j,k)=0d0
               lscc(i,j,k)=0d0               
               do l=1,3
                  lssss(i,j,k,l)=0d0
                  lsspp(i,j,k,l)=0d0
                  lsshh(i,j,k,l)=0d0
                  lssaa(i,j,k,l)=0d0
                  lsscc(i,j,k,l)=0d0
               enddo
            enddo
         enddo
      enddo

c     quartic neutral couplings
 
      lssss(1,1,1,1) = gb2/16d0

      lssss(2,2,2,2) = gb2/16d0

      lssss(3,3,3,3) = kk**2/4d0

      lssss(1,1,2,2) = (2*ll**2-gb2)/48d0
      lssss(1,2,1,2) = lssss(1,1,2,2) ! symmetrize the couplings
      lssss(1,2,2,1) = lssss(1,1,2,2) ! (got a better proposal? ;-)
      lssss(2,1,1,2) = lssss(1,1,2,2)
      lssss(2,1,2,1) = lssss(1,1,2,2)
      lssss(2,2,1,1) = lssss(1,1,2,2)

      lssss(1,2,3,3) = -ll*kk/24d0
      lssss(2,1,3,3) = lssss(1,2,3,3)
      lssss(2,3,1,3) = lssss(1,2,3,3)
      lssss(1,3,2,3) = lssss(1,2,3,3)
      lssss(3,1,2,3) = lssss(1,2,3,3)
      lssss(3,2,1,3) = lssss(1,2,3,3)
      lssss(3,1,3,2) = lssss(1,2,3,3)
      lssss(3,2,3,1) = lssss(1,2,3,3)
      lssss(3,3,1,2) = lssss(1,2,3,3)
      lssss(3,3,2,1) = lssss(1,2,3,3)
      lssss(1,3,3,2) = lssss(1,2,3,3)
      lssss(2,3,3,1) = lssss(1,2,3,3)      
      
      lssss(1,1,3,3) = ll**2/24d0
      lssss(1,3,1,3) = lssss(1,1,3,3)
      lssss(1,3,3,1) = lssss(1,1,3,3)
      lssss(3,1,1,3) = lssss(1,1,3,3)
      lssss(3,1,3,1) = lssss(1,1,3,3)
      lssss(3,3,1,1) = lssss(1,1,3,3)
      lssss(2,2,3,3) = lssss(1,1,3,3)
      lssss(2,3,2,3) = lssss(1,1,3,3)
      lssss(2,3,3,2) = lssss(1,1,3,3)
      lssss(3,2,2,3) = lssss(1,1,3,3)
      lssss(3,2,3,2) = lssss(1,1,3,3)
      lssss(3,3,2,2) = lssss(1,1,3,3)

      lsspp(1,1,1,1) = gb2/8d0
      lsspp(2,2,2,2) = lsspp(1,1,1,1)

      lsspp(1,1,2,2) = (2*ll**2-gb2)/8d0
      lsspp(2,2,1,1) = lsspp(1,1,2,2)

      lsspp(1,1,3,3) = ll**2/4d0
      lsspp(2,2,3,3) = lsspp(1,1,3,3)
      lsspp(3,3,1,1) = lsspp(1,1,3,3)
      lsspp(3,3,2,2) = lsspp(1,1,3,3)

      lsspp(1,2,3,3) = ll*kk/4d0
      lsspp(2,1,3,3) = lsspp(1,2,3,3)
      lsspp(3,3,1,2) = lsspp(1,2,3,3)
      lsspp(3,3,2,1) = lsspp(1,2,3,3)
      
      lsspp(1,3,2,3) = -ll*kk/4d0
      lsspp(3,1,2,3) = lsspp(1,3,2,3)
      lsspp(1,3,3,2) = lsspp(1,3,2,3)
      lsspp(3,1,3,2) = lsspp(1,3,2,3)
      lsspp(2,3,1,3) = lsspp(1,3,2,3)
      lsspp(3,2,1,3) = lsspp(1,3,2,3)
      lsspp(2,3,3,1) = lsspp(1,3,2,3)
      lsspp(3,2,3,1) = lsspp(1,3,2,3)
      
      lsspp(3,3,3,3) = kk**2/2d0

c     rotate the quartics

      do i=1,3              
         do j=1,3
            do k=1,3
               do l=1,3
                  do a=1,3
                     do b=1,3
                        lsshh(i,j,k,l) = lsshh(i,j,k,l)
     $                       +6*RS(k,a)*RS(l,b)*lssss(i,j,a,b)
                        lssaa(i,j,k,l) = lssaa(i,j,k,l)
     $                       +  RP(k,a)*RP(l,b)*lsspp(i,j,a,b)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      
c     trilinear neutral couplings

      lsss(1,1,1) = gb2*v1/2d0/sq2

      lsss(2,2,2) = gb2*v2/2d0/sq2
      
      lsss(1,2,2) = v1/6d0/sq2*(2*ll**2-gb2)
      lsss(2,1,2) = lsss(1,2,2)
      lsss(2,2,1) = lsss(1,2,2)
      
      lsss(1,1,2) = v2/6d0/sq2*(2*ll**2-gb2)
      lsss(1,2,1) = lsss(1,1,2)
      lsss(2,1,1) = lsss(1,1,2)
      
      lsss(3,1,1) = ll**2*xx/3d0/sq2
      lsss(1,1,3) = lsss(3,1,1)
      lsss(1,3,1) = lsss(3,1,1)      
      lsss(3,2,2) = lsss(3,1,1)
      lsss(2,2,3) = lsss(3,1,1)
      lsss(2,3,2) = lsss(3,1,1)

      lsss(3,3,3) = kk*Ak/3d0/sq2 + sq2*kk**2*xx
      
      lsss(1,3,3) = ll/3d0/sq2*(ll*v1-kk*v2)
      lsss(3,1,3) = lsss(1,3,3)
      lsss(3,3,1) = lsss(1,3,3)

      lsss(2,3,3) = ll/3d0/sq2*(ll*v2-kk*v1)
      lsss(3,2,3) = lsss(2,3,3)
      lsss(3,3,2) = lsss(2,3,3)

      lsss(1,2,3) = -ll*Al/6d0/sq2 - ll*kk*xx/3d0/sq2
      lsss(1,3,2) = lsss(1,2,3)
      lsss(2,1,3) = lsss(1,2,3)
      lsss(2,3,1) = lsss(1,2,3)
      lsss(3,1,2) = lsss(1,2,3)
      lsss(3,2,1) = lsss(1,2,3)


      lspp(1,1,1) = gb2*v1/2d0/sq2

      lspp(2,2,2) = gb2*v2/2d0/sq2
      
      lspp(1,2,2) = v1/2d0/sq2*(2*ll**2-gb2)

      lspp(2,1,1) = v2/2d0/sq2*(2*ll**2-gb2)

      lspp(3,1,1) = ll**2*xx/sq2
      lspp(3,2,2) = lspp(3,1,1)

      lspp(3,3,3) = -kk*Ak/sq2 + sq2*kk**2*xx
      
      lspp(1,3,3) = ll/sq2*(ll*v1+kk*v2)

      lspp(2,3,3) = ll/sq2*(ll*v2+kk*v1)

      lspp(3,1,3) = -ll*kk*v2/sq2
      lspp(3,3,1) = lspp(3,1,3)

      lspp(3,2,3) = -ll*kk*v1/sq2
      lspp(3,3,2) = lspp(3,2,3)

      lspp(1,2,3) = ll*Al/2d0/sq2-ll*kk*xx/sq2
      lspp(1,3,2) = lspp(1,2,3)
      lspp(2,1,3) = lspp(1,2,3)
      lspp(2,3,1) = lspp(1,2,3)      

      lspp(3,1,2) = ll*Al/2d0/sq2+ll*kk*xx/sq2
      lspp(3,2,1) = lspp(3,1,2)

c     rotate the trilinears

      do i=1,3              
         do k=1,3
            do l=1,3
               do a=1,3
                  do b=1,3
                     lshh(i,k,l) = lshh(i,k,l)
     $                    +3*RS(k,a)*RS(l,b)*lsss(i,a,b)
                     lsaa(i,k,l) = lsaa(i,k,l)
     $                    +  RP(k,a)*RP(l,b)*lspp(i,a,b)
                  enddo
               enddo
            enddo
         enddo
      enddo

c     quartic charged couplings

      c2b = (v1**2-v2**2)/(v1**2+v2**2)
      s2b = 2*v1*v2/(v1**2+v2**2)

      lsscc(1,1,1,1) = (g**2+gp**2*c2b)/8d0
      lsscc(2,2,2,2) = lsscc(1,1,1,1)     

      lsscc(1,1,2,2) = (g**2-gp**2*c2b)/8d0
      lsscc(2,2,1,1) = lsscc(1,1,2,2)     
      
      lsscc(1,2,1,1) = (2*ll**2-g**2)*s2b/8d0
      lsscc(2,1,1,1) = lsscc(1,2,1,1)     
      
      lsscc(1,2,2,2) = -(2*ll**2-g**2)*s2b/8d0
      lsscc(2,1,2,2) = lsscc(1,2,2,2)     

      lsscc(1,1,1,2) = -gp**2/8d0*s2b
      lsscc(1,1,2,1) = lsscc(1,1,1,2)

      lsscc(2,2,1,2) = gp**2/8d0*s2b
      lsscc(2,2,2,1) = lsscc(2,2,1,2)

      lsscc(1,2,1,2) = (2*ll**2-g**2)*c2b/8d0
      lsscc(1,2,2,1) = lsscc(1,2,1,2)
      lsscc(2,1,1,2) = lsscc(1,2,1,2)
      lsscc(2,1,2,1) = lsscc(1,2,1,2)

      lsscc(3,3,1,2) = -kk*ll/2d0*c2b
      lsscc(3,3,2,1) = lsscc(3,3,1,2)

      lsscc(3,3,1,1) = ll/2d0*(ll-kk*s2b)

      lsscc(3,3,2,2) = ll/2d0*(ll+kk*s2b)
      
c     trilinear charged couplings
      
      lscc(1,1,1) = (v1*(g**2+gp**2*c2b)+v2*(2*ll**2-g**2)*s2b)/2d0/sq2

      lscc(1,2,2) = (v1*(g**2-gp**2*c2b)-v2*(2*ll**2-g**2)*s2b)/2d0/sq2

      lscc(1,1,2) = (-v1*gp**2*s2b+v2*(2*ll**2-g**2)*c2b)/2d0/sq2
      lscc(1,2,1) = lscc(1,1,2)

      lscc(2,1,1) = (v2*(g**2-gp**2*c2b)+v1*(2*ll**2-g**2)*s2b)/2d0/sq2

      lscc(2,2,2) = (v2*(g**2+gp**2*c2b)-v1*(2*ll**2-g**2)*s2b)/2d0/sq2

      lscc(2,1,2) = (v2*gp**2*s2b+v1*(2*ll**2-g**2)*c2b)/2d0/sq2
      lscc(2,2,1) = lscc(2,1,2)

      lscc(3,1,1) = ll/sq2*(2*ll*xx-(Al+2*kk*xx)*s2b)

      lscc(3,2,2) = ll/sq2*(2*ll*xx+(Al+2*kk*xx)*s2b)

      lscc(3,1,2) = -ll/sq2*(Al+2*kk*xx)*c2b
      lscc(3,2,1) = lscc(3,1,2)

      return
      end

c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     

      subroutine coupl_s_ino(g,gp,ll,kk,NN,UU,VV,lsnene,lschch)
      
      implicit none

      double precision g,gp,ll,kk,NN(5,5),UU(2,2),VV(2,2),
     $     lsnene(3,5,5),lschch(3,2,2)

      double precision ls00(3,5,5),lspm(3,2,2),sq2

      integer i,k,l,a,b

      sq2 = sqrt(2d0)

c     neutralino couplings

      do i=1,3
         do k=1,5
            do l=1,5
               ls00(i,k,l) = 0d0
               lsnene(i,k,l) = 0d0
            enddo
         enddo
      enddo

      ls00(1,1,3) = - gp/4d0
      ls00(1,3,1) = ls00(1,1,3)      

      ls00(2,1,4) = gp/4d0
      ls00(2,4,1) = ls00(2,1,4)      

      ls00(1,2,3) = g/4d0
      ls00(1,3,2) = ls00(1,2,3)      

      ls00(2,2,4) = -g/4d0
      ls00(2,4,2) = ls00(2,2,4)      

      ls00(3,5,5) = kk/sq2

      ls00(1,4,5) = -ll/2d0/sq2
      ls00(1,5,4) = ls00(1,4,5)
      ls00(2,3,5) = ls00(1,4,5)
      ls00(2,5,3) = ls00(1,4,5)
      ls00(3,3,4) = ls00(1,4,5)
      ls00(3,4,3) = ls00(1,4,5)

      do i=1,3
         do k=1,5
            do l=1,5
               do a=1,5
                  do b=1,5
                     lsnene(i,k,l) = lsnene(i,k,l)
     $                    +NN(k,a)*NN(l,b)*ls00(i,a,b)
                  enddo
               enddo
            enddo
         enddo
      enddo
         
c     chargino couplings

      do i=1,3
         do k=1,2
            do l=1,2
               lspm(i,k,l) = 0d0
               lschch(i,k,l) = 0d0
            enddo
         enddo
      enddo

      lspm(1,1,2) = g/sq2

      lspm(2,2,1) = g/sq2

      lspm(3,2,2) = ll/sq2

      do i=1,3
         do k=1,2
            do l=1,2
               do a=1,2
                  do b=1,2
                     lschch(i,k,l) = lschch(i,k,l)
     $                    +VV(k,a)*UU(l,b)*lspm(i,a,b)
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
      end
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     

      subroutine coupl_p_sf(g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     $     Rt,Rb,Rtau,lpptt,lppbb,lpptata,lppntnt,lppuu,lppdd,
     $     lppee,lppnn,lptt,lpbb,lptata,lpntnt,lpuu,lpdd,lpee,lpnn)

      implicit none

      double precision g,gp,ll,ht,hb,htau,v1,v2,xx,At,Ab,Atau,
     $     Rt(2,2),Rb(2,2),Rtau(2,2),lpptt(3,3,2,2),lppbb(3,3,2,2),
     $     lpptata(3,3,2,2),lppntnt(3,3),lppuu(3,3,2,2),lppdd(3,3,2,2),
     $     lppee(3,3,2,2),lppnn(3,3),lptt(3,2,2),lpbb(3,2,2),
     $     lptata(3,2,2),lpntnt(3),lpuu(3,2,2),lpdd(3,2,2),lpee(3,2,2),
     $     lpnn(3),lpptt_int(3,3,2,2),lppbb_int(3,3,2,2),
     $     lpptata_int(3,3,2,2),lptt_int(3,2,2),lpbb_int(3,2,2),
     $     lptata_int(3,2,2)

      double precision sq2,gb2,sw2,guL,guR,gdL,gdR,geL,geR,gnu

      integer i,j,k,l,a,b

      sq2 = sqrt(2d0)
      gb2 = (g**2+gp**2)/2d0
      sw2 = gp**2/2d0/gb2

      guL =  .5d0 - 2d0/3d0*sw2
      gdL = -.5d0 + 1d0/3d0*sw2
      gnu =  .5d0
      geL = -.5d0 + sw2
      guR =  2d0/3d0*sw2
      gdR = -1d0/3d0*sw2
      geR = -sw2

c     FIRST TWO GENERATIONS

c     quartic, up squarks

      do i=1,3
         do j=1,3
            do k=1,2
               do l=1,2
                  lppuu(i,j,k,l) = 0d0
               enddo
            enddo
         enddo
      enddo

      lppuu(1,1,1,1) = gb2*guL/2d0
      lppuu(1,1,2,2) = gb2*guR/2d0
      lppuu(2,2,1,1) = -gb2*guL/2d0
      lppuu(2,2,2,2) = -gb2*guR/2d0
      
c     trilinear, up squarks

      do i=1,3
         do k=1,2
            do l=1,2
               lpuu(i,k,l) = 0d0
            enddo
         enddo
      enddo

c     quartic, down squarks

      do i=1,3
         do j=1,3
            do k=1,2
               do l=1,2
                  lppdd(i,j,k,l) = 0d0
               enddo
            enddo
         enddo
      enddo

      lppdd(1,1,1,1) = gb2*gdL/2d0
      lppdd(1,1,2,2) = gb2*gdR/2d0
      lppdd(2,2,1,1) = -gb2*gdL/2d0
      lppdd(2,2,2,2) = -gb2*gdR/2d0
      
c     trilinear, down squarks

      do i=1,3
         do k=1,2
            do l=1,2
               lpdd(i,k,l) = 0d0
            enddo
         enddo
      enddo

c     quartic, charged sleptons

      do i=1,3
         do j=1,3
            do k=1,2
               do l=1,2
                  lppee(i,j,k,l) = 0d0
               enddo
            enddo
         enddo
      enddo

      lppee(1,1,1,1) = gb2*geL/2d0
      lppee(1,1,2,2) = gb2*geR/2d0
      lppee(2,2,1,1) = -gb2*geL/2d0
      lppee(2,2,2,2) = -gb2*geR/2d0

c     trilinear, charged sleptons

      do i=1,3
         do k=1,2
            do l=1,2
               lpee(i,k,l) = 0d0
            enddo
         enddo
      enddo

c     quartic, sneutrinos

      do i=1,3
         do j=1,3
            lppnn(i,j) = 0d0
         enddo
      enddo

      lppnn(1,1) = gb2*gnu/2d0
      lppnn(2,2) = -gb2*gnu/2d0
      
c     trilinear, sneutrinos

      lpnn(1) = 0d0
      lpnn(2) = 0d0
      lpnn(3) = 0d0

c     THIRD GENERATION

c     quartic, top squarks

      do i=1,3
         do j=1,3
            do k=1,2
               do l=1,2
                  lpptt_int(i,j,k,l) = lppuu(i,j,k,l)
               enddo
            enddo
         enddo
      enddo

      lpptt_int(2,2,1,1) = lpptt_int(2,2,1,1) + ht**2/2d0
      lpptt_int(2,2,2,2) = lpptt_int(2,2,2,2) + ht**2/2d0

      lpptt_int(1,3,1,2) = lpptt_int(1,3,1,2) + ht*ll/4d0
      lpptt_int(1,3,2,1) = lpptt_int(1,3,1,2)
      lpptt_int(3,1,1,2) = lpptt_int(1,3,1,2)
      lpptt_int(3,1,2,1) = lpptt_int(1,3,1,2)

c     trilinear, top squarks

      do i=1,3
         do k=1,2
            do l=1,2
               lptt_int(i,k,l) = 0d0
            enddo
         enddo
      enddo

      lptt_int(1,1,2) = ht*ll*xx/sq2
      lptt_int(2,1,2) = ht*At/sq2
      lptt_int(3,1,2) = ht*ll*v1/sq2

      lptt_int(1,2,1) = -lptt_int(1,1,2) 
      lptt_int(2,2,1) = -lptt_int(2,1,2) 
      lptt_int(3,2,1) = -lptt_int(3,1,2) 

c     quartic, bottom squarks

      do i=1,3
         do j=1,3
            do k=1,2
               do l=1,2
                  lppbb_int(i,j,k,l) = lppdd(i,j,k,l)
               enddo
            enddo
         enddo
      enddo

      lppbb_int(1,1,1,1) = lppbb_int(1,1,1,1) + hb**2/2d0
      lppbb_int(1,1,2,2) = lppbb_int(1,1,2,2) + hb**2/2d0

      lppbb_int(2,3,1,2) = lppbb_int(2,3,1,2) + hb*ll/4d0
      lppbb_int(2,3,2,1) = lppbb_int(2,3,1,2)
      lppbb_int(3,2,1,2) = lppbb_int(2,3,1,2)
      lppbb_int(3,2,2,1) = lppbb_int(2,3,1,2)

c     trilinear, bottom squarks

      do i=1,3
         do k=1,2
            do l=1,2
               lpbb_int(i,k,l) = 0d0
            enddo
         enddo
      enddo

      lpbb_int(1,1,2) = hb*Ab/sq2
      lpbb_int(2,1,2) = hb*ll*xx/sq2
      lpbb_int(3,1,2) = hb*ll*v2/sq2

      lpbb_int(1,2,1) = -lpbb_int(1,1,2) 
      lpbb_int(2,2,1) = -lpbb_int(2,1,2) 
      lpbb_int(3,2,1) = -lpbb_int(3,1,2) 

c     quartic, staus

      do i=1,3
         do j=1,3
            do k=1,2
               do l=1,2
                   lpptata_int(i,j,k,l) = lppee(i,j,k,l)
               enddo
            enddo
         enddo
      enddo

      lpptata_int(1,1,1,1) = lpptata_int(1,1,1,1) + htau**2/2d0
      lpptata_int(1,1,2,2) = lpptata_int(1,1,2,2) + htau**2/2d0

      lpptata_int(2,3,1,2) = lpptata_int(2,3,1,2) + htau*ll/4d0
      lpptata_int(2,3,2,1) = lpptata_int(2,3,1,2)
      lpptata_int(3,2,1,2) = lpptata_int(2,3,1,2)
      lpptata_int(3,2,2,1) = lpptata_int(2,3,1,2)
      
c     trilinear, staus

      do i=1,3
         do k=1,2
            do l=1,2
               lptata_int(i,k,l) = 0d0
            enddo
         enddo
      enddo

      lptata_int(1,1,2) = htau*Atau/sq2
      lptata_int(2,1,2) = htau*ll*xx/sq2
      lptata_int(3,1,2) = htau*ll*v2/sq2

      lptata_int(1,2,1) = -lptata_int(1,1,2) 
      lptata_int(2,2,1) = -lptata_int(2,1,2) 
      lptata_int(3,2,1) = -lptata_int(3,1,2) 

c     sneutrinos (same as first two generations)

      do i=1,3
         lpntnt(i) = lpnn(i)
         do j=1,3
            lppntnt(i,j) = lppnn(i,j)
         enddo
      enddo

c     now rotate the third-generation couplings

c     quartic
      do i = 1,3
         do j = 1,3
            do k=1,2
               do l=1,2
                  lpptt(i,j,k,l) = 0d0
                  lppbb(i,j,k,l) = 0d0
                  lpptata(i,j,k,l) = 0d0
                  do a=1,2
                     do b=1,2
                        lpptt(i,j,k,l) = lpptt(i,j,k,l)
     $                       +Rt(k,a)*Rt(l,b)*lpptt_int(i,j,a,b)
                        lppbb(i,j,k,l) = lppbb(i,j,k,l)
     $                       +Rb(k,a)*Rb(l,b)*lppbb_int(i,j,a,b)
                        lpptata(i,j,k,l) = lpptata(i,j,k,l)
     $                       +Rtau(k,a)*Rtau(l,b)*lpptata_int(i,j,a,b)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
           
c     trilinear
      do i = 1,3         
         do k=1,2
            do l=1,2
               lptt(i,k,l) = 0d0
               lpbb(i,k,l) = 0d0
               lptata(i,k,l) = 0d0
               do a=1,2
                  do b=1,2
                     lptt(i,k,l) = lptt(i,k,l)
     $                    +Rt(k,a)*Rt(l,b)*lptt_int(i,a,b)
                     lpbb(i,k,l) = lpbb(i,k,l)
     $                    +Rb(k,a)*Rb(l,b)*lpbb_int(i,a,b)
                     lptata(i,k,l) = lptata(i,k,l)
     $                    +Rtau(k,a)*Rtau(l,b)*lptata_int(i,a,b)
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
      end


c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     

      subroutine coupl_p_hh(g,gp,ll,kk,v1,v2,xx,Al,Ak,RS,RP,
     $     lpphh,lppaa,lpah,lppcc,lpcc)

      implicit none

      double precision g,gp,ll,kk,v1,v2,xx,Al,Ak,RS(3,3),RP(3,3),
     $     lpphh(3,3,3,3),lppaa(3,3,3,3),lpah(3,3,3),
     $     lppcc(3,3,2,2),lpcc(3,2,2)

      double precision lpppp(3,3,3,3),lsspp(3,3,3,3),
     $     lspp(3,3,3),gb2,sq2,c2b,s2b

      integer i,j,k,l,a,b
      
      gb2 = (g**2+gp**2)/2d0
      sq2 = sqrt(2d0)
      
      do i=1,3                  ! initialize
         do j=1,3
            do k=1,3
               lspp(i,j,k)=0d0
               lpah(i,j,k)=0d0
               lpcc(i,j,k)=0d0               
               do l=1,3
                  lpppp(i,j,k,l)=0d0
                  lsspp(i,j,k,l)=0d0
                  lpphh(i,j,k,l)=0d0
                  lppaa(i,j,k,l)=0d0
                  lppcc(i,j,k,l)=0d0
               enddo
            enddo
         enddo
      enddo

c     quartic neutral couplings
 
      lpppp(1,1,1,1) = gb2/16d0

      lpppp(2,2,2,2) = gb2/16d0

      lpppp(3,3,3,3) = kk**2/4d0

      lpppp(1,1,2,2) = (2*ll**2-gb2)/48d0
      lpppp(1,2,1,2) = lpppp(1,1,2,2) ! symmetrize the couplings
      lpppp(1,2,2,1) = lpppp(1,1,2,2) ! (got a better proposal? ;-)
      lpppp(2,1,1,2) = lpppp(1,1,2,2)
      lpppp(2,1,2,1) = lpppp(1,1,2,2)
      lpppp(2,2,1,1) = lpppp(1,1,2,2)

      lpppp(1,2,3,3) = -ll*kk/24d0
      lpppp(2,1,3,3) = lpppp(1,2,3,3)
      lpppp(2,3,1,3) = lpppp(1,2,3,3)
      lpppp(1,3,2,3) = lpppp(1,2,3,3)
      lpppp(3,1,2,3) = lpppp(1,2,3,3)
      lpppp(3,2,1,3) = lpppp(1,2,3,3)
      lpppp(3,1,3,2) = lpppp(1,2,3,3)
      lpppp(3,2,3,1) = lpppp(1,2,3,3)
      lpppp(3,3,1,2) = lpppp(1,2,3,3)
      lpppp(3,3,2,1) = lpppp(1,2,3,3)
      lpppp(1,3,3,2) = lpppp(1,2,3,3)
      lpppp(2,3,3,1) = lpppp(1,2,3,3)      
      
      lpppp(1,1,3,3) = ll**2/24d0
      lpppp(1,3,1,3) = lpppp(1,1,3,3)
      lpppp(1,3,3,1) = lpppp(1,1,3,3)
      lpppp(3,1,1,3) = lpppp(1,1,3,3)
      lpppp(3,1,3,1) = lpppp(1,1,3,3)
      lpppp(3,3,1,1) = lpppp(1,1,3,3)
      lpppp(2,2,3,3) = lpppp(1,1,3,3)
      lpppp(2,3,2,3) = lpppp(1,1,3,3)
      lpppp(2,3,3,2) = lpppp(1,1,3,3)
      lpppp(3,2,2,3) = lpppp(1,1,3,3)
      lpppp(3,2,3,2) = lpppp(1,1,3,3)
      lpppp(3,3,2,2) = lpppp(1,1,3,3)

      lsspp(1,1,1,1) = gb2/8d0
      lsspp(2,2,2,2) = lsspp(1,1,1,1)

      lsspp(1,1,2,2) = (2*ll**2-gb2)/8d0
      lsspp(2,2,1,1) = lsspp(1,1,2,2)

      lsspp(1,1,3,3) = ll**2/4d0
      lsspp(2,2,3,3) = lsspp(1,1,3,3)
      lsspp(3,3,1,1) = lsspp(1,1,3,3)
      lsspp(3,3,2,2) = lsspp(1,1,3,3)

      lsspp(1,2,3,3) = ll*kk/4d0
      lsspp(2,1,3,3) = lsspp(1,2,3,3)
      lsspp(3,3,1,2) = lsspp(1,2,3,3)
      lsspp(3,3,2,1) = lsspp(1,2,3,3)
      
      lsspp(1,3,2,3) = -ll*kk/4d0
      lsspp(3,1,2,3) = lsspp(1,3,2,3)
      lsspp(1,3,3,2) = lsspp(1,3,2,3)
      lsspp(3,1,3,2) = lsspp(1,3,2,3)
      lsspp(2,3,1,3) = lsspp(1,3,2,3)
      lsspp(3,2,1,3) = lsspp(1,3,2,3)
      lsspp(2,3,3,1) = lsspp(1,3,2,3)
      lsspp(3,2,3,1) = lsspp(1,3,2,3)
      
      lsspp(3,3,3,3) = kk**2/2d0

c     rotate the quartics

      do i=1,3              
         do j=1,3
            do k=1,3
               do l=1,3
                  do a=1,3
                     do b=1,3
                        lpphh(i,j,k,l) = lpphh(i,j,k,l)
     $                       + RS(k,a)*RS(l,b)*lsspp(a,b,i,j)
                        lppaa(i,j,k,l) = lppaa(i,j,k,l)
     $                       +  6*RP(k,a)*RP(l,b)*lpppp(i,j,a,b)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      
c     trilinear neutral couplings

      lspp(1,1,1) = gb2*v1/2d0/sq2

      lspp(2,2,2) = gb2*v2/2d0/sq2
      
      lspp(1,2,2) = v1/2d0/sq2*(2*ll**2-gb2)

      lspp(2,1,1) = v2/2d0/sq2*(2*ll**2-gb2)

      lspp(3,1,1) = ll**2*xx/sq2
      lspp(3,2,2) = lspp(3,1,1)

      lspp(3,3,3) = -kk*Ak/sq2 + sq2*kk**2*xx
      
      lspp(1,3,3) = ll/sq2*(ll*v1+kk*v2)

      lspp(2,3,3) = ll/sq2*(ll*v2+kk*v1)

      lspp(3,1,3) = -ll*kk*v2/sq2
      lspp(3,3,1) = lspp(3,1,3)

      lspp(3,2,3) = -ll*kk*v1/sq2
      lspp(3,3,2) = lspp(3,2,3)

      lspp(1,2,3) = ll*Al/2d0/sq2-ll*kk*xx/sq2
      lspp(1,3,2) = lspp(1,2,3)
      lspp(2,1,3) = lspp(1,2,3)
      lspp(2,3,1) = lspp(1,2,3)      

      lspp(3,1,2) = ll*Al/2d0/sq2+ll*kk*xx/sq2
      lspp(3,2,1) = lspp(3,1,2)

c     rotate the trilinears

      do i=1,3              
         do k=1,3
            do l=1,3
               do a=1,3
                  do b=1,3
                     lpah(i,k,l) = lpah(i,k,l)
     $                    +  2*RS(l,a)*RP(k,b)*lspp(a,b,i)
                  enddo
               enddo
            enddo
         enddo
      enddo

c     quartic charged couplings

      c2b = (v1**2-v2**2)/(v1**2+v2**2)
      s2b = 2*v1*v2/(v1**2+v2**2)

      lppcc(1,1,1,1) = (g**2+gp**2*c2b)/8d0
      lppcc(2,2,2,2) = lppcc(1,1,1,1)     

      lppcc(1,1,2,2) = (g**2-gp**2*c2b)/8d0
      lppcc(2,2,1,1) = lppcc(1,1,2,2)     
      
      lppcc(1,2,1,1) = -(2*ll**2-g**2)*s2b/8d0
      lppcc(2,1,1,1) = lppcc(1,2,1,1)     
      
      lppcc(1,2,2,2) = (2*ll**2-g**2)*s2b/8d0
      lppcc(2,1,2,2) = lppcc(1,2,2,2)     

      lppcc(1,1,1,2) = -gp**2/8d0*s2b
      lppcc(1,1,2,1) = lppcc(1,1,1,2)

      lppcc(2,2,1,2) = gp**2/8d0*s2b
      lppcc(2,2,2,1) = lppcc(2,2,1,2)

      lppcc(1,2,1,2) = -(2*ll**2-g**2)*s2b/8d0
      lppcc(1,2,2,1) = lppcc(1,2,1,2)
      lppcc(2,1,1,2) = lppcc(1,2,1,2)
      lppcc(2,1,2,1) = lppcc(1,2,1,2)

      lppcc(3,3,1,2) = kk*ll/2d0*c2b
      lppcc(3,3,2,1) = lppcc(3,3,1,2)

      lppcc(3,3,1,1) = ll/2d0*(ll+kk*s2b)

      lppcc(3,3,2,2) = ll/2d0*(ll-kk*s2b)
      
c     trilinear charged couplings
      
      lpcc(1,1,2) = v2*(2*ll**2-g**2)/2d0/sq2
      lpcc(1,2,1) = -lpcc(1,1,2)

      lpcc(2,1,2) = v1*(2*ll**2-g**2)/2d0/sq2
      lpcc(2,2,1) = -lpcc(2,1,2)

      lpcc(3,1,2) = ll/sq2*(Al-2*kk*xx)
      lpcc(3,2,1) = -lpcc(3,1,2)

      return
      end

c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     

      subroutine coupl_p_ino(g,gp,ll,kk,NN,UU,VV,lpnene,lpchch)
      
      implicit none

      double precision g,gp,ll,kk,NN(5,5),UU(2,2),VV(2,2),
     $     lpnene(3,5,5),lpchch(3,2,2)

      double precision lp00(3,5,5),lppm(3,2,2),sq2

      integer i,k,l,a,b

      sq2 = sqrt(2d0)

c     neutralino couplings

      do i=1,3
         do k=1,5
            do l=1,5
               lp00(i,k,l) = 0d0
               lpnene(i,k,l) = 0d0
            enddo
         enddo
      enddo

      lp00(1,1,3) = gp/4d0
      lp00(1,3,1) = lp00(1,1,3)      

      lp00(2,1,4) = -gp/4d0
      lp00(2,4,1) = lp00(2,1,4)      

      lp00(1,2,3) = -g/4d0
      lp00(1,3,2) = lp00(1,2,3)      

      lp00(2,2,4) = g/4d0
      lp00(2,4,2) = lp00(2,2,4)      

      lp00(3,5,5) = kk/sq2

      lp00(1,4,5) = -ll/2d0/sq2
      lp00(1,5,4) = lp00(1,4,5)
      lp00(2,3,5) = lp00(1,4,5)
      lp00(2,5,3) = lp00(1,4,5)
      lp00(3,3,4) = lp00(1,4,5)
      lp00(3,4,3) = lp00(1,4,5)

      do i=1,3
         do k=1,5
            do l=1,5
               do a=1,5
                  do b=1,5
                     lpnene(i,k,l) = lpnene(i,k,l)
     $                    +NN(k,a)*NN(l,b)*lp00(i,a,b)
                  enddo
               enddo
            enddo
         enddo
      enddo
         
c     chargino couplings

      do i=1,3
         do k=1,2
            do l=1,2
               lppm(i,k,l) = 0d0
               lpchch(i,k,l) = 0d0
            enddo
         enddo
      enddo

      lppm(1,1,2) = -g/sq2

      lppm(2,2,1) = -g/sq2

      lppm(3,2,2) = ll/sq2

      do i=1,3
         do k=1,2
            do l=1,2
               do a=1,2
                  do b=1,2
                     lpchch(i,k,l) = lpchch(i,k,l)
     $                    +VV(k,a)*UU(l,b)*lppm(i,a,b)
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
      end
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     

      subroutine coupl_Z_ino(g,gp,NN,UU,VV,lznene,azchch,bzchch)
      
      implicit none

      double precision g,gp,NN(5,5),UU(2,2),VV(2,2),
     $     lznene(5,5),azchch(2,2),bzchch(2,2)

      double precision lz00(5,5),lzpm(2,2),sq2,gb2,c2w

      integer i,j,k,l

      sq2 = sqrt(2d0)
      gb2 = (g**2+gp**2)/2d0
      c2w = (g**2-gp**2)/(g**2+gp**2)

c     neutralino couplings

      do i=1,5
         do j=1,5
            lz00(i,j) = 0d0
            lznene(i,j) = 0d0
         enddo
      enddo

      lz00(3,3) = sqrt(gb2)/2d0/sq2
      
      lz00(4,4) = -lz00(3,3)
      
      do i=1,5
         do j=1,5
            do k=1,5
               do l=1,5
                  lznene(i,j) = lznene(i,j)
     $                 +NN(i,k)*NN(j,l)*lz00(k,l)
               enddo
            enddo
         enddo
      enddo
         
c     chargino couplings

      do i=1,2
         do j=1,2
            lzpm(i,j) = 0d0
            azchch(i,j) = 0d0
            bzchch(i,j) = 0d0
         enddo
      enddo

      lzpm(1,1) = g**2/sqrt(gb2)/sq2
      
      lzpm(2,2) = sqrt(gb2)/sq2*c2w

      do i=1,2
         do j=1,2
            do k=1,2
               do l=1,2
                  azchch(i,j) = azchch(i,j)
     $                 +VV(i,k)*VV(j,l)*lzpm(k,l)
                  bzchch(i,j) = bzchch(i,j)
     $                 +UU(i,k)*UU(j,l)*lzpm(k,l)
               enddo
            enddo
         enddo
      enddo

      return
      end
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     

      subroutine coupl_W_ino(g,NN,UU,VV,awnech,bwnech)
      
      implicit none

      double precision g,NN(5,5),UU(2,2),VV(2,2),
     $     awnech(5,2),bwnech(5,2)

      double precision aw0p(5,2),bw0p(5,2),sq2

      integer i,j,k,l

      sq2 = sqrt(2d0)

      do i=1,5
         do j=1,2
            aw0p(i,j) = 0d0
            bw0p(i,j) = 0d0
            awnech(i,j) = 0d0
            bwnech(i,j) = 0d0
         enddo
      enddo

      aw0p(2,1) = -g
      aw0p(4,2) = g/sq2
      
      bw0p(2,1) = -g
      bw0p(3,2) = -g/sq2
      
      do i=1,5
         do j=1,2
            do k=1,5
               do l=1,2
                  awnech(i,j) = awnech(i,j)
     $                 +NN(i,k)*VV(j,l)*aw0p(k,l)
                  bwnech(i,j) = bwnech(i,j)
     $                 +NN(i,k)*UU(j,l)*bw0p(k,l)
               enddo
            enddo
         enddo
      enddo
         
      return
      end

*
***********************************************************************
*

      double precision function myA0(m,q)      
      double precision m,q

      if(m.ne.0d0) then
         myA0 = m*(1d0-Log(m/q))
      else
         myA0 = 0d0
      endif

      return
      end

*
***********************************************************************
*


      double precision function myB0NM(q,m1,m2,mu2) 

c     from Degrassi and Sirlin, Phys. Rev. D46 (1992) 3104.
      
      double precision q,m1,m2,OmegaNM,mu2

      if(q.eq.0d0) then

         if(m1.eq.0d0.and.m2.ne.0d0) then
            myB0NM = 1d0-Log(m2/mu2)
         elseif(m1.ne.0d0.and.m2.eq.0d0) then
            myB0NM = 1d0-Log(m1/mu2)
         elseif(abs(m1-m2).le.1d-8) then
            myB0NM = -Log(m1/mu2)
         else
            myB0NM = 1d0 - Log(m2/mu2) + m1/(m1-m2)*Log(m2/m1)
         endif
         
      else

         if(m1.eq.0d0.and.m2.ne.0d0) then
            
            if(m2.ne.q) then
               myB0NM = -(Log(m2/mu2)-2-(m2/q-1d0)*Log(abs(1d0-q/m2))) 
            else 
               myB0NM = -(Log(m2/mu2) - 2)
            endif
            
         elseif(m2.eq.0d0.and.m1.ne.0d0) then
            
            if(m1.ne.q) then
               myB0NM = -(Log(m1/mu2)-2-(m1/q-1d0)*Log(abs(1d0-q/m1))) 
            else
               myB0NM = -(Log(m1/mu2) - 2)
            endif
            
         elseif(m2.eq.0d0.and.m1.eq.0d0) then
            
            myB0NM = -(Log(q/mu2) - 2) ! cut the imaginary part (I Pi)
            
         else
c$$$            print *, "in general case for myB0NM."
            myB0NM = -( log(q/mu2)-2.d0 + 
     1           1.d0/2.d0*( 1.d0 + (m1/q-m2/q))*log(m1/q) +
     2           1.d0/2.d0*( 1.d0 - (m1/q-m2/q))*log(m2/q) +
     3           2.d0*OmegaNM(m1/q,m2/q))
c$$$            print *, log(q/mu2)-2.d0
c$$$            print *, 1.d0/2.d0*( 1.d0 + (m1/q-m2/q))*log(m1/q)
c$$$            print *, 1.d0/2.d0*( 1.d0 - (m1/q-m2/q))*log(m2/q)
c$$$            print *, 2.d0*OmegaNM(m1/q,m2/q)
c$$$            print *, "m1/q = ", m1/q
c$$$            print *, "m2/q = ", m2/q
         endif
         
      endif

      return
      end
      
c     function OmegaNM(a,b) contained in myB0NM
      double precision function OmegaNM(a,b)
      double precision a,b,cbig
      Cbig = (a+b)/2.d0 - (a-b)**2.d0/4.d0 -1.d0/4.d0
      if(Cbig.gt.0.d0) then
         OmegaNM = sqrt(Cbig)*
     1        (atan((1.d0 + a - b)/(2.d0*sqrt(Cbig))) +
     2        atan((1.d0 - a + b)/(2.d0*sqrt(Cbig))) )
      elseif(Cbig.lt.0d0) then
         Cbig = - Cbig
         OmegaNM = 1.d0/2.d0*sqrt(Cbig)*
     1        log((a/2.d0 +b/2.d0 -1.d0/2.d0 -sqrt(Cbig))/
     2        (a/2.d0 + b/2.d0 -1.d0/2.d0 + sqrt(Cbig)))
      else
         OmegaNM = 0         
      endif

      return
      end

*
**********************************************************************
*
      
      double precision function myB1NM(p,m1,m2,q)

      implicit none

      double precision p,m1,m2,q
      double precision myA0,myB0NM
      
      if(p.eq.0d0) then
         if(abs(m1-m2).le.1d-8) then
            myB1NM = -Log(m1/q)/2d0
         else
            if(m1.eq.0d0) then
               myB1NM = (1d0-2d0*Log(m2/q))/4d0
            elseif(m2.eq.0d0) then
               myB1NM = (3d0-2d0*Log(m1/q))/4d0
            else
               myB1NM = (1d0-Log(m2/q)+m1**2/(m1-m2)**2*Log(m2/m1)
     $              +(m1+m2)/(m1-m2)/2d0)/2d0
            endif
         endif
      else
         myB1NM = (myA0(m2,q)-myA0(m1,q)
     $        +(p+m1-m2)*myB0NM(p,m1,m2,q))/2d0/p
      endif

      return
      end

*
**********************************************************************
*

      double precision function myF(q,m1,m2,mu2) 
      
      implicit none
      double precision q,m1,m2,mu2,myA0,myB0NM

      myF = myA0(m1,mu2)-2*myA0(m2,mu2)
     $     -(2*q+2*m1-m2)*myB0NM(q,m1,m2,mu2)

      return
      end

*
***********************************************************************
*

      double precision function myG(q,m1,m2,mu2) 
      
      implicit none
      double precision q,m1,m2,mu2,myA0,myB0NM

      if(q.eq.0d0.and.m1.eq.0d0.and.m2.eq.0d0) then
         myG = 0d0
      else
         myG = (q-m1-m2)*myB0NM(q,m1,m2,mu2) - myA0(m1,mu2) 
     $        - myA0(m2,mu2)
      endif

      return
      end

*
***********************************************************************
*

      double precision function myB22(q,m1,m2,mu2) 
      
      implicit none
      double precision q,m1,m2,mu2,myA0,myB0NM,myB1NM

      if(q.eq.0d0.and.m1.eq.0d0.and.m2.eq.0d0) then
         myB22 = 0d0
c$$$         print *, "IN myB22 q,m1 or m2 is zero"
      else
         myB22 = ((myA0(m1,mu2)+myA0(m2,mu2))/2d0
     $        +(m1+m2-q/2d0)*myB0NM(q,m1,m2,mu2)
     $        +(m2-m1)*(myB1NM(q,m1,m2,mu2)-myB0NM(q,m1,m2,mu2)/2d0)
     $        +m1+m2-q/3d0)/6d0
c$$$         print *, " myA0(m1,mu2) = ", myA0(m1,mu2)
c$$$         print *, " myA0(m2,mu2) = ", myA0(m2,mu2)
c$$$         print *, " myB0NM(q,m1,m2,mu2) = ", myB0NM(q,m1,m2,mu2)
c$$$         print *, " myB1NM(q,m1,m2,mu2) = ", myB1NM(q,m1,m2,mu2)
c$$$         print *, "m1 = ", m1 
c$$$         print *, "m2 = ", m2
c$$$         print *, "q = ", q 
         
      endif
      
      return
      end

*
***********************************************************************
*

      double precision function myB22T(q,m1,m2,mu2) 
      
      implicit none
      double precision q,m1,m2,mu2,myA0,myB22

      myB22T = myB22(q,m1,m2,mu2) - myA0(m1,mu2)/4d0 - myA0(m2,mu2)/4d0
c$$$      print *, "myB22(q,m1,m2,mu2) = ", myB22(q,m1,m2,mu2)
c$$$      print *, "myA0(m1,mu2)/4d0 = ", myA0(m1,mu2)/4d0
c$$$      print *, "myA0(m2,mu2)/4d0 = ", myA0(m2,mu2)/4d0
      

      return
      end

*
**********************************************************************
*

      double precision function myH(q,m1,m2,mu2) 
      
      implicit none
      double precision q,m1,m2,mu2,myG,myB22

      myH = 4*myB22(q,m1,m2,mu2) + myG(q,m1,m2,mu2)

      return
      end

c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

      SUBROUTINE jacobi(a,n,np,d,v,nrot)

c     just calls Tomas Hahn's diagonalization routine

      INTEGER n,np,nrot,i,j
      double precision a(np,np),d(np),v(np,np)      
      double complex M(np,np),U(np,np)

      do i=1,np                  ! turn to complex
         do j=1,np
            M(i,j) = DCMPLX(a(i,j))
         enddo
      enddo

      call HEigensystem(np, M, np, d, U, np, 0)
      
      do i=1,np                  ! back to real
         do j=1,np
            v(i,j) = DBLE(U(j,i)) ! the other jacobi had this convention
         enddo
      enddo

      return
      end

* diagonalization of a Hermitian n-by-n matrix using the Jacobi algorithm
* code adapted from the "Handbook" routines for complex A
* (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
* this file is part of the Diag library
* last modified 27 Sep 07 th
************************************************************************
** HEigensystem diagonalizes a Hermitian n-by-n matrix.
** Input: n, A = n-by-n matrix, Hermitian
** (only the upper triangle of A needs to be filled).
** Output: d = vector of eigenvalues, U = transformation matrix
** these fulfill diag(d) = U A U^+ = U A U^-1 with U unitary.

	subroutine HEigensystem(n, A,ldA, d, U,ldU, sort)
	implicit none
	integer n, ldA, ldU, sort
	double complex A(ldA,*), U(ldU,*)
	double precision d(*)

	integer p, q, j
	double precision red, off, thresh
	double precision delta, t, invc, s
	double complex x, y, Apq
	double precision ev(2,16)

	integer sweep
	common /nsweeps/ sweep

	double precision sq
	double complex c
	sq(c) = DBLE(c*DCONJG(c))

	if( n .gt. 16 ) then
	  print *, "Dimension too large"
	  d(1) = -999
	  return
	endif

	do p = 1, n
	  ev(1,p) = 0
	  ev(2,p) = DBLE(A(p,p))
	  d(p) = ev(2,p)
	enddo

	do p = 1, n
	  do q = 1, n
	    U(q,p) = 0
	  enddo
	  U(p,p) = 1
	enddo

	red = .04D0/n**4

	do sweep = 1, 50
	  off = 0
	  do q = 2, n
	    do p = 1, q - 1
	      off = off + sq(A(p,q))
	    enddo
	  enddo
	  if( off .lt. 2D0**(-103) ) goto 1

	  thresh = 0
	  if( sweep .lt. 4 ) thresh = off*red

	  do q = 2, n
	    do p = 1, q - 1
	      off = sq(A(p,q))
	      if( sweep .gt. 4 .and. off .lt.
     &              2D0**(-103)*max(ev(2,p)**2, ev(2,q)**2) ) then
	        A(p,q) = 0
	      else
	        if( off .gt. thresh ) then
	          t = .5D0*(ev(2,p) - ev(2,q))
	          t = 1/(t + sign(sqrt(t**2 + off), t))

	          delta = t*off
	          ev(1,p) = ev(1,p) + delta
	          ev(2,p) = d(p) + ev(1,p)
	          ev(1,q) = ev(1,q) - delta
	          ev(2,q) = d(q) + ev(1,q)

	          invc = sqrt(delta*t + 1)
	          s = t/invc
	          t = delta/(invc + 1)

	          Apq = A(p,q)

	          do j = 1, p - 1
	            x = A(j,p)
	            y = A(j,q)
	            A(j,p) = x + s*(DCONJG(Apq)*y - t*x)
	            A(j,q) = y - s*(Apq*x + t*y)
	          enddo

	          do j = p + 1, q - 1
	            x = A(p,j)
	            y = A(j,q)
	            A(p,j) = x + s*(Apq*DCONJG(y) - t*x)
	            A(j,q) = y - s*(Apq*DCONJG(x) + t*y)
	          enddo

	          do j = q + 1, n
	            x = A(p,j)
	            y = A(q,j)
	            A(p,j) = x + s*(Apq*y - t*x)
	            A(q,j) = y - s*(DCONJG(Apq)*x + t*y)
	          enddo

	          A(p,q) = 0

	          do j = 1, n
	            x = U(p,j)
	            y = U(q,j)
	            U(p,j) = x + s*(Apq*y - t*x)
	            U(q,j) = y - s*(DCONJG(Apq)*x + t*y)
	          enddo
	        endif
	      endif
	    enddo
	  enddo

	  do p = 1, n
	    ev(1,p) = 0
	    d(p) = ev(2,p)
	  enddo
	enddo

	print *, "Bad convergence in HEigensystem"

1	if( sort .eq. 0 ) return

* sort the eigenvalues

	do p = 1, n - 1
	  j = p
	  t = d(p)
	  do q = p + 1, n
	    if( sort*(t - d(q)) .gt. 0 ) then
	      j = q
	      t = d(q)
	    endif
	  enddo

	  if( j .ne. p ) then
	    d(j) = d(p)
	    d(p) = t
	    do q = 1, n
	      x = U(p,q)
	      U(p,q) = U(j,q)
	      U(j,q) = x
	    enddo
	  endif
	enddo
	end
