# 1 "cool1d_multi_bdf.src"
# 1 "<built-in>"
# 1 "<command line>"
# 1 "cool1d_multi_bdf.src"

# 1 "fortran.def" 1




















# 2 "cool1d_multi_bdf.src" 2
c=======================================================================
c//////////////////////  SUBROUTINE COOL1D_MULTI  \\\\\\\\\\\\\\\\\\\\\c

      subroutine cool1d_multi_bdf(
     &                d, e, ge, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, idual, imethod,
     &                iexpand, ispecies, imetal, idim,
     &                is, ie, j, k, ih2co, ipiht, iter,
     &                aye, temstart, temend,
     &                utem, uxyz, uaye, urho, utim,
     &                eta1, eta2, gamma,
     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &                reHeII2a, reHeIIIa, brema, compa, 
     &                comp_xraya, comp_temp,
     &                piHI, piHeI, piHeII, comp1, comp2,
     &                HM, H2I, H2II, DI, DII, HDI, metal,
     &                hyd01ka, h2k01a, vibha, rotha, rotla,
     &                hyd01k, h2k01, vibh, roth, rotl,
     &                gpldla, gphdla, gpldl, gphdl,
     &                hdltea, hdlowa, hdlte, hdlow, 
     &                hdcoola, hdcool, ciecoa, cieco, 
     &                ceHI, ceHeI, ceHeII, ciHI, ciHeI, ciHeIS, ciHeII,
     &                reHII, reHeII1, reHeII2, reHeIII, brem,
     &                indixe, t1, t2, logtem, tdef, edotplus,
     &                edotminus, tgas, tgasold, p2d,
     &                inutot, iradtype, nfreq, imetalregen,
     &                iradshield, avgsighp, avgsighep, 
     &                avgsighe2p, itmask, iciecool,ih2optical,ifedtgas)
c
c  SOLVE RADIATIVE COOLING/HEATING EQUATIONS
c
c  written by: Yu Zhang, Peter Anninos and Tom Abel
c  date:       
c  modified1: January, 1996 by Greg Bryan; adapted to KRONOS
c  modified2: October, 1996 by GB; moved to AMR
c
c  PURPOSE:
c    Solve the energy cooling equations.
c
c  INPUTS:
c    is,ie   - start and end indicies of active region (zero-based!)
c
c  PARAMETERS:
c
c-----------------------------------------------------------------------
c
      implicit NONE
c
c  Arguments
c
      integer in, jn, kn, is, ie, j, k, nratec, imethod, idim,
     &        idual, iexpand, ih2co, ipiht, ispecies, imetal, 
     &        nfreq, iradshield, iradtype, imetalregen, iciecool,
     &        ih2optical, ifedtgas

      real    aye, temstart, temend,
     &        utem, uxyz, uaye, urho, utim,
     &        eta1, eta2, gamma
      real    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &       de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &      HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      real
     &        metal(in,jn,kn)
      real    hyd01ka(nratec), h2k01a(nratec), vibha(nratec), 
     &        rotha(nratec), rotla(nratec), gpldla(nratec),
     &        gphdla(nratec), hdltea(nratec), hdlowa(nratec),
     $        hdcoola(nratec, 5), ciecoa(nratec)
      real    ceHIa(nratec), ceHeIa(nratec), ceHeIIa(nratec),
     &        ciHIa(nratec), ciHeIa(nratec), ciHeISa(nratec), 
     &        ciHeIIa(nratec), reHIIa(nratec), reHeII1a(nratec), 
     &        reHeII2a(nratec), reHeIIIa(nratec), brema(nratec)
      real    compa, piHI, piHeI, piHeII, comp_xraya, comp_temp,
     &        inutot(nfreq), avgsighp, avgsighep, avgsighe2p,
     &        ciecont, mciecont, medot
c
c  Parameters
c
      double precision mh
      real    ZSOLAR
      parameter (mh = 1.67d-24, ZSOLAR = 0.02041)
c
c  Locals
c
      integer i, m, j1, iter, iradfield, mciei
      real dom, qq, vibl, logtem0, logtem9, dlogtem, energy, zr
      real dt2, ttmin, comp1, comp2, scoef, acoef, HIdot,
     &     hdlte1, hdlow1, gamma2(in), x, nH2, nother, fudge, fH2,
     &     gphdl1, factor, ciefudge, tau, comp_edot, comp_xray_edot
      double precision coolunit, dbase1, tbase1, xbase1
c
c  Slice locals
c 
      integer indixe(in)
      real t1(in), t2(in), logtem(in), tdef(in), p2d(in),
     &     tgas(in), tgasold(in)
      double precision edotplus(in), edotminus(in)
c
c  Cooling/heating slice locals
c
      double precision ceHI(in), ceHeI(in), ceHeII(in),
     &     ciHI(in), ciHeI(in), ciHeIS(in), ciHeII(in),
     &     reHII(in), reHeII1(in), reHeII2(in), reHeIII(in),
     &     brem(in)
      real hyd01k(in), h2k01(in), vibh(in), roth(in), rotl(in),
     &     gpldl(in), gphdl(in), hdlte(in), hdlow(in), 
     &     hdcool(in, 5), hdcoolval(in), cieco(in)
c

      integer nti, ndi, cmgenerate, NIT, NID, NIB
      real    TEMMIN, DELT, DENMIN, DELD, FREQDEL, FREQMIN
      PARAMETER(NIT=200,TEMMIN=3.0,DELT=0.03)
      PARAMETER(NID=300,DENMIN=-12.0,DELD=0.05)
      PARAMETER(NIB=400,FREQDEL=0.02,FREQMIN=1.0)
      real    metal_cool, metal_heat, xi
      real    cbovcool(NIT,NID), cbovheat(NIT,NID), denwk(NID), 
     &        radt(NIB), eb(NIB)
      common  /cen_metal_com/ cbovcool, cbovheat, denwk


c     Iteration mask
      
      logical itmask(in)

c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
c     Set log values of start and end of lookup tables
c
      mciecont = 0
      logtem0 = log(temstart)
      logtem9 = log(temend)
      dlogtem= (log(temend) - log(temstart))/real(nratec-1)
c
c     Set units
c
      dom      = urho*(aye**3)/mh
      tbase1   = utim
      xbase1   = uxyz/(aye*uaye)    ! uxyz is [x]*a      = [x]*[a]*a'        '
      dbase1   = urho*(aye*uaye)**3 ! urho is [dens]/a^3 = [dens]/([a]*a')^3 '
      coolunit = (uaye**5 * xbase1**2 * mh**2) / (tbase1**3 * dbase1)
c      write(6,*) "coolunit: ", coolunit
c      write(6,*) "dom:      ", dom
      zr       = 1.0/(aye*uaye) - 1.0
      fudge    = 1.0
      ciefudge = 1.0

      do i = is+1, ie+1
        edotminus(i) = 0.
        edotplus(i) = 0.
      enddo

c
c     Set compton cooling coefficients (and temperature)
c
      if (iexpand .eq. 1) then
         comp1 = compa * (1.0 + zr)**4
         comp2 = 2.73  * (1.0 + zr)
      else
         comp1 = 1.0e-30
         comp2 = 1.0e-30
      endif
c
c     Compute Pressure
c
      if (ifedtgas .eq. 0) then
      if (imethod .eq. 2) then
c
c        Zeus - e() is really gas energy
c
         do i = is+1, ie+1
            if (itmask(i)) then
               p2d(i) = (gamma - 1.0)*d(i,j,k)*e(i,j,k)
            endif
         enddo
      else
         if (idual .eq. 1) then
c
c           PPM with dual energy -- use gas energy
c
            do i = is+1, ie+1
               if (itmask(i)) then
                  p2d(i) = (gamma - 1.0)*d(i,j,k)*ge(i,j,k)
               endif
            enddo
         else
c
c           PPM without dual energy -- use total energy
c
            do i = is+1, ie+1
               if (itmask(i)) then
                  p2d(i) = e(i,j,k) - 0.5*u(i,j,k)**2
                  if (idim .gt. 1) p2d(i) = p2d(i) - 0.5*v(i,j,k)**2
                  if (idim .gt. 2) p2d(i) = p2d(i) - 0.5*w(i,j,k)**2
                  p2d(i) = max((gamma - 1.0)*d(i,j,k)*p2d(i), 1.0e-30)
               endif
            enddo
         endif
      endif
c
c     Compute temperature
c
      do i = is+1, ie+1
         if (itmask(i)) then
            tgas(i) = 
     &           (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))/4.0 +
     &           HI(i,j,k) + HII(i,j,k) + de(i,j,k)
         endif
      enddo
c
c          (include molecular hydrogen, but ignore deuterium)
c
      if (ispecies .gt. 1) then
         do i = is+1, ie+1
            if (itmask(i)) then
               tgas(i) = tgas(i) +
     &              HM(i,j,k) + (H2I(i,j,k) + H2II(i,j,k))/2.0
            endif
         enddo
      endif
c
        do i = is+1, ie+1
          if (itmask(i)) then
c             write(6,*), 'tgas ', tgas(i)
             tgas(i) = max(p2d(i)*(utem/tgas(i)), temstart)
             if(tgas(i).ne.tgas(i))then
                write(6,*) 'cool1d_multi tgas0a:',p2d(i),tgas(i)
                write(6,*) 'cool1d_multi tgas0b:',HI(i,j,k),HII(i,j,k),HM(i,j,k)
                write(6,*) 'cool1d_multi tgas0c:',H2I(i,j,k)
             endif
          endif
        enddo
c   Omukai
c
c     Correct temperature for gamma from H2
c
      if (ispecies .gt. 1) then
         do i = is+1, ie+1
            if (itmask(i)) then
               nH2 = 0.5*(H2I(i,j,k) + H2II(i,j,k))
               nother = (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))/4.0 +
     &              HI(i,j,k) + HII(i,j,k) + de(i,j,k) + HM(i,j,k)
               if (nH2/nother .gt. 1.0e-3) then
                  x = 6100/tgas(i) ! not quite self-consistent
                  if (x .gt. 10.0) then
                     gamma2(i) = 0.5*5.0
                  else
                    gamma2(i) = 0.5*(5.0 + 2.0*x**2 * exp(x)
     &                        /(exp(x)-1)**2)
                  endif
               else
                  gamma2(i) = 2.5
               endif
               gamma2(i) = 1.0 + (nH2 + nother)/
     &              (nH2*gamma2(i) + nother/(gamma-1.0))
               if (gamma2(i).ne.gamma2(i))then
                 write(6,*)'c1dm_bdf: g20a', gamma2(i),nH2
                 write(6,*)'c1dm_bdf: g20b', nother,tgas(i),ifedtgas
                 write(6,*)'c1dm_bdf: g20c', H2I(i,j,k), H2II(i,j,k)
                 write(6,*)'c1dm_bdf: g20d', HeI(i,j,k), HeII(i,j,k),
     &                                       HeIII(i,j,k)
                 write(6,*)'c1dm_bdf: g20e', HI(i,j,k), HII(i,j,k),
     &                                       de(i,j,k)
                 write(6,*)'c1dm_bdf: g20f', HM(i,j,k)
                 write(6,*)'c1dm_bdf: g20g', i,j,k,iter
               endif
               tgas(i) = tgas(i) * (gamma2(i) - 1.0)/(gamma - 1.0)
               if(tgas(i).ne.tgas(i))then
                  write(6,*) 'cool1d_multi tgas0a:',p2d(i),tgas(i)
                  write(6,*) 'cool1d_multi tgas0b:',HI(i,j,k),HII(i,j,k),HM(i,j,k)
                  write(6,*) 'cool1d_multi tgas0c:',H2I(i,j,k)
               endif
            endif
         enddo
      endif
c
c     If this is the first time through, just set tgasold to tgas
c
      if (iter .eq. 1) then
         do i = is+1, ie+1
            if (itmask(i)) then 
c               write(6,*) i, tgas(i)
               tgasold(i) = tgas(i)
            endif
         enddo
      endif
      endif
c
c     --- 6 species cooling ---
c
      do i = is+1, ie+1

         if (itmask(i)) then

c
c        Compute log temperature and truncate if above/below table max/min
c
         logtem(i) = log(0.5*(tgas(i)+tgasold(i)))
         logtem(i) = max(logtem(i), logtem0)
         logtem(i) = min(logtem(i), logtem9)
c
c        Compute index into the table and precompute parts of linear interp
c
         indixe(i) = min(nratec-1,
     .                  max(1,int((logtem(i)-logtem0)/dlogtem)+1))
         t1(i) = (logtem0 + (indixe(i) - 1)*dlogtem)
         t2(i) = (logtem0 + (indixe(i)    )*dlogtem)
         tdef(i) = t2(i) - t1(i)
c
c        Lookup cooling values and do a linear temperature in log(T)
c
         ceHI(i) = ceHIa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ceHIa(indixe(i)+1) -ceHIa(indixe(i)))/tdef(i)
         ceHeI(i) = ceHeIa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ceHeIa(indixe(i)+1) -ceHeIa(indixe(i)))/tdef(i)
         ceHeII(i) = ceHeIIa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ceHeIIa(indixe(i)+1) -ceHeIIa(indixe(i)))/tdef(i)
         ciHI(i) = ciHIa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ciHIa(indixe(i)+1) -ciHIa(indixe(i)))/tdef(i)
         ciHeI(i) = ciHeIa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ciHeIa(indixe(i)+1) -ciHeIa(indixe(i)))/tdef(i)
         ciHeIS(i) = ciHeISa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ciHeISa(indixe(i)+1) -ciHeISa(indixe(i)))/tdef(i)
         ciHeII(i) = ciHeIIa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ciHeIIa(indixe(i)+1) -ciHeIIa(indixe(i)))/tdef(i)
         reHII(i) = reHIIa(indixe(i)) + (logtem(i) - t1(i))
     .         *(reHIIa(indixe(i)+1) -reHIIa(indixe(i)))/tdef(i)
         reHeII1(i)=reHeII1a(indixe(i)) + (logtem(i) - t1(i))
     .        *(reHeII1a(indixe(i)+1)-reHeII1a(indixe(i)))/tdef(i)
         reHeII2(i)=reHeII2a(indixe(i)) + (logtem(i) - t1(i))
     .        *(reHeII2a(indixe(i)+1)-reHeII2a(indixe(i)))/tdef(i)
         reHeIII(i)=reHeIIIa(indixe(i)) + (logtem(i) - t1(i))
     .        *(reHeIIIa(indixe(i)+1)-reHeIIIa(indixe(i)))/tdef(i)
         brem(i) = brema(indixe(i)) + (logtem(i) - t1(i))
     .         *(brema(indixe(i)+1) -brema(indixe(i)))/tdef(i)

      endif
      enddo
c
c     Compute the cooling function
c
      do i = is+1, ie+1
         if (itmask(i)) then
         
         edotminus(i) = abs(
c
c                    Collisional excitations
c
     .             - (ceHI  (i)*HI  (i,j,k))*de(i,j,k)           ! ce of HI
     .             - ((ceHeI (i)*HeII(i,j,k))*de(i,j,k))
     .                  *de(i,j,k)*dom/4.0                       ! ce of HeI
     .             - (ceHeII(i)*HeII(i,j,k))*de(i,j,k)/4.0         ! ce of HeII
c
c                    Collisional ionizations
c
     .             - (ciHI  (i)*HI  (i,j,k))*de(i,j,k)             ! ci of HI
     .             - (ciHeI (i)*HeI (i,j,k))*de(i,j,k)/4.0         ! ci of HeI
     .             - (ciHeII(i)*HeII(i,j,k))*de(i,j,k)/4.0         ! ci of HeII
     .             - ((ciHeIS(i)*HeII(i,j,k))*de(i,j,k))
     .                  *de(i,j,k)*dom/4.0  ! ci of HeIS
c
c                    Recombinations
c
     .             - (reHII  (i)*HII  (i,j,k))*de(i,j,k)          ! re of HII
     .             - (reHeII1(i)*HeII (i,j,k))*de(i,j,k)/4.0      ! re of HeII
     .             - (reHeII2(i)*HeII (i,j,k))*de(i,j,k)/4.0      ! re of HeII
     .             - (reHeIII(i)*HeIII(i,j,k))*de(i,j,k)/4.0      ! re of HeIII
c
c                    Bremsstrahlung
c
     .             - (brem(i)*(HII(i,j,k)+HeII(i,j,k)/4.0+HeIII(i,j,k)))
     .                        *de(i,j,k)
     .                 ) 
          
c
c                    Compton cooling or heating
c
           comp_edot = - comp1*(tgas(i)-comp2)*de(i,j,k)/dom
           if (comp_edot.lt.0.0d0) then
             edotminus(i) = edotminus(i) +
     .            abs(comp_edot)
           else
             edotplus(i) = edotplus(i) +
     .            abs(comp_edot)
           endif
c
c                    X-ray compton heating
c
           comp_xray_edot=-comp_xraya*(tgas(i)-comp_temp)*de(i,j,k)/dom
           if (comp_xray_edot.lt.0.0d0) then
             edotminus(i) = edotminus(i) +
     .               abs(comp_xray_edot)
           else
             edotplus(i) = edotplus(i) +
     .               abs(comp_xray_edot)
           endif
         endif
      enddo
c     
c     --- H2 cooling ---
c
      if (ispecies .gt. 1) then
c


c        Use the Galli and Palla (1999) cooling rates for molecular H.
c
# 453 "cool1d_multi_bdf.src"
cc
         do i = is+1, ie+1
            if (itmask(i)) then
            hyd01k(i) = hyd01ka(indixe(i)) + (logtem(i) - t1(i))
     .         *(hyd01ka(indixe(i)+1)-hyd01ka(indixe(i)))/tdef(i)
            h2k01(i) = h2k01a(indixe(i)) + (logtem(i) - t1(i))
     .         *(h2k01a(indixe(i)+1) - h2k01a(indixe(i)))/tdef(i)
            vibh(i) = vibha(indixe(i)) + (logtem(i) - t1(i))
     .         *(vibha(indixe(i)+1) - vibha(indixe(i)))/tdef(i)
            roth(i) = rotha(indixe(i)) + (logtem(i) - t1(i))
     .         *(rotha(indixe(i)+1) - rotha(indixe(i)))/tdef(i)
            rotl(i) = rotla(indixe(i)) + (logtem(i) - t1(i))
     .         *(rotla(indixe(i)+1) - rotla(indixe(i)))/tdef(i)
            cieco(i) = ciecoa(indixe(i)) + (logtem(i) - t1(i))
     .         *(ciecoa(indixe(i)+1) - ciecoa(indixe(i)))/tdef(i)
            endif
         enddo
c
         do i = is+1, ie+1
            if (itmask(i)) then
            qq   = 1.2*(HI(i,j,k)*dom)**0.77 + 
     .                (H2I(i,j,k)*dom/2.)**0.77
            vibl = (HI(i,j,k)*hyd01k(i) + 
     .             H2I(i,j,k)/2.*h2k01(i))
     .             *dom*8.18e-13
c
            if (ih2optical .eq. 1) then
c            nH2 = 0.5*H2I(i,j,k)
c            nother = (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))/4.0 +
c     .                HI(i,j,k) + HII(i,j,k) + de(i,j,k)
c            fH2 = nH2/(nH2 + nother)
c            fudge = sqrt((40.0 * 10**(4.8 * 
c     .               sqrt(max(log10(tgas(i)),2.0)-2.0)) / fH2**2)/
c     .               ((nH2 + nother)*dom) )
              fudge = (0.76*d(i,j,k)*dom/8.0e9)**(-0.45)
              fudge = min(fudge, 1.0)
            endif
c
            edotminus(i)=edotminus(i) + float(ih2co)*fudge*H2I(i,j,k)*(
     .          vibh(i)/(1.0+vibh(i)/max(   vibl     ,1.0e-30)) +
     .          roth(i)/(1.0+roth(i)/max(qq*rotl(i),1.0e-30))     
     .                                                      )/2.0/dom
            endif
         enddo
c

c

c     CIE
c     cooling from H2-H2 and He-H2 collisional induced emission comes
C     With its own radiative transfer correction as discussed in
C     Ripamonti & Abel 2003
         if (iciecool .eq. 1) then
            do i = is+1, ie+1
            if (itmask(i)) then
c     Only calculate if H2I(i,j,k) is a substantial fraction
              if ((H2I(i,j,k).gt.0.2*d(i,j,k)).and.
     &            (d(i,j,k)*dom.gt.1e10)) then
                ciefudge = 1.
                tau = ((d(i,j,k)/2e16)*dom)**2.8  ! 2e16 is in units of cm^-3
                tau = max(tau, 1.e-5)
                ciefudge = min((1.-exp(-tau))/tau,1.0)
c                ciecont = H2I(i,j,k)*(d(i,j,k)*cieco(i))*ciefudge
                edotminus(i) = edotminus(i) + H2I(i,j,k)*
     &                     (d(i,j,k)*cieco(i))*ciefudge
              endif
            endif
            enddo
         endif

      endif                     !  ispecies = 2
c
c     --- Cooling from HD ---
c
      if (ispecies .gt. 2) then
         do i = is+1, ie+1
            if (itmask(i)) then
c            hdlte(i) = hdltea(indixe(i)) + (logtem(i) - t1(i))
c     .         *(hdltea(indixe(i)+1) - hdltea(indixe(i)))/tdef(i)
c            hdlow(i) = hdlowa(indixe(i)) + (logtem(i) - t1(i))
c     .         *(hdlowa(indixe(i)+1) - hdlowa(indixe(i)))/tdef(i)
            do m = 1, 5
              hdcool(i,m) = hdcoola(indixe(i),m) + (logtem(i) - t1(i))
     .           *(hdcoola(indixe(i)+1,m) - hdcoola(indixe(i),m))/tdef(i)
            enddo
            endif
         enddo
c
         do i = is+1, ie+1
            if (itmask(i)) then
c               hdlte1 = hdlte(i)/(HDI(i,j,k)*dom/2.0)
c               hdlow1 = max(hdlow(i), 1.0e-30)
c               edotminus(i) = edotminus(i) + HDI(i,j,k)*
c     .                     (hdlte1/(1.0 + hdlte1/hdlow1)/(2.0*dom))
c   Do we want to consider HII as well?
c   Note that, as per Flower 2000, we do not include H2 collisions
                hdcoolval(i) = 0.0e0
                do m = 1,5
                  hdcoolval(i) = (hdcool(i,m)**(m-1))
     .                          *((HI(i,j,k)+HII(i,j,k))*dom)
                enddo
                edotminus(i) = edotminus(i)
c     .                       + HDI(i,j,k)*(hdcoolval(i)/coolunit)
            endif
         enddo
      endif
c
c     --- Compute (external) radiative heating terms ---
c
c                       Photoionization heating
c

      if (iradshield .eq. 0) then
c
c        regular version
c
         if (iradtype .eq. 8) then
c
c           1) heating assuming high energy photons produces secondary
c              electrons which do the heating (Shull & Steenberg, 1985).
c
            do i = is+1, ie+1
               if (itmask(i)) then
               x = max(HII(i,j,k)/(HI(i,j,k)+HII(i,j,k)), 1.0e-4)
               factor = 0.9971*(1.0-(1.0-x**0.2663)**1.3163)
               edotplus(i) = edotplus(i) + float(ipiht)*factor*(
     .                + piHI  *HI  (i,j,k)         ! pi of HI
     .                + piHeI *HeI (i,j,k)*0.25     ! pi of HeI
     .                + piHeII*HeII(i,j,k)*0.25     ! pi of HeII
     .              )/dom
               endif
            enddo
         else
c
c           2) standard heating
c
            do i = is+1, ie+1
               if (itmask(i)) then
               edotplus(i) = edotplus(i) + float(ipiht)*(
     .                + piHI  *HI  (i,j,k)         ! pi of HI
     .                + piHeI *HeI (i,j,k)*0.25     ! pi of HeI
     .                + piHeII*HeII(i,j,k)*0.25     ! pi of HeII
     .              )/dom
               endif
            enddo
         endif
      else
c
c        version with approximate self-shielding
c
         do i = is+1, ie+1
            if (itmask(i)) then
            edotplus(i) = edotplus(i) + float(ipiht)*(
     .                + piHI  *HI  (i,j,k)*
     .                   exp(-avgsighp*HI(i,j,k)*dom)
     .                + piHeI *HeI (i,j,k)*0.25*
     .                   exp(-avgsighep*0.25*HeI(i,j,k)*dom)
     .                + piHeII*HeII(i,j,k)*0.25*
     .                   exp(-avgsighe2p*0.25*HeII(i,j,k)*dom)
     .           )/dom
            endif
         enddo
      endif
c

c

c
c     Set tgasold
c
      if (ifedtgas.eq.0) then
        do i=is+1, ie+1
           if (itmask(i)) then
              tgasold(i) = tgas(i)
           endif
        enddo
      endif
c
      return
      end

