# 1 "solve_rate_cool.src"
# 1 "<built-in>"
# 1 "<command line>"
# 1 "solve_rate_cool.src"

# 1 "fortran.def" 1




















# 2 "solve_rate_cool.src" 2


c=======================================================================
c/////////////////////  SUBROUTINE SOLVE_RATE  \\\\\\\\\\\\\\\\\\\\\\\\c

      subroutine solve_rate_cool(d, e, ge, u, v, w, de,
     &                HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, iexpand, imethod,
     &                idual, ispecies, imetal, idim,
     &                is, js, ks, ie, je, ke, ih2co, ipiht,
     &                dt, aye, temstart, temend, 
     &                utem, uxyz, uaye, urho, utim,
     &                eta1, eta2, gamma, fh, dtoh,
     &                k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &                k11a, k12a, k13a, k13dda, k14a, k15a,
     &                k16a, k17a, k18a, k19a, k22a, k23a,
     &                k24, k25, k26, k27, k28, k29, k30, k31,
     &                k50a, k51a, k52a, k53a, k54a, k55a, k56a,
     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &                reHeII2a, reHeIIIa, brema, compa,
     &                comp_xraya, comp_temp, piHI, piHeI, piHeII,
     &                HM, H2I, H2II, DI, DII, HDI, metal,
     &                hyd01ka, h2k01a, vibha, rotha, rotla, 
     &                gpldla, gphdla, hdltea, hdlowa, hdcoola, ciecoa, 
     &                inutot, iradtype, nfreq, imetalregen,
     &                iradshield, avgsighp, avgsighep, avgsighe2p,
     &                iciecool, ih2optical, errcode, omaskflag )
c
c  SOLVE MULTI-SPECIES RATE EQUATIONS AND RADIATIVE COOLING
c
c  written by: Yu Zhang, Peter Anninos and Tom Abel
c  date:       
c  modified1:  January, 1996 by Greg Bryan; converted to KRONOS
c  modified2:  October, 1996 by GB; adapted to AMR
c  modified3:  May,     1999 by GB and Tom Abel, 3bodyH2, solver, HD
c  modified4:  June,    2005 by GB to solve rate & cool at same time
c  modified5:  May,     2006 by Matt Turk to integrate temperature, not energy
c
c  PURPOSE:
c    Solve the multi-species rate and cool equations.
c
c  INPUTS:
c    in,jn,kn - dimensions of 3D fields
c
c    d        - total density field
c    de       - electron density field
c    HI,HII   - H density fields (neutral & ionized)
c    HeI/II/III - He density fields
c    DI/II    - D density fields (neutral & ionized)
c    HDI      - neutral HD molecule density field
c    HM       - H- density field
c    H2I      - H_2 (molecular H) density field
c    H2II     - H_2+ density field
c
c    is,ie    - start and end indices of active region (zero based)
c    idual    - dual energy formalism flag (0 = off, 1 = on)
c    iexpand  - comoving coordinates flag (0 = off, 1 = on)
c    idim     - dimensionality (rank) of problem
c    ispecies - chemistry module (1 - H/He only, 2 - molecular H, 3 - D) 
c    iradshield - flag for crude radiative shielding correction (rad type 12)
c    iradtype - type of radiative field (only used if = 8)
c    imetal   - flag if metal field is active (0 = no, 1 = yes)
c    imethod  - Hydro method (0 = PPMDE, 2 = ZEUS-type)
c    ih2co    - flag to include H2 cooling (1 = on, 0 = off)
c    ipiht    - flag to include photoionization heating (1 = on, 0 = off)
c    iciecool - flag to include CIE cooling (1 = on, 0 = off)
c    ih2optical - flag to include H2 optical depth fudge (1 = on, 0 = off)
c
c    fh       - Hydrogen mass fraction (typically 0.76)
c    dtoh     - Deuterium to H mass ratio
c    dt       - timestep to integrate over
c    aye      - expansion factor (in code units)
c
c    utim     - time units (i.e. code units to CGS conversion factor)
c    uaye     - expansion factor conversion factor (uaye = 1/(1+zinit))
c    urho     - density units
c    uxyz     - length units
c    utem     - temperature(-like) units
c
c    temstart, temend - start and end of temperature range for rate table
c    nratec   - dimensions of chemical rate arrays (functions of temperature)
c
c  OUTPUTS:
c    update chemical rate densities (HI, HII, etc)
c
c  PARAMETERS:
c    itmax   - maximum allowed sub-cycle iterations
c    mh      - H mass in cgs units
c
c-----------------------------------------------------------------------
c
      implicit NONE
c
c  General Arguments
c
      integer in, jn, kn, is, js, ks, ie, je, ke, nratec, imethod,
     &        idual, iexpand, ih2co, ipiht, ispecies, imetal, idim,
     &        iradtype, nfreq, imetalregen, iradshield, iciecool, 
     &        ih2optical, errcode, omaskflag
      real    dt, aye, temstart, temend, eta1, eta2, gamma,
     &        utim, uxyz, uaye, urho, utem, fh, dtoh
c
c  Density, energy and velocity fields fields
c
      real    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      real    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      real    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &        metal(in,jn,kn)
c
c  Cooling tables (coolings rates as a function of temperature)
c
      real    hyd01ka(nratec), h2k01a(nratec), vibha(nratec), 
     &        rotha(nratec), rotla(nratec), gpldla(nratec),
     &        gphdla(nratec), hdltea(nratec), hdlowa(nratec),
     $        hdcoola(nratec, 5), ciecoa(nratec)
      real    ceHIa(nratec), ceHeIa(nratec), ceHeIIa(nratec),
     &        ciHIa(nratec), ciHeIa(nratec), ciHeISa(nratec), 
     &        ciHeIIa(nratec), reHIIa(nratec), reHeII1a(nratec), 
     &        reHeII2a(nratec), reHeIIIa(nratec), brema(nratec)
      real    compa, piHI, piHeI, piHeII, comp_xraya, comp_temp,
     &        inutot(nfreq), avgsighp, avgsighep, avgsighe2p
c
c  Chemistry tables (rates as a function of temperature)
c
      real k1a (nratec), k2a (nratec), k3a (nratec), k4a (nratec), 
     &     k5a (nratec), k6a (nratec), k7a (nratec), k8a (nratec), 
     &     k9a (nratec), k10a(nratec), k11a(nratec), k12a(nratec), 
     &     k13a(nratec), k14a(nratec), k15a(nratec), k16a(nratec), 
     &     k17a(nratec), k18a(nratec), k19a(nratec), k22a(nratec),
     &     k23a(nratec), k50a(nratec), k51a(nratec), k52a(nratec), 
     &     k53a(nratec), k54a(nratec), k55a(nratec), k56a(nratec),
     &     k13dda(nratec, 7),
     &     k24, k25, k26, k27, k28, k29, k30, k31
c
c  Parameters
c
      integer itmax, ijk
      parameter (itmax = 3000, ijk = 1031)
      double precision mh
      parameter (mh = 1.67d-24)
c
c  Locals
c
      integer i, j, k, iter, imin, ifedtgas
      real ttmin, dom, energy, comp1, comp2
      double precision coolunit, dbase1, tbase1, xbase1, kunit,
     &                 kunit_3bdy, chunit
c
c  row temporaries
c 
      integer indixe(ijk)
      real t1(ijk), t2(ijk), logtem(ijk), tdef(ijk),
     &     dtit(ijk), ttot(ijk), p2d(ijk), tgas(ijk), tgasold(ijk),
     &     tdotplus(ijk), tdotminus(ijk), tdot(ijk), olddtit(ijk),
     &     tgasolder(ijk), olderdtit(ijk), usedtit(ijk)
c
c  Rate equation row temporaries
c
      real*8 totalN(ijk), mu(ijk), totalMass(ijk)
      real gamma2(ijk), temp, tto
      real*8 dedot_prev(ijk), HIdot_prev(ijk), H2Idot_prev(ijk),
     &       HDIdot_prev(ijk), DIdot_prev(ijk),
     &     HIdot(ijk), dedot(ijk), H2Idot(ijk), HDIdot(ijk), DIdot(ijk),
     &     mudot(ijk), chdot(ijk),
     &     dedotplus(ijk), dedotminus(ijk), 
     &     HIdotplus(ijk), HIdotminus(ijk), 
     &     H2Idotplus(ijk), H2Idotminus(ijk),
     &     HDIdotplus(ijk), HDIdotminus(ijk),
     &     DIdotplus(ijk), DIdotminus(ijk)
      real k24shield(ijk), k25shield(ijk), k26shield(ijk)
      real*8 HIp(ijk), HIIp(ijk), HeIp(ijk), HeIIp(ijk), HeIIIp(ijk),
     &     HMp(ijk), H2Ip(ijk), H2IIp(ijk), dep(ijk),
     &     DIp(ijk), DIIp(ijk), HDIp(ijk)
      real k1 (ijk), k2 (ijk), k3 (ijk), k4 (ijk), k5 (ijk),
     &     k6 (ijk), k7 (ijk), k8 (ijk), k9 (ijk), k10(ijk),
     &     k11(ijk), k12(ijk), k13(ijk), k14(ijk), k15(ijk),
     &     k16(ijk), k17(ijk), k18(ijk), k19(ijk), k22(ijk),
     &     k23(ijk), k50(ijk), k51(ijk), k52(ijk), k53(ijk), 
     &     k54(ijk), k55(ijk), k56(ijk), k13dd(ijk, 7)
c
c  Cooling/heating row locals
c
      double precision ceHI(ijk), ceHeI(ijk), ceHeII(ijk),
     &     ciHI(ijk), ciHeI(ijk), ciHeIS(ijk), ciHeII(ijk),
     &     reHII(ijk), reHeII1(ijk), reHeII2(ijk), reHeIII(ijk),
     &     brem(ijk), edot(ijk), 
     &     edotplus(ijk), edotminus(ijk)
      real hyd01k(ijk), h2k01(ijk), vibh(ijk), roth(ijk), rotl(ijk),
     &     gpldl(ijk), gphdl(ijk), hdlte(ijk), hdlow(ijk), 
     &     hdcool(ijk, 5), cieco(ijk)

      real toler
      parameter (toler = 1e-20 )

c
c     Iteration mask, output mask, and output flag
c     
      logical itmask(ijk), omask(ijk), eqmask(ijk), timemask(ijk),
     &        omaskflagl

      double precision everg, e24, e26
      parameter(everg = 1.60184d-12, e24 = 13.6d0, e26 = 24.6d0)

c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
c     Set units
c
      dom      = urho*(aye**3)/mh
      tbase1   = utim
      xbase1   = uxyz/(aye*uaye)    ! uxyz is [x]*a      = [x]*[a]*a'        '
      dbase1   = urho*(aye*uaye)**3 ! urho is [dens]/a^3 = [dens]/([a]*a')^3 '
      coolunit = (uaye**5 * xbase1**2 * mh**2) / (tbase1**3 * dbase1)
      kunit   = (uaye**3 * mh) / (dbase1 * tbase1)
      kunit_3bdy  = kunit * (uaye**3 * mh) / dbase1
      chunit = (coolunit * utim/dom) / everg  ! divide eV by chunit for units

c      write(6,*) 'kunit = ', kunit

c
c  Convert densities from comoving to proper
c
      call scale_fields(d, de, HI, HII, HeI, HeII, HeIII,
     &                  HM, H2I, H2II, DI, DII, HDI, metal,
     &                  is, ie, js, je, ks, ke,
     &                  in, jn, kn, ispecies, imetal, aye**(-3))
c
c  Loop over zones, and do an entire i-column in one go
c
c      tgas(ijk+1) = 0
      ifedtgas = 1
      if (omaskflag.eq.1) then
        write(6,*) 'setting omask to true'
        omaskflagl = .true.
      endif
      do k = ks+1, ke+1
      do j = js+1, je+1
c
c        Initialize variables
c
         do i = is+1, ie+1
            dtit(i) = 0.0
            ttot(i) = 0.0
            itmask(i) = .true.
            omask(i) = omaskflagl
            eqmask(i) = .false.
            timemask(i) = .false.
            tgas(i) = 0.
            tgasold(i) = 0.
            tdotplus(i) = 0.
            tdotminus(i)=0.
            tdot(i)=1.0e-30
            H2Idot_prev(i) = 0.
            HIdot_prev(i) = 0.
            dedot_prev(i) = 0.
            H2Idot(i) = 0.
            HDIdot(i) = 0.
            DIdot(i) = 0.
            HIdot(i) = 0.
            dedot(i) = 0.
            mudot(i) = 0.
            iter = 0
            gamma2(i) = 0.
            dedot_prev(i) = 0.
            HIdot_prev(i) = 0.
            H2Idot_prev(i) = 0.
            HDIdot_prev(i) = 0.
            DIdot_prev(i) = 0.
            HIdot(i) = 0.
            dedot(i) = 0.
            H2Idot(i) = 0.
            mudot(i) = 0.
            dedotplus(i) = 0.
            dedotminus(i) = 0.
            HIdotplus(i) = 0.
            HIdotminus(i) = 0.
            H2Idotplus(i) = 0.
            H2Idotminus(i) = 0.
            HDIdotplus(i) = 0.
            HDIdotminus(i) = 0.
            DIdotplus(i) = 0.
            DIdotminus(i) = 0.
            totalN(i) = 0.
            totalMass(i) = 0.
         enddo
c
c        We calculate an initial mu, to save time
c
         call calculate_mu(
     &                d, e, ge, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, is, ie, j, k,
     &                HM, H2I, H2II, DI, DII, HDI, ispecies,
     &                totalMass, totalN, mu, mudot, dtit, itmask )
c
c        We calculate an initial temperature, as enzo still globally tracks
c        energy, not temperature.
c
         call calculate_tgas(
     &                d, e, ge, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, idual, imethod, ispecies, idim,
     &                is, ie, j, k, aye, temstart, temend,
     &                utem, uxyz, uaye, urho, utim, gamma,
     &                HM, H2I, H2II, DI, DII, HDI, 
     &                tgas, tgasold, p2d, gamma2, totalN, itmask )

          do i = is+1, ie+1
            tgasolder(i) = tgas(i)
          enddo

c
c        ------------------ Loop over subcycles ----------------
c
         do iter = 1, itmax
c
c           Compute the cooling rate for this row
c
            call cool1d_multi_bdf(
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
c        Look-up rates as a function of temperature for 1D set of zones
c
            call lookup_cool_rates1d(temstart, temend, nratec, j, k,
     &                is, ie, ijk, iradtype, iradshield, in, jn, kn,
     &                ispecies, tgas, tgasold, HI, HII, HeI, HeII,
     &                k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &                k11a, k12a, k13a, k13dda, k14a, k15a, k16a,
     &                k17a, k18a, k19a, k22a, k23a,
     &                k50a, k51a, k52a, k53a, k54a, k55a, k56a,
     &                avgsighp, avgsighep, avgsighe2p, piHI, piHeI,
     &                k1, k2, k3, k4, k5, k6, k7, k8, k9, k10,
     &                k11, k12, k13, k14, k15, k16, k17, k18,
     &                k19, k22, k24, k25, k26,
     &                k50, k51, k52, k53, k54, k55,
     &                k56, k13dd, k24shield, k25shield, k26shield,
     &                t1, t2, tdef, logtem, indixe, 
     &                dom, coolunit, tbase1, itmask )
c
c           Compute dedot, HIdot and H2Idot, the rates of change of 
c           de, HI and HII
c
            call rate_timestep(dedot, HIdot, H2Idot, HDIdot, DIdot,
     &                     ispecies,
     &                     de, HI, HII, HeI, HeII, HeIII, d,
     &                     HM, H2I, H2II, DI, DII, HDI,
     &                     in, jn, kn, is, ie, j, k, 
     &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
     &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
     &                     k23, k24, k25, k26, k27, k28, k29, k30, k31,
     &                     k50, k51, k52, k53, k54, k55,
     &                     k56, k24shield, k25shield, k26shield,
     &                     dedotplus,  dedotminus,
     &                     HIdotplus,  HIdotminus,
     &                     H2Idotplus, H2Idotminus,
     &                     HDIdotplus, HDIdotminus, 
     &                     DIdotplus, DIdotminus, 
     &                     dtit, iter,
     &                     dedot_prev, HIdot_prev, H2Idot_prev,
     &                     HDIdot_prev, DIdot_prev,
     &                     itmask )
c
c           Find timestep that keeps relative chemical changes below 10%
c
            do i = is+1, ie+1
            if (itmask(i).eqv..true.) then

c
c              Bound from below to prevent numerical errors
c
	       if (abs(dedot(i)) .lt. 1.0e-30*1d-10) 
     &             dedot(i) = min(1.0e-30*1d-10,de(i,j,k))
	       if (abs(HIdot(i)) .lt. 1.0e-30*1d-10)
     &             HIdot(i) = min(1.0e-30*1d-10,HI(i,j,k))
	       if (abs(H2Idot(i)) .lt. 1.0e-30*1d-10)
     &             H2Idot(i) = min(1.0e-30*1d-10,H2I(i,j,k))
c
c              If the net HI ionization rate is almost perfectly balanced then set
c                  it to zero (since it is zero to available precision)
c
               if (min(abs(k1(i)* de(i,j,k)*HI(i,j,k)),
     &                 abs(k2(i)*HII(i,j,k)*de(i,j,k)))/
     &             max(abs(dedot(i)),abs(HIdot(i)))
     &               .gt. 1.0e6) then
                  dedot(i) = 1.0e-30
                  HIdot(i) = 1.0e-30
               endif
               if (ispecies.gt.1) then
                 if (min(
     &                abs(2.*DBLE(k22(i))*DBLE(HI(i,j,k))*HI(i,j,k)**2),
     &                abs(2.*k13(i)*HI(i,j,k)*H2I(i,j,k)/2.))/
     &               max(abs(H2Idot(i)),abs(HIdot(i)))
     &                 .gt. 1.0e6) then
                    H2Idot(i) = 1.0e-30
                    HIdot(i) = 1.0e-30
                 endif
               endif
               if (timemask(i).eqv..true.) then
                 H2Idot(i) = H2Idot_prev(i)
               endif
c
c              compute minimum rate timestep

               if (iter.eq.1) olddtit(i)=1e30
               dtit(i) = min(0.1*max(abs(de(i,j,k)/dedot(i)),
     &                               abs(HI(i,j,k)/HIdot(i))),
     &                       0.1*abs(H2I(i,j,k)/H2Idot(i)),
     &                       dt-ttot(i), 0.5*dt, 10*olddtit(i))
               if (ispecies.gt.2) then
                 dtit(i) = min(dtit(i), 0.1*abs(HDI(i,j,k)/HDIdot(i)))
c     &                         0.1*abs(DI(i,j,k)/DIdot(i)))
               endif
c              if dtit(i) is ridiculously small, set into equilibrium
               if ((dtit(i)/dt.lt.1d-8.and.d(i,j,k)*dom.gt.1e10
     &             .and.tgas(i).gt.1000) 
     &             .or. timemask(i).eqv..true.) then
c                 write(6,*) 'putting into equilibrium, bitches', i,j,k
c                 eqmask(i) = .true.
c                 timemask(i) = .true.
c                 dtit(i) = 1.0e+30
               else
c                 write(6,*) 'taking out of equilibrium, bitches', i,j,k
c                 eqmask(i) = .false.
c                 timemask(i) = .false.
               endif
               olderdtit(i) = olddtit(i)
               olddtit(i) = dtit(i)
c
c              Output some debugging information if required
c
               if (dtit(i)/dt .lt. 1.0e-2 .and. iter .gt. 300 .and.
     &             abs((dt-ttot(i))/dt) .gt. 1.0e-3) then
                  write(4,1000) iter,i,j,k,dtit(i),
     &              ttot(i),dt,de(i,j,k),dedot(i),HI(i,j,k),HIdot(i),
     &              H2Idot(i),tgas(i), dedot_prev(i)
                  write(4,1100) HI(i,j,k),HII(i,j,k),
     &              HeI(i,j,k),HeII(i,j,k),HeIII(i,j,k),
     &              HM(i,j,k),H2I(i,j,k),H2II(i,j,k),de(i,j,k)
                  write(4,1100)
     &               -    k1(i) *de(i,j,k)    *HI(i,j,k)  ,
     &               -    k7(i) *de(i,j,k)    *HI(i,j,k),
     &               -    k8(i) *HM(i,j,k)    *HI(i,j,k),
     &               -    k9(i) *HII(i,j,k)   *HI(i,j,k),
     &               -    k10(i)*H2II(i,j,k)  *HI(i,j,k)/2.,
     &               - 2.*k22(i)*HI(i,j,k)**2 *HI(i,j,k),
     &               +    k2(i) *HII(i,j,k)   *de(i,j,k) ,
     &               + 2.*k13(i)*HI(i,j,k)    *H2I(i,j,k)/2.,
     &               +    k11(i)*HII(i,j,k)   *H2I(i,j,k)/2.,
     &               + 2.*k12(i)*de(i,j,k)    *H2I(i,j,k)/2.,
     &               +    k14(i)*HM(i,j,k)    *de(i,j,k),
     &               +    k15(i)*HM(i,j,k)    *HI(i,j,k),
     &               + 2.*k16(i)*HM(i,j,k)    *HII(i,j,k),
     &               + 2.*k18(i)*H2II(i,j,k)  *de(i,j,k)/2.,
     &               +    k19(i)*H2II(i,j,k)  *HM(i,j,k)/2.
               endif
 1000          format(i5,3(i3,1x),1p,9(e11.3))
 1100          format(20(e11.3))
c
            endif
            enddo   ! end loop over i
c
c           Compute maximum timestep for cooling/heating
c
            do i = is+1, ie+1
            if (itmask(i).eqv..true.) then
c
c              Set energy per unit volume of this cell based in the pressure
c              (the gamma used here is the right one even for H2 since p2d 
c               is calculated with this gamma).
c
               if (imethod .eq. 2) then
c
c        Zeus - e() is really gas energy
c
                  energy = d(i,j,k)*e(i,j,k)
               else
                  if (idual .eq. 1) then
c
c           PPM with dual energy -- use gas energy
c
                    energy = d(i,j,k)*ge(i,j,k)
                  else
c
c           PPM without dual energy -- use total energy
c
                    energy = e(i,j,k) - 0.5*u(i,j,k)**2
                    if (idim .gt. 1) energy = energy - 0.5*v(i,j,k)**2
                    if (idim .gt. 2) energy = energy - 0.5*w(i,j,k)**2
                    energy = d(i,j,k)*energy
                  endif
               endif
               energy = max(energy, 1.0e-30)
               edot(i) = edotplus(i) - edotminus(i)
c
c              If the temperature is at the bottom of the temperature look-up 
c              table and edot < 0, then shut off the cooling.
c
               if (tgas(i) .le. temstart .and. edot(i) .lt. 0.0) then
                    edot(i) = 1.0e-30
                    edotminus(i) = 1.0e-30
                    edotplus(i) = 1.0e-30
                endif
                    
	       if (abs(edot(i)) .lt. 1.0e-30) edot(i) = 1.0e-30
c
c              Compute timestep for 10% change
c
c              Add a timestep based on chemical heating
                  if(ispecies.gt.1.and.d(i,j,k)*dom.gt.1e7) then
                    chdot(i) = H2Idot(i)*4.48/chunit
                  else
                    chdot(i) = 0
                  endif
c                  if (abs(tgas(i)-tgasold(i))/tgas(i).lt.1e-3) then
c                    tdot(i) = 1.0e-30*1d-10
c                  endif
c
c                 Even if eqmask is set, we timestep on cooling
c
                  dtit(i) = min(
     &                        abs(0.1*energy/edot(i)), 
     &                        dt-ttot(i), dtit(i),
     &                        abs(0.001*tgas(i)/tdot(i)),
     &                        abs(0.01*energy/chdot(i)))
               if (mod(iter,10**4).eq.0.and.itmask(i).eqv..true.) then
                 write(6,*), i,j,k,iter,d(i,j,k),H2I(i,j,k),tgas(i),
     &                       dt,dtit(i),(dt-ttot(i))/dtit(i),
     &                       0.1*H2I(i,j,k)/H2Idot(i),
     &                       0.1*H2I(i,j,k)/H2Idot_prev(i),
     &                       H2Idot_prev(i), H2Idot(i)
c     &                       0.1*energy/edot(i),
c     &                       0.1*HI(i,j,k)/HIdot(i),
c     &                       0.1*HI(i,j,k)/HIdot_prev(i),
c     &                       0.01*energy/chdot(i),
c     &                       0.1*tgas(i)/tdot(i),
c     &                       10*olderdtit(i),
c     &                       0.1*de(i,j,k)/dedot(i),
c     &                       0.1*de(i,j,k)/dedot_prev(i)
               endif
               if (i.eq.13) then
c                 write(6,*) 'H2I:', H2I(i,j,k), H2Idot(i), H2Idot_prev(i)
               endif
c              endif
c
c


               if (ge(i,j,k) .le. 0.0 .and. idual .eq. 1)
     &           write(6,*) 'solve_rate_cool: ge0a:',
     &             i,j,k,ge(i,j,k),energy,d(i,j,k),e(i,j,k),iter,tgas(i)
c               if (idual .eq. 1 .and.
c     &           ge(i,j,k)+edot(i)/d(i,j,k)*dtit(i) .le. 0.0)
c     &         write(6,*) 'solve_rate_cool: ge1a:',
c     &              i,j,k,iter,ge(i,j,k),edot(i),tgas(i),
c     &              energy,de(i,j,k),ttot(i),d(i,j,k),e(i,j,k)

c
c              If the timestep is too small, then output some debugging info
c
c               if (((dtit(i)/dt .lt. 1.0e-2 .and. iter .gt. 100) 
c     &               .or. iter .gt. itmax-100) .and.
c     &              abs((dt-ttot(i))/dt) .gt. 1.0e-3) then
                if (iter.gt.2e4) then
c                 omask(i) = .true.
c                 write(6,*) 'solve_rate_cool: omask = 1', i,j,k, iter
               endif
c
            endif
              H2Idot_prev(i) = 0.0d0
              HIdot_prev(i) = 0.0d0
              dedot_prev(i) = 0.0d0
            enddo   ! end loop over i
c   For RK-2 we cut dt in half
           do i=is+1,ie+1
             usedtit(i) = 0.5*dtit(i)
           enddo
c
c
c           Solve rate equations with backward differencing ---
            call step_rate(de, HI, HII, HeI, HeII, HeIII, d,
     &                     HM, H2I, H2II, DI, DII, HDI, usedtit,
     &                     in, jn, kn, is, ie, j, k, ispecies,
     &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
     &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
     &                     k23, k24, k25, k26, k27, k28, k29, k30, k31,
     &                     k50, k51, k52, k53, k54, k55,
     &                     k56, k24shield, k25shield, k26shield,
     &                     HIp, HIIp, HeIp, HeIIp, HeIIIp, dep,
     &                    HMp, H2Ip, H2IIp, DIp, DIIp, HDIp, dedot_prev,
     &                 HIdot_prev, H2Idot_prev, HDIdot_prev, DIdot_prev,
     &                     fh, dtoh, tgas, tgasold, itmask, eqmask )
c
c We calculate mu here to get a mudot
c
            call calculate_mu(
     &                d, e, ge, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, is, ie, j, k,
     &                HM, H2I, H2II, DI, DII, HDI, ispecies,
     &                totalMass, totalN, mu, mudot, dtit, itmask )
            call cool1d_multi_bdf(
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
c           Update temp, total and gas energy one timestep
c
            do i = is+1, ie+1
            if (itmask(i).eqv..true.) then
c
c           First chemical heating from creation of H2
c           Assume that during the three-body reaction, the H2 molecule
c           is in a fully excited state, and that it decays immediately.
c           Only use at densities where three-body is important.
c
c            write(6,*) 'hello!', i,k
              if (ispecies.gt.1.and.d(i,j,k)*dom.gt.1e7) then
                if(H2Idot_prev(i).gt.0.0d0) then
                  edotplus(i) = edotplus(i)
     &                        + abs(H2Idot_prev(i))*4.48/chunit
                elseif(H2Idot_prev(i).lt.0.0d0) then
                  edotminus(i) = edotminus(i)
     &                        +  abs(H2Idot_prev(i))*4.48/chunit
                endif
              endif
c
c
             if (imethod .eq. 2) then
c
c        Zeus - e() is really gas energy
c
                energy = d(i,j,k)*e(i,j,k)
             else
                if (idual .eq. 1) then
c
c           PPM with dual energy -- use gas energy
c
                  energy = d(i,j,k)*ge(i,j,k)
                else
c
c           PPM without dual energy -- use total energy
c
                  energy = e(i,j,k) - 0.5*u(i,j,k)**2
                  if (idim .gt. 1) energy = energy - 0.5*v(i,j,k)**2
                  if (idim .gt. 2) energy = energy - 0.5*w(i,j,k)**2
                  energy = d(i,j,k)*energy
                endif
              endif
              energy = max(energy, 1.0e-30)
              tdotplus(i) = tgas(i)*edotplus(i)/energy
              tdotminus(i) = edotminus(i)/energy
              edot(i) = edotplus(i) - edotminus(i)
c             Take into account changes in mu
              if (mudot(i) .gt. 0.0) then
                tdotplus(i) = tdotplus(i) 
     &                        + abs(tgas(i)*mudot(i)/mu(i))
              else
                tdotminus(i) = tdotminus(i)
     &                        + abs(mudot(i)/mu(i))
              endif
              if (abs(tdotplus(i)-tgas(i)*tdotminus(i)) / 
     &            (tdotplus(i)+tgas(i)*tdotminus(i)).lt.toler) then
c                write(6,*) 'zeroing tdot:',i,j,tgas(i),tdot
c                tdotminus(i) = 0.0d0
c                tdotplus(i) = 0.0d0
              endif
              tdot(i) = tdotplus(i) - tgas(i)*tdotminus(i)
              if (tdot(i).ne.tdot(i)) then
                write(6,*) 'solve_rate_cool: tdot0a',tdot(i),
     &                     tdotplus(i),tdotminus(i)
                write(6,*) 'solve_rate_cool: tdot0b',mudot(i),mu(i)
                write(6,*) 'solve_rate_cool: tdot0c',e(i,j,k),ge(i,j,k)
                write(6,*) 'solve_rate_cool: tdot0d',totalMass(i),totalN(i)
                write(6,*) 'solve_rate_cool: tdot0e',temp,dtit(i)
              endif
              tgasold(i) = tgas(i)
              tgas(i) = (tdotplus(i)*dtit(i) + tgas(i)) 
     &                   / (1.0e0 + (tdotminus(i)*dtit(i)))
               if ((abs((tgas(i)-tgasold(i))/tgas(i)).gt.0.5).and.
     &              (tgas(i).gt.100))then
c                   write(6,*) 'HUGE CHANGE: a', e(i,j,k), temp, iter
c                   write(6,*) 'HUGE CHANGE: b', tgas(i), tgasold(i), gamma2(i)
c                   write(6,*) 'HUGE CHANGE: c', totalN(i), d(i,j,k), totalMass(i)
c                   write(6,*) 'HUGE CHANGE: d', mudot(i), mu(i), mudot(i)/mu(i)
c                   write(6,*) 'HUGE CHANGE: e', tdotplus(i),tdotminus(i),dtit(i)
c                   write(6,*) 'HUGE CHANGE: f', edotplus(i),edotminus(i)
c                   write(6,*) 'HUGE CHANGE: g', edot(i)*dtit(i)
c                   write(6,*) 'HUGE CHANGE: h', H2I(i,j,k)/d(i,j,k), mudot(i)*dtit(i)
c                   write(6,*) 'HUGE CHANGE: i', i,j,k
c                   write(6,*) 'HUGE CHANGE: j', H2Idot_prev(i)*dtit(i), H2Idot_prev(i), H2I(i,j,k)
c                   write(6,*) 'HUGE CHANGE: k', (energy/edot(i))/dtit(i), (tgas(i)/tdot(i))/dtit(i)
c                   write(6,*) 'HUGE CHANGE:'
c                   omask(i)=.true.
                   write(6,*) 'TGAS NOTE: a', i,j,k
                   write(6,*) 'TGAS NOTE: b', d(i,j,k),tgas(i),tgasold(i)
                   write(6,*) 'TGAS NOTE: c', dtit(i)*edot(i)/energy, dtit(i)*tdot(i)/tgas(i)
                   write(6,*) 'TGAS NOTE: d', dtit(i)*chdot(i)/energy, H2Idot_prev(i)*dtit(i)
                   write(6,*) 'TGAS NOTE: e', iter,dtit(i),dt
                   write(6,*) 'TGAS NOTE: f', mudot(i)*dtit(i),mu(i),mudot(i)
                   write(6,*) 'TGAS NOTE: g', mudot(i)/mu(i), edot(i)/energy, 
     &                tgas(i)*dtit(i)*(mudot(i)/mu(i)+edot(i)/energy)
                   write(6,*) 'TGAS NOTE: h', chdot(i)/energy*tgas(i)*dtit(i),
     & (H2Idot_prev(i)/energy)*4.48/chunit*tgas(i)*dtit(i)
                   write(6,*) 'TGAS NOTE: i', H2Idot_prev(i), H2Idot(i), H2Idot_prev(i)/H2Idot(i)
                   write(6,*) 'TGAS NOTE: j', H2Ip(i), H2I(i,j,k), H2Ip(i)/H2I(i,j,k)
                   write(6,*) 'TGAS NOTE: k', HIp(i), HI(i,j,k), HIp(i)/HI(i,j,k)
                   write(6,*) 'TGAS NOTE: l', HIIp(i), HII(i,j,k), HIIp(i)/HII(i,j,k)
               endif
               if ((tgas(i).ne.tgas(i)).or.
     &             (tgas(i).lt.0.0d0)) then
                 write(6,*) 'solve_rate_cool: tgas1a',tgas(i),
     &                      tdotplus(i),tdotminus(i)
                 write(6,*) 'solve_rate_cool: tgas1b',dtit(i),i,j,k,iter
                 write(6,*) 'solve_rate_cool: tgas1c',d(i,j,k),mudot(i),mu(i)
                 write(6,*) 'solve_rate_cool: tgas1d',tgasold(i)
               endif
            endif
            enddo
c
c Now we recalculate the energy, based on the new temperature
c
            if (ispecies.gt.1) then
                call calculate_gamma2(
     &                       d, e, ge, de, HI, HII, HeI, HeII, HeIII,
     &                       in, jn, kn, is, ie, j, k,
     &                       HM, H2I, H2II, DI, DII, HDI, tgas,
     &                       gamma, gamma2, totalN, itmask )
            else
              do i = is+1, ie+1
                gamma2(i) = gamma
              enddo
            endif
            do i = is+1, ie+1
            if (itmask(i).eqv..true.) then
              if (gamma2(i).gt.1.8d0.or.gamma2(i).lt.1.0d0) then
                write(6,*) 'solve_rate_cool: gamma2', gamma2(i), H2I(i,j,k)
              endif
              if (imethod.eq.2) then ! Zeus
                e(i,j,k) = totalN(i) * tgas(i)
     &                   / ( d(i,j,k) * (gamma2(i)-1.0)
     &                      * utem )
              else
                e(i,j,k) = totalN(i) * tgas(i)
     &                   / ( d(i,j,k) * (gamma2(i)-1.0)
     &                       * utem )
                if (idual.eq.1) ge(i,j,k) = e(i,j,k)
                e(i,j,k) = e(i,j,k) + 0.5*u(i,j,k)**2
                if (idim .gt. 1) e(i,j,k) = e(i,j,k) + 0.5*v(i,j,k)**2
                if (idim .gt. 2) e(i,j,k) = e(i,j,k) + 0.5*w(i,j,k)**2
               endif
               if((ge(i,j,k).ne.ge(i,j,k)).or.
     &            (ge(i,j,k).lt.0.0d0))then
                 write(6,*)'solve_rate_cool: ge2a', ge(i,j,k),i,j,k,iter
                 write(6,*)'solve_rate_cool: ge2b', totalN(i),gamma2(i)
                 write(6,*)'solve_rate_cool: ge2c', energy,d(i,j,k),tgas(i)
               endif
            endif
c            write(6,*) e(i,j,k)-ge(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k)
            enddo
            call lookup_cool_rates1d(temstart, temend, nratec, j, k,
     &                is, ie, ijk, iradtype, iradshield, in, jn, kn,
     &                ispecies, tgas, tgasold, HI, HII, HeI, HeII,
     &                k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &                k11a, k12a, k13a, k13dda, k14a, k15a, k16a,
     &                k17a, k18a, k19a, k22a, k23a,
     &                k50a, k51a, k52a, k53a, k54a, k55a, k56a,
     &                avgsighp, avgsighep, avgsighe2p, piHI, piHeI,
     &                k1, k2, k3, k4, k5, k6, k7, k8, k9, k10,
     &                k11, k12, k13, k14, k15, k16, k17, k18,
     &                k19, k22, k24, k25, k26,
     &                k50, k51, k52, k53, k54, k55,
     &                k56, k13dd, k24shield, k25shield, k26shield,
     &                t1, t2, tdef, logtem, indixe, 
     &                dom, coolunit, tbase1, itmask )
            call step_rate(de, HI, HII, HeI, HeII, HeIII, d,
     &                     HM, H2I, H2II, DI, DII, HDI, usedtit,
     &                     in, jn, kn, is, ie, j, k, ispecies,
     &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
     &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
     &                     k23, k24, k25, k26, k27, k28, k29, k30, k31,
     &                     k50, k51, k52, k53, k54, k55,
     &                     k56, k24shield, k25shield, k26shield,
     &                     HIp, HIIp, HeIp, HeIIp, HeIIIp, dep,
     &                    HMp, H2Ip, H2IIp, DIp, DIIp, HDIp, dedot_prev,
     &                 HIdot_prev, H2Idot_prev, HDIdot_prev, DIdot_prev,
     &                     fh, dtoh, tgas, tgasold, itmask, eqmask )
            call calculate_mu(
     &                d, e, ge, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, is, ie, j, k,
     &                HM, H2I, H2II, DI, DII, HDI, ispecies,
     &                totalMass, totalN, mu, mudot, dtit, itmask )
c
c           Add the timestep to the elapsed time for each cell and find
c            minimum elapsed time step in this row
c
            ttmin = 1.0e+30
            do i = is+1, ie+1
            if (itmask(i).eqv..true.) then
               ttot(i) = ttot(i) + dtit(i)
               if (ttot(i).lt.ttmin) then
                 ttmin = ttot(i)
                 imin = i
               endif
               if (omask(i) .eqv. .true.) then
                 tto = max(ttot(i),dtit(i))
                 write(3,2000),i,j,k,iter,ge(i,j,k),edot(i),tgas(i),
     &             energy,de(i,j,k),tto,d(i,j,k),e(i,j,k),dtit(i),
     &             HI(i,j,k), HII(i,j,k), HeI(i,j,k),HeII(i,j,k), 
     &             HeIII(i,j,k),H2I(i,j,k),H2II(i,j,k),HM(i,j,k),
     &             DI(i,j,k), DII(i,j,k), HDI(i,j,k),
     &             dedot(i),HIdot(i),H2Idot(i),energy,dedot_prev(i),
     &             HIdot_prev(i), H2Idot_prev(i),mudot(i),mu(i),dep(i),tdot(i),
     &             tdotplus(i),tdotminus(i),H2Idotplus(i),H2Idotminus(i),
     &             k1(i),k2(i),k3(i),k4(i),k5(i),k6(i),k7(i),k8(i),
     &             k9(i),k10(i),k11(i),k12(i),k13(i),k14(i),k15(i),
     &             k16(i),k17(i),k18(i),k19(i),k22(i),k23(i),
     &             k50(i),k51(i),k52(i),k53(i),k54(i),k55(i),
     &             HIp(i),HIIp(i),HeIp(i),HeIIp(i),HeIIIp(i),
     &             H2Ip(i),H2IIp(i),HMp(i),p2d(i)

 2000          format(4(i7,1x),1p,85(e20.9))
               endif
               if (abs(dt-ttmin).lt.0.001*dt) then
c                 write(6,*) 'setting masks', i,j,k,iter
                 omask(i) = .false.
                 itmask(i) = .false.
               endif
            endif
            enddo
c     
c           If all cells are done (on this slice), then exit
c     
            if (abs(dt-ttmin).lt. 0.001*dt) go to 9999
c     
c           Next subcycle iteration
c     
         enddo
c     
 9999    continue
c     
c        Check to see if we exceed the maximum iterations
c     
         if (iter .gt. itmax ) then
            write(6,*) 'solve_rate_cool: RATE iter exceeds ',itmax,
     &                 ' at i,j,k =',imin,j,k
          write(6,*) 'solve_rate_cool 1:', d(imin,j,k), de(imin,j,k)
          write(6,*) 'solve_rate_cool 2:', HI(imin,j,k), HII(imin,j,k)
          write(6,*) 'solve_rate_cool 3:', HeI(imin,j,k), HeII(imin,j,k)
          write(6,*) 'solve_rate_cool 4:', HeIII(imin,j,k), HM(imin,j,k)
          write(6,*) 'solve_rate_cool 5:', H2I(imin,j,k), H2II(imin,j,k)
          write(6,*) 'solve_rate_cool 6:', HIdot(imin), dedot(imin)
          write(6,*) 'solve_rate_cool 7:', H2Idot(imin), dtit(imin)
          write(6,*) 'solve_rate_cool 8:',0.1*de(imin,j,k)/dedot(imin),
     & 0.1*HI(imin,j,k)/HIdot(imin)
          write(6,*) 'solve_rate_cool 9:',0.1*H2I(imin,j,k)/H2Idot(imin)
          write(6,*) 'solve_rate_cool 0:',abs(dt-ttmin),0.001*dt
          write(6,*) 'solve_rate_cool A:',urho,aye,uxyz
          write(6,*) 'solve_rate_cool B:',utim,uaye
          write(6,*) 'solve_rate_cool C:',tgas(imin)
            errcode = 1
c            return
          do i=is+1, ie+1
            write(6,*) 'solve_rate_cool icheck 0:', i, ttot(i)/dt, 
     &                  d(i,j,k), H2I(i,j,k)/d(i,j,k), tgas(i)
          enddo
         endif

c We will now check our calculated tgas versus our old tgas
      
c      do i = is+1, ie+1
c        tgasolder(i)=tgas(i)
c      enddo
c         call calculate_tgas(
c     &                d, e, ge, u, v, w, de, HI, HII, HeI, HeII, HeIII,
c     &                in, jn, kn, idual, imethod, ispecies, idim,
c     &                is, ie, j, k, aye, temstart, temend,
c     &                utem, uxyz, uaye, urho, utim, gamma,
c     &                HM, H2I, H2II, DI, DII, HDI, 
c     &                tgas, tgasold, p2d, gamma2, totalN, itmask )
c
c      do i = is+1, ie+1
c        if (abs(log10(tgas(i)) - log10(tgasolder(i))).gt.0.05) then
c            write(6,*),'s_r_c: TGAS PROBLEM a:', tgas(i),tgasolder(i)
c            write(6,*),'s_r_c: TGAS PROBLEM b:', iter, gamma2(i)
c            write(6,*),'s_r_c: TGAS PROBLEM c:', i,j,k
c            write(6,*),'s_r_c: TGAS PROBLEM d:'
c        endif
c      enddo

c     
c     Next j,k
c     
c      write(6,*) 'j = ',j,'iter = ',iter

c      do i = is+1, ie
c        if (abs(log10(tgas(i))-log10(tgasolder(i))).gt.0.2
c     &      .and.d(i,j,k)*dom.gt.1e10) then
c            write(6,*) 'src: tgd0a:', tgas(i),tgasolder(i)
c            write(6,*) 'src: tgd0b:', d(i,j,k), H2I(i,j,k)
c            write(6,*) 'src: tgd0c:', ge(i,j,k), e(i,j,k)
c        endif
c        if (abs(log10(tgas(i))-log10(tgas(i+1))).gt.0.2
c     &      .and.d(i,j,k)*dom.gt.1e10) then
c            write(6,*) 'src: tgd1a:', tgas(i),tgas(i+1)
c            write(6,*) 'src: tgd1b:', d(i,j,k), H2I(i,j,k)
c            write(6,*) 'src: tgd1c:', d(i+1,j,k), H2I(i+1,j,k)
c            write(6,*) 'src: tgd1d:', ge(i,j,k), ge(i+1,j,k)
c            write(6,*) 'src: tgd1e:', e(i,j,k), e(i+1,j,k)
c            if (j.le.je) then
c              write(6,*) 'src: tgd1f:', d(i,j+1,k), H2I(i,j+1,k)
c              write(6,*) 'src: tgd1g:', ge(i,j+1,k), ge(i,j+1,k)
c              write(6,*) 'src: tgd1h:', e(i,j+1,k), e(i,j+1,k)
c            endif
c            if (k.le.ke) then
c              write(6,*) 'src: tgd1i:', d(i,j,k+1), H2I(i,j,k+1)
c              write(6,*) 'src: tgd1j:', ge(i,j,k+1), ge(i,j,k+1)
c              write(6,*) 'src: tgd1k:', e(i,j,k+1), e(i,j,k+1)
c            endif
c        endif
c      enddo

      enddo
      enddo
c
c     Convert densities back to comoving from proper
c
      call scale_fields(d, de, HI, HII, HeI, HeII, HeIII,
     &                  HM, H2I, H2II, DI, DII, HDI, metal,
     &                  is, ie, js, je, ks, ke,
     &                  in, jn, kn, ispecies, imetal, aye**3)
c
c     Correct the species to ensure consistency (i.e. type conservation)
c
      call make_consistent(de, HI, HII, HeI, HeII, HeIII,
     &                     HM, H2I, H2II, DI, DII, HDI, d,
     &                     is, ie, js, je, ks, ke,
     &                     in, jn, kn, ispecies, fh, dtoh, iter )

      return
      end
c
c -----------------------------------------------------------
c   This routine scales the density fields from comoving to
c     proper densities (and back again).
c
      subroutine scale_fields(d, de, HI, HII, HeI, HeII, HeIII,
     &                        HM, H2I, H2II, DI, DII, HDI, metal,
     &                        is, ie, js, je, ks, ke,
     &                        in, jn, kn, ispecies, imetal, factor)
c -------------------------------------------------------------------
c
      implicit NONE
c
c     Arguments
c
      integer in, jn, kn, is, ie, js, je, ks, ke, ispecies, imetal
      real    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      real    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      real     d(in,jn,kn),metal(in,jn,kn)
      real    factor
c
c     locals
c
      integer i, j, k
c
c     Multiply all fields by factor (1/a^3 or a^3)
c
      do k = ks+1, ke+1
         do j = js+1, je+1
            do i = is+1, ie+1
               d(i,j,k)     = d(i,j,k)*factor
               de(i,j,k)    = de(i,j,k)*factor
               HI(i,j,k)    = HI(i,j,k)*factor
               HII(i,j,k)   = HII(i,j,k)*factor
               HeI(i,j,k)   = HeI(i,j,k)*factor
               HeII(i,j,k)  = HeII(i,j,k)*factor
               HeIII(i,j,k) = HeIII(i,j,k)*factor
            enddo
            if (ispecies .gt. 1) then
               do i = is+1, ie+1
                  HM(i,j,k)   = HM(i,j,k)*factor
                  H2I(i,j,k)  = H2I(i,j,k)*factor
                  H2II(i,j,k) = H2II(i,j,k)*factor
               enddo
            endif
            if (ispecies .gt. 2) then
               do i = is+1, ie+1
                  DI(i,j,k)  = DI(i,j,k)*factor
                  DII(i,j,k) = DII(i,j,k)*factor
                  HDI(i,j,k) = HDI(i,j,k)*factor
               enddo
            endif
            if (imetal .eq. 1) then
               do i = is+1, ie+1
                  metal(i,j,k) = metal(i,j,k)*factor
               enddo
            endif
         enddo
      enddo
c
      return
      end

c
c -----------------------------------------------------------
c This routine uses the temperature to look up the chemical
c   rates which are tabulated in a log table as a function
c   of temperature.
c
      subroutine lookup_cool_rates1d(temstart, temend, nratec, j, k,
     &                is, ie, ijk, iradtype, iradshield, in, jn, kn,
     &                ispecies, tgas1d, tgasold, HI, HII, HeI, HeII,
     &                k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &                k11a, k12a, k13a, k13dda, k14a, k15a, k16a,
     &                k17a, k18a, k19a, k22a, k23a,
     &                k50a, k51a, k52a, k53a, k54a, k55a, k56a,
     &                avgsighp, avgsighep, avgsighe2p, piHI, piHeI,
     &                k1, k2, k3, k4, k5, k6, k7, k8, k9, k10,
     &                k11, k12, k13, k14, k15, k16, k17, k18,
     &                k19, k22, k24, k25, k26,
     &                k50, k51, k52, k53, k54, k55,
     &                k56, k13dd, k24shield, k25shield, k26shield,
     &                t1, t2, tdef, logtem, indixe, 
     &                dom, coolunit, tbase1, itmask )
c -------------------------------------------------------------------
c
      implicit NONE
c
c     Arguments
c
      integer is, ie, ijk, iradtype, iradshield, nratec, 
     &        in, jn, kn, ispecies, j, k
      real temstart, temend, tgas1d(in), tgasold(in), dom
      double precision coolunit, tbase1
      logical itmask(in)
c
c     Chemistry rates as a function of temperature
c
      real k1a (nratec), k2a (nratec), k3a (nratec), k4a (nratec), 
     &     k5a (nratec), k6a (nratec), k7a (nratec), k8a (nratec), 
     &     k9a (nratec), k10a(nratec), k11a(nratec), k12a(nratec), 
     &     k13a(nratec), k14a(nratec), k15a(nratec), k16a(nratec), 
     &     k17a(nratec), k18a(nratec), k19a(nratec), k22a(nratec),
     &     k23a(nratec), k50a(nratec), k51a(nratec), k52a(nratec),
     &     k53a(nratec), k54a(nratec), k55a(nratec), k56a(nratec),
     &     k13dda(nratec, 7),
     &     k24, k25, k26,
     &     avgsighp, avgsighep, avgsighe2p, piHI, piHeI
c
c     Density fields
c
      real    HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn)
c
c     Returned rate values
c
      real k1 (ijk), k2 (ijk), k3 (ijk), k4 (ijk), k5 (ijk),
     &     k6 (ijk), k7 (ijk), k8 (ijk), k9 (ijk), k10(ijk),
     &     k11(ijk), k12(ijk), k13(ijk), k14(ijk), k15(ijk),
     &     k16(ijk), k17(ijk), k18(ijk), k19(ijk), k22(ijk),
     &     k23(ijk), k50(ijk), k51(ijk), k52(ijk), k53(ijk),
     &     k54(ijk), k55(ijk), k56(ijk), k13dd(ijk, 7),
     &     k24shield(ijk), k25shield(ijk), k26shield(ijk)
c
c     1D temporaries (passed in)
c 
      integer indixe(ijk)
      real t1(ijk), t2(ijk), logtem(ijk), tdef(ijk)
c
c     Parameters
c
      double precision everg, e24, e26
      parameter(everg = 1.60184d-12, e24 = 13.6d0, e26 = 24.6d0)
c
c     locals
c
      integer i, n1
      real factor, x, logtem0, logtem9, dlogtem, nh
c
c     Set log values of start and end of lookup tables
c
      logtem0 = log(temstart)
      logtem9 = log(temend)
      dlogtem = (log(temend) - log(temstart))/real(nratec-1)
c
      do i = is+1, ie+1
      if (itmask(i).eqv..true.) then
c
c        Compute temp-centered temperature (and log)
c
         logtem(i) = log(tgas1d(i))
c         logtem(i) = log(0.5*(tgas1d(i)+tgasold(i)))
         logtem(i) = max(logtem(i), logtem0)
         logtem(i) = min(logtem(i), logtem9)
c
c        Find index into tble and precompute interpolation values
c
         indixe(i) = min(nratec-1,
     &                max(1,int((logtem(i)-logtem0)/dlogtem)+1))
         t1(i) = (logtem0 + (indixe(i) - 1)*dlogtem)
         t2(i) = (logtem0 + (indixe(i)    )*dlogtem)
         tdef(i) = t2(i) - t1(i)
c
c        Do linear table lookup (in log temperature)
c
         k1(i) = k1a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k1a(indixe(i)+1) -k1a(indixe(i)))/tdef(i)
         k2(i) = k2a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k2a(indixe(i)+1) -k2a(indixe(i)))/tdef(i)
         k3(i) = k3a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k3a(indixe(i)+1) -k3a(indixe(i)))/tdef(i)
         k4(i) = k4a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k4a(indixe(i)+1) -k4a(indixe(i)))/tdef(i)
         k5(i) = k5a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k5a(indixe(i)+1) -k5a(indixe(i)))/tdef(i)
         k6(i) = k6a(indixe(i)) + (logtem(i) - t1(i))
     &           *(k6a(indixe(i)+1) -k6a(indixe(i)))/tdef(i)
c
      endif
      enddo
c
c     Look-up for 9-species model
c
      if (ispecies .gt. 1) then
         do i = is+1, ie+1
      if (itmask(i).eqv..true.) then
            k7(i) = k7a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k7a(indixe(i)+1) -k7a(indixe(i)))/tdef(i)
            k8(i) = k8a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k8a(indixe(i)+1) -k8a(indixe(i)))/tdef(i)
            k9(i) = k9a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k9a(indixe(i)+1) -k9a(indixe(i)))/tdef(i)
            k10(i) = k10a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k10a(indixe(i)+1) -k10a(indixe(i)))/tdef(i)
            k11(i) = k11a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k11a(indixe(i)+1) -k11a(indixe(i)))/tdef(i)
            k12(i) = k12a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k12a(indixe(i)+1) -k12a(indixe(i)))/tdef(i)
            k13(i) = k13a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k13a(indixe(i)+1) -k13a(indixe(i)))/tdef(i)
            k14(i) = k14a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k14a(indixe(i)+1) -k14a(indixe(i)))/tdef(i)
            k15(i) = k15a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k15a(indixe(i)+1) -k15a(indixe(i)))/tdef(i)
            k16(i) = k16a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k16a(indixe(i)+1) -k16a(indixe(i)))/tdef(i)
            k17(i) = k17a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k17a(indixe(i)+1) -k17a(indixe(i)))/tdef(i)
            k18(i) = k18a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k18a(indixe(i)+1) -k18a(indixe(i)))/tdef(i)
            k19(i) = k19a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k19a(indixe(i)+1) -k19a(indixe(i)))/tdef(i)
            k22(i) = k22a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k22a(indixe(i)+1) -k22a(indixe(i)))/tdef(i)
            k23(i) = k23a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k23a(indixe(i)+1) -k23a(indixe(i)))/tdef(i)
         endif
         enddo
c     
         do n1 = 1, 7
            do i = is+1, ie+1
            if (itmask(i).eqv..true.) then
               k13dd(i,n1) = k13dda(indixe(i),n1) + (logtem(i) - t1(i))
     &             *(k13dda(indixe(i)+1,n1) - 
     &               k13dda(indixe(i)  ,n1) )/tdef(i)
            endif
            enddo
         enddo
      endif
c
c     Look-up for 12-species model
c
      if (ispecies .gt. 2) then
         do i = is+1, ie+1
         if (itmask(i).eqv..true.) then
            k50(i) = k50a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k50a(indixe(i)+1) -k50a(indixe(i)))/tdef(i)
            k51(i) = k51a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k51a(indixe(i)+1) -k51a(indixe(i)))/tdef(i)
            k52(i) = k52a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k52a(indixe(i)+1) -k52a(indixe(i)))/tdef(i)
            k53(i) = k53a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k53a(indixe(i)+1) -k53a(indixe(i)))/tdef(i)
            k54(i) = k54a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k54a(indixe(i)+1) -k54a(indixe(i)))/tdef(i)
            k55(i) = k55a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k55a(indixe(i)+1) -k55a(indixe(i)))/tdef(i)
            k56(i) = k56a(indixe(i)) + (logtem(i) - t1(i))
     &            *(k56a(indixe(i)+1) -k56a(indixe(i)))/tdef(i)
         endif
         enddo
      endif
c
c        Include approximate self-shielding factors if requested
c

      do i = is+1, ie+1
      if (itmask(i).eqv..true.) then
         k24shield(i) = k24
         k25shield(i) = k25
         k26shield(i) = k26
      endif
      enddo
      if (iradshield .eq. 1) then
         do i = is+1, ie+1
         if (itmask(i).eqv..true.) then
            k24shield(i) = k24shield(i)*exp(-HI(i,j,k)*avgsighp*dom)
            k25shield(i) = k25shield(i)*exp(-HeII(i,j,k)*avgsighe2p*dom)
            k26shield(i) = k26shield(i)*exp(-HeI(i,j,k)*avgsighep*dom)
         endif
         enddo
      endif
c
c        If using a high-energy radiation field, then account for
c          effects of secondary elections (Shull * Steenberg 1985)
c          (see calc_rate.src)
c
      if (iradtype .eq. 8) then
         do i = is+1, ie+1
         if (itmask(i).eqv..true.) then
            x = max(HII(i,j,k)/(HI(i,j,k)+HII(i,j,k)), 1.0e-4)
            factor = 0.3908*(1.0 - x**0.4092)**1.7592
            k24shield(i) = k24shield(i) + 
     &         factor*(piHI + 0.08*piHeI)/(e24*everg) * coolunit*tbase1
            factor = 0.0554*(1.0 - x**0.4614)**1.6660
            k26shield(i) = k26shield(i) + 
     &         factor*(piHI/0.08 + piHeI)/(e26*everg) * coolunit*tbase1
         endif
         enddo
      endif
c

c
c           If using H2, and using the density-dependent collisional
c             H2 dissociation rate, then replace the the density-independant
c                k13 rate with the new one.
c         May/00: there appears to be a problem with the density-dependent
c             collisional rates.  Currently turned off until further notice.
c


            if (ispecies .gt. 1) then
               do i = is+1, ie+1
               if (itmask(i).eqv..true.) then
                  nh = min(HI(i,j,k)*dom, 1.0e9)
                  k13(i) = 1.0e-30
                  if (tgas1d(i) .ge. 500.0 .and.
     &                tgas1d(i) .lt. 1.0e6) then
                     k13(i) = k13dd(i,1)-k13dd(i,2)/
     &                          (1.0+(nh/k13dd(i,5))**k13dd(i,7))
     &                     + k13dd(i,3)-k13dd(i,4)/
     &                          (1.0+(nh/k13dd(i,6))**k13dd(i,7))
                     k13(i) = max(10.0**k13(i), 1.0e-30)
                  endif
               endif
               enddo
            endif

c
      return
      end

c -------------------------------------------------------------------
c  This routine calculates the electron and HI rates of change in
c    order to determine the maximum permitted timestep
c
      subroutine rate_timestep(dedot, HIdot, H2Idot, HDIdot, DIdot,
     &                     ispecies,
     &                     de, HI, HII, HeI, HeII, HeIII, d,
     &                     HM, H2I, H2II, DI, DII, HDI,
     &                     in, jn, kn, is, ie, j, k, 
     &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
     &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
     &                     k23, k24, k25, k26, k27, k28, k29, k30, k31,
     &                     k50, k51, k52, k53, k54, k55,
     &                     k56, k24shield, k25shield, k26shield,
     &                     dedotplus,  dedotminus,
     &                     HIdotplus,  HIdotminus,
     &                     H2Idotplus, H2Idotminus,
     &                     HDIdotplus, HDIdotminus, 
     &                     DIdotplus, DIdotminus, 
     &                     dtit, iter,
     &                     dedot_prev, HIdot_prev, H2Idot_prev,
     &                     HDIdot_prev, DIdot_prev,
     &                     itmask )
c -------------------------------------------------------------------
c
      implicit NONE
c
c     arguments
c
      integer ispecies, is, ie, j, k, in, jn, kn
      real*8 dedot(in), HIdot(in), H2Idot(in), HDIdot(in),
     &       DIdot(in),
     &       dedotplus(in), dedotminus(in),
     &       HIdotplus(in), HIdotminus(in),
     &       H2Idotplus(in), H2Idotminus(in),
     &       HDIdotplus(in), HDIdotminus(in),
     &       DIdotplus(in), DIdotminus(in),
     &       dedot_prev(in), HIdot_prev(in), H2Idot_prev(in),
     &       HDIdot_prev(in), DIdot_prev(in)
      logical itmask(in)
c
c     Density fields
c
      real    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &         d(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      real    DI(in,jn,kn),  DII(in,jn,kn),  HDI(in,jn,kn)
      real    dtit(in)
      integer    iter
c
c     Rate values
c
      real k1 (in), k2 (in), k3 (in), k4 (in), k5 (in),
     &     k6 (in), k7 (in), k8 (in), k9 (in), k10(in),
     &     k11(in), k12(in), k13(in), k14(in), k15(in),
     &     k16(in), k17(in), k18(in), k19(in), k22(in),
     &     k23(in), k50(in), k51(in), k52(in), k53(in), 
     &     k54(in), k55(in), k56(in), 
     &     k24shield(in), k25shield(in), k26shield(in),
     &     k24, k25, k26, k27, k28, k29, k30, k31
c
c     locals
c
      integer i
      real toler
      parameter (toler = 1e-20 )
c
      if (ispecies .eq. 1) then
c     
         do i = is+1, ie+1
         if (itmask(i).eqv..true.) then
c
c     Compute the electron density rate-of-change
c
            dedot(i) = 
     &               + DBLE(k1(i))*DBLE(HI(i,j,k)*de(i,j,k))
     &               + DBLE(k3(i))*DBLE(HeI(i,j,k)*de(i,j,k)/4.0)
     &               + DBLE(k5(i))*DBLE(HeII(i,j,k)*de(i,j,k)/4.0)
     &               - DBLE(k2(i))*DBLE(HII(i,j,k)*de(i,j,k))
     &               - DBLE(k4(i))*DBLE(HeII(i,j,k)*de(i,j,k)/4.0)
     &               - DBLE(k6(i))*DBLE(HeIII(i,j,k)*de(i,j,k)/4.0)

     &               +      ( k24shield(i)*HI(i,j,k)
     &               + k25shield(i)*HeII(i,j,k)/4.0
     &               + k26shield(i)*HeI(i,j,k)/4.0)

c
c     Compute the HI density rate-of-change
c     
            HIdot(i) =
     &               - DBLE(k1(i))*DBLE(HI(i,j,k)*de(i,j,k))
     &               + DBLE(k2(i))*DBLE(HII(i,j,k)*de(i,j,k))

     &               -      k24shield(i)*HI(i,j,k)

c     

            H2Idot(i) = 0.0d0

         endif
         enddo
      else
c
c         Include molecular hydrogen rates for HIdot
c     
c  Reordered to reduce roundoff-error
         do i = is+1, ie+1
         if (itmask(i).eqv..true.) then
            HIdotminus(i) = abs(
     &               -(2.*DBLE(k22(i))*DBLE(HI(i,j,k))  *HI(i,j,k)**2)
     &               -    DBLE(k1(i)) *DBLE(de(i,j,k)   *HI(i,j,k)  )
     &               -    DBLE(k7(i)) *DBLE(de(i,j,k)   *HI(i,j,k))
     &               -    DBLE(k8(i)) *DBLE(HM(i,j,k)   *HI(i,j,k))
     &               -    DBLE(k9(i)) *DBLE(HII(i,j,k)  *HI(i,j,k))
     &               -    DBLE(k10(i))*DBLE(H2II(i,j,k) *HI(i,j,k)/2. )

     &               -      k24shield(i)*HI(i,j,k)

     &               )
            HIdotplus(i) =
     &               + 2.*DBLE(k13(i))*DBLE(HI(i,j,k)   *H2I(i,j,k)/2. )
     &               +    DBLE(k15(i))*DBLE(HM(i,j,k)   *HI(i,j,k))
     &               +    DBLE(k2(i)) *DBLE(HII(i,j,k)  *de(i,j,k) )
     &               +    DBLE(k11(i))*DBLE(HII(i,j,k)  *H2I(i,j,k)/2.)
     &               + 2.*DBLE(k12(i))*DBLE(de(i,j,k)   *H2I(i,j,k)/2.)
     &               +    DBLE(k14(i))*DBLE(HM(i,j,k)   *de(i,j,k))
     &               + 2.*DBLE(k16(i))*DBLE(HM(i,j,k)   *HII(i,j,k))
     &               + 2.*DBLE(k18(i))*DBLE(H2II(i,j,k) *de(i,j,k)/2.)
     &               +    DBLE(k19(i))*DBLE(H2II(i,j,k) *HM(i,j,k)/2.)
     &               + 2.*DBLE(k23(i))*DBLE(H2I(i,j,k)  *H2I(i,j,k)/4.)
            if (abs(HIdotplus(i)-HIdotminus(i)) / 
     &          abs(HIdotplus(i)+HIdotminus(i)).lt.toler) then
               HIdot(i) = HIdot_prev(i)
            else
               HIdot(i) = HIdotplus(i) - HIdotminus(i)
            endif

            if ((HIdot(i).ne.HIdot(i)).or.(HIdot(i).gt.1d50))then
                write(6,*),'solve_rate_cool: HIdot0a:',i,j,k
                write(6,*),'solve_rate_cool: HIdot0b:',
     &              -(2.*(k22(i)*HI(i,j,k)) ) *HI(i,j,k)**2, HIdot(i)
                write(6,*),'solve_rate_cool: HIdot0c:',HI(i,j,k),k22(i),
     &                    d(i,j,k)
                write(6,*),'solve_rate_cool: HIdot0d:',HIdotplus(i),HIdotminus(i)
                write(6,*),'solve_rate_cool: HIdot0e:',dtit(i),iter
            endif
c
c     Compute the electron density rate-of-change
c
c   Re-ordered to reduce roundoff-error
            dedotplus(i) = 
     &               +  k8(i) * DBLE(HM(i,j,k)    * HI(i,j,k)      )
     &               +  k15(i)* DBLE(HM(i,j,k)    * HI(i,j,k)      )
     &               +  k1(i) * DBLE(HI(i,j,k)    * de(i,j,k)      )
     &               +  k3(i) * DBLE(HeI(i,j,k)   * de(i,j,k)/4.0d0)
     &               +  k5(i) * DBLE(HeII(i,j,k)  * de(i,j,k)/4.0d0)
     &               +  k17(i)* DBLE(HM(i,j,k)    * HII(i,j,k)     )
     &               +  k14(i)* DBLE(HM(i,j,k)    * de(i,j,k)      )

     &               + (k24shield(i)*HI(i,j,k)
     &               +  k25shield(i)*HeII(i,j,k)/4.0d0
     &               +  k26shield(i)*HeI(i,j,k)/4.0d0)

            dedotminus(i) = abs(
     &               -  k7(i) * DBLE(HI(i,j,k)    * de(i,j,k)      )
     &               -  k2(i) * DBLE(HII(i,j,k)   * de(i,j,k)      )
     &               -  k4(i) * DBLE(HeII(i,j,k)  * de(i,j,k)/4.0d0)
     &               -  k6(i) * DBLE(HeIII(i,j,k) * de(i,j,k)/4.0d0)
     &               -  k18(i)* DBLE(H2II(i,j,k)  * de(i,j,k)/2.0d0)
     &              )
c
            if (abs(dedotplus(i)-dedotminus(i)) / 
     &          abs(dedotplus(i)+dedotminus(i)).lt.toler) then
               dedot(i) = dedot_prev(i)
            else
               dedot(i) = dedotplus(i) - dedotminus(i)
            endif

            if (dedot(i).ne.dedot(i)) then
                write(6,*)'solve_rate_cool: dedot1a',i,j,k
                write(6,*)'solve_rate_cool: dedot1b',de(i,j,k),HI(i,j,k)
                write(6,*)'solve_rate_cool: dedot1c',k1(i),k3(i),k5(i)
                write(6,*)'solve_rate_cool: dedot1d',k8(i),k15(i),k17(i)
                write(6,*)'solve_rate_cool: dedot1e',k14(i),k2(i),k4(i)
                write(6,*)'solve_rate_cool: dedot1f',k6(i),k7(i),k18(i)
                write(6,*)'solve_rate_cool: dedot1g',HeI(i,j,k),HeII(i,j,k)
                write(6,*)'solve_rate_cool: dedot1h',dedot(i)
            endif
c
c
c     Compute the H2 rate of change
c
            H2Idotplus(i) = 
     &            +    (k8(i) * DBLE(HM(i,j,k)  * HI(i,j,k))
     &            +    k10(i) * DBLE(H2II(i,j,k)* HI(i,j,k)/2.)
     &            +    k19(i) * DBLE(H2II(i,j,k)* HM(i,j,k)/2.)
     &          + DBLE(k22(i) *      HI(i,j,k) )* DBLE(HI(i,j,k)**2))
            H2Idotminus(i) =
     &            +    k13(i) * DBLE(HI(i,j,k)  * H2I(i,j,k)/2.)
     &            +    k11(i) * DBLE(HII(i,j,k) * H2I(i,j,k)/2.)
     &            +    k12(i) * DBLE(de(i,j,k)  * H2I(i,j,k)/2.)
     &            +    k23(i) * DBLE(H2I(i,j,k) * H2I(i,j,k)/4.)

     &            +     k29 * H2I(i,j,k)/2.
     &            +     k31 * H2I(i,j,k)/2.

            if (abs(H2Idotplus(i)-H2Idotminus(i)) / 
     &          abs(H2Idotplus(i)+H2Idotminus(i)).lt.toler) then
                H2Idot(i) = H2Idot_prev(i)
            else
               H2Idot(i) = H2Idotplus(i) - H2Idotminus(i)
            endif
            if (abs(H2Idot(i)).gt.1e300) then
               write(6,*) 'solve_rate_cool: H2Idot0a', i,j,k
               write(6,*) 'solve_rate_cool: H2Idot0b', H2Idotplus(i),H2Idotminus(i)
               write(6,*) 'solve_rate_cool: H2Idot0c', (H2Idotplus(i)-H2Idotminus(i))/(H2Idotplus(i)+H2Idotminus(i))
               write(6,*) 'solve_rate_cool: H2Idot0d', H2I(i,j,k), dtit(i), iter
            endif
c

         endif
         enddo
         if (ispecies.gt.2) then
           do i = is+1, ie+1
           if (itmask(i).eqv..true.) then
             HDIdotplus(i) = 3.0d0*( k52(i) * DII(i,j,k)* H2I(i,j,k)/4.0d0
     &                      +        k54(i) * DI(i,j,k) * H2I(i,j,k)/4.0d0
     &                      )
             HDIdotminus(i) = 3.0d0*(k53(i) * HII(i,j,k) * HDI(i,j,k)/3.0d0
     &                      +        k55(i) * HI(i,j,k) * HDI(i,j,k)/3.0d0
     &                      )
             if (abs(HDIdotplus(i)-HDIdotminus(i)) / 
     &           abs(HDIdotplus(i)+HDIdotminus(i)).lt.toler) then
                 HDIdot(i) = HDIdot_prev(i)
             else
                HDIdot(i) = HDIdotminus(i) - HDIdotplus(i)
             endif
             DIdotplus(i) = 2.0d0*(
     &                        k2(i) * DII(i,j,k) * de(i,j,k)/2.0d0
     &                    +  k51(i) * DII(i,j,k) * HI(i,j,k)/2.0d0
     &                    +  k55(i) * HDI(i,j,k) * HI(i,j,k)/3.0d0
     &                     )
             DIdotminus(i) = 2.0d0*(
     &                         k1(i) * de(i,j,k)  * DI(i,j,k)/2.0d0
     &                     +  k50(i) * HII(i,j,k) * DI(i,j,k)/2.0d0
     &                     +  k54(i) * H2I(i,j,k) * DI(i,j,k)/4.0d0 
     &                      )
             if (abs(DIdotplus(i)-DIdotminus(i)) / 
     &           abs(DIdotplus(i)+DIdotminus(i)).lt.toler) then
                 DIdot(i) = DIdot_prev(i)
             else
                DIdot(i) = DIdotplus(i) - DIdotminus(i) 
             endif
           endif
           enddo
         endif
      endif
c
      return
      end


c -----------------------------------------------------------
c  This routine uses a backward-finite difference time integrator
c   to advance the rate equations by one (sub-)cycle (dtit).
c
      subroutine step_rate(de, HI, HII, HeI, HeII, HeIII, d,
     &                     HM, H2I, H2II, DI, DII, HDI, dtit,
     &                     in, jn, kn, is, ie, j, k, ispecies,
     &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
     &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
     &                     k23, k24, k25, k26, k27, k28, k29, k30, k31,
     &                     k50, k51, k52, k53, k54, k55,
     &                     k56, k24shield, k25shield, k26shield,
     &                     HIp, HIIp, HeIp, HeIIp, HeIIIp, dep,
     &                    HMp, H2Ip, H2IIp, DIp, DIIp, HDIp, dedot_prev,
     &                 HIdot_prev, H2Idot_prev, HDIdot_prev, DIdot_prev,
     &                     fh, dtoh, tgas, tgasold, itmask, eqmask )
c -------------------------------------------------------------------
c
      implicit NONE
c
c     arguments
c
      integer ispecies, in, jn, kn, is, ie, j, k
      real    dtit(in)
      real*8  dedot_prev(in),HIdot_prev(in),H2Idot_prev(in),
     &        HDIdot_prev(in),DIdot_prev(in)
      logical itmask(in), eqmask(in)
c
c     Density fields
c
      real    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &         d(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      real    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
c
c     Rate values
c
      real k1 (in), k2 (in), k3 (in), k4 (in), k5 (in),
     &     k6 (in), k7 (in), k8 (in), k9 (in), k10(in),
     &     k11(in), k12(in), k13(in), k14(in), k15(in),
     &     k16(in), k17(in), k18(in), k19(in), k22(in),
     &     k23(in), k50(in), k51(in), k52(in), k53(in), 
     &     k54(in), k55(in), k56(in), 
     &     k24shield(in), k25shield(in), k26shield(in),
     &     k24, k25, k26, k27, k28, k29, k30, k31
c
c     temporaries (passed in)
c
      real*8 HIp(in), HIIp(in), HeIp(in), HeIIp(in), HeIIIp(in),
     &     HMp(in), H2Ip(in), H2IIp(in), dep(in),
     &     DIp(in), DIIp(in), HDIp(in)
      real fh, dtoh, tgas(in), tgasold(in)
c
c     locals
c
      integer i, dostop(in)
      real*8 temp
      real*8 scoef, acoef 
      real*8 totalH, totalHe, totalD, totalDen,
     &       correctH, correctHe, correctD
      real tempde
      real toler, Dtoler
      parameter (toler = 1e-20 )
      parameter (Dtoler = 1.0e-30)

      do i = is+1, ie+1
        dostop(i) = 0
      enddo
c     scoef is creation
c     acoef is destruction
c
c   A) the 6-species integrator
c      
      if (ispecies .eq. 1) then
c
         do i = is+1, ie+1
         if (itmask(i).eqv..true.) then
c
c        1) HI
c
            scoef  = k2(i)*HII(i,j,k)*de(i,j,k)
            acoef  = k1(i)*de(i,j,k)

     &             + k24shield(i)

            HIp(i)  = (scoef*DBLE(dtit(i)) + HI(i,j,k))
     &                / (1.0d0 + acoef*DBLE(dtit(i)))
c
c        2) HII
c 
            scoef  = (k1(i)*HIp(i))*de(i,j,k)

     &                   + k24shield(i)*HIp(i)
c     &                   + k24shield(i)*HI(i,j,k)

            acoef  = k2(i)*de (i,j,k)
            HIIp(i) = (scoef*dtit(i) + HII(i,j,k))
     &                / (1.0d0 + acoef*dtit(i))
c
c                 3) Electron density
c
            scoef = 0.0d0

     &                 + k24shield(i)*HI(i,j,k)
     &                 + k25shield(i)*HeII(i,j,k)/4.
     &                 + k26shield(i)*HeI(i,j,k)/4.

            acoef = -(k1(i)*HI(i,j,k)      - k2(i)*HII(i,j,k)
     &              + k3(i)*HeI(i,j,k)/4.  - k6(i)*HeIII(i,j,k)/4.0d0
     &              + k5(i)*HeII(i,j,k)/4. - k4(i)*HeII(i,j,k)/4.0d0)
            dep(i)   = (scoef*dtit(i) + de(i,j,k))
     &                     / (1.0d0 + acoef*dtit(i))
c
c  --- (B) Do helium chemistry in any case: (for ispecies = 1) ---
c
c        4) HeI
c 
            scoef  = k4(i)*HeII(i,j,k)*de(i,j,k)
            acoef  = k3(i)*de(i,j,k)

     &                + k26shield(i)

            HeIp(i)   = ( scoef*dtit(i) + HeI(i,j,k) ) 
     &                 / ( 1.0d0 + acoef*dtit(i) )
c
c        5) HeII
c
            scoef  = (k3(i)*HeIp(i))*de(i,j,k)
     &             + (k6(i)*HeIII(i,j,k))*de(i,j,k)

     &             + k26shield(i)*HeIp(i)

            acoef  = k4(i)*de(i,j,k) + k5(i)*de(i,j,k)

     &             + k25shield(i)

            scoef=max(scoef,1.0e-30)
            HeIIp(i)  = ( scoef*dtit(i) + HeII(i,j,k) )
     &                 / ( 1.0d0 + acoef*dtit(i) )
c
c        6) HeIII
c
            scoef   = k5(i)*HeIIp(i)*de(i,j,k)

     &              + k25shield(i)*HeIIp(i)

            acoef   = k6(i)*de(i,j,k)
            scoef=max(scoef,1.0e-30)
            HeIIIp(i)  = ( scoef*dtit(i) + HeIII(i,j,k) )
     &                   / ( 1.0d0 + acoef*dtit(i) )
         endif
         enddo
      endif
c
c --- (C) Now do extra 3-species for molecular hydrogen ---
c
      if (ispecies .gt. 1) then
c
c        First, do HI/HII with molecular hydrogen terms
c
         do i = is+1, ie+1
         if (itmask(i).eqv..true.) then
c
c        1) HI
c     
            scoef  =       DBLE(k2(i)) * DBLE(HII(i,j,k) * de(i,j,k))
     &             + 2.0d0*DBLE(k13(i))* DBLE(HI(i,j,k)  * H2I(i,j,k))/2.0d0
     &             +       DBLE(k11(i))* DBLE(HII(i,j,k) * H2I(i,j,k))/2.0d0
     &             + 2.0d0*DBLE(k12(i))* DBLE(de(i,j,k)  * H2I(i,j,k))/2.0d0
     &             +       DBLE(k14(i))* DBLE(HM(i,j,k)  * de(i,j,k))
     &             +       DBLE(k15(i))* DBLE(HM(i,j,k)  * HI(i,j,k))
     &             + 2.0d0*DBLE(k16(i))* DBLE(HM(i,j,k)  * HII(i,j,k))
     &             + 2.0d0*DBLE(k18(i))* DBLE(H2II(i,j,k)* de(i,j,k))/2.0d0
     &             +       DBLE(k19(i))* DBLE(H2II(i,j,k)* HM(i,j,k))/2.0d0
     &             + 2.0d0*DBLE(k23(i))* DBLE(H2I(i,j,k) * H2I(i,j,k))/4.0d0

     &             + 2.0d0*(k31   * H2I(i,j,k)/2.0d0)

            acoef  =         k1(i) * de(i,j,k)
     &             +         k7(i) * de(i,j,k)  
     &             +         k8(i) * HM(i,j,k)
     &             +         k9(i) *HII(i,j,k)
     &             +        k10(i)*H2II(i,j,k)/2.0d0
c     &             +        k13(i) * H2I(i,j,k)/2.0d0
     &             +2.0d0*(DBLE(k22(i)* HI(i,j,k))*(HI(i,j,k)))
c Changed above to prevent overflow

     &             + k24shield(i)


            if (ispecies.gt.4) then
              scoef = scoef + (
     &                          k50(i) * HII(i,j,k) * DI(i,j,k) / 2.0d0
     &                        + k54(i) * H2I(i,j,k) * DI(i,j,k) / 4.0d0
     &                        )
              acoef = acoef + (
     &                          k51(i) * DII(i,j,k) / 2.0d0
     &                        + k55(i) * HDI(i,j,k) / 3.0d0
     &                        )
            endif
c            HIp(i)  = ( scoef*dtit(i) + HI(i,j,k) ) / 
c     &                      ( 1.0d0 + acoef*dtit(i) )
            if (abs(scoef-HI(i,j,k)*acoef)/(scoef+HI(i,j,k)*acoef)
     &          .lt.toler .or. eqmask(i).eqv..true.)then
              HIp(i) = scoef/acoef
            else
c              HIp(i) = ((scoef*dtit(i))/(1.0d0 + acoef*dtit(i)))
c     &               + (      HI(i,j,k)/(1.0d0 + acoef*dtit(i)))
            HIp(i)  = ( scoef*dtit(i) + HI(i,j,k) ) / 
     &                      ( 1.0d0 + acoef*dtit(i) )
            endif
            if ((HIp(i).ne.HIp(i)).or.(HIp(i).gt.1.0e+30))then
              write(6,*) 'solve_rate_cool: HIp0a:', i,j,k
              write(6,*) 'solve_rate_cool: HIp0b:', HIp(i),HI(i,j,k)
              write(6,*) 'solve_rate_cool: HIp0c:', scoef,acoef, d(i,j,k)
              write(6,*) 'solve_rate_cool: HIp0d:', H2I(i,j,k), HII(i,j,k)
              write(6,*) 'solve_rate_cool: HIp0e:', HM(i,j,k), de(i,j,k)
              write(6,*) 'solve_rate_cool: HIp0f:', H2II(i,j,k),tgas(i)
              write(6,*) 'solve_rate_cool: HIp0g:', k2(i),k13(i),k11(i)
              write(6,*) 'solve_rate_cool: HIp0h:', k12(i),k14(i),k15(i)
              write(6,*) 'solve_rate_cool: HIp0i:', k16(i),k18(i),k19(i)
              write(6,*) 'solve_rate_cool: HIp0j:', DBLE(k12(i))*DBLE(de(i,j,k)*H2I(i,j,k))
            endif
c
c          2) HII
c 
            scoef  =    (k1(i)  * HIp(i)) * de(i,j,k)
     &             +    (k10(i) * H2II(i,j,k))*HIp(i)/2.0d0

     &             + k24shield(i)*HIp(i)

            acoef  =     k2(i)  * de(i,j,k)
     &             +     k9(i)  * HIp(i)
     &             +     k11(i) * H2I(i,j,k)/2.0d0
     &             +     k16(i) * HM(i,j,k)
     &             +     k17(i) * HM(i,j,k)

            if (ispecies.gt.4) then
              scoef = scoef + (
     &                          k51(i) * HIp(i)     * DII(i,j,k) / 2.0d0
     &                        + k52(i) * H2I(i,j,k) * DII(i,j,k) / 4.0d0
     &                        )
              acoef = acoef + ( 
     &                          k50(i) * DI(i,j,k) / 2.0d0
     &                        + k53(i) * HDI(i,j,k) / 3.0d0
     &                        )
            endif

            if (abs(scoef-HII(i,j,k)*acoef)/(scoef+HII(i,j,k)*acoef)
     &          .lt.toler .or. eqmask(i).eqv..true.)then
              HIIp(i) = scoef/acoef
            else
            HIIp(i)   = ( scoef*dtit(i) + HII(i,j,k) )
     &                      / ( 1.0d0 + acoef*dtit(i) )
            endif
c
c          3) HeI
c 
            scoef  = k4(i)*DBLE(HeII(i,j,k)*de(i,j,k))
            acoef  = k3(i)*de(i,j,k)

     &                   + k26shield(i)

            if (abs(scoef-HeI(i,j,k)*acoef)/(scoef+HeI(i,j,k)*acoef)
     &          .lt.toler .or. eqmask(i).eqv..true.)then
              HeIp(i) = scoef/acoef
            else
            HeIp(i)   = ( scoef*dtit(i) + HeI(i,j,k) ) 
     &                 / ( 1.0d0 + acoef*dtit(i) )
            endif
c
c          4) HeII
c
            scoef  = k3(i)*DBLE(HeIp(i)*de(i,j,k))
     &             + k6(i)*DBLE(HeIII(i,j,k)*de(i,j,k))

     &             + k26shield(i)*HeIp(i)

            acoef  = k4(i)*de(i,j,k) + k5(i)*de(i,j,k)

     &             + k25shield(i)

            scoef=max(scoef,1.0e-30)
            if (abs(scoef-HeII(i,j,k)*acoef)/(scoef+HeII(i,j,k)*acoef)
     &          .lt.toler .or. eqmask(i).eqv..true.)then
              HeIIp(i) = scoef/acoef
            else
              HeIIp(i)  = ( scoef*dtit(i) + HeII(i,j,k) )
     &                   / ( 1.0d0 + acoef*dtit(i) )
            endif
c
c          5) HeIII
c
            scoef   = k5(i)*HeIIp(i)*de(i,j,k)

     &              + k25shield(i)*HeIIp(i)

            acoef   = k6(i)*de(i,j,k)
            scoef=max(scoef,1.0e-30)
            HeIIIp(i)  = ( scoef*dtit(i) + HeIII(i,j,k) )
     &                  / ( 1.0d0 + acoef*dtit(i) )
c     
c          6) electrons:
c
            scoef =   k8(i) * DBLE(HM(i,j,k)      * HIp(i)    )
     &             +  k15(i)* DBLE(HM(i,j,k)      * HIp(i)    )
     &             +  k17(i)* DBLE(HM(i,j,k)      * HII(i,j,k))
     &             +  k1(i) * DBLE(HIp(i)         * de(i,j,k) )
     &             +  k3(i) * DBLE(HeIp(i) /4.0d0 * de(i,j,k) )
     &             +  k5(i) * DBLE(HeIIp(i)/4.0d0 * de(i,j,k) )
     &             +  k14(i)* DBLE(HM(i,j,k)      * de(i,j,k) )
c                  

     &             + k24shield(i)*HIp(i)
     &             + k25shield(i)*HeIIp(i)/4.0d0
     &             + k26shield(i)*HeIp(i)/4.0d0

            acoef =   k2(i)   *HII(i,j,k)
     &              + k6(i)   *HeIIIp(i)/4.0d0
     &              + k4(i)   *HeIIp(i)/4.0d0
     &              + k7(i)   *HIp(i)
     &              + k18(i)  *H2II(i,j,k)/2.0d0
c           acoef = - (k1(i) *HIp(i)        - k2(i)*HII(i,j,k)
c    &              +  k3(i) *HeIp(i)/4.0d0 - k6(i)*HeIIIp(i)/4.0d0
c    &              +  k5(i) *HeIIp(i)/4.0d0- k4(i)*HeIIp(i)/4.0d0
c    &              +  k14(i)*HM(i,j,k)
c    &              -  k7(i) *HIp(i)
c    &              -  k18(i)*H2II(i,j,k)/2.0d0)
            if (ispecies.gt.2) then
              scoef = scoef + k1(i) * DI(i,j,k) * de(i,j,k) / 2.0d0
              acoef = acoef + k2(i) * DII(i,j,k) / 2.0d0
            endif
            if (abs(scoef-de(i,j,k)*acoef)/(scoef+de(i,j,k)*acoef)
     &          .lt.toler .or. eqmask(i).eqv..true.)then
              dep(i) = scoef/acoef
            else
            dep(i)  = ( scoef*dtit(i) + de(i,j,k) )
     &                / ( 1.0d0 + acoef*dtit(i) )
            endif
            if ((dep(i).lt.0.0d0).or.(dep(i).ne.dep(i))) then
              write(6,*) 's_r_c: deplt00a:',k8(i),k15(i)
              write(6,*) 's_r_c: deplt00b:',k17(i),k1(i)
              write(6,*) 's_r_c: deplt00c:',k1(i),k3(i)
              write(6,*) 's_r_c: deplt00d:',k5(i),k14(i)
              write(6,*) 's_r_c: deplt00e:',k2(i),k6(i)
              write(6,*) 's_r_c: deplt00f:',k4(i),k7(i)
              write(6,*) 's_r_c: deplt00g:',k18(i),HM(i,j,k)
              write(6,*) 's_r_c: deplt00h:',HIp(i),HeIp(i)
              write(6,*) 's_r_c: deplt00i:',HeIIp(i),HeIIIp(i)
              write(6,*) 's_r_c: deplt00j:',HIIp(i),H2II(i,j,k)
              write(6,*) 's_r_c: deplt00k:',d(i,j,k),scoef,acoef
              write(6,*) 's_r_c: deplt00l:',dtit(i),dep(i),de(i,j,k)
            endif
c
c           7) H-
c
            HMp(i) = ( (k7(i)*HIp(i))*dep(i) )
     &             / ( (k8(i)+k15(i))*HIp(i)
     &             + ( k16(i)+k17(i))*HIIp(i)+k14(i)*dep(i)
c            HMp(i) = ( (k7(i) *HIp(i))*dep(i) )
c     &             / ( (k8(i) *HIp(i) + k15(i)*HIp(i)  )
c     &             +   (k16(i)*HIIp(i)+ k17(i)*HIIp(i) )
c     &             +   (k14(i)*dep(i)                  ) 

     &             + k27

     &                 )
c
c           8) H2+
c
            H2IIp(i) = 2.0d0*( (k9 (i)*HIp(i))*HIIp(i)
     &                       + (k11(i)*H2I(i,j,k)/2.0d0)*HIIp(i)
     &                       + (k17(i)*HMp(i))*HIIp(i)

     &                    + k29*H2I(i,j,k)

     &                    )
     &                 /  ( k10(i)*HIp(i) + k18(i)*dep(i)
     &                    + k19(i)*HMp(i)

     &                    + (k28+k30)

     &                    )
c
c           9) H2
c
            scoef = 2.0d0*((k8(i) * DBLE(HMp(i))   * HIp(i))
     &            +       (k10(i) * DBLE(H2IIp(i)) * HIp(i))/2.0d0
     &            +       (k19(i) * DBLE(H2IIp(i)) * HMp(i))/2.0d0
     &            +   DBLE(k22(i) * HIp(i)  ) * DBLE(HIp(i)**2.0d0))
            acoef = ( k13(i)*HIp(i) + k11(i)*HIIp(i)
     &              + k12(i)*dep(i) + k23(i)*H2I(i,j,k)/2.0d0)

     &              + k29 + k31

c
            if (ispecies.gt.4) then
              scoef = scoef + 2.0d0 * (
     &                          k53(i) * HDI(i,j,k) * HIIp(i) / 3.0d0
     &                        + k55(i) * HDI(i,j,k) * HIp(i)  / 3.0d0
     &                        )
              acoef = acoef + (
     &                          k52(i) * DII(i,j,k) / 2.0d0
     &                        + k54(i) * DI(i,j,k) / 2.0d0
     &                        )
            endif
            if (abs(scoef-H2I(i,j,k)*acoef)/(scoef+H2I(i,j,k)*acoef)
     &          .lt.toler .or. eqmask(i).eqv..true.)then
              H2Ip(i) = scoef/acoef
c              write(6,*) H2I(i,j,k), d(i,j,k), HIp(i), HI(i,j,k), HMp(i), H2Ip(i), scoef, acoef, tgas(i)
            else
            H2Ip(i) = ( scoef*DBLE(dtit(i)) + H2I(i,j,k) )
     &                / ( 1.0d0 + acoef*DBLE(dtit(i)))
            endif
c            H2Ip(i) = (scoef*dtit(i))/(1.0d0 + acoef*dtit(i))
c     &              + (    H2I(i,j,k)/(1.0d0 + acoef*dtit(i)))

         endif
         enddo
      endif                     ! H2
c
c  --- (D) Now do extra 3-species for molecular HD --- c     
      if (ispecies .gt. 2) then
         do i = is+1, ie+1
         if (itmask(i).eqv..true.) then
c
c            2a) DII
c 
            if (ispecies.gt.2) then
              scoef = 2.0d0*(
     &                   k1(i) * dep(i)   * DI(i,j,k)/2.0d0
     &              +  (k50(i) * HIIp(i)) * DI(i,j,k)/2.0d0
     &              +  (k53(i) * HIIp(i)) * HDI(i,j,k)/3.0d0

     &              + k24shield(i)*DI(i,j,k)

     &              )
              acoef = 1.0d0*(
     &                 (k2(i) * dep(i))
     &              +  (k51(i) * HIp(i))
     &              +  (k52(i) * H2Ip(i))/2.0d0
     &                      )
c Set into equilibrium, and then use mass conservation to get DI
c
c              if (abs(scoef-DII(i,j,k)*acoef)/(scoef+DII(i,j,k)*acoef).lt.Dtoler)then
                DIIp(i) = scoef/acoef
c              else
c              DIIp(i) = ( scoef*DBLE(dtit(i)) + DII(i,j,k) )
c     &                  / ( 1.0d0 + acoef*DBLE(dtit(i)))
c              endif
            endif
c     
c             1a) DI
c     
            if (ispecies.gt.2) then
              scoef = 2.0d0*(
     &                  k2(i) * DII(i,j,k) * dep(i)/2.0d0
     &              +  k51(i) * DII(i,j,k) * HIp(i)/2.0d0
     &              +  k55(i) * HDI(i,j,k) * HIp(i)/3.0d0
     &                   )
              acoef = 1.0d0*(
     &                  k1(i) * dep(i)
     &              +  k50(i) * HIIp(i)
     &              +  k54(i) * H2Ip(i)/2.0d0

     &               + k24shield(i)

     &                    )
c              if ((abs(scoef-DI(i,j,k)*acoef)/(scoef+DI(i,j,k)*acoef)).lt.Dtoler)then
                DIp(i) = scoef/acoef
c              else
c                write(6,*) 'DI', (abs(scoef-DI(i,j,k)*acoef)/(scoef+DI(i,j,k)*acoef))
c              DIp(i) = ( scoef*DBLE(dtit(i)) + DI(i,j,k) )
c     &                  / ( 1.0d0 + acoef*DBLE(dtit(i)))
c              endif
c   Actually, let's do this via D-conservation
              DIp(i) = fh*dtoh*d(i,j,k)-(DIIp(i)+(2./3.)*HDI(i,j,k))
              if (DIp(i) .lt. 0.0d0) then
c                write(6,*) 'solve_rate_cool: dilt0a', i,j,k
c                write(6,*) 'solve_rate_cool: dilt0b', DIp(i),DIIp(i),HDI(i,j,k)
c                write(6,*) 'solve_rate_cool: dilt0c', fh*dtoh*d(i,j,k)
                DIp(i) = scoef/acoef
                DIIp(i) = max(fh*dtoh*d(i,j,k)-(DIp(i) + (2./3.)*HDI(i,j,k)),1.0e-30)
c                write(6,*) 'solve_rate_cool: dilt0d', DIp(i),DIIp(i)
              endif
            endif
c
c            9a) HDI
c 
            if (ispecies.gt.2) then
              scoef = 3.0d0*(
     &                       k52(i) * DIIp(i)* H2Ip(i)/4.0d0
     &               +       k54(i) * DIp(i) * H2Ip(i)/4.0d0
     &                   )
              acoef = 1.0d0*(
     &                       k53(i) * HIIp(i)
     &                  +    k55(i) * HIp(i)
     &                    )
c
              if (abs(scoef-HDI(i,j,k)*acoef)/(scoef+HDI(i,j,k)*acoef)
     &          .lt.toler .or. eqmask(i).eqv..true.)then
                HDIp(i) = scoef/acoef
              else
c                write(6,*) 'HDI', (abs(scoef-HDI(i,j,k)*acoef)/(scoef+HDI(i,j,k)*acoef))
                HDIp(i) = ( scoef*DBLE(dtit(i)) + HDI(i,j,k) )
     &                    / ( 1.0d0 + acoef*DBLE(dtit(i)))
              endif
            endif
         endif
         enddo
      endif
c
c   --- (E) Set densities from 1D temps to 3D fields ---
c
      do i = is+1, ie+1
      if (itmask(i).eqv..true.) then
c  
         tempde = HII(i,j,k) + HeII(i,j,k)/4. + HeIII(i,j,k)/2.0d0
         if (ispecies .gt. 1) then
           tempde = tempde - HM(i,j,k) + H2II(i,j,k)/2.0d0

           if((HIp(i).lt.0.0d0).or.   (HIp(i).ne.HIp(i)).or.
     &        (HIIp(i).lt.0.0d0).or.  (HIIp(i).ne.HIIp(i)).or.
     &        (HeIp(i).lt.0.0d0).or.  (HeIp(i).ne.HeIp(i)).or.
     &        (HeIIp(i).lt.0.0d0).or. (HeIIp(i).ne.HeIIp(i)).or.
     &        (HeIIIp(i).lt.0.0d0).or.(HeIIIp(i).ne.HeIIIp(i)).or.
     &        (H2Ip(i).lt.0.0d0).or.  (H2Ip(i).ne.H2Ip(i)).or.
     &        (H2IIp(i).lt.0.0d0).or. (H2IIp(i).ne.H2IIp(i)).or.
     &        (HMp(i).lt.0.0d0).or.   (HMp(i).ne.HMp(i)).or.
     &        (dep(i).lt.0.0d0).or.   (dep(i).ne.dep(i)))then
c     &        (tempde.lt.0.0d0).or.   (tempde.ne.tempde))then
              write(6,*) 'step_rate: lt0a', HIp(i),HI(i,j,k)
              write(6,*) 'step_rate: lt0b', HIIp(i),HII(i,j,k)
              write(6,*) 'step_rate: lt0c', HeIp(i),HeI(i,j,k)
              write(6,*) 'step_rate: lt0d', HeIIp(i),HeII(i,j,k)
              write(6,*) 'step_rate: lt0e', HeIIIp(i),HeIII(i,j,k)
              write(6,*) 'step_rate: lt0f', H2Ip(i),H2I(i,j,k)
              write(6,*) 'step_rate: lt0g', H2IIp(i),H2II(i,j,k)
              write(6,*) 'step_rate: lt0h', HMp(i),HM(i,j,k)
              write(6,*) 'step_rate: lt0i', tempde,de(i,j,k)
              write(6,*) 'step_rate: lt0j', dep(i),tgas(i)
              write(6,*) 'step_rate: lt0k', i,j,k
              write(6,*) 'step_rate: lt0l', d(i,j,k), temp
              write(6,*) 'step_rate: lt0n', DIp(i), DI(i,j,k)
              write(6,*) 'step_rate: lt0o', DIIp(i), DII(i,j,k)
              write(6,*) 'step_rate: lt0p', HDIp(i), HDI(i,j,k)
              if (tempde.gt.0.0d0)
     &            dostop(i) = 1
           endif
           if ((ispecies.gt.2).and.
     &        ((DIp(i).lt.0.0d0).or.
     &         (DIIp(i).lt.0.0d0).or.
     &         (HDIp(i).lt.0.0d0)))then
              write(6,*) 'step_rate: lt1a', DIp(i), DI(i,j,k)
              write(6,*) 'step_rate: lt1b', DIIp(i), DII(i,j,k)
              write(6,*) 'step_rate: lt1c', HDIp(i), HDI(i,j,k)
            endif
         endif

         totalH = HIp(i) + HIIp(i)
         totalHe = HeIp(i) + HeIIp(i) + HeIIIp(i)
         totalD = 0

         if (ispecies.gt.1) then
            totalH = totalH + HMp(i) + H2Ip(i) + H2IIp(i)
         endif
         if (ispecies .gt. 2) then
            totalH = totalH + (1.0/3.0)*HDIp(i)
            totalD = DIp(i) + DIIp(i) + (2.0/3.0)*HDIp(i)
         endif

         totalDen = totalH+totalHe+totalD

c        Do consistency correction here, rather than in make_consistent

         correctH = fh*(totalDen/totalH)
         correctHe = (1.0 - fh)*(totalDen/totalHe)
         correctH = fh*(d(i,j,k)/totalH)
         correctHe = (1.0-fh)*(d(i,j,k)/totalHe)
c         if (abs(correctH-1.0).gt.5e-1) then
c            write(6,*) 'HEAVY CORRECT a', correctH, correctHe
c            write(6,*) 'HEAVY CORRECT b', H2I(i,j,k), H2Ip(i)
c            write(6,*) 'HEAVY CORRECT c', HI(i,j,k), HIp(i)
c            write(6,*) 'HEAVY CORRECT d', HII(i,j,k), HIIp(i)
c            write(6,*) 'HEAVY CORRECT e', tgas(i)
c         endif
         correctH = 1
         correctHe = 1

         HIdot_prev(i) = HI(i,j,k)
         HI(i,j,k)    = HIp(i)*correctH
         if (abs(HI(i,j,k)-HIdot_prev(i))/HI(i,j,k).lt.toler) then
           HIdot_prev(i) = 1.0e-30
         else
           HIdot_prev(i) = DBLE(HI(i,j,k)-HIdot_prev(i))/
     &                     max(abs(DBLE(dtit(i))),1.0e-30)
         endif
         
         HII(i,j,k)   = HIIp(i)*correctH
         HeI(i,j,k)   = HeIp(i)*correctHe
         HeII(i,j,k)  = HeIIp(i)*correctHe
         HeIII(i,j,k) = max(HeIIIp(i)*correctHe, 1.0e-25)
c         HeIII(i,j,k) = max(HeIIIp(i), 1.0e-30)
c
         if (ispecies .gt. 1) then
            H2Idot_prev(i) = H2I(i,j,k)
            HM(i,j,k)    = HMp(i)*correctH
            H2I(i,j,k)   = H2Ip(i)*correctH
c            if (abs(H2I(i,j,k)-H2Idot_prev(i))/H2I(i,j,k).lt.toler) then
c              H2Idot_prev(i) = 1.0e-30
c            else
              H2Idot_prev(i) = DBLE(H2I(i,j,k)-H2Idot_prev(i))/
     &                             max(abs(DBLE(dtit(i))),1.0e-30)
c            endif
            H2II(i,j,k)  = H2IIp(i)*correctH
         endif
c
         if (ispecies .gt. 2) then
            correctD = fh*dtoh*(d(i,j,k)/totalD)
            DIdot_prev(i) = DI(i,j,k)
            DI(i,j,k)    = DIp(i)*correctD
            if (abs(DI(i,j,k)-DIdot_prev(i))/DI(i,j,k).lt.toler) then
              DIdot_prev(i) = 1.0e-30
            else
              DIdot_prev(i) = DBLE(DI(i,j,k)-DIdot_prev(i))/
     &                        max(abs(DBLE(dtit(i))),1.0e-30)
            endif
            DII(i,j,k)   = DIIp(i)*correctD
            HDI(i,j,k)   = HDIp(i)*correctD
            HDIdot_prev(i) = HDI(i,j,k)
            HDI(i,j,k)    = HDIp(i)*correctD
            if (abs(HDI(i,j,k)-HDIdot_prev(i))/HDI(i,j,k).lt.toler) then
              HDIdot_prev(i) = 1.0e-30
            else
              HDIdot_prev(i) = DBLE(HDI(i,j,k)-HDIdot_prev(i))/
     &                        max(abs(DBLE(dtit(i))),1.0e-30)
            endif
         endif
c
c        de(i,j,k)    = dep(i)
c
c        Use charge conservation to determine electron fraction
c 
         dedot_prev(i) = de(i,j,k)
         temp = dedot_prev(i)
         de(i,j,k) = HII(i,j,k) + HeII(i,j,k)/4. + HeIII(i,j,k)/2.0d0
         if (ispecies .gt. 1) 
     &      de(i,j,k) = de(i,j,k) - HM(i,j,k) + H2II(i,j,k)/2.0d0
         if (ispecies .gt. 2)
     &      de(i,j,k) = de(i,j,k) + DII(i,j,k)/2.0d0
         if (abs(de(i,j,k)-dedot_prev(i))/de(i,j,k).lt.toler) then
            dedot_prev(i) = 0.0d0
         else
            dedot_prev(i) = DBLE(de(i,j,k)-dedot_prev(i))/
     &                         max(abs(DBLE(dtit(i))),1.0e-30)
         endif
         if (i.eq.24) then
c            write(6,*) 'correctH:', correctH, H2I(i,j,k)/d(i,j,k), tgas(i), H2Idot_prev(i)
         endif
         if (de(i,j,k).lt.0.0d0) then
c           write(6,*) 'step_rate: delt0a:', k7(i),HIp(i),dep(i)
c           write(6,*) 'step_rate: delt0b:', k8(i),k15(i),HIp(i)
c           write(6,*) 'step_rate: delt0c:', k16(i),k16(i),HIIp(i)
c           write(6,*) 'step_rate: delt0d:', k14(i),de(i,j,k),d(i,j,k)
c           If we are at zero to approximately machine precision, then
c           truncate.
           if (abs(de(i,j,k)/d(i,j,k)).lt.1e-6) then
c             write(6,*) 'solve_rate_cool: de within machine precision'
c             write(6,*) 'solve_rate_cool: resetting to tiny',i,j,k
             de(i,j,k) = 1.0e-30
             dedot_prev(i) = 1.0e-30
             dostop(i)=0
           else
             write(6,*) 'solve_rate_cool: de problematic'
             write(6,*) 'solve_rate_cool: ', i,j,k
             write(6,*) 'solve_rate_cool: ', de(i,j,k), d(i,j,k)
             dostop(i)=1
           endif
         endif

      endif
      enddo                     ! end loop over i

      do i = is+1, ie+1
        if (dostop(i).eq.1)then
          write(6,*) 'Stopping as a result of dostop(i) = 1'
          write(6,*), i,j,k
          stop
        endif
      enddo
c
      return
      end

c ------------------------------------------------------------------
c   This routine correct the highest abundence species to
c     ensure conservation of particle number and charge.
c
      subroutine make_consistent(de, HI, HII, HeI, HeII, HeIII,
     &                        HM, H2I, H2II, DI, DII, HDI, d,
     &                        is, ie, js, je, ks, ke,
     &                        in, jn, kn, ispecies, fh, dtoh, iter)
c -------------------------------------------------------------------
c
      implicit NONE
c
c     Arguments
c
      integer in, jn, kn, is, ie, js, je, ks, ke, ispecies
      real    de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &       HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &         d(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)
      real    DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)  
      real    fh, dtoh
c
c     Parameters
c
      integer ijk, iter
      parameter (ijk = 1031)
c
c     locals
c
      integer i, j, k
      real totalDen, totalH(ijk), totalHe(ijk), 
     &     totalD, correctH, correctHe, correctD
c
c     Loop over all zones
c
      do k = ks+1, ke+1
      do j = js+1, je+1
c
c     Compute total densities of H and He
c         (ensure non-negativity)
c
      do i = is+1, ie+1
         HI   (i,j,k) = abs(HI   (i,j,k))
         HII  (i,j,k) = abs(HII  (i,j,k))
         HeI  (i,j,k) = abs(HeI  (i,j,k))
         HeII (i,j,k) = abs(HeII (i,j,k))
         HeIII(i,j,k) = abs(HeIII(i,j,k))
         totalH(i) = HI(i,j,k) + HII(i,j,k)
         totalHe(i) = HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k)
         totalDen = totalH(i) + totalHe(i)
      enddo
c
c     include molecular hydrogen
c
      if (ispecies .gt. 1) then
         do i = is+1, ie+1
            HM   (i,j,k) = abs(HM   (i,j,k))
            H2II (i,j,k) = abs(H2II (i,j,k))
            H2I  (i,j,k) = abs(H2I  (i,j,k))
            totalH(i) = totalH(i) + HM(i,j,k) + H2I(i,j,k) + H2II(i,j,k)
            totalDen = totalDen + HM(i,j,k) + H2I(i,j,k)
     &                  + H2II(i,j,k)
         enddo
      endif
c
c     Correct densities by keeping fractions the same
c
      do i = is+1, ie+1
c         if (totalDen(i).ne.d(i,j,k))then
c            write(6,*) 'solve_rate_cool: totalDen0a', i,j,k
c            write(6,*) 'solve_rate_cool: totalDen0b', d(i,j,k), 
c     &                 totalDen
c         endif
         correctH = fh*(d(i,j,k)/totalH(i))
         HI(i,j,k)  = HI(i,j,k)*correctH
         HII(i,j,k) = HII(i,j,k)*correctH
c
         correctHe = (1.0 - fh)*(d(i,j,k)/totalHe(i))
         HeI(i,j,k)   = HeI(i,j,k)*correctHe
         HeII(i,j,k)  = HeII(i,j,k)*correctHe
         HeIII(i,j,k) = HeIII(i,j,k)*correctHe
c
c     Correct molecular hydrogen-related fractions
c
         if (ispecies .gt. 1) then
            HM   (i,j,k) = HM(i,j,k)*correctH
            H2II (i,j,k) = H2II(i,j,k)*correctH
            H2I  (i,j,k) = H2I(i,j,k)*correctH
         endif
      enddo
c
c     Do the same thing for deuterium (ignore HD) Assumes dtoh is small
c
      if (ispecies .gt. 2) then
         do i = is+1, ie+1
            DI  (i,j,k) = abs(DI  (i,j,k))
            DII (i,j,k) = abs(DII (i,j,k))
            HDI (i,j,k) = abs(HDI (i,j,k))
            totalD = DI(i,j,k) + DII(i,j,k) + (2.0/3.0)*HDI(i,j,k)
            correctD = fh*dtoh*(d(i,j,k)/totalD)
            DI  (i,j,k) = DI (i,j,k)*correctD
            DII (i,j,k) = DII(i,j,k)*correctD
c            HDI (i,j,k) = HDI(i,j,k)*(3.0/2.0)*correctD
            HDI (i,j,k) = HDI(i,j,k)*correctD
         enddo
      endif
c
c       Set the electron density and check total density
c
      do i = is+1, ie+1
         totalDen = HI(i,j,k) + HII(i,j,k)
         totalDen = totalDen + HeI(i,j,k) + HeII(i,j,k)
     &               + HeIII(i,j,k)
         de (i,j,k) = HII(i,j,k) + HeII(i,j,k)/4. + HeIII(i,j,k)/2.
         if (ispecies .gt. 1) then 
            de(i,j,k) = de(i,j,k) - HM(i,j,k) + H2II(i,j,k)/2.
            totalDen = totalDen + HM(i,j,k) + H2I(i,j,k) 
     &                  + H2II(i,j,k)
         endif
         if (abs(log10(totalDen)-log10(d(i,j,k))).gt.0.01) then
            write(6,*) 'solve_rate_cool: totalDen0a', i,j,k
            write(6,*) 'solve_rate_cool: totalDen0b', d(i,j,k),
     &                 totalDen, iter
         endif
      enddo
c
      enddo  ! end loop over j
      enddo  ! end loop over k
c
      return
      end

c ------------------------------------------------------------------
c   This routine takes an input of fields and energies, and returns
c     the temperature of the gas.
c
      subroutine calculate_tgas(
     &                d, e, ge, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, idual, imethod, ispecies, idim,
     &                is, ie, j, k, aye, temstart, temend,
     &                utem, uxyz, uaye, urho, utim, gamma,
     &                HM, H2I, H2II, DI, DII, HDI, 
     &                tgas, tgasold, p2d, gamma2, totalN, itmask )
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
      integer in, jn, kn, is, ie, j, k, imethod, idim,
     &        idual, ispecies

      real    aye, temstart, temend,
     &        utem, uxyz, uaye, urho, utim,
     &        gamma
      real    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &       de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &      HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      logical itmask(in)
c
c  Parameters
c
      real*8 mh
      real    mintemp
      parameter (mh = 1.67d-24, mintemp=0.001)
c
c  Locals
c
      integer i, j1 
      real dom, logtem0, logtem9, energy,
     &     gamma2(in) 
      real*8 coolunit, dbase1, tbase1, xbase1, totalN(in)
c
c  Slice locals
c 
      integer indixe(in)
      real t1(in), t2(in), logtem(in), tdef(in), p2d(in),
     &     tgas(in), tgasold(in)

c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
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
c
c     Compute Pressure
c

      if (imethod .eq. 2) then
c
c        Zeus - e() is really gas energy
c
         do i = is+1, ie+1
             p2d(i) = (gamma - 1.0)*d(i,j,k)*e(i,j,k)
         enddo
      else
         if (idual .eq. 1) then
c
c           PPM with dual energy -- use gas energy
c
            do i = is+1, ie+1
                p2d(i) = (gamma - 1.0)*d(i,j,k)*ge(i,j,k)
            enddo
         else
c
c           PPM without dual energy -- use total energy
c
            do i = is+1, ie+1
               p2d(i) = e(i,j,k) - 0.5*u(i,j,k)**2
               if (idim .gt. 1) p2d(i) = p2d(i) - 0.5*v(i,j,k)**2
               if (idim .gt. 2) p2d(i) = p2d(i) - 0.5*w(i,j,k)**2
               p2d(i) = max((gamma - 1.0)*d(i,j,k)*p2d(i), 1.0e-30)
            enddo
         endif
      endif

cc Assume totalN is already calculated (true if calculate_mu precedes)
c
c      do i = is+1, ie+1
c          totalN(i) = 
c     &         (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))/4.0 +
c     &         HI(i,j,k) + HII(i,j,k) + de(i,j,k)
c      enddo
cc
cc          (include molecular hydrogen, but ignore deuterium)
cc
c      if (ispecies .gt. 1) then
c         do i = is+1, ie+1
c             totalN(i) = totalN(i) +
c     &            HM(i,j,k) + (H2I(i,j,k) + H2II(i,j,k))/2.0
c         enddo
c      endif
c
      do i = is+1, ie+1
        tgas(i) = max(p2d(i)*(utem/totalN(i)), temstart)
        if(tgas(i).ne.tgas(i))then
           write(6,*) 'solve_rate_cool tgas1a:',p2d(i),totalN(i)
           write(6,*) 'solve_rate_cool tgas1b:',HI(i,j,k),HII(i,j,k),HM(i,j,k)
           write(6,*) 'solve_rate_cool tgas1c:',H2I(i,j,k)
        endif
      enddo
c   Omukai
c
c     Correct temperature for gamma from H2
c
      if (ispecies .gt. 1) then
         call calculate_gamma2(
     &                d, e, ge, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, is, ie, j, k,
     &                HM, H2I, H2II, DI, DII, HDI, tgas,
     &                gamma, gamma2, totalN, itmask )
         do i = is+1, ie+1
           tgas(i) = tgas(i) * (gamma2(i) - 1.0)/(gamma - 1.0)
           if(tgas(i).ne.tgas(i))then
              write(6,*) 'solve_rate_cool tgas2a:',p2d(i),totalN(i)
              write(6,*) 'solve_rate_cool tgas2b:',HI(i,j,k),HII(i,j,k),HM(i,j,k)
              write(6,*) 'solve_rate_cool tgas2c:',H2I(i,j,k)
           endif
         enddo
      endif

      do i = is+1, ie+1
        if (tgas(i).lt.mintemp)then
          tgas(i)=mintemp
        endif
      enddo
      do i = is+1, ie+1
        tgasold(i) = tgas(i)
      enddo
      return
      end

c ------------------------------------------------------------------
c   This routine takes an input of fields and energies, and returns
c     the gamma2 correction.  Should *only* be called when 
c     ispecies.gt.1 and so makes no checks as to that case.
c
      subroutine calculate_gamma2(
     &                d, e, ge, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, is, ie, j, k,
     &                HM, H2I, H2II, DI, DII, HDI, tgas,
     &                gamma, gamma2, totalN, itmask )
c
      integer in, jn, kn, is, ie, j, k

      real    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &       de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &      HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      logical itmask(in)
c
c  Locals
c
      integer i, j1 
      real gamma2(in), x, nH2, nother, tgas(in), gamma
      real*8 totalN(in)

c     Correct temperature for gamma from H2
c
      do i = is+1, ie+1
         nH2 = 0.5*(H2I(i,j,k) + H2II(i,j,k))
c       totalN should be pre-calculated, thus including D species if needed
         nother = totalN(i) - nH2
         if (nH2/nother .gt. 1.0e-3) then
            x = 6100/tgas(i) ! not quite self-consistent
            if (x .gt. 10.0) then
               gamma2(i) = 0.5*5.0
            else
              gamma2(i) = 0.5*(5.0 + 2.0*x**2 * exp(x)
     &                  /(exp(x)-1)**2)
            endif
         else
            gamma2(i) = 2.5
         endif
         gamma2(i) = 1.0 + (nH2 + nother)/
     &        (nH2*gamma2(i) + nother/(gamma - 1.0))
         if ((gamma2(i).ne.gamma2(i)).or.
     &       (gamma2(i).lt.0.0d0))then
           write(6,*)'solve_rate_cool: g20a', gamma2(i),nH2
           write(6,*)'solve_rate_cool: g20b', nother,tgas(i),gamma
           write(6,*)'solve_rate_cool: g20c', H2I(i,j,k), H2II(i,j,k)
           write(6,*)'solve_rate_cool: g20d', HeI(i,j,k), HeII(i,j,k),
     &                                 HeIII(i,j,k)
           write(6,*)'solve_rate_cool: g20e', HI(i,j,k), HII(i,j,k),
     &                                 de(i,j,k)
           write(6,*)'solve_rate_cool: g20f', HM(i,j,k)
           write(6,*)'solve_rate_cool: g20g', i,j,k
         endif
         gamma2(i) = gamma
      enddo
      return
      end


c ------------------------------------------------------------------
c   This routine takes an input of fields and energies, and returns
c     the gamma2 correction
c
      subroutine calculate_mu(
     &                d, e, ge, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, is, ie, j, k,
     &                HM, H2I, H2II, DI, DII, HDI, ispecies,
     &                totalMass, totalN, mu, mudot, dtit, itmask )
c
      integer in, jn, kn, is, ie, j, k

      real    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &       de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &      HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      logical itmask(in)
c
c  Locals
c
      integer i, ispecies
      real*8 totalMass(in), totalN(in), mu(in), mudot(in)
      real dtit(in), toler
      parameter (toler = 1e-20 )
c
c  Perform the first half of our mudot calculation
c  Note that in some extreme circumstances, totalMass != d
c  mainly as a result of difficult integration of the ODEs
c
      do i = is+1, ie+1
        totalMass(i) = 0.0d0
        totalN(i) = 0.0d0
      if (itmask(i).eqv..true.) then
        totalN(i) = ((HI(i,j,k)+HII(i,j,k))
     &              +(HeI(i,j,k)+HeII(i,j,k)+HeIII(i,j,k))/4.0
     &              +de(i,j,k))
        totalMass(i) = HI(i,j,k)+HII(i,j,k)
     &                 +HeI(i,j,k)+HeII(i,j,k)+HeIII(i,j,k)
        if (ispecies.gt.1) then
           totalN(i)=totalN(i)+(H2I(i,j,k)+H2II(i,j,k))/2.0
     &                        +HM(i,j,k)
           totalMass(i)=totalMass(i)+(H2I(i,j,k)+H2II(i,j,k))
     &                              +HM(i,j,k)
        endif
        if (ispecies.gt.2) then
           totalN(i)=totalN(i)+(DI(i,j,k)+DII(i,j,k))
     &                        +(HDI(i,j,k)/3.0)
           totalMass(i)=totalMass(i)+DI(i,j,k)+DII(i,j,k)
     &                              +HDI(i,j,k)
        endif
        totalMass(i) = d(i,j,k)
        if (dtit(i).gt.0.0) then
          mudot(i) = (totalMass(i)/totalN(i) - mu(i))/dtit(i)
          if (abs(dtit(i)*mudot(i)/mu(i)).lt.toler)then
            mudot(i) = 0.0d0
          endif
        else
          mudot(i) = 0.0d0
        endif
        mu(i) = (totalMass(i)/totalN(i))
        if((mudot(i).ne.mudot(i)).or.(abs(mudot(i)).gt.DBLE(1.0e+30)))then
          write(6,*) 'solve_rate_cool: mudot0a',mudot(i),dtit(i),1.0e+30
          write(6,*) 'solve_rate_cool: mudot0b',mu(i),totalMass(i)/totalN(i)
        endif
      endif
      enddo
      return
      end
