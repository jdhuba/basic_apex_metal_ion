    module apex_module

    use parameter_apex_mod

    implicit none
!      include 'params_apex.inc'


    integer,parameter   :: &
    nal=nalt,  & ! 85
    nlo=nlon,  & ! 180
    nla=nlat  ! 90

    real :: gpalt(nalt)    ! 124
    real :: gplat(nlatp1)  ! 91
    real :: gplon(nlonp2)  ! 182
    integer,parameter :: &
    lwk = nlatp1*nlonp2*nalt*5 + nlatp1+nlonp2+nalt
    real :: wk(lwk)

    contains
!-----------------------------------------------------------------------
    subroutine apxparm(date)
    use apex,only: apex_mka,apex_mall,apex_q2g
!      use parameter_apex_mod

! Original file for this subroutine
! ~richmond/prog/tgcmtst/modsrc/apxparm.mk was copied on 2/25/00.

! 5/02 B.Foster: adapted for tiegcm1.

! Args:
    real,intent(in) :: date

! Local:

! msgun = Fortran unit number for diagnostics
! ngrf  = Number of epochs in the current DGRF/IGRF; see COFRM in
!         file magfld.f

    integer,parameter :: &
    msgun=6, &
    ngrf=8

! Local:
    integer :: j,i,ii,ist,jjm,jjg
    integer :: iplot=0
    real :: rekm,h0km,alt,hr,ror03,glat,glon,dellat,bmag,alon,xlatm, &
    vmp,w,d,be3,sim,xlatqd,f,xlonmi,qdlon,qdlat,gdlat,xlongi,frki, &
    frkj,dellon,si,gdlon,date1(1),re
    real :: h0,hs,pi,rtd,dtr,h00
    character(len=80) :: title

! Local:
    real,parameter :: e=1.e-6, r1=1.06e7, alfa=1.668
    real :: &
    tanth0(nmlat), &
    tanths(nmlat), &
    theta0(nmlat), &
    rcos0s(nmlat), &
    hamh0(nmlat)
    real :: &
    wt(4,nmlonp1,nmlat), &
    dim(nlonp1,0:nlatp1), &
    djm(nlonp1,0:nlatp1)
    real :: dtheta,tanths2
    real :: fac,rmin,rmax,rmag
    real :: dlatg,dlong,dlatm,dlonm,dmagphrlon,r0
    real(8) :: pi_dyn
    real :: mageq_lat(72),mageq_lon(72)


! Specify grid values
! Center min, max altitudes about 130 km
! nalt = 124

!     gpalt = (/87.0000,  94.4443,  103.558,  113.209,  123.492,
!    |  134.419,  145.968,  158.080,  170.645,  183.130,  194.626,
!    |  205.053,  215.141,  225.241,  235.341,  245.441,  255.540,
!    |  265.639,  275.740,  285.840,  295.940,  306.040,  316.141,
!    |  326.240,  336.340,  346.441,  356.540,  366.640,  376.741,
!    |  386.840,  396.939,  407.039,  417.139,  427.239,  437.339,
!    |  447.439,  457.539,  467.638,  477.737,  487.837,  497.938,
!    |  508.037,  518.136,  528.236,  538.335,  548.436,  558.535,
!    |  568.635,  578.735,  588.834,  598.936,  609.034,  619.133,
!    |  629.232,  639.332,  649.431,  659.531,  669.631,  679.729,
!    |  689.830,  699.929,  710.029,  720.128,  730.227,  740.327,
!    |  750.427,  760.525,  770.624,  780.724,  790.823,  800.923,
!    |  811.021,  821.122,  831.220,  841.320,  851.419,  861.518,
!    |  871.617,  881.716,  891.814,  901.914,  912.013,  922.112,
!    |  932.212,  942.311,  952.411,  962.510,  972.608,  982.707,
!    |  992.806,  1002.90,  1013.00,  1023.10,  1033.20,  1043.30,
!    |  1053.40,  1063.50,  1073.60,  1083.70,  1093.79,  1103.89,
!    |  1113.99,  1124.09,  1134.19,  1144.29,  1154.39,  1164.48,
!    |  1174.58,  1184.68,  1199.77,  1265.92,  1463.17,  1806.23,
!    |  2246.80,  2767.43,  3383.29,  4113.77,  4984.00,  6026.80,
!    |  7285.54,  8818.13,  10702.4,  13043.4,  15000. /)

! nalt = 85
    gpalt = (/10., 87.00, 90.00,   92.16,  94.37,  96.67,  99.10, &
    &    101.69,  104.48,  107.50,  110.80, 114.44, 118.50, 123.10, &
    &    128.40,  134.58,  141.89,  150.55, 160.65, 172.06, 184.34, &
    &    196.83,  208.81,  219.80,  229.62, 238.45, 246.64, 254.55, &
    &    262.45,  270.49,  278.76,  287.29, 296.10, 305.22, 314.68, &
    &    324.48,  334.67,  345.28,  356.33, 367.87, 379.94, 392.59, &
    &    405.89,  419.90,  434.70,  450.38, 467.05, 484.85, 503.93, &
    &    524.49,  546.75,  571.00,  597.60, 626.96, 659.62, 696.40, &
    &    738.92,  790.87,  860.07,  960.78, 1114.21, 1345.61, 1677.69, &
    &    2123.44, 2682.74, 3344.79, 4094.94, 4921.55, 5819.94, 6792.88, &
    &    7849.50, 9004.43, 10278.11, 11698.49, 13303.99, 15147.06, &
    &    17296.33, 29833.68, 42840.96, 86373.68, 130426.70, 334906.46, &
    &    848000.64, 4.e7,  4.e10/)

    pi_dyn=3.14159265358979312
    dmagphrlon = 360./float(nmagphrlon)
    re = 6.378165e8   ! earth radius for apex
!     h0 = 9.7e6        ! reference height
!    h0 = 5.7e6        ! reference height
    h00 = h0 * 1.e5        ! reference height set as km in parameter_apex_mod
    hs = 1.3e7
    r0 =re+h00
    pi = 4.*atan(1.)
    rtd = 180./pi                   ! radians to degrees
    dtr = pi/180.                   ! degrees to radians
! Set grid deltas:
    dlatg = pi/float(nlat)
    dlong = 2.*pi/float(nlon)
    dlatm = pi_dyn/float(nmlat-1) ! note use of pi_dyn
    dlonm = 2.*pi_dyn/float(nmlon)
    dmagphrlon = 360./float(nmagphrlon)

    dellat = 180./float(nlatp1-1)
    do j=1,nlatp1
        gplat(j) = (j-1)*dellat - 90.
    enddo
    dellon = 360./float(nlonp2-1)
    do i=1,nlonp2
        gplon(i) = (float(i)-1.0)*dellon - 180.
    enddo

!$$$      write(6,"('apxparm: nlatp1=',i4,' gplat=',/,(8f9.2))")nlatp1,gplat
!$$$      write(6,"('apxparm: nlonp2=',i4,' gplon=',/,(8f9.2))")nlonp2,gplon
!$$$      write(6,"('apxparm: nalt  =',i4,' gpalt=',/,(8g9.1))") nalt,gpalt

!  Initialize interpolation arrays, but do not write them

    date1(1) = date
!     write(6,"('apxparm call apxmka: date = ',f7.0)") date

!     SUBROUTINE APXMKA (MSGUN,EPOCH,nepoch,GPLAT,GPLON,GPALT,NLAT,NLON,
!    |                   NALT,WK,LWK, IST)
!     call apxmka (msgun, date1,1,gplat,gplon,gpalt,nlatp1,nlonp2,nalt,
!    |             wk,lwk,ist)

! subroutine apex_mka(date,gplat,gplon,gpalt,nlat,nlon,nalt,ier)
!!!!      write(6,"('apxparm call apex_mka: date = ',f7.0)") date

    call apex_mka(date,gplat,gplon,gpalt,nlatp1,nlonp2,nalt,ist)

!!!!      write(6,"('apxparm after apex_mka: ist = ',i3)") ist
    if (ist /= 0) call shutdown('apxmka')

! Caculate the geolat and geolon for magnetic equator
!     do i=1,72
!       qdlat = 0.
!       qdlon = 5*i
!       call apex_q2g(qdlat,qdlon,alt,gdlat,gdlon,ist)
!       mageq_lat(i) = gdlat
!       mageq_lon(i) = gdlon
!     enddo
!     write(6,"('mageq_lat(deg) = ',/,(8f8.2))") mageq_lat(:)
!     write(6,"('mageq_lon(deg) = ',/,(8f8.2))") mageq_lon(:)

    end subroutine apxparm
    end module apex_module
