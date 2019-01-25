!     namelist data

module namelist_mod

  use parameter_mod

    logical :: hall,restart
    logical :: lmadala,lcr,lvs,lweimer,lhwm93,lhwm14

    integer :: psmooth,nion1,nion2
    integer :: maxstep,mmass

    real :: snn(nneut)
    real :: hrmax, dthr, hrpr, dt0, &
            rmin, altmin, fbar, f10p7, ap, &                            
            year, day, hrinit, tvn0, tvexb0, ver, veh, vw,&
            gamss, alt_crit, cqe, alt_crit_avg
    real :: storm_ti, storm_tf, vexb_max, &
            decay_time, pcrit, anu_drag0, &
            blat_max4, stn,denmin,blat_min

    real :: alt_metal,del_alt_metal,deni_mg0,deni_fe0

end module namelist_mod
