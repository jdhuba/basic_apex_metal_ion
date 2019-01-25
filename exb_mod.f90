module exb_mod

  use parameter_mod

    real :: vexbs(nzp1,nf,nl,nion)
    real :: vexbp(nz,nfp1,nl,nion)
    real :: vexbh(nz,nf,nlp1,nion)

    real :: vexbs_phi(nzp1,nf,nl)
    real :: vexbp_phi(nz,nfp1,nl)
    real :: vexbh_phi(nz,nf,nlp1)

    real :: u1p(nz,nf,nl),u2s(nz,nf,nl),u3h(nz,nf,nl)

    real :: eps(nzp1,nf,nl),eph(nz,nf,nlp1)
    real :: ehs(nzp1,nf,nl),ehp(nz,nfp1,nl)

    real :: ehpx(nz,nfp1,nl),ehpy(nz,nfp1,nl),ehpz(nz,nfp1,nl)
    real :: ephx(nz,nf,nlp1),ephy(nz,nf,nlp1),ephz(nz,nf,nlp1)

    real :: jp(nz,nf,nl),jphi(nz,nf,nl)

    real :: nuinoci(nzp1,nfp1,nlp1,nion)
    real :: gpoci(nzp1,nfp1,nlp1,nion)
    real :: gsoci(nzp1,nfp1,nlp1,nion)

end module exb_mod
