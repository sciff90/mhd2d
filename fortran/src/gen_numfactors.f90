subroutine gen_numfactors()
  use global
  implicit none

  !Local Arrays
  double precision,dimension(:,:),allocatable ::sigp_eps_arr,sigh_eps_arr &
    esig_arr,esig2_arr,csig_arr,csig2_arr,ssig_arr,ssig2_arr


  !Allocate Local Arrays
  allocate(sigp_eps_arr(lbound(sigp_arr,dim=1):ubound(sigp_arr,dim=1),&
    lbound(sigp_arr,dim=2):ubound(sigp_arr,dim=2)))

  allocate(sigh_eps_arr(lbound(sigp_arr,dim=1):ubound(sigp_arr,dim=1),&
    lbound(sigp_arr,dim=2):ubound(sigp_arr,dim=2)))

  allocate(esig_arr(lbound(sigp_arr,dim=1):ubound(sigp_arr,dim=1),&
    lbound(sigp_arr,dim=2):ubound(sigp_arr,dim=2)))

  allocate(esig2_arr(lbound(sigp_arr,dim=1):ubound(sigp_arr,dim=1),&
    lbound(sigp_arr,dim=2):ubound(sigp_arr,dim=2)))

  allocate(csig_arr(lbound(sigp_arr,dim=1):ubound(sigp_arr,dim=1),&
    lbound(sigp_arr,dim=2):ubound(sigp_arr,dim=2)))

  allocate(csig2_arr(lbound(sigp_arr,dim=1):ubound(sigp_arr,dim=1),&
    lbound(sigp_arr,dim=2):ubound(sigp_arr,dim=2)))

  allocate(ssig_arr(lbound(sigp_arr,dim=1):ubound(sigp_arr,dim=1),&
    lbound(sigp_arr,dim=2):ubound(sigp_arr,dim=2)))

  allocate(ssig2_arr(lbound(sigp_arr,dim=1):ubound(sigp_arr,dim=1),&
    lbound(sigp_arr,dim=2):ubound(sigp_arr,dim=2)))

  !Allocate Global Factor arrays





end subroutine gen_numfactors
