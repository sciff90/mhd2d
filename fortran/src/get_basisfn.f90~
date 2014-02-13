subroutine get_basisfn()
    use global
    use vector_operations
    implicit none

    double precision,parameter :: scaleup = 4
    integer,parameter :: NBsets = 1
    integer,parameter :: nmodes = int(num_u1*modefrac)
    integer,parameter :: nm = nmodes-1

    !arrays
    double precision,dimension(0:num_u1-1,0:num_u1-3) :: evecs
    double precision,dimension(0:num_u1-1,0:nmodes-1,0:nbsets-1)::ev_temp
    double precision,dimension(0:nmodes-1,0:nbsets-1) :: evl_temp
    double precision,dimension(:),allocatable :: u1_sc,d2co,d1co,d0co,dufac
    double precision,dimension(:,:),allocatable :: a,a2
    integer,dimension(:),allocatable :: indx,index_array

    !counters
    integer :: ii,jj,kk,sets

    !variables
    integer :: ixbc_0,ixbc_n,Npts,Ipts,max_index
    double precision :: du1_sc,du1_2,du1_sq,p0,pn,normfac
    integer :: Mmax,Nmax,LWORK

    !Spepherical Harmonics
    integer :: n ,lda,ldvl,ldvr
    integer :: info,IFAIL,ERROR
    double precision,dimension(:),allocatable :: rwork,wi,wr
    double precision,dimension(:,:),allocatable :: vl,vr,evec
    double precision,dimension(:),allocatable :: w,eval,temp_array
    double precision,dimension(:),allocatable :: work
    double precision :: ABNRM
    double precision,dimension(:),allocatable :: rconde,rcondv,sc
    integer :: ilo,ihi
    double precision,dimension(:,:),allocatable :: Ev,eve,evt,eveinv;
    double precision,dimension(:),allocatable :: Evl,Ev2,nulm;
    integer,dimension(:),allocatable :: IPIV,iwork,Eidx,Oidx

    !Coefficient Arrays
    double precision,dimension(:,:),allocatable :: b3coef,b3b3,psifaci,&
      psifacg,psiicoef,psigcoef
    double precision,dimension(:),allocatable :: nufacg,nufaci,nufac0,temp_trans

    do sets = 0, nbsets-1

        if (sets .eq. 0) then
                ixbc_0 = 1      !Low L Boundary Condition switch for basis functions ixbc eq 1 means derivative =0 BC else zero
                ixbc_n = 1      !High L Boundary Condition switch for basis functions ixbc eq 1 means derivative =0 BC else zero
        end if

        if (sets .eq. 1) then
                ixbc_0 = 0
                ixbc_n = 0
        end if

        Npts = Scaleup*(num_u1-1)+1
        Ipts = Npts-2

        du1_sc = (u1(num_u1-1)-u1(0))/(Npts-1)
        if(sets .eq. 0) then
          Allocate(u1_sc(0:Npts-1))
          Allocate(indx(0:int(npts/scaleup)))
        end if

        do ii = 0, Npts-1
                u1_sc(ii) = ii*du1_sc+u1(0)
                end do

        do ii = 0, int(Npts/scaleup)
                indx(ii) = ii*scaleup
        end do

        du1_2 = 2*du1_sc
        du1_sq = du1_sc**2
        if(sets.eq.0) then
          Allocate(a(0:Ipts-1,0:Ipts-1))          !Finite difference matrix A
          Allocate(a2(0:Ipts-1,0:IPts-1))
          Allocate(d2co(0:Npts-1))
          Allocate(d1co(0:Npts-1))
          Allocate(d0co(0:Npts-1))
        end if
        a = 0;


        d2co    = 4.0*u1_sc*(1+u1_sc)/du1_sq
        d1co    = (2.0+3.0*u1_sc)/du1_sc
        d0co    = m*m/u1_sc

        !Populate inner coeffs into matrix A

        do ii = 2, Ipts-1
                a(ii-2,ii-1) = d2co(ii)-d1co(ii)
                a(ii-1,ii-1) = -(d0co(ii)+2*d2co(ii))
                a(ii,ii-1) = d2co(ii)+d1co(ii)
        end do

        if (ixbc_0 .eq.1) then
                p0 = d2co(1)-d1co(1)
        else
                p0 = 0
        end if

        if(ixbc_n .eq. 1) then
                pn = d2co(npts-2) + d1co(npts-2)
        else
                pn = 0.0
        end if

        a(0,0) = -2.0d0*d2co(1)-d0co(1)+4.d0*p0/3.d0
        a(1,0) = d2co(1)+d1co(1)-p0/3.d0
        a(Npts-3,Npts-3) = -2.0*d2co(Npts-2)-d0co(Npts-2)+4.0*pn/3.0
        a(Npts-4,Npts-3) = d2co(Npts-2)-d1co(Npts-2)-pn/3.0

        lwork = -1
        N = Ipts
        lda = N
        ldvl = N
        ldvr = N

        if(sets.eq.0) then
          allocate(vl(0:ldvl-1,0:n-1))
          allocate(evec(0:Npts-1,0:Npts-3))
          allocate(vr(0:ldvl-1,0:n-1))
          allocate(wr(0:n-1))
          allocate(wi(0:n-1))
          allocate(eval(0:n-1))
          allocate(temp_array(0:n-1))
          allocate(index_array(0:n-1))
          allocate(work(0:100))
          allocate(rwork(0:(2*N)-1))
          allocate(iwork(0:2*N-1))
          allocate(sc(0:N-1))
          allocate(rconde(0:N-1))
          allocate(rcondv(0:N-1))
        end if


        CALL dgeevx('B','V', 'V','N', N, a, LDA, Wr,wi, VL,LDVL,VR, LDVR,ilo,ihi,sc,abnrm,rconde,rcondv,work, LWORK, iwork, INFO)
        LWORK = int(work(0))
        deallocate(work)
        allocate(work(0:lwork-1))
        CALL dgeevx('B','V', 'V','N', N, a, LDA, Wr,wi, VL,LDVL,VR, LDVR,ilo,ihi,sc,abnrm,rconde,rcondv,work, LWORK, iwork, INFO)
        IF( INFO.GT.0 ) THEN
           WRITE(*,*)'The algorithm failed to compute eigenvalues.'
        END IF
        eval = wr
        do ii=1,Npts-2
          do jj = 0,N-1
            evec(ii,jj) = vl(ii-1,jj)
          end do
        end do

        if (ixbc_0 .eq. 1) then
                evec(0,:) = (4.0*evec(1,:) -evec(2,:))/3.0
        else
                evec(0,:) = 0.0
        end if

        if (ixbc_n .eq. 1) then
                evec(Npts-1,:) = (4.0*evec(Npts-2,:)-evec(Npts-3,:))/3.0
        else
                evec(Npts-1,:) = 0.0
        end if

        !Sort Eigenvectors
        !Array index sort function
        do kk = 0,Ipts-1
                max_index = maxloc(real(eval),1)-1
                index_array(Ipts-1-kk) = max_index
                temp_array(Ipts-1-kk) = eval(max_index)
                eval(max_index) = -1
        end do
        eval = temp_array

        do kk = 0,num_u1-3
                !Sort and resample to smaller grid
                do jj = 0,int(npts/scaleup)
                        evecs(jj,kk) = evec(indx(jj),index_array(kk))
                end do
        end do
        ev_temp(:,0:nm,sets) = evecs(:,0:nm)
        evl_temp(0:nm,sets) = eval(0:nm)

     end do

     allocate(evl(0:NBsets*nmodes-1))
     allocate(ev(0:num_u1-1,0:NBsets*nmodes-1))
     allocate(ev2(0:num_u1-1))


     do ii = 0, NBsets*nmodes-1
        if (ii<size(evl_temp,1)) then
          evl(ii) = evl_temp(ii,0)
          ev(:,ii) = ev_temp(:,ii,0)
        else
          evl(ii) = evl_temp(ii-size(evl_temp,1),1)
          ev(:,ii) = ev_temp(:,ii-size(evl_temp,1),1)
        end if
    end do

    write(*,*)"number of basis functions generated = ",nmodes*nbsets-1


!    deallocate(vl)
    !deallocate(evec)
    !deallocate(vr)
    !deallocate(w)
    !deallocate(eval)
    !deallocate(temp_array)
    !deallocate(index_array)
    !deallocate(work)
    !deallocate(rwork)
    !deallocate(sc)
    !deallocate(rconde)
    !deallocate(rcondv)

    !Normalize the basis functions
    allocate(dufac(0:num_u1-1))
    dufac = del_u1/(2.0*sqrt(1.0+u1))
    dufac(0) = dufac(0)*0.5;
    dufac(num_u1-1) = dufac(num_u1-1)*0.5

    do ii = 0,NBsets*nmodes-1
      ev2 = ev(:,ii)*ev(:,ii)*dufac
      normfac = sqrt(sum(ev2))
      ev(:,ii) = ev(:,ii)/normfac
    end do

    allocate(evt(0:NBsets*nmodes-1,0:num_u1-1))
    allocate(eve(0:NBsets*nmodes-1,0:NBsets*nmodes-1))
    allocate(eveinv(0:NBsets*nmodes-1,0:NBsets*nmodes-1))
    allocate(nulm(0:NBsets*nmodes-1))


    nulm = (sqrt(1.0+4.0*evl)-1.0)*0.5;
    evt = transpose(ev)
    eve = matmul(evt,ev)
    eveinv = eve

    !Calculate Inverse

    Mmax = nmodes*nbsets
    Nmax = nmodes*nbsets
    LDA = Mmax
    LWORK = 64*Nmax

    deallocate(WORK)
    allocate(IPIV(0:Nmax-1))
    allocate(WORK(0:LWORK-1))

    write(*,*) 'Calculating Inverse...'
    call dgetrf(Mmax,Nmax,eveinv,LDA,IPIV,INFO)

    if (info.eq.0) then
      call dgetri(Nmax,eveinv,LDA,IPIV,WORK,LWORK,INFO)
    else
      write(*,*) 'The Factor U is Singular'
    end if

    write(*,*) 'Matrix Inversion Completed...'

    !For fitting at only Br locations
    Allocate(Eidx(0:num_u1/2))
    Allocate(Oidx(0:num_u1/2-1))

    !Array index of odd and even gridpoints
    do ii = 0, num_u1/2
      Eidx(ii) = ii*2
      if(ii*2<num_u1-1)then
        Oidx(ii) = ii*2+1
      end if
    end do

    !Allocate Coefficient Arrays
    Allocate(b3coef(0:NBsets*nmodes-1,0:num_u1/2))
    Allocate(b3b3(0:num_u1-1,0:num_u1/2))
    Allocate(nufac0(0:NBsets*nmodes-1))
    Allocate(nufaci(0:NBsets*nmodes-1))
    Allocate(nufacg(0:NBsets*nmodes-1))
    Allocate(psifaci(0:NBsets*nmodes-1,0:num_u1/2))
    Allocate(psifacg(0:NBsets*nmodes-1,0:num_u1/2))
    Allocate(psiicoef(0:NBsets*nmodes-1,0:num_u1/2))
    Allocate(psigcoef(0:NBsets*nmodes-1,0:num_u1/2))
    Allocate(psiib3_S(0:num_u1-1,0:num_u1/2))
    Allocate(psigb3_S(0:num_u1-1,0:num_u1/2))
    Allocate(psiib3_N(0:num_u1-1,0:num_u1/2))
    Allocate(psigb3_N(0:num_u1-1,0:num_u1/2))
    Allocate(temp_trans(0:num_u1/2))


    b3coef = matmul(eveinv,evt(:,Eidx))
    b3b3 = matmul(ev,b3coef)

    nufac0 = nulm*r_iono*(1.d0 - 1.d0/r_iono**(2*nulm+1))
    nufaci = (1.d0 + nulm/(nulm+1)/r_iono**(2*nulm+1))/nufac0
    nufacg = (2*nulm+1)/(nulm+1)/r_iono**nulm/nufac0

    !For Southern Ionosphere/Atmopshere
    psifaci = matmul_GDS(hratm_S(Eidx),nufaci)
    psifacg = matmul_GDS(hratm_S(Eidx),nufacg)
    psiicoef = b3coef*psifaci
    psigcoef = b3coef*psifacg
    psiib3_S = matmul(ev,psiicoef)   !now b3b3 should be nxp1 x nxp1 unit matrix
    psigb3_S = matmul(ev,psigcoef)

    !!For Northern Ionosphere/Atmopshere
    psifaci = matmul_GDS(hratm_N(Eidx),nufaci)
    psifacg = matmul_GDS(hratm_N(Eidx),nufacg)
    psiicoef = b3coef*psifaci
    psigcoef = b3coef*psifacg
    psiib3_N = matmul(ev,psiicoef)  !now b3b3 should be nxp1 x nxp1 unit matrix
    psigb3_N = matmul(ev,psigcoef)

    write(*,*)"done"
end subroutine get_basisfn




