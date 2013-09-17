subroutine get_basisfn()
    use global
    implicit none
    
    double precision,parameter :: scaleup = 4
    integer,parameter :: NBsets = 1
    integer,parameter :: nmodes = int(num_u1*modefrac)
    integer,parameter :: nm = nmodes-1

    !arrays
    complex,dimension(0:num_u1-1,0:num_u1-3) :: evecs
    complex,dimension(0:num_u1-1,0:nmodes-1,0:nbsets-1)::ev_temp
    complex,dimension(0:nmodes-1,0:nbsets-1) :: evl_temp
    double precision,dimension(:),allocatable :: u1_sc,d2co,d1co,d0co
    complex,dimension(:,:),allocatable :: a,a2
    integer,dimension(:),allocatable :: indx,index_array

    !counters
    integer :: ii,jj,kk,sets

    !variables
    integer :: ixbc_0,ixbc_n,Npts,Ipts,max_index
    double precision :: du1_sc,du1_2,du1_sq,p0,pn

    !Spepherical Harmonics
    integer :: n ,lda,ldvl,ldvr
    integer :: info,lwork
    complex,dimension(:),allocatable :: rwork(:)
    complex,dimension(:,:),allocatable :: vl,vr,evec    
    complex,dimension(:),allocatable :: w,eval,temp_array
    complex,dimension(:),allocatable :: work
    double precision :: ABNRM
    double precision,dimension(:),allocatable :: rconde,rcondv,sc
    integer :: ilo,ihi

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
        Allocate(u1_sc(0:Npts-1))
        Allocate(indx(0:int(npts/scaleup)))

        do ii = 0, Npts-1
                u1_sc(ii) = ii*du1_sc+u1(0)
        end do       

        do ii = 0, int(Npts/scaleup)
                indx(ii) = ii*scaleup
        end do

        du1_2 = 2*du1_sc
        du1_sq = du1_sc**2

        Allocate(a(0:Ipts-1,0:Ipts-1))          !Finite difference matrix A
        Allocate(a2(0:Ipts-1,0:IPts-1))
        Allocate(d2co(0:Npts-1))
        Allocate(d1co(0:Npts-1))
        Allocate(d0co(0:Npts-1))


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
        
        allocate(vl(0:ldvl-1,0:n-1))
        allocate(evec(0:Npts-1,0:Npts-3))
        allocate(vr(0:ldvl-1,0:n-1))
        allocate(w(0:n-1))    
        allocate(eval(0:n-1))
        allocate(temp_array(0:n-1))
        allocate(index_array(0:n-1))
        allocate(work(0:100))
        allocate(rwork(0:(2*N)-1))
        allocate(sc(0:N-1))
        allocate(rconde(0:N-1))
        allocate(rcondv(0:N-1))


        CALL cgeevx('B','V', 'V','N', N, a, LDA, W, VL,LDVL,VR, LDVR,ilo,ihi,sc,abnrm,rconde,rcondv,work, LWORK, RWORK, INFO)
        LWORK = int(work(0))
        deallocate(work)
        allocate(work(0:lwork-1))
        CALL cgeevx('B','V', 'V','N', N, a, LDA, W, VL,LDVL,VR, LDVR,ilo,ihi,sc,abnrm,rconde,rcondv,work, LWORK, RWORK, INFO)
        IF( INFO.GT.0 ) THEN
           WRITE(*,*)'The algorithm failed to compute eigenvalues.'           
        END IF
        eval = w
        evec = vl(1:Npts-2,:)

        !open (unit=10, file='../evec.dat' ,form = 'unformatted')
        !write(10)evec
        !close(10)

        if (ixbc_0 .eq. 1) then
                evec(:,0) = (4.0*evec(1,:) -evec(2,:))/3.0
        else
                evec(:,0) = 0.0
        end if

        if (ixbc_n .eq. 1) then
                evec(:,Npts-1) = (4.0*evec(Npts-2,:)-evec(Npts-3,:))/3.0
        else
                evec(:,Npts-1) = 0.0
        end if

        !Sort Eigenvectors
        !Array index sort function
        do kk = 0,Ipts-1
                max_index = maxloc(real(eval),1) - 1
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
        

        open (unit=10, file='../evecs.dat' ,form = 'unformatted')
        write(10)evecs
        close(10)
    end do

    write(*,*) "done"






end subroutine get_basisfn
