subroutine get_basisfn()
    use global
    implicit none
    
    real,parameter :: scaleup = 4
    integer,parameter :: NBsets = 1
    integer,parameter :: nmodes = int(num_u1*modefrac)
    integer,parameter :: nm = nmodes-1

    !arrays
    complex,dimension(0:num_u1-1,0:num_u1-3) :: evecs
    complex,dimension(0:num_u1-1,0:nmodes-1,0:nbsets-1)::ev_temp
    complex,dimension(0:nmodes-1,0:nbsets-1) :: evl_temp
    real,dimension(:),allocatable :: u1_sc,d2co,d1co,d0co
    real,dimension(:,:),allocatable :: a
    integer,dimension(:),allocatable :: indx

    !counters
    integer :: ii,jj,kk,sets

    !variables
    integer :: ixbc_0,ixbc_n,Npts,Ipts
    real :: du1_sc,du1_2,du1_sq,p0,pn

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

    end do

end subroutine get_basisfn
