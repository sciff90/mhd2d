subroutine get_Va()
        use interpolation
        use global
        implicit none
        
        !arrays
        double precision,dimension(0:num_u1-1,0:num_u3-1) :: rho,Va,No_E,No_F1,No_F2,No_Mag,LVal_arr
        double precision,dimension(0:num_u1-1,0:num_u3-1) :: Lat,b0,nH_eq,rho_mag,No_tot,Rho_a      
        double precision,dimension(0:num_u1-1,0:num_u3-1) :: omega_p,omega_e
        double precision,dimension(0:num_u1-1):: xd,tanhm,tterm,va_sp,b_eq,nH
        double precision,dimension(0:num_u1-1,0:num_u3-1) :: Ne_cc,Nu_ei,Nu_en,Nu_in,Nu_e,Nu_i,lambda
        double precision,dimension(0:num_u1-1) :: Thick_Sd,Thick_Sh,Thick_Sp 
        double precision,dimension(0:num_u1-1,0:num_u3-1) :: V_arr,Coll_Mod,Eperp,Wave_sp,Va_arr
        double precision,dimension(0:num_u1-1,0:Num_u3-1) :: courx,courz
        integer,dimension(0:1) :: i_loc
        !double precision,dimension(0:num_u1-1,0:7) :: flr_a

        double precision,dimension(0:num_u1-1) :: cosalpha_N,sin2alpha_N,cosalpha_S,sin2alpha_S
        double precision,dimension(0:num_u1-1) :: musigp,musigh,musig0,musigzz
        double precision,dimension(0:num_u1-1) :: sigtt,sigtp,sigpt,sigpp
        double precision,dimension(0:num_u1-1) :: sig11,sig12,sig21,sig22,sigden
        double precision,dimension(0:num_u1-1) :: e1b1atm_N,e1b2atm_N,e2b1atm_N,e2b2atm_N
        double precision,dimension(0:num_u1-1) :: e1b1atm_S,e1b2atm_S,e2b1atm_S,e2b2atm_S


        !Ionospheric Data File Arrays
        double precision,dimension(0:num_u1-1,0:num_u3-1) :: No_N,Mav,Temp
        complex*16,dimension(0:num_u1-1) :: sigpatm_N,sighatm_N,sig0atm_N
        complex*16,dimension(0:num_u1-1) :: sigpatm_S,sighatm_S,sig0atm_S


        !Data File Variables
        integer,parameter :: pt = 29  
        double precision,dimension(0:pt+5,0:3):: neutral

        double precision::SclHt

        !Parameters
         double precision,parameter:: Eregion = 100e3
         double precision,parameter:: Mol_H = 1.0e-3
         double precision,parameter:: Av_n = 6.023e23
         double precision,parameter:: Mp = 1.67e-27
         double precision,parameter:: Me = 9.109e-31
         double precision,parameter:: q_i = 1.602e-19
         double precision,parameter:: q_e = -1.602e-19
         double precision,parameter:: k_0 = 8.0e15

         !Counters
         integer::ii,jj

         !Variables
         double precision :: Plaw,Ne_E,Alt_E,Esht_t,Esht_b,Ne_F1,Alt_F1,F1sht_t,F1sht_b
         double precision :: Ne_F2,Alt_F2,F2sht_t,F2sht_b,X
         double precision :: Scale_Sd,Scale_Sh,Scale_Sp,Tol,dtx,dtz,dty 
         double precision :: nt,dtxmax,dtzmax,cmax,cmaz,musigfac
         integer :: i1x,i1z,i3x,i3z

         !Allocations
         allocate(eps_arr(0:num_u1-1,0:Num_u3-1))
         allocate(va2_arr(0:num_u1-1,0:Num_u3-1))
         
         !Magnetic Field Strength
         do ii = 0, num_u3-1
                 Lval_arr(:,ii) = x_arr(:,num_u3_half)
         end do

         Lat = pi/2 - acos(costh)
         b0 = sqrt(4.0-3.0*cos(Lat)**2)/(cos(Lat)**6) * k_0/(LVal_arr*(Re))**3        
         
         !Using colin's FLR Va Functions
         xd = Lval_arr(:,num_u3_half)

         tanhm = 4.0*(xd-5.6)   !Plasma pause at L = 5.6
         tterm = (157.0/xd**2.6)*(tanh(tanhm)+1.0)/1.6
         va_sp = (2.09/xd**1.1+tterm)*990.0-50.0
         va_sp = va_sp*1000.0

         b_eq = b0(:,num_u3_half)
         nH = b_eq**2/(u0_u*va_sp**2)   !Equatorial plasma density
         do ii = 0, num_u3-1
                 nH_eq(:,ii) = nH
         end do

         Plaw = 4.0     !Power law of density profile
         rho_mag = nH_eq*(Lval_arr/r_arr)**Plaw
         No_Mag = rho_mag/mp !Number of electrons and protons in magnetosphere
         
         do ii = 0, num_u1-1
                do jj = 0, num_u3-1
                        !E Region
                        Ne_E = 1.5e11   !Numbeer of electron at E Region peek
                        Alt_E = Eregion !Altitude of e region max
                        Esht_t = 300.0e3!Scale Height above E region
                        Esht_b = 20.0e3 !Scale Height below E region

                        if (z_arr(ii,jj) .gt. Alt_E) then
                                X = (z_arr(ii,jj)-Alt_E)/Esht_t
                        end if

                        if (z_arr(ii,jj) .le. Alt_E) then
                                X = (z_arr(ii,jj)-Alt_E)/Esht_b
                        end if
                        No_E(ii,jj) = Ne_E*exp(1.0-X-exp(-X))

                        !F1 Region
                        Ne_F1 = 2.5e11  !Number of electron at F1 peek (per m^3)
                        Alt_F1 = 200.0e3!Altitude of F1 region (in m)
                        F1sht_t = 200.0e3!Scale Height above F1 region (in m)
                        F1sht_b = 40.0e3!Scale Height below F1 region (in m)

                        if (z_arr(ii,jj) .gt. Alt_F1) then
                                X = (z_arr(ii,jj)-Alt_F1)/F1sht_t
                        end if

                        if (z_arr(ii,jj) .le. Alt_F1) then
                                X = (z_arr(ii,jj)-Alt_F1)/F1sht_b
                        end if
                        No_F1(ii,jj) = Ne_F1*exp(1.0-X-exp(-X))
                        
                        !F2 Region
                        Ne_F2 = 2.0e12  !Number of electron at F1 peek (per m^3)
                        Alt_F2 = 350.0e3!Altitude of F1 region (in m)
                        F2sht_t = 175.0e3!Scale Height above F1 region (in m)
                        F2sht_b = 75.0e3!Scale Height below F1 region (in m)

                        if (z_arr(ii,jj) .gt. Alt_F2) then
                                X = (z_arr(ii,jj)-Alt_F2)/F2sht_t
                        end if

                        if (z_arr(ii,jj) .le. Alt_F2) then
                                X = (z_arr(ii,jj)-Alt_F2)/F2sht_b
                        end if
                        No_F2(ii,jj) = Ne_F2*exp(1.0-X-exp(-X))

                end do
           end do
           No_tot = No_E+No_F1+No_F2+No_mag     !Total number of electrons and protons
           
           !Total mass Density (Assuming Hydrogen Plasma)
           Rho_a = (No_tot/Av_N)*Mol_H  !Mass density in kg/m^3

           !Read ionospheric data from file

           !Gyro Frequency 
           omega_p = q_i*b0/(Mp)        !Gyrofrequency of Protons
           omega_e = q_e*b0/(Me)        !Gyro Frequency of electrons

           !Neutral atmosphere for collision frequency (From kelly)
           SclHt = 2673.0
           open(1,file = 'iono_data/NeutralAtm_Min_Mod.txt',status='old',access='sequential',form='formatted',action='read')
                read(1,*)
                read(1,*)
                read(1,*)
                do ii = 1, pt
                        read(1,*) neutral(ii,0),neutral(ii,1),neutral(ii,2),neutral(ii,3)
                end do                
           close(1)

           do jj = 0, 3
                   neutral(0,jj) = neutral(1,jj)                   
           end do
           
           neutral(ii,0) = 3000
           neutral(ii+1,0) = 4000
           neutral(ii+2,0) = 5000
           neutral(ii+3,0) = 6000
           neutral(ii+4,0) = 100000

           do jj = 0,4 
                   neutral(jj+ii,1) = neutral(ii-1,1)
                   neutral(jj+ii,2) = neutral(ii-1,2)
                   neutral(jj+ii,3)= neutral(ii-1,3)*exp(-(neutral(ii-1+jj,0)-neutral(pt,0))/SclHt)
           end do

           neutral(0,0) = 0
          
           do ii = 0, Num_u1-1

                call quad_inter(neutral(:,3),neutral(:,0)*1.0e3,z_arr(ii,:),No_N(ii,:)) !quadratic interpolation routine
                call quad_inter(neutral(:,2),neutral(:,0)*1.0e3,z_arr(ii,:),Mav(ii,:))
                call quad_inter(neutral(:,1),neutral(:,0)*1.0e3,z_arr(ii,:),Temp(ii,:))               

           end do


           !Collision Frequencies

           where (mav <0) 
               mav = 1
           end where

           Ne_cc = No_tot/1.0e6         !# electrons (per cc)
           
           !from Kelley (possible mistake)
           !Nu_ei	= [34.0+4.18*ALog((Temp^3)/Ne_cc)]*Ne_cc*Temp^(-1.5)
           !Nu_en	= 5.4e-10*No_N*sqrt(Temp)
           !Nu_in	= 2.6e-9*(No_N+Ne_cc)/sqrt(Mav)

           !From my 1D code (Zhang and Cole)
           Nu_ei = (34.0+4.18*ALog((real(Temp)**3))/Ne_cc)*Ne_cc*real(Temp)**(-1.5)
           Nu_en = 5.4e-10*No_N*sqrt(Temp(:,:))
           Nu_in = 2.6e-9*(No_N+Ne_cc)*sqrt(Mav)

           Nu_e = Nu_en + Nu_ei
           Nu_i = Nu_in !+ Nu_ei                !Leaving Nu_ei in Nu_i leads to small negatives in Hall but gives reasonable profile

           lambda = sqrt(Me/(4.0*pi*1.0e-7*No_tot*q_e**2))

           !Conductivities from 1D code
           allocate(sig0_arr(0:num_u1-1,0:Num_u3-1))
           allocate(sigp_arr(0:num_u1-1,0:Num_u3-1))
           allocate(sigh_arr(0:num_u1-1,0:Num_u3-1))
           allocate(eta_arr(0:num_u1-1,0:Num_u3-1))

           sig0_arr = No_tot*(q_e**2/(Me*Nu_e) + q_i**2/(Mp*Nu_i))        !from 1D code
           sigp_arr = (No_tot*q_i**2/(Mp))* (Nu_i)/((Nu_i)**2+omega_p**2) + (No_tot*q_e**2/(Me))* (Nu_e)/((Nu_e)**2+omega_e**2)
           sigh_arr = -( (No_tot*q_i**2/(Mp))* (omega_p)/((Nu_i)**2+omega_p**2) + (No_tot*q_e**2/(Me))* (omega_e)/((Nu_e)**2+omega_e**2) )
           eta_arr = (1.0/(u0_u*sig0_arr))/Re**2                           ! Re^2 comes from Scaling? (it works!)

           !Conductances is thin sheet 
           Thick_Sd = z_arr(:,0)
           Scale_Sd = 1.0/0.2e3         !1/ Direct scale Height (in m) 
           sig0atm_N = sig0_arr(:,0)*((1.0-Exp(-Scale_sd*Thick_sd))/Scale_sd)
           sig0atm_S = sig0_arr(:,Num_u3-1)*((1.0-Exp(-Scale_sd*Thick_sd))/Scale_sd)
           
           Thick_Sp = z_arr(:,0)           
           Scale_Sp = 1.0/0.2e3         !1/ Pederson scale Height (in m)
           sigpatm_N = sigp_arr(:,0)       *((1.0-Exp(-Scale_sp*Thick_sp))/Scale_sp)
           sigpatm_S = sigp_arr(:,Num_u3-1)*((1.0-Exp(-Scale_sp*Thick_sp))/Scale_sp)
           
           Thick_Sh = z_arr(:,0)           
           Scale_Sh = 1.0/0.2e3         !1/ Hall scale Height (in m)
           sighatm_N = sigh_arr(:,0)       *((1.0-Exp(-Scale_sh*Thick_sh))/Scale_sh)
           sighatm_S = sigh_arr(:,Num_u3-1)*((1.0-Exp(-Scale_sh*Thick_sh))/Scale_sh)
           
           Tol = 0.0e-0
           where (sigh_arr .le. Tol) 
               sigh_arr = Tol
           end where

           Tol = minval(eta_arr)
           !zero = where (eta_arr(0.3*Num_u1:0.7*Num_u1,*) GE Tol)+0.3*(Num_u1-1)
           !eta_arr(*,0.3*Num_u3:0.7*Num_u3) = Tol
           
           !!!	Atmospheric thin sheet conductances
           !sigpatm_N(*) = 1.05
           !sigpatm_S(*) = 1.05
           !sighatm_N(*) = 1.05
           !sighatm_S(*) = 1.05
           !sig0atm_N(*) = 500.0
           !sig0atm_S(*) = 500.0
           
           
           !!Height Distributed Conductivity Arrays (set to zero if thin sheet is only required)
           !sigp_arr[*,*]= 0.0
           !sigh_arr[*,*]= 0.0
           !sigh_arr[0,*]= 0.0          !Sets inner L shells Sigma_h to zero -> no mode conversion
           !sigh_arr[Num_u1-1,*]= 0.0   !Sets inner L shells Sigma_h to zero -> no mode conversion
           
           
           !!Set Corner Conductance in thin sheet to ?
           !sigpatm_N(0) = 0.0
           !sighatm_N(0) = 0.0
           !sigpatm_S(0) = 0.0
           !sighatm_S(0) = 0.0

           !Calculate Va and Wave Speed

           V_arr = B0/sqrt(u0_u*rho_a)                  !Alfven Speed in m/s ()
           Coll_Mod = omega_p**2/(omega_p**2 + Nu_i**2)    !see Lysak '99 eq 7: Note ignore electron term since me << Mav
           Eperp = e0_u*(1+(c_u**2/V_arr**2)*Coll_mod)    !E perp in SI (including collision modification)
           eps_arr = Eperp!/Re
           Wave_sp = c_u/sqrt(Eperp/e0_u)               !The actually wave speed including collisions modification
           Va_arr = Wave_sp/Re
           !Va_arr = V_arr/Re                           !No Modification by collisions to wave speed
           va2_arr = va_arr**2

           !Determine time steps

           dtz = cour*minval(abs(d32*h3)/Va_arr)
           dtx = cour*minval(abs(2*del_u1*h1)/Va_arr)

           if (m .ne. 0) then
                   dty = cour*minval(h2/m/Va_arr)
                   dt=1.0/sqrt(3.0)/(sqrt(1.0/dtx**2+1.0/dty**2+1.0/dtz**2))
                   write(*,*)"dtx = ",dtx
                   write(*,*)"dty= ",dty
                   write(*,*)"dtz = ",dtz
           else
                   dt=1.0/sqrt(2.0)/sqrt(1.0/dtx**2+1.0/dtz**2)
                   write(*,*)"dtx = ",dtx
                   write(*,*)"dtz = ",dtz

           end if

           write(*,*)"min(Va_arr) = ",minval(Va_arr)
           write(*,*)"max(Va_arr) = ",maxval(Va_arr)

           nt = tmax/dt
           dtzmax = cour*maxval(abs(d32*h3)/Va_arr)
           dtxmax = cour*maxval(abs(2*del_u1*h1)/Va_arr)

           write(*,*) "Min spacing in u1 = ",minval(abs(2.0*del_u1*h1))*Re,"m"
           write(*,*) "Min spacing in u3 = ",minval(abs(d32*h3))*Re,"m"
           write(*,*) "Max Va = ",maxval(abs(Va_arr))*Re,"m/s"

           write(*,*) "******************************************"
           write(*,*) "Total Simulation Time,tmax = ",tmax
           write(*,*) "Time steps = ",dt
           write(*,*) "Number of time loops = ",nt

           courx = Va_arr*dt/abs(del_u1*h1)
           i_loc = maxloc(courx)
           cmax = maxval(courx)
           i1z = int(i_loc(0)/(num_u1))
           i1x = i_loc(0)-i1z*(num_u1)
           write(*,*)"Max Cour_len u1 = ",cmax," at i1x,i1z = ",i1x,i1z


           courz = Va_arr*dt/(abs(h3))
           i_loc = maxloc(courz)
           cmaz = maxval(courz)
           i3z = int(i_loc(0)/(num_u1))
           i3x = i_loc(0)-i3z*(num_u1)
           write(*,*)"Max Cour_len u3 = ",cmaz," at i3x,i3z = ",i3x,i3z

           write(*,*)"num_u1 x nzp1 x nt = ",num_u1*num_u3*nt
           write(*,*)"Max travel time between grid points u1,u3 directions = ",dtxmax,dtzmax
           write(*,*) "****************************************************"    


           !factors involving conductances in the numerics of thin sheet current sheet
           !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           
           !musigfac gives muSig in s/Re for 1 mho
           musigfac = u0_u*Re
           
           !For Northern Hemipshere
           cosalpha_N = -2.0*costh(:,0)/bsqrt(:,0) !set as postive for now, costh0 only goes 0< x < 90 deg
           sin2alpha_N = 1.0 - cosalpha_N*cosalpha_N
           musigp  = sigpatm_N*musigfac
           musigh  = sighatm_N*musigfac
           musig0  = sig0atm_N*musigfac
           
           musigzz = musig0 + (musigp-musig0)*sin2alpha_N
           sigtt   = musig0*musigp/musigzz
           sigtp   =-musig0*musigh*cosalpha_N/musigzz
           sigpt   = musig0*musigh*cosalpha_N/musigzz
           sigpp   = musigp + musigh*musigh*sin2alpha_N/musigzz
           
           sig11   = sigtt*hphiatm_N/hthatm_N
           sig12   = sigtp
           sig21   = sigpt
           sig22   = sigpp*hthatm_N/hphiatm_N
           sigden  = sig11*sig22-sig12*sig21

           e1b1atm_N =-sig12/sigden
           e1b2atm_N =-sig22/sigden
           e2b1atm_N = sig11/sigden
           e2b2atm_N = sig21/sigden

           !For southern Hemipshere
           cosalpha_S  = -2.0*costh(:,Num_u3-1)/bsqrt(:,Num_u3-1)
           sin2alpha_S  = 1.0 - cosalpha_S*cosalpha_S

           musigp = sigpatm_S*musigfac
           musigh = sighatm_S*musigfac
           musig0 = sig0atm_S*musigfac

           musigzz = musig0 + (musigp-musig0)*sin2alpha_S
           sigtt = musig0*musigp/musigzz
           sigtp =-musig0*musigh*cosalpha_S/musigzz
           sigpt = musig0*musigh*cosalpha_S/musigzz
           sigpp = musigp + musigh*musigh*sin2alpha_S/musigzz

           sig11 = sigtt*hphiatm_S/hthatm_S
           sig12 = sigtp
           sig21 = sigpt
           sig22 = sigpp*hthatm_S/hphiatm_S
           sigden = sig11*sig22-sig12*sig21

           e1b1atm_S =-sig12/sigden
           e1b2atm_S =-sig22/sigden
           e2b1atm_S = sig11/sigden
           e2b2atm_S = sig21/sigden

           write(*,*) 'Get Va Complete...'

end subroutine get_Va
