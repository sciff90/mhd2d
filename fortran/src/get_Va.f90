subroutine get_Va()
        
        use global
        implicit none
        
        !arrays
        real,dimension(0:num_u1-1,0:num_u3-1) :: rho,Va,No_E,No_F1,No_F2,No_Mag,LVal_arr
        real,dimension(0:num_u1-1,num_u3-1) :: Lat,b0,nH_eq,rho_mag,No_tot,Rho_a      
        real,dimension(0:num_u1-1,0:num_u3-1) :: omega_p,omega_e
        real,dimension(0:num_u1-1):: xd,tanhm,tterm,va_sp,b_eq,nH

        !Ionospheric Data File Arrays
        real,dimension(0:num_u1-1,0:num_u3-1) :: sigp_arr,sigh_arr,sig0_arr,eta_arr
        real,dimension(0:num_u1-1,0:num_u3-1) :: No_N,Mav,Temp
        complex,dimension(0:num_u1-1) :: sigpatm_N,sighatm_N,sig0atm_N
        complex,dimension(0:num_u1-1) :: sigpatm_S,sighatm_S,sig0atm_S


        !Data File Variables
        integer,parameter :: pt = 29  
        real,dimension(0:pt+5,0:3):: neutral
        real::SclHt

        !Parameters
         real,parameter:: Eregion = 100e3
         real,parameter:: Mol_H = 1.0e-3
         real,parameter:: Av_n = 6.023e23
         real,parameter:: Mp = 1.67e-27
         real,parameter:: Me = 9.109e-31
         real,parameter:: q_i = 1.602e-19
         real,parameter:: q_e = -1.602e-19
         real,parameter:: k_0 = 8.0e15

         !Counters
         integer::ii,jj

         !Variables
         real:: Plaw,Ne_E,Alt_E,Esht_t,Esht_b,Ne_F1,Alt_F1,F1sht_t,F1sht_b
         real :: Ne_F2,Alt_F2,F2sht_t,F2sht_b,X
         
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
           write(*,*) 'ii = ',ii

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
           write(*,*) 'output = ' ,neutral(:,0)
           write(*,*) 'output = ' ,neutral(:,1)                     
           write(*,*) 'output = ' ,neutral(:,2)
           write(*,*) 'output = ' ,neutral(:,3)
end subroutine get_Va
