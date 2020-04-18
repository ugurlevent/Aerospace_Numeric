c-----------------------------------------------------------------------
c.........AE305 HOMEWORK 1 SOLUTION CODE IN FORTRAN - GROUP 27..........
c-----------------------------------------------------------------------

c.___________________________INTRODUCTION_______________________________

        program main
        implicit none

c...Defining the variables and their types..............................
        real :: error_e, error_mid, error_rk2, error_all
        real :: h_e, h_mid, h_rk2, act
        common  /error/ error_e, error_mid, error_rk2, error_all
        common  /stepsize/ h_e, h_mid, h_rk2, act
        integer :: w = 1, u=1, q=1
        integer :: choice1, choice2
        character :: choice3
        real :: altitude, pressure, rho, V_ground, V_air, V_r
        real :: T_0 = 288.15, g = 9.81, b=0.0065
        real :: pressure_0 = 101325, R=8.314 , M = 0.02896

c...Choosing altitude...................................................
        print*, 'Please choose the location or enter altitude: '
        print*, '(Enter the chosen number)'
        do while (u .eq. 1)
           print*, '1- Sea Level'
           print*, '2- Ankara Esenboga Airport (Altitude=952 km)'
           print*, '3- Istanbul Ataturk Airport(Altitude=50 km)'
           print*, '4- El Alto International Airport(Altitude=4061.5km)'
           print*, '5- Enter the altitude manually'
           read(*,*) choice1
           if (choice1 .eq. 1) then
              altitude = 0
              u = 0
           else if (choice1 .eq. 2) then
                altitude = 952
                u = 0
           else if (choice1 .eq. 3) then
                altitude = 50
                u = 0
           else if (choice1 .eq. 4) then
               altitude = 4061.5
               u = 0
           else if (choice1 .eq. 5) then
                print*, 'Enter the altitude: '
                read (*,*) altitude
                u = 0
           else
               print*, 'Please enter 1, 2, 3, 4 or 5'
           end if
        end do


c...Calculating pressure and density of the air.........................
        pressure = pressure_0*((T_0-b*altitude)/T_0)**(M*g/(R*b))
        rho = M*pressure/R/(T_0-b*altitude)


c...Setting atmospheric wind and calculating relative velocity..........
        print*, 'Enter the atmospheric wind'
        print*, '--->+ PLANE -<---'
        print*, '(Use positive for tail winds)'
        print*, '(Use negative for head winds)'
        read(*,*) V_air
        
        
        do while (q .eq. 1)
        print*, 'Activate adaptive step size mode? (y or n)'
        read(*,*) choice3
        if (choice3 .eq. 'y') then
           act = 1
           q = 0
           print *, 'Adaptive step size activated! Enter allowed error:'
           read(*,*) error_all
        else if (choice3 .eq. 'n') then
             act = 0
             q = 0
        else
            print*, 'Please enter y or n'
        end if
        end do


c....Choosing the method................................................
        print*, 'Please choose a method (Enter its number)'
        do while (w .eq. 1)
           print*, '1- Euler Method'
           print*, '2- Midpoint Method'
           print*, '3- Second Order Runge-Kutta with p=0.8'
           print*, '4- Decreasing Euler stepsize for better accuracy'
           read(*,*) choice2
           if (choice2 .eq. 1) then
            call euler(pressure, rho, V_air)
            w = 0
            h_e = 0
             else if (choice2 .eq. 2) then
                  call midpoint(pressure, rho, V_air)
                  h_mid = 0
                  w = 0
             else if (choice2 .eq. 3)   then
                  call rk2(pressure, rho, V_air)
                  h_rk2 = 0
                  w = 0
             else if (choice2 .eq. 4) then
                  call decrease_h_euler(pressure, rho, V_air)
                  w=0
             else
                  print*, 'Please enter 1, 2, 3 or 4'
             end if
        end do
        

        pause
        stop
        end
c_______________________________________________________________________






c_______________________CODES OF THE METHODS____________________________


c-------------------------MIDPOINT METHOD-------------------------------

        subroutine midpoint(pressure, rho, V_air)

c....Defining variables.................................................
        character*40 fname
        real :: error_e, error_mid, error_rk2, error_all
        real :: h_e, h_mid, h_rk2, act
        common  /error/ error_e, error_mid, error_rk2, error_all
        common  /stepsize/ h_e, h_mid, h_rk2, act
        real :: rho, k1, k2, V, V_new, t, pressure, x
        real :: m = 15000.0
        real :: g = 9.81
        real :: p = 0.5

c....Set stepsize.......................................................
        if (h_mid .eq. 0) then
           print*, '  Enter stepsize of Midpoint Method: '
           read(*,*) h_mid
        end if

c....Openning outout file...............................................
        print*,'  Enter the output file name for Midpoint Method:'
        read(*,'(a)') fname
        if( fname .eq. ' ') fname = 'velocity.dat'
        open(1,file=fname,form='formatted')

c...Setting initial conditions..........................................
        t = 0.0
        x = 0.0
        V=0.0
        V_r = V - V_air
        write(1,'(2f12.5)') t, V

c...Solution loop.......................................................
        do while(xL(V_r,rho) .le. m*g)
           V_r = V - V_air
           k1 = ODE(t, V_r, rho, pressure)
10         k2 = ODE(t+p*h_mid, V_r+p*h_mid*k1, rho, pressure)
           V_new = V + k2*h_mid
           error_mid = abs(V_new-V)*100/V_new
           if (error_mid .gt. error_all .and. act. eq. 1) then
              h_mid = h_mid*abs(error_all/error_mid)**0.2
              go to 10
           end if
           t = t + h_mid
           x = x + (V_new+V)*0.5*h_mid
           write(1,'(4f12.5)') t, V_new, x, error_mid
           V = V_new
        end do

c...Printing results....................................................
        print*, 'Lift off distance: ', x
        print*, 'Lift off time: ', t

        close(1)
        return
        end
c-----------------------------------------------------------------------




c---------------------------EULERS METHOD-------------------------------

        subroutine euler(pressure, rho, V_air)

c....Defining variables.................................................
        character*40 fname
        real :: error_e, error_mid, error_rk2, error_all
        real :: h_e, h_mid, h_rk2, act
        common  /error/ error_e, error_mid, error_rk2, error_all
        common  /stepsize/ h_e, h_mid, h_rk2, act
        real :: rho, k, V, V_new, t, pressure, x
        real :: m = 15000.0
        real :: g = 9.81

c....Set stepsize.......................................................
        if (h_e .eq. 0) then
           print*, 'Enter stepsize of Euler Method: '
           read(*,*) h_e
        end if

c....Openning outout file...............................................
        print*,'  Enter the output file name for Euler Method:'
        read(*,'(a)') fname
        if( fname .eq. ' ') fname = 'velocity.dat'
        open(1,file=fname,form='formatted')

c...Setting initial conditions..........................................
        t = 0.0
        x = 0.0
        V=0.0
        V_r = V - V_air
        write(1,'(2f12.5)') t, V

c...Solution loop.......................................................
        do while(xL(V_r,rho) .le. m*g)
           V_r = V - V_air
20         k = ODE(t, V_r, rho, pressure)
           V_new = V_r + k*h_e
           error_e = abs(V_new-V)*100/V_new
           if (error_e .gt. error_all .and. act. eq. 1) then  !Adaptive
              h_e = h_e*abs(error_all/error_e)**0.2
              go to 20
           end if
           t = t + h_e
           x = x + (V_new+V)*0.5*h_e
           write(1,'(4f12.5)') t, V_new, x, error_e
           V = V_new
        end do

c...Printing results....................................................
        print*, 'Lift off distance: ', x
        print*, 'Lift off time: ', t

        close(1)
        return
        end
c-----------------------------------------------------------------------



c-------------------------RK2 METHOD (with p=0.8)-----------------------

        subroutine rk2(pressure, rho, V_air)

c....Defining variables.................................................
        character*40 fname
        real :: error_e, error_mid, error_rk2, error_all
        real :: h_e, h_mid, h_rk2, act
        common  /error/ error_e, error_mid, error_rk2, error_all
        common  /stepsize/ h_e, h_mid, h_rk2, act
        real :: rho, k, k1, k2, V, V_new, t, pressure, x
        real :: m = 15000.0
        real :: g = 9.81
        real :: p = 0.8

c....Set stepsize.......................................................
        if (h_rk2 .eq. 0) then
           print*, 'Enter stepsize of RK-2: '
           read(*,*) h_rk2
        end if

c...Openning outout file................................................
        print*,'  Enter the output file name for RK2 Method:'
        read(*,'(a)') fname
        if( fname .eq. ' ') fname = 'velocity.dat'
        open(1,file=fname,form='formatted')

c...Setting initial conditions..........................................
        t = 0.0
        x = 0.0
        V=0.0
        V_r = V - V_air
        write(1,'(2f12.5)') t, V

c...Solution loop.......................................................
        do while(xL(V_r,rho) .le. m*g)
           V_r = V - V_air
           k1 = ODE(t, V_r, rho, pressure)
30         k2 = ODE(t+p*h_rk2, V_r+p*h_rk2*k1, rho, pressure)
           k = 0.375*k1 + 0.625*k2
           V_new = V + k*h_rk2
           error_rk2 = abs(V_new-V)*100/V_new
           if (error_rk2 .gt. error_all .and. act. eq. 1) then  !Adaptive
              h_rk2 = h_rk2*abs(error_all/error_rk2)**0.2
              go to 30
           end if
           t = t + h_rk2
           x = x + (V_new+V)*0.5*h_rk2
           write(1,'(4f12.5)') t, V_new, x, error_rk2
           V = V_new
        end do

c...Printing results....................................................
        print*, 'Lift off distance: ', x
        print*, 'Lift off time: ', t

        close(1)
        return
        end
c-----------------------------------------------------------------------

c-------------------DECREASING STEPSIZE FOR EULER-----------------------

        subroutine decrease_h_euler(pressure, rho, V_air)

c....Defining variables.................................................
                real :: error_e, error_mid, error_rk2, error_all
                real :: h_e, h_mid, h_rk2, act
                common  /error/ error_e, error_mid, error_rk2, error_all
                common  /stepsize/ h_e, h_mid, h_rk2, act

c....Setting stepsize...................................................
                print*, 'Enter the stepsize of the midpoint method: '
                read (*,*) h_mid
                h_e = h_mid

c....Checking errors....................................................
                call midpoint(pressure, rho, V_air)
                call euler(pressure, rho, V_air)
                print*, 'Error of midpoint method: ', error_mid, '%'
                print*, 'Error of euler method: ', error_e, '%'
                print*, 'Stepsize of euler method: ', h_e

c....Decreasing stepsize(By halfing)....................................
                do while (error_e .gt. error_mid)
                   h_e = h_e/2.0
                   call euler(pressure, rho, V_air)
                end do

c....Printing the result................................................
                print*, 'Error of midpoint method: ', error_mid, '%'
                print*, 'Error of euler method: ', error_e, '%'
                print*, 'Stepsize of euler method: ', h_e

        pause
        return
        end
c-----------------------------------------------------------------------






c_________________________FUNCTIONS_____________________________________


c....Ordinary Differential Equation.....................................
        real function ODE(t, V, rho, p)
             real :: m = 15000.0
             ODE = (T_max(V, rho, p) - D(V, rho) - F_fric(V, rho))/m
        return
        end

        real function D(V, rho)
             real :: c_D = 0.14
             real :: S = 50.0
             D = c_D * 0.5 * rho * S * V**2
        return
        end


c.....Lift..............................................................
        real function xL(V, rho)
             real :: c_L
             real :: S = 50.0
             if (V .le. 0) then
                c_L = 0
             else
                c_L = 1.5
             end if
             xL = c_L * 0.5 * rho * S * V**2
        return
        end


c...Friction Force......................................................
        real function F_fric(V, rho)
             real :: k_f = 0.02
             real :: m = 15000.0
             real :: g = 9.81
             F_fric = k_f * (m*g - xL(V, rho))
        return
        end


c...Maximum Thrust......................................................
        real function T_max(V, rho, p)
             real :: T_max_sl = 40000.0
             real :: rho_sl = 1.225
             real :: gama = 1.4
             real :: a
             a = sqrt(gama*p/rho)
             T_max = T_max_sl*(1+0.2*V/a)*((rho/rho_sl)**0.8)
        end


