c------------------------------------------------------------------------
c  RK4 SOLVER for a system of ODEs                                      |
c  Course:  AE305                                                       |
c------------------------------------------------------------------------
      program sysRK4
      parameter (neq=4)
      real y(neq)
      character*40 fname

c..read the input data
      print*, 'Enter the step size and the end point:'
      read(*,*) dt,tmax
c..open the output file
      print*, ' Enter the output file name :'
      read(*,'(a)') fname
      if( fname .eq. ' ') fname = 'solution.dat'
      open(1,file=fname,form='formatted')

c..set the initial conditions
      time = 0.
      y(1) = 0.
      y(2) = 0.
      y(3) = 5.*3.14159/180.
      y(4) = 0.
      write(1,'(5E14.5)') time,(y(n),n=1,neq)

c..solution loop
      DO WHILE (time.lt.tmax)
        call SRK4(time,dt,y)
        time = time + dt
        write(1,'(5E14.5)') time,(y(n),n=1,neq)
      ENDDO

      close(1)
      pause
      stop
      end

c------------------------------------------------------------------------
       subroutine SRK4(time,dt,y)
       parameter (neq=4)
       real y(neq),ytmp(neq),k1(neq),k2(neq),k3(neq),k4(neq)

c..obtaining slopes
       dt2 = 0.5*dt
       call ODEs(time,y,k1)
       
       do n = 1,neq
         ytmp(n)  = y(n) + k1(n)*dt2
       enddo
       call ODEs(time+dt2,ytmp,k2)
       
       do n = 1,neq
          ytmp(n) = y(n) + k2(n)*dt2
       end do
       call ODEs(time+dt2,ytmp,k3)
       
       do n = 1,neq
          ytmp(n) = y(n) + k3(n)*dt
       enddo
       call ODEs(time+dt,ytmp,k4)


c..obtain the solution at t+dt and update y for the next step
       do n = 1,neq
          phi  = (k1(n) + 2.*(k2(n)+k3(n)) + k4(n))/6.
          y(n) = y(n) + phi*dt
       enddo

       return
       end

c------------------------------------------------------------------------
       subroutine ODEs(time,y,f)
       parameter (neq=4)
       real :: y(neq),f(neq), c_L
c..importing data
       real :: I = 50.
       real :: U=40.
       real :: c_a = 150.
       real :: a = 0.15
       real :: k_z = 20000.
       real :: S = 15.
       real :: c_z = 300.
       real :: M = 150.
       real :: k_a = 15000.
       real :: rho = 1.225

       
c..setting c_L value due to the angle of attack
       if (y(3)*180/3.14159 .lt. 12)  then
          c_L = 2*3.14159*y(3)
       else
           c_L = 2*3.14159*12*3.14159/180
        end if


c..define the ODE's & return the slopes in the "f" array
       f(1) = y(2)
       f(2) =(0.5*rho*U*S*c_L*sqrt(U**2+y(2)**2.)-c_z*y(2)-k_z*y(1))/M
       f(3) = y(4)
       f(4) =(0.5*rho*U*S*c_L*a*sqrt(U**2+y(2)**2.)-c_a*y(4)-k_a*y(3))/I

       return
       end
