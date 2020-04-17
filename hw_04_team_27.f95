		Module vars
 		integer, parameter :: imax = 301
 		integer :: ntmax,ntout
 		real :: dnum, sigma
 		real, dimension(imax) :: x,qn,qnp
		End module

        

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
		program CONVECTION_DIFFUSION
   		Use vars

		character :: choice                              !...Select the method and problem
        print*, 'Choose the equation and the method: '
        print*, 'a- Conduction-FTCS '
        print*, 'b- Conduction-FTBS'
        print*, 'c- Linear Convection-FTCS'
        print*, 'd- Linear Convection-FTBS'
        print*, 'e- Conduction-BTCS (Implicit Solution)'
        print*, 'f- Linear Convection-BTCS (Implicit Solution)'
        print*, 'g- Conduction-FTCS with insulated wall at x=0'
        print*, 'h- Diffusion Convection-FTCS (Explicit Solution)'
        print*, 'i- Diffusion Convection-BTCS (Implicit Solution)'
        read(*,*) choice

   		call INIT(choice)             !..Read the input data, and initialize the problem
   		DO nt = 2,ntmax               !..Start the solution loop 
	!..Solve for q^n+1
    	 if(choice.eq.'a' .or. choice.eq.'b' .or. choice.eq.'c' .or. choice.eq.'d' .or. choice.eq.'g' .or. choice.eq.'h')then  !...Explicit Parts
    	 call EXPLICIT(choice)
         else if (choice.eq.'e' .or. choice.eq.'f' .or. choice.eq.'i')then                                                     !...Implicit parts
		 call fIMPLICIT(choice)
         end if
    	 qn(:) = qnp(:)       !..Update the solution
  !..Output intermediate/final solutions
    	 if( MOD(nt,ntout) .eq. 0 .or. nt .eq. ntmax ) call QOUT(nt)
   		ENDDO                        
   		stop 'DONE'
         End     

!------------------------------------------------------------------------
		subroutine INIT(choice)
   		use vars
        character choice
        if (choice.eq.'a' .or. choice.eq.'b' .or. choice.eq.'e' .or. choice.eq.'g')then   !...Read dnum
   			write(*,'(a)',advance='no') 'Enter dnum,ntmax & ntout : '
  			read(*,*) dnum,ntmax,ntout
        else if (choice.eq.'c' .or. choice.eq.'d' .or. choice.eq.'f') then                !...Read sigma
          	write(*,'(a)',advance='no') 'Enter sigma,ntmax & ntout : '
  			read(*,*) sigma,ntmax,ntout
        else if (choice.eq.'h' .or. choice.eq.'i') then                                   !...Read sigma and dnum
          	write(*,'(a)',advance='no') 'Enter dnum,sigma,ntmax & ntout : '
  			read(*,*) dnum,sigma,ntmax,ntout
        endif
        
        if (choice.eq.'a' .or. choice.eq.'b' .or. choice.eq.'e' .or. choice.eq.'g')then   !...x values for Part a
   			x(1)    = 0.
   			x(imax) = 10.
   			dx      = (x(imax)-x(1))/(imax-1)                                             !...x values for Part b and c
        else if (choice.eq.'c' .or. choice.eq.'d' .or. choice.eq.'f' .or. choice.eq.'h' .or. choice.eq.'i') then
          	x(1)    = -10.
   			x(imax) = 20.
   			dx      = (x(imax)-x(1))/(imax-1)
        endif
!..discretize 1D space
   		do i = 2,imax-1
     	 x(i)  = x(i-1)+dx
  		 enddo
!..Apply IC/BC
		if (choice.eq.'a' .or. choice.eq.'b' .or. choice.eq.'e')then              !...BCs for conduction
  			qn       = 20.
   			qn(1)    = 10.
   			qn(imax) = 100.
   			qnp(1)    = 10.
   			qnp(imax) = 100.
   		call QOUT(1)
        else if (choice.eq.'c' .or. choice.eq.'d' .or. choice.eq.'f' .or. choice.eq.'h' .or. choice.eq.'i') then   !...BCs for convection
          	i = 1
          	do while (x(i).lt.(-3.1416))
        	qn(i) = -2.0
            i = i+1
            enddo
            do while (x(i).le.3.1416)
        	qn(i) = 2*COS(x(i))
            i = i+1
            enddo
            do while (x(i).le.20)
        	qn(i) = 2.0
            i = i+1
            enddo
            qnp(imax) = 2.0
   			qn(1)    = 0.
   			qnp(1)   = 0.
            if (choice.eq.'h' .or. choice.eq.'i')then    !...Additional BC for Part c
              qn(imax)  =0.
              qnp(imax) =0.
            endif
   		call QOUT(1)
        else if (choice.eq.'g') then                     !...BCs for insulated case in part a
			qn       = 20.
            qn(imax) = 100.
            qnp(imax) = 100.
        endif
        
   		return 
		End

!------------------------------------------------------------------------
		subroutine EXPLICIT(choice)
  		Use vars
        character :: choice
        
        if (choice.eq.'a')then              !...Conduction FTCS
          do i=2,imax-1
    		qnp(i) = dnum*qn(i+1) + (1-2*dnum)*qn(i) + dnum*qn(i-1)
          enddo
        else if (choice.eq.'b')then              !...Conduction FTBS
          qnp(2) = dnum*qn(3) + (1-2*dnum)*qn(2) + dnum*qn(1)
          do i=3,imax-1
    		qnp(i) = (1+dnum)*qn(i) - 2*dnum*qn(i-1) + dnum*qn(i-2)
          enddo
         else if (choice.eq.'c')then              !...Convection FTCS
          do i=2,imax-1
    		qnp(i) = 0.5*sigma*qn(i-1) + qn(i) - 0.5*sigma*qn(i+1)
          enddo
          qnp(imax) = qn(imax)
         else if (choice.eq.'d')then              !...Convection FTBS
          do i=2,imax-1
    		qnp(i) = (1-sigma)*qn(i) + sigma*qn(i-1)
          enddo
          qnp(imax) = qn(imax)
        else if (choice.eq.'g')then              !...Insulated case
          qnp(1) = 2*dnum*qn(2) + (1-2*dnum)*qn(1)
          do i=2,imax-1
    		qnp(i) = dnum*qn(i+1) + (1-2*dnum)*qn(i) + dnum*qn(i-1)
          enddo
        else if (choice.eq.'h')then              !...Diffusion Convection FTCS
          do i=2,imax-1
    		qnp(i) = (dnum+sigma*0.5)*qn(i-1) + (1-2*dnum)*qn(i) + (dnum-sigma*0.5)*qn(i+1)
          enddo
        endif
  		return 
		End

!------------------------------------------------------------------------
		subroutine fIMPLICIT(choice)
   		Use vars
        character :: choice
   		real, dimension(imax) :: cl,cm,cu,rhs

        
		if (choice.eq.'e')then                      !...Conduction BTCS
!..Fill the 3 diagonal and rhs arrays
   		do i = 2,imax-1
			cl(i)  = -dnum                          !..lower diagonal
			cm(i)  = 1+2*dnum                       !..main diagonal
			cu(i)  = -dnum                          !..upper diagonal
			rhs(i) = qn(i)
   		enddo
!..apply BC
		rhs(2)      = rhs(2) + dnum*qn(1)
		rhs(imax-1) = rhs(imax-1) + dnum*qn(imax)
   		call THOMAS(2,imax-1, cl,cm,cu,rhs)   !..Solve the system of equations
   		qnp(2:imax-1) = rhs(2:imax-1)         !..Extract the solution

        else if (choice.eq.'f')then                     !...Convection BTCS
!..Fill the 3 diagonal and rhs arrays
   		do i = 2,imax-1
			cl(i)  = -0.5*sigma                         !..lower diagonal
			cm(i)  = 1                                  !..main diagonal
			cu(i)  = sigma*0.5                          !..upper diagonal
			rhs(i) = qn(i)
   		enddo
!..apply BC
		rhs(2)      = rhs(2) + sigma*0.5*qn(1)
		rhs(imax-1) = rhs(imax-1) -sigma*0.5*qn(imax)
   		call THOMAS(2,imax-1, cl,cm,cu,rhs)   !..Solve the system of equations
   		qnp(2:imax-1) = rhs(2:imax-1)         !..Extract the solution

        else if (choice.eq.'i')then                   !...Diffusion Convection BTCS
!..Fill the 3 diagonal and rhs arrays
   		do i = 2,imax-1
			cl(i)  = -dnum-sigma*0.5                          !..lower diagonal
			cm(i)  = 1+2*dnum                                 !..main diagonal
			cu(i)  = sigma*0.5-dnum                           !..upper diagonal
			rhs(i) = qn(i)
   		enddo
!..apply BC
		rhs(2)      = rhs(2) + (dnum+sigma*0.5)*qn(1)
		rhs(imax-1) = rhs(imax-1) + (dnum-sigma*0.5)*qn(imax)
   		call THOMAS(2,imax-1, cl,cm,cu,rhs)   !..Solve the system of equations
   		qnp(2:imax-1) = rhs(2:imax-1)         !..Extract the solution

        endif
 		return 
		End

!-------------------------------------------------------------------
		subroutine QOUT(nt)
 		 Use vars
  		character fname*32,string*7,ext*5
  		write(string,'(f7.5)') float(nt)/100000
  		read(string,'(2x,a5)') ext
  		fname = 'q-'//ext//'.dat' 
  		open(1,file=fname,form='formatted')
  		do i=1,imax
   		  write(1,'(2e14.6)') x(i),qn(i)
 		 enddo
  		close(1)
		 return 
		End

!-------------------------------------------------------------------
		subroutine THOMAS(il,iu, A,B,C,F)
!............................................................
! Solution of a tridiagonal system of equations of the form
!  A(i)*x(i-1) + B(i)*x(i) + C(i)*x(i+1) = F(i)  for k=il,iu
!  the solution X(i) is stored in F(i)
!  A(il-1) and C(iu+1) are not used.
!  A,B,C,F are arrays to be filled by the caller program
!............................................................
 		 Use vars
 		 real, dimension(imax) :: a,b,c,f,tmp
 		 tmp(il)=c(il)/b(il)
 		 f(il)=f(il)/b(il)
 		 ilp1 = il+1
 		 do i=ilp1,iu
    		 z=1./(b(i)-a(i)*tmp(i-1))
     		tmp(i)=c(i)*z
     		f(i)=(f(i)-a(i)*f(i-1))*z
 		 enddo
  		iupil=iu+il
  		do ii=ilp1,iu
     		i=iupil-ii
     		f(i)=f(i)-tmp(i)*f(i+1)
  		enddo
  		return
		End
!------------------------------------END----------------------------