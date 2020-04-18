	Module vars
  	integer, parameter :: imax=51, jmax=26
  	integer, parameter :: kmax=100000, kres=10, ksol=1000
  	integer :: iws,iwe,jws,jwe,iws1,iwe1,jws1,jwe1
  	real, parameter :: xl=10., yl=5., omega=1.2, rl2allow=-5.
  	real, dimension(imax,jmax) :: phi_k,phi_kp1,u,v
    real, dimension(kmax) :: rl2                              !..Defining Rl2 array
    real,dimension((imax-2)*(jmax-2),(imax-2)*(jmax-2)) :: a
    real,dimension((imax-2)*(jmax-2)) :: b
  	real :: dx,dy,beta2,dx2i,dy2i,uinf,vinf, x(imax), y(jmax)
  	logical, dimension(imax,jmax) :: body,bc_o_t,bc_o_r,bc_o_b,bc_o_l,bc_i_t,bc_i_r,bc_i_b,bc_i_l,c_1,c_2,c_3,c_4,a_1,a_2,a_3
    logical, dimension(imax,jmax) :: a_4,a_5,a_6,a_7,a_8
    character :: choice, activate
	End module

!------------------------------------------------------------------------------|
!..A POINT/LINE ITERATIVE SOLVER FOR ELLIPTIC PDEs                             |
!  Course:  AE305                                                              |
!------------------------------------------------------------------------------|
 	program ELLIPTIC
 	use vars
!..Read the input data, generate the grid data and initialize the solution
  	call INIT()
!..Start the iterative solution loop 
  	k = 0
 	 rl2(1) = 0.
  	DO WHILE (k .eq. 0 .or. (k .lt. kmax .and. rl2(k) .gt. rl2allow ))
     	k = k+1
!..Apply BCs 
    	 call BC()
!..Point iterative solutions
		if (choice.eq.'a') call JACOBI()   	
		if (choice.eq.'b') call GAUSS_SEIDEL()
		if (choice.eq.'c') call SOR()
!..or Direct/Line iterative solutions
		if (choice.eq.'d') call DIRECT()                       !...The code script is not working properly
!..Update phi_kp1 array and evaluate the L2 norm of the residual
     	rl2(k) = SQRT( SUM( (phi_kp1-phi_k)**2 ) )
     	phi_k = phi_kp1
    	 if(k .eq. 1) rl2_1=rl2(k)
     	 rl2(k) = LOG(rl2(k)/rl2_1)
    	 if( MOD(k,kres).eq.0 .or. k .eq. 1) print*, 'Residual @ k =',k,rl2(k)
         
!..Output intermediate solutions
     	if( MOD(k,ksol).eq.0 .and. k. ne. kmax) call QOUT(k)
  	ENDDO
  	print*, 'Residual @ k =',k,rl2(k)  
  	call QOUT(k)
    open(2,file='residual.dat',form='formatted')
	do k=1,kmax
    	write(2,*) k,rl2(k)     !..Writing Rl2 values to a dat file
    enddo
    close(2)

 	stop
	end

!------------------------------------------------------------------------
	subroutine INIT()
 	use vars

	print*, 'Please choose the method. Enter its letter.  '
    print*, 'a- Point Jacobi Iteration'
    print*, 'b- Gauss-Seidel Iteration'
    print*, 'c- Succesive Over Relaxation (SOR)'
!    print*, 'd- Direct Solution (BONUS)'         Since it is not working truely, it has been commanded out :(
    read(*,*) choice

    print*, 'Would you like to add a second rectangle? y or n'
    read (*,*) activate

  	write(*,'(a)',advance='no') 'Enter the angle of attack: '
  	read(*,*) aoa

  	dx = xl/(imax-1)
  	dy = yl/(jmax-1)
  	beta2 = (dx/dy)**2
  	dx2i = 1./(2*dx)
  	dy2i = 1./(2*dy)
  	aoa  = aoa/180.*3.14159
  	uinf = COS(aoa)
  	vinf = SIN(aoa) 

!..Set wall/object boundary locations and body grid points
	body = .false.
  	iws = 4/dx + 1
  	iwe = 6/dx + 1
  	jws = 2/dy + 1
  	jwe = 3/dy + 1
    body( iws:iwe, jws:jwe ) = .true.
    
    if(activate.eq.'y') then       !..Boundaries for second rectangle
   	 	iws1 = 7/dx + 1
  		iwe1 = 8/dx + 1
  		jws1 = 1.2/dy + 1
  		jwe1 = 3.8/dy + 1
        body( iws1:iwe1, jws1:jwe1 ) = .true.
    endif
  	
    

!..Grid generation 
  	do i=1,imax
     	x(i)= dx*(i-1) 
  	enddo
  	do j=1,jmax
     	y(j)= dy*(j-1) 
  	enddo

 	 phi_k   = 0.       !..Initial guess, a bad one
  	 phi_kp1 = phi_k

 	return 
	end

!..Implement, Point Jacobi, Gauss-Seidel and SOR methods  and solve for phi^k+1
	subroutine JACOBI()
 	use vars

  	do j = 2,jmax-1
  	do i = 2,imax-1
!..Check for the solution domain
	if (.not.body(i,j)) then
     phi_kp1(i,j) = (phi_k(i-1,j) + phi_k(i+1,j) + beta2*phi_k(i,j-1) + beta2*phi_k(i,j+1))/(2.0*(1.0+beta2))  
    endif       
  	enddo
  	enddo

 	return 
	end

    subroutine GAUSS_SEIDEL()
    use vars

    do j = 2,jmax-1
    do i = 2,imax-1
      if (.not.body(i,j)) then
        phi_kp1(i,j) = (phi_kp1(i-1,j) + phi_k(i+1,j) + beta2*phi_kp1(i,j-1) + beta2*phi_k(i,j+1))/(2.0*(1.0+beta2))
      endif
    enddo
    enddo

    return 
    end

    subroutine SOR()
    use vars
    real::gs

    do j = 2,jmax-1
    do i = 2,imax-1
      if (.not.body(i,j)) then
        gs = (phi_kp1(i-1,j) + phi_k(i+1,j) + beta2*phi_kp1(i,j-1) + beta2*phi_k(i,j+1))/(2.0*(1.0+beta2))
        phi_kp1(i,j) = phi_k(i,j) + omega * (gs-phi_k(i,j))
      endif
    enddo
    enddo

    return 
    end

    subroutine DIRECT()
    use vars
    integer :: n,m
	n = (imax-2)*(jmax-2)
	a = 0.
    b = 0.

    bc_o_t( 3:imax-2,jmax-1 ) = .true.
    bc_o_r( imax-1,3:jmax-2 ) = .true.
    bc_o_b( 3:imax-2, 2 ) = .true.
    bc_o_l( 2,3:jmax-2 ) = .true.
    bc_i_t(iws:iwe, jwe+1) = .true.
    bc_i_r(iwe+1, jws:jwe) = .true.
    bc_i_b(iws:iwe, jws-1) = .true.
    bc_i_l(iws-1, jws:jwe) = .true.
    c_1(2,2) = .true.
    c_2(2,jmax-1) = .true.
    c_3(imax-1,jmax-1) = .true.
    c_4(imax-1,2) = .true.
    a_1(iws-1,jws) = .true.
    a_2(iws,jws-1) = .true.
    a_3(iws-1,jwe) = .true.
    a_4(iws,jwe+1) = .true.
    a_5(iwe+1,jwe) = .true.
    a_6(iwe,jwe+1) = .true.
    a_7(iwe+1,jws) = .true.
    a_8(iwe,jws-1) = .true.

    do j=2,jmax-1
    do i=2,imax-1
        m = (j-2)*(imax-2)+(i-2)+1
        if(body(i,j)) then
          a(m,(j-2)*(imax-2)+(i-2)+1) = 1.0
        else if(c_1(i,j)) then
          a(m,(j-2)*(imax-2)+(i-2)+1) = 1.0 + beta2 - 2.0*(1.0+beta2)
          a(m,(j-1)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-2)*(imax-2)+(i-1)+1) = 1.0
        else if(c_2(i,j)) then
          a(m,(j-2)*(imax-2)+(i-2)+1) = 1.0 + beta2 - 2.0*(1.0+beta2)
          a(m,(j-3)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-2)*(imax-2)+(i-1)+1) = 1.0
        else if(c_3(i,j)) then
          a(m,(j-2)*(imax-2)+(i-2)+1) = 1.0 + beta2 - 2.0*(1.0+beta2)
          a(m,(j-3)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-2)*(imax-2)+(i-3)+1) = 1.0
        else if(c_4(i,j)) then
          a(m,(j-2)*(imax-2)+(i-2)+1) = 1.0 + beta2 - 2.0*(1.0+beta2)
          a(m,(j-1)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-2)*(imax-2)+(i-3)+1) = 1.0
        else if(a_1(i,j)) then
          a(m,(j-2)*(imax-2)+(i-2)+1) = 0.5 - 2.0*(1.0+beta2)
          a(m,(j-3)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-1)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-2)*(imax-2)+(i-3)+1) = 1.0
          a(m,(j-3)*(imax-2)+(i-1)+1) = 0.5
        else if(a_2(i,j)) then
          a(m,(j-2)*(imax-2)+(i-2)+1) = 0.5*beta2 - 2.0*(1.0+beta2)
          a(m,(j-3)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-2)*(imax-2)+(i-3)+1) = 1.0
          a(m,(j-2)*(imax-2)+(i-1)+1) = 1.0
          a(m,(j-1)*(imax-2)+(i-3)+1) = 0.5*beta2
        else if(a_3(i,j)) then
          a(m,(j-2)*(imax-2)+(i-2)+1) = 0.5 - 2.0*(1.0+beta2)
          a(m,(j-3)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-1)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-2)*(imax-2)+(i-3)+1) = 1.0
          a(m,(j-1)*(imax-2)+(i-1)+1) = 0.5
        else if(a_4(i,j)) then
          a(m,(j-2)*(imax-2)+(i-2)+1) = 0.5*beta2 - 2.0*(1.0+beta2)
          a(m,(j-1)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-2)*(imax-2)+(i-3)+1) = 1.0
          a(m,(j-2)*(imax-2)+(i-1)+1) = 1.0
          a(m,(j-3)*(imax-2)+(i-3)+1) = 0.5*beta2
        else if(a_5(i,j)) then
          a(m,(j-2)*(imax-2)+(i-2)+1) = 0.5 - 2.0*(1.0+beta2)
          a(m,(j-3)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-1)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-2)*(imax-2)+(i-1)+1) = 1.0
          a(m,(j-1)*(imax-2)+(i-3)+1) = 0.5
        else if(a_6(i,j)) then
          a(m,(j-2)*(imax-2)+(i-2)+1) = 0.5*beta2 - 2.0*(1.0+beta2)
          a(m,(j-1)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-2)*(imax-2)+(i-3)+1) = 1.0
          a(m,(j-2)*(imax-2)+(i-1)+1) = 1.0
          a(m,(j-3)*(imax-2)+(i-1)+1) = 0.5*beta2
        else if(a_7(i,j)) then
          a(m,(j-2)*(imax-2)+(i-2)+1) = 0.5 - 2.0*(1.0+beta2)
          a(m,(j-3)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-1)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-2)*(imax-2)+(i-1)+1) = 1.0
          a(m,(j-3)*(imax-2)+(i-3)+1) = 0.5
        else if(a_8(i,j)) then
          a(m,(j-2)*(imax-2)+(i-2)+1) = 0.5*beta2 - 2.0*(1.0+beta2)
          a(m,(j-3)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-2)*(imax-2)+(i-3)+1) = 1.0
          a(m,(j-2)*(imax-2)+(i-1)+1) = 1.0
          a(m,(j-1)*(imax-2)+(i-1)+1) = 0.5*beta2
        else if(bc_o_t(i,j)) then
          a(m,(j-2)*(imax-2)+(i-2)+1) = beta2 - 2.0*(1.0+beta2)
          a(m,(j-3)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-2)*(imax-2)+(i-3)+1) = 1.0
          a(m,(j-2)*(imax-2)+(i-1)+1) = 1.0
        else if(bc_o_r(i,j)) then
          a(m,(j-2)*(imax-2)+(i-2)+1) = 1.0 - 2.0*(1.0+beta2)
          a(m,(j-3)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-1)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-2)*(imax-2)+(i-3)+1) = 1.0
        else if(bc_o_b(i,j)) then
          a(m,(j-2)*(imax-2)+(i-2)+1) = beta2 - 2.0*(1.0+beta2)
          a(m,(j-1)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-2)*(imax-2)+(i-3)+1) = 1.0
          a(m,(j-2)*(imax-2)+(i-1)+1) = 1.0
        else if(bc_o_l(i,j)) then
          a(m,(j-2)*(imax-2)+(i-2)+1) = 1.0 - 2.0*(1.0+beta2)
          a(m,(j-3)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-1)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-2)*(imax-2)+(i-1)+1) = 1.0
        else if(bc_i_t(i,j)) then
          a(m,(j-2)*(imax-2)+(i-2)+1) = beta2 - 2.0*(1.0+beta2)
          a(m,(j-1)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-2)*(imax-2)+(i-3)+1) = 1.0
          a(m,(j-2)*(imax-2)+(i-1)+1) = 1.0
        else if(bc_i_r(i,j)) then
          a(m,(j-2)*(imax-2)+(i-2)+1) = 1.0 - 2.0*(1.0+beta2)
          a(m,(j-3)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-1)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-2)*(imax-2)+(i-1)+1) = 1.0
        else if(bc_i_b(i,j)) then
          a(m,(j-2)*(imax-2)+(i-2)+1) = beta2 - 2.0*(1.0+beta2)
          a(m,(j-3)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-2)*(imax-2)+(i-3)+1) = 1.0
          a(m,(j-2)*(imax-2)+(i-1)+1) = 1.0
        else if(bc_i_l(i,j)) then
          a(m,(j-2)*(imax-2)+(i-2)+1) = 1.0 - 2.0*(1.0+beta2)
          a(m,(j-3)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-1)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-2)*(imax-2)+(i-3)+1) = 1.0
        else 
          a(m,(j-2)*(imax-2)+(i-2)+1) = - 2.0*(1.0+beta2)
          a(m,(j-3)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-1)*(imax-2)+(i-2)+1) = beta2
          a(m,(j-2)*(imax-2)+(i-3)+1) = 1.0
          a(m,(j-2)*(imax-2)+(i-1)+1) = 1.0
        endif
        b((j-2)*(imax-2)+(i-2)+1) = 0.
      enddo
      enddo

        do i=2,imax-1
      		j=2
      		b((j-2)*(imax-2)+(i-2)+1) = b((j-2)*(imax-2)+(i-2)+1) + beta2*vinf
      		j = jmax-1
      		b((j-2)*(imax-2)+(i-2)+1) = b((j-2)*(imax-2)+(i-2)+1) - beta2*vinf
    	enddo


    	do j=2,jmax-1
        	i=2
			b((j-2)*(imax-2)+(i-2)+1) = b((j-2)*(imax-2)+(i-2)+1) + uinf
        	i = imax-1
			b((j-2)*(imax-2)+(i-2)+1) = b((j-2)*(imax-2)+(i-2)+1) - uinf
    	enddo


    call GAUSS(n,a,b)
    
	do j=2,jmax-1
    do i=2,imax-1
      m = (j-2)*(imax-2)+(i-2)+1
      phi_k(i,j) = b(m)
    enddo
    enddo

    do j = 2,jmax-1     
     phi_k  (1,j)    = phi_k(2,j)      - uinf*dx
     phi_k  (imax,j) = phi_k(imax-1,j) + uinf*dx
  	enddo

  	do i = 2,imax-1    
     	phi_k(i,1)      = phi_k(i,2)      - vinf*dy
     	phi_k(i,jmax)   = phi_k(i,jmax-1) + vinf*dy
  	enddo

    phi_k(1,1) = 0.5*(phi_k(1,2) + phi_k(2,1))
    phi_k(1,jmax) = 0.5*(phi_k(2,jmax) + phi_k(1,jmax-1))
    phi_k(imax,jmax) = 0.5*(phi_k(imax-1,jmax)+phi_k(imax,jmax-1))
    phi_k(imax,1) = 0.5*(phi_k(imax-1,2) + phi_k(imax,2))

    do i = iws,iwe
 	   phi_k(i,jws)  = phi_k(i,jws-1)
	   phi_k(i,jwe)  = phi_k(i,jwe+1)
  	enddo
    do j = jws,jwe
       phi_k(iws,j)  = phi_k(iws-1,j)
	   phi_k(iwe,j)  = phi_k(iwe+1,j)
    enddo
   

    phi_k(iws,jws) = 0.5*(phi_k(iws-1,jws) + phi_k(iws,jws-1))
    phi_k(iwe,jws) = 0.5*(phi_k(iwe+1,jws) + phi_k(iwe,jws-1))
    phi_k(iwe,jwe) = 0.5*(phi_k(iwe+1,jwe) + phi_k(iwe,jwe+1))
    phi_k(iws,jwe) = 0.5*(phi_k(iws-1,jwe) + phi_k(iws,jwe+1))


    call QOUT(k)

	stop
    return 
    end
    
    
	

!-------------------------------------------------------------------
	subroutine BC()
 	use vars
  	do j = 2,jmax-1      !..Set left/right farfield BC
     phi_k  (1,j)    = phi_k(2,j)      - uinf*dx
     phi_k  (imax,j) = phi_k(imax-1,j) + uinf*dx
     phi_kp1(1,j)    = phi_k(1,j)
     phi_kp1(imax,j) = phi_k(imax,j)
  	enddo

  	do i = 2,imax-1      !..Set bottom/top farfield BC
     	phi_k(i,1)      = phi_k(i,2)      - vinf*dy
     	phi_k(i,jmax)   = phi_k(i,jmax-1) + vinf*dy
     	phi_kp1(i,1)    = phi_k(i,1)
     	phi_kp1(i,jmax) = phi_k(i,jmax)
  	enddo

!..Set wall BCs around the object
  	do i = iws,iwe
 	   phi_k(i,jws)  = phi_k(i,jws-1)
	   phi_k(i,jwe)  = phi_k(i,jwe+1)
       phi_kp1(i,jws) = phi_k(i,jws)
       phi_kp1(i,jwe) = phi_k(i,jwe)
  	enddo
    do j = jws,jwe
       phi_k(iws,j)  = phi_k(iws-1,j)
	   phi_k(iwe,j)  = phi_k(iwe+1,j)
       phi_kp1(iws,j) = phi_k(iws,j)
       phi_kp1(iwe,j) = phi_k(iwe,j)
    enddo

!	BCs for corners of the object
    phi_k(iws,jws) = 0.5*(phi_k(iws-1,jws) + phi_k(iws,jws-1))
    phi_k(iwe,jws) = 0.5*(phi_k(iwe+1,jws) + phi_k(iwe,jws-1))
    phi_k(iwe,jwe) = 0.5*(phi_k(iwe+1,jwe) + phi_k(iwe,jwe+1))
    phi_k(iws,jwe) = 0.5*(phi_k(iws-1,jwe) + phi_k(iws,jwe+1))
    phi_kp1(iws,jws) = phi_k(iws,jws)
    phi_kp1(iwe,jws) = phi_k(iwe,jws)
    phi_kp1(iwe,jwe) = phi_k(iwe,jwe)
    phi_kp1(iws,jwe) = phi_k(iws,jwe)

    if(activate.eq.'y') then         !...BCs for second rectangle
    do i = iws1,iwe1
 	   phi_k(i,jws1)  = phi_k(i,jws1-1)
	   phi_k(i,jwe1)  = phi_k(i,jwe1+1)
       phi_kp1(i,jws1) = phi_k(i,jws1)
       phi_kp1(i,jwe1) = phi_k(i,jwe1)
  	enddo
    do j = jws1,jwe1
       phi_k(iws1,j)  = phi_k(iws1-1,j)
	   phi_k(iwe1,j)  = phi_k(iwe1+1,j)
       phi_kp1(iws1,j) = phi_k(iws1,j)
       phi_kp1(iwe1,j) = phi_k(iwe1,j)
    enddo
    !	BCs for coeners of the object
    phi_k(iws1,jws1) = 0.5*(phi_k(iws1-1,jws1) + phi_k(iws1,jws1-1))
    phi_k(iwe1,jws1) = 0.5*(phi_k(iwe1+1,jws1) + phi_k(iwe1,jws1-1))
    phi_k(iwe1,jwe1) = 0.5*(phi_k(iwe1+1,jwe1) + phi_k(iwe1,jwe1+1))
    phi_k(iws1,jwe1) = 0.5*(phi_k(iws1-1,jwe1) + phi_k(iws1,jwe1+1))
    phi_kp1(iws1,jws1) = phi_k(iws1,jws1)
    phi_kp1(iwe1,jws1) = phi_k(iwe1,jws1)
    phi_kp1(iwe1,jwe1) = phi_k(iwe1,jwe1)
    phi_kp1(iws1,jwe1) = phi_k(iws1,jwe1)
    endif
    
 	return 
	end

!-------------------------------------------------------------------
	subroutine QOUT(k)          !..Output the solution in tecplot format
 	use vars
   		character fname*32,string*9,ext*6
   		write(string,'(f9.6)') float(k)/1000000
   		read(string,'(3x,a6)') ext
   		fname = 'q-'//ext//'.tec'
   		open(1,file=fname,form='formatted')
   		write(1,*) ' variables="x","y","phi","u","v", '
   		write(1,*) ' zone i=',imax, ', j=',jmax
   		call VELOCITY()
   		do j=1,jmax
   		do i=1,imax
      		write(1,*) x(i),y(j),phi_k(i,j),u(i,j),v(i,j)
   		enddo
   		enddo
   		close(1)
 		return
	end

!-------------------------------------------------------------------
		subroutine VELOCITY()
 		use vars
!..Compute u and v velocity components & enforce the FS/wall velocities
 		u = uinf
  		v = vinf
  		do j = 2,jmax-1
  		do i = 2,imax-1
     		if( body(i,j) ) then
       			u(i,j) = 0.
       			v(i,j) = 0.
     		else
       			u(i,j) = (phi_k(i+1,j)- phi_k(i-1,j))*dx2i
       			v(i,j) = (phi_k(i,j+1)- phi_k(i,j-1))*dy2i
     		endif
  		enddo
  		enddo
 		return 
		end

!-------------------------------------------------------------------
	subroutine THOMAS(is,ie, a,b,c,f)
  	use vars, only : imax
  	real, dimension(imax) ::  a,b,c,f,x
!
!  Solution of a tridiagonal system of n equations of the form
!  A(i)*x(i-1) + B(i)*x(i) + C(i)*x(i+1) = F(i)  for k=is,ie
!  the solution X(i) is stored in F(i)
!  A(is-1) and C(ie+1) are not used.
!  Solution is returned in array F
!
   	x(is)=c(is)/b(is)
   	f(is)=f(is)/b(is)
   	isp1 = is+1
   	do i=isp1,ie
      	z   =1./(b(i)-a(i)*x(i-1))
      	x(i)=z*c(i)
      	f(i)=(f(i)-a(i)*f(i-1))*z
   	enddo
   	iepis=ie+is
   	do ii=isp1,ie
      	i=iepis-ii
      	f(i)=f(i)-x(i)*f(i+1)
   	enddo
     
 	return
	end

!-------------------------------------------------------------------
	subroutine GAUSS(n,a,b)
   	real a(n,n),b(n)
!..Convert to upper triangular form
   	do k = 1,N-1
   	if (ABS(a(K,K)).gt.1.E-6) THEN
      do i = k+1, n
      x = a(i,k)/a(k,k)
         do j = k+1, n
            a(i,j) = a(i,j) -a(k,j)*x
         enddo
      b(i) = b(i) - b(k)*x
      enddo
   	else
      write (6,*) 'zero pivot found in line:', k
      stop
   	endif
   	enddo
!..Back substitution
   	do i = n,1,-1
     sum = b(i)
     if (i.lt.n) then
       do j= i+1,n
         sum = sum - a(i,j)*b(j)
       enddo
     endif
     b(i) = sum/a(i,i)
   	enddo
  	return
	end
