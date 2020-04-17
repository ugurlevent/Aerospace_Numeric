		Module vars
        integer :: ncell, nnode
        integer, dimension(:,:), allocatable :: node, neigh
        real, dimension(:,:), allocatable :: xy
        real, dimension(:), allocatable :: Tcell, dTdx, dTdy,area,Tnode
        real :: Tbc(10)
        End module

        !------------------------------------------------------------------------------|
        !..A 2-D FINITE VOLUME SOLVER FOR THE HEAT TRANSFER EQUATION                   |
        !                      ON UNSTRUCTURED GRIDS                                   |
        !  Course: AE305                                                               |
        !------------------------------------------------------------------------------|
         program FiniteVolume
         Use vars
         integer :: ntmax=10000, ntio=100
         real    :: dt=1.E-5, delTallow = 1.E-3, heat_in=0., heat_out=0.
         character :: choice
         print*, 'a- Without Convection'
         print*, 'b- With Convection'
         print*, 'Please enter a or b'
         read(*,*) choice

        !..Read the input data and initialize the solution
       call INIT()
       !..Start the solution loop
       DO nt=1,ntmax
       if(choice.eq.'a') then !..Evaluate temperature gradients for cells
       call GRADIENTa
       elseif (choice .eq. 'b') then
       call GRADIENTb
       endif
       delTmax = 0.
       do n = 1,ncell     !..Sweep all the cells and solve for T^n+1
       if(choice.eq.'a') then
       delT     = - dt/area(n) * TFLUXa(n,qin,qout)
       elseif (choice .eq. 'b') then
       delT     = - dt/area(n) * TFLUXb(n,qin,qout)
       endif
       heat_in = heat_in + qin*dt                     !..summing heat values
       heat_out= heat_out +qout*dt
       Tcell(n) = Tcell(n) + delT
       delTmax  = MAX(ABS(delT),delTmax)
       
       enddo
       time = time + dt
       print*, ' nt, Time, DelTmax :',nt,time,delTmax
       !..Output the intermediate solutions
       if(delTmax .lt. delTallow) exit
       if( MOD(nt,ntio) .eq. 0 .and. nt .ne. ntmax ) call TECout(nt)
       ENDDO
       call TECout(nt-1)       !..Output the final solution
       print*, 'Heat extracted from holes: ', heat_out, ' J' !....Print heat etracted
       print*, 'Heatdelivered to the blade: ', heat_in, ' J'
       print*, 'Total heat transfer: ', heat_out+heat_in, ' J'
       stop  'DONE'
       end

!------------------------------------------------------------------------
      subroutine INIT()
      Use vars
      character :: fn*16
      logical ok

      !..Read the grid data
      write(*,'(/(a))',advance='no')'  Enter the grid file name [grid.dat]: '
      read(*,'(a)') fn
      if( fn .eq. ' ') fn = 'grid.dat'
      inquire(FILE=fn,EXIST=ok)
      if( .not. ok ) then
       print*, '  ', fn, ' does not exist! \n\n'
       stop
       endif
       open(1,file=fn,form='formatted')
       read(1,*) ncell,nnode
       allocate( node(3,ncell),neigh(3,ncell),xy(2,nnode),area(ncell), &
             Tcell(ncell),Tnode(nnode),dTdx(ncell),dTdy(ncell))
             read(1,*) (no,(xy(i,n),i=1,2),n=1,nnode)
             read(1,*) (no,(node(i,n),i=1,3),(neigh(i,n),i=1,3),n=1,ncell)
             close(1)
             print*, ' # of cells :',ncell
             print*, ' # of nodes :',nnode

             !..Compute cell areas
             do n = 1,ncell
      n1 = node(1,n)
      n2 = node(2,n)
      n3 = node(3,n)
      area(n) = 0.5*( (xy(1,n2)-xy(1,n1))*(xy(2,n3)-xy(2,n1)) &
                     -(xy(2,n2)-xy(2,n1))*(xy(1,n3)-xy(1,n1))  )
                      enddo

!..Set Initial and Boundary Conditions
       Tbc(1)= 1200.
       Tbc(2)= 200.
       Tbc(3)=100.
       Tbc(4)=100.


       !..Initialize the solution
        Tcell = 100.
        
        call TECout(0)

        return
        end

!-------------------------------------------------------------------
         subroutine  GRADIENTa()
         Use vars

         dTdx = 0.
         dTdy = 0.
         DO n = 1,ncell
         do ns = 1,3
        n1 = node(ns,n)
        ni2 = MOD(ns,3) + 1
        n2 = node(ni2,n)
        dx = xy(1,n2)-xy(1,n1)
        dy = xy(2,n2)-xy(2,n1)
        ne = neigh(ns,n)
        if(ne .gt. 0) then       !..averaging Tface for real neighbor
           Tface = 0.5*(Tcell(n) +Tcell(ne))
        else                     !..set Tface to boundary condition
           Tface = Tbc(ne*(-1))
        endif                    !...Calculate gradient components
        dTdx(n) = dTdx(n) + Tface*dy
        dTdy(n) = dTdy(n) - Tface*dx
         enddo
         dTdx(n) = dTdx(n)/area(n)
         dTdy(n) = dTdy(n)/area(n)
         ENDDO

         return
         end


         
         subroutine  GRADIENTb()
         Use vars

         dTdx = 0.
         dTdy = 0.
         DO n = 1,ncell
         do ns = 1,3
        n1 = node(ns,n)
        ni2 = MOD(ns,3) + 1
        n2 = node(ni2,n)
        dx = xy(1,n2)-xy(1,n1)
        dy = xy(2,n2)-xy(2,n1)
        ne = neigh(ns,n)
        if(ne .gt. 0) then       !..averaging Tface for real neighbor
           Tface = 0.5*(Tcell(n) +Tcell(ne))
        else                     !..Set Tface by averaging Tcell and Tfluid values
           Tface = 0.5*(Tcell(n) +Tbc(-1*ne))
        endif                    !..Calculate gradient 
        dTdx(n) = dTdx(n) + Tface*dy
        dTdy(n) = dTdy(n) - Tface*dx
         enddo
         dTdx(n) = dTdx(n)/area(n)
         dTdy(n) = dTdy(n)/area(n)
         ENDDO

         return
         end

!------------------------------------------------------------------------

        function TFLUXa(n,qin,qout)
        Use vars
        real, parameter :: alpha = 9.3E-6 , Nul = 2500 , Nud = 500  
        real, parameter :: d = 0.007 , kcoolair = 0.0018 , khotgas = 0.0075 
        real, parameter :: L = 0.1 , ktit = 22  
        tfluxa = 0.
        qin=0.
        qout=0.
        do ns = 1,3               !..Add the surface fluxes
        n1 = node(ns,n)
        ni2 = MOD(ns,3) + 1
        n2 = node(ni2,n)
        dx = xy(1,n2)-xy(1,n1)
        dy = xy(2,n2)-xy(2,n1)
        !..Compute face fluxes
        ne = neigh(ns,n)
        if( ne .gt. 0 ) then        !..real neighbor...
         	f = alpha*(dTdx(n)+dTdx(ne))*0.5
         	g = alpha*(dTdy(n)+dTdy(ne))*0.5
         elseif( ne .eq. -1 ) then   !..surface 1
         	f = alpha*dTdx(n)
         	g = alpha*dTdy(n)
         elseif( ne .eq. -2 ) then   !..surface 2
         	f = alpha*dTdx(n)			
         	g = alpha*dTdy(n)
         elseif( ne .eq. -3 ) then   !..surface 3
         	f = alpha*dTdx(n)			
         	g = alpha*dTdy(n)
         elseif( ne .eq. -4) then   !..surface 4
         	f = alpha*dTdx(n)		
         	g = alpha*dTdy(n)
         endif
            tfluxa = tfluxa - (f*dy - g*dx)
         if (ne .lt. 0) then				!..Heat rate calculation
         if ( ne .eq. -1 ) then         
         	qin = qin + ktit/alpha*(f*dy-g*dx)
         else 
           qout = qout + ktit/alpha*(f*dy-g*dx)
          endif
          endif
         enddo

         return
         end


         
        function TFLUXb(n,qin,qout)
        Use vars
        real, parameter :: alpha = 9.3E-6 , Nul = 2500. , Nud = 500.  
        real, parameter :: d = 0.007 , kcoolair = 0.018 , khotgas = 0.075 
        real, parameter :: L = 0.1 , ktit = 22., d2 = 0.002, d3 = 0.0025
		qin=0.
        qout=0.
        tfluxb = 0.
        do ns = 1,3               !..Add the surface fluxes
        n1 = node(ns,n)
        ni2 = MOD(ns,3) + 1
        n2 = node(ni2,n)
        dx = xy(1,n2)-xy(1,n1)
        dy = xy(2,n2)-xy(2,n1)
        !..computing the face fluxes
        ne = neigh(ns,n)
        if( ne .gt. 0 ) then        !..real neighbor...
         	f = alpha*(dTdx(n)+dTdx(ne))/2.
         	g = alpha*(dTdy(n)+dTdy(ne))/2.
         elseif( ne .eq. -1 ) then   !..surface 1
         	f = (((Nul/L*khotgas*alpha/ktit)*((Tcell(n)+Tbc(1))/2.-Tbc(1))*dy/(sqrt(dx**2+dy**2)))+alpha*dTdx(n))
         	g = (((Nul/L*khotgas*alpha/ktit)*((Tcell(n)+Tbc(1))/2.-Tbc(1))*dx/(sqrt(dx**2+dy**2)))+alpha*dTdy(n))
         elseif( ne .eq. -2 ) then   !..surface 2
         	f = (((Nud/d*kcoolair*alpha/ktit)*((Tcell(n)+Tbc(2))/2.-Tbc(2))*dy/(sqrt(dx**2+dy**2)))+alpha*dTdx(n))
         	g = (((Nud/d*kcoolair*alpha/ktit)*((Tcell(n)+Tbc(2))/2.-Tbc(2))*dx/(sqrt(dx**2+dy**2)))+alpha*dTdy(n))
         elseif( ne .eq. -3 ) then   !..surface 3
         	f = (((Nud/d2*kcoolair*alpha/ktit)*((Tcell(n)+Tbc(3))/2.-Tbc(3))*dy/(sqrt(dx**2+dy**2)))+alpha*dTdx(n))
         	g = (((Nud/d2*kcoolair*alpha/ktit)*((Tcell(n)+Tbc(3))/2.-Tbc(3))*dx/(sqrt(dx**2+dy**2)))+alpha*dTdy(n))
         elseif( ne .eq. -4 ) then   !..surface 4
         	f = (((Nud/d3*kcoolair*alpha/ktit)*((Tcell(n)+Tbc(4))/2.-Tbc(4))*dy/(sqrt(dx**2+dy**2)))+alpha*dTdx(n))
         	g = (((Nud/d3*kcoolair*alpha/ktit)*((Tcell(n)+Tbc(4))/2.-Tbc(4))*dx/(sqrt(dx**2+dy**2)))+alpha*dTdy(n))
         endif
         tfluxb = tfluxb - (f*dy - g*dx)
          if (ne .lt. 0) then				!..Heat rate calculation
         if ( ne .eq. -1 ) then         
         	qin = qin + ktit/alpha*(f*dy-g*dx)
         else 
           qout = qout + ktit/alpha*(f*dy-g*dx)
          endif
          endif
         enddo

         return
         end

!-------------------------------------------------------------------
       subroutine  TECout(nstep)
!..Output the solution and the grid in TECPLOT format
       Use vars
       character :: fname*32, string*8, ext*5
       call QNODE()             !..Evaluate average temperatures at nodes
!..Set the output file name
        write(string,'(f8.4)') float(nstep)/10000000
        read(string,'(3x,a5)') ext
       fname = 'temp-'//ext//'.tec'
       open(1,file=fname, form='formatted')
        write(1,100) nnode,ncell
       write(1,101) (xy(1,n),xy(2,n),Tnode(n),n=1,nnode)
       write(1,102) (node(1,n),node(2,n),node(3,n),n=1,ncell)
       close(1)
        100 format (' VARIABLES= "X", "Y", "TEMPERATURE"'/, &
              ' ZONE N=', I6,' E=', I6,' F=FEPOINT ',' ET=triangle'  )

        101 format (3(1x,e12.5))
        102 format (3(1x,i6))
        return
          end
!-------------------------------------------------------------------
        subroutine  QNODE()
!..Evaluate node temperatures by averaging the cell temperatures
         Use vars
         integer, dimension(:) :: npass(nnode)
         Tnode = 0.
         npass = 0
!..Find the contribution of cells to the node temperatures
        do n=1,ncell
        do nf=1,3
        nn = node(nf,n)
         Tnode(nn)=Tnode(nn)+Tcell(n)
         npass(nn)=npass(nn)+1
         enddo
        enddo
!..Average the total node temperature with # of contributing cells
           do n=1,nnode
           Tnode(n)=Tnode(n)/npass(n)
           enddo
           return
           end
!------------------------------END----------------------------------