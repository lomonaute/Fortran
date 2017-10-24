module fea

    !! This module contains procedures that are common to FEM-analysis
    
    implicit none
    save
    private
    public :: displ, initial, buildload, buildstiff, enforce, recover

contains

!
!--------------------------------------------------------------------------------------------------
!
    subroutine initial
        
        !! This subroutine is mainly used to allocate vectors and matrices
        
        use fedata
        use link1
        use plane42rect
               
! Hint for continuum elements:nt!eger, parameter :: mdim = 8 
        !integer, dimension(mdim) :: edof
        integer :: e
        ! This subroutine computes the number of global equation,
        ! half bandwidth, etc and allocates global arrays.

        ! Calculate number of equations
        neqn = 2*nn

        if (.not. banded) then
            allocate (kmat(neqn, neqn))
        else
            bw = 0 
			do e=1,ne
            	bw = max( bw, 2*(maxval(element(e)%ix) - minval(element(e)%ix) +1))
                
            end do
            allocate (kmat(bw,neqn))
            
        end if
        allocate (p(neqn), d(neqn))
        allocate (strain(ne, 3), stress(ne, 3))
     
        ! Initial stress and strain
        strain = 0
        stress = 0        
    end subroutine initial
!
!--------------------------------------------------------------------------------------------------
!
    subroutine displ

        !! This subroutine calculates displacements

        use fedata
        use numeth
        use processor

        integer :: e
        real(wp), dimension(:), allocatable :: plotval
		real(wp) :: von_mises,principal_stress(ne, 3) , bcos, bsin
		
        ! Build load-vector
        call buildload

        ! Build stiffness matrix
        call buildstiff

        ! Remove rigid body modes
        call enforce

        if (.not. banded) then
            ! Factor stiffness matrix
                       
            call factor(kmat)
            ! Solve for displacement vector
            call solve(kmat, p)
        else
          call bfactor(kmat)
          call bsolve(kmat, p)
        end if

        ! Transfer results
        d(1:neqn) = p(1:neqn)
        
        ! Recover stress
        call recover
                
        ! Output results
        call output
		!call stopwatch('stop')
        ! Plot deformed shape
        call plot( undeformed + deformed )

        ! Plot element values
        allocate (plotval(ne))
        !print*, stress
        do e = 1, ne
            if (element(e)%id == 1) then
                plotval(e) = stress(e,1)
            else if (element(e)%id == 2) then
                
                plotval(e) = stress(e,1)

            end if
        end do
        !print *, 'Displacement last node', d(neqn - 1),d(neqn)
        call plot( elements, eval=plotval, title="Stress", legend=.true. )
        
		! Von mises stress
        do e= 1, ne
          
                plotval(e) = sqrt(stress(e,1)**2 + stress(e,2)**2 - stress(e,1)*stress(e,2) + 3*stress(e,3)**2)
        
        end do
        call plot( elements, eval=plotval, title="Von Mises Stress", legend=.true. )
        ! Compute principal stress and direction
        do e=1,ne
          principal_stress(e,1) = 0.5*(stress(e,1)+stress(e,2))* sqrt( ( 0.5*(stress(e,1)+stress(e,2)))**2 + (0.5*stress(e,3))**2) 
          principal_stress(e,2) = 0.5*(stress(e,1)-stress(e,2))* sqrt( ( 0.5*(stress(e,1)+stress(e,2)))**2 + (0.5*stress(e,3))**2) 
          bcos = (stress(e,1)-stress(e,2))/(principal_stress(e,1)-principal_stress(e,2))
           
	      principal_stress(e,3) = 0.5 * atan2( bcos , - 2*stress(e,3)/(principal_stress(e,1)-principal_stress(e,2)))
		  
		  plotval(e) = principal_stress(e,3)
        end do
        call plot( elements, eval=plotval, title="Principal Stress 3", legend=.true. )
        

    end subroutine displ
!
!--------------------------------------------------------------------------------------------------
!
    subroutine buildload

        !! This subroutine builds the global load vector

        use fedata
        use plane42rect

        integer :: i, e,j, eface
! Hint for continuum elements:
        integer, dimension(mdim) :: edof
        real(wp), dimension(mdim) :: xe
        real(wp), dimension(mdim) :: re
        real(wp) ::  thk, fe

        ! Build load vector
        p(1:neqn) = 0
        do i = 1, np            
            select case(int(loads(i, 1)))
            case( 1 )
            	! Build nodal load contribution
                p(int(2*loads(i,2)-2+loads(i,3))) = loads(i,4)
                              
            case( 2 )
            	! Build uniformly distributed surface (pressure) load contribution
                !print *, 'ERROR in fea/buildload'
                !print *, 'Distributed loads not defined -- you need to add your own code here'
                !stop
                e     = loads(i, 2)
                eface = loads(i, 3)
                fe	  = loads(i, 4)
				thk   = mprop(element(e)%mat)%thk
                xe    = 0
                edof  = 0 
            do j = 1, element(i)%numnode
                 xe(2*j-1) = x(element(e)%ix(j),1)
                 xe(2*j  ) = x(element(e)%ix(j),2)
                 edof(2*j-1) = 2 * element(e)%ix(j) - 1
                 edof(2*j)   = 2 * element(e)%ix(j)
            end do
                
            call plane42rect_re( xe, eface, fe, thk, re)
                
			do j=1,8 !allocate the loads in {re} onto the {p} vector
            	p(edof(j)) = p(edof(j)) + re(j);
            end do
            
            case default
                print *, 'ERROR in fea/buildload'
                print *, 'Load type not known'
                stop
            end select
        end do
    end subroutine buildload
!
!--------------------------------------------------------------------------------------------------
!
    subroutine buildstiff

        !! This subroutine builds the global stiffness matrix from
        !! the local element stiffness matrices

        use fedata
        use link1
        use plane42rect

        integer :: e, i, j
        integer :: nen 
! Hint for system matrix in band form:
!        integer :: irow, icol
        integer, dimension(mdim) :: edof
        real(wp), dimension(mdim) :: xe
        real(wp), dimension(mdim, mdim) :: ke 
! Hint for modal analysis:
!        real(wp), dimension(mdim, mdim) :: me 
        real(wp) :: young, area 
! Hint for modal analysis and continuum elements:
        real(wp) :: nu, dens, thk

        ! Reset stiffness matrix
        if (.not. banded) then
            kmat = 0
        else
            !print*,'ERROR in fea/buildstiff'
            !print*,'Band form not implemented -- you need to add your own code here'
            !stop

			kmat = 0
            
        end if

        do e = 1, ne
          
            ! Find coordinates and degrees of freedom
            nen = element(e)%numnode
            do i = 1, nen
                 xe(2*i-1) = x(element(e)%ix(i),1)
                 xe(2*i  ) = x(element(e)%ix(i),2)
                 edof(2*i-1) = 2 * element(e)%ix(i) - 1  
                 edof(2*i)   = 2 * element(e)%ix(i)
            end do

            ! Gather material properties and find element stiffness matrix
            select case( element(e)%id )
            case( 1 )
                 young = mprop(element(e)%mat)%young
                 area  = mprop(element(e)%mat)%area
                 call link1_ke(xe, young, area, ke)
            case( 2 )
				 young = mprop(element(e)%mat)%young
                 nu = mprop(element(e)%mat)%nu
                 thk = mprop(element(e)%mat)%thk
            	 call plane42rect_ke(xe, young, nu, thk, ke)
!$$$$$$                  print *, 'ERROR in fea/buildstiff:'
!$$$$$$                  print *, 'Stiffness matrix for plane42rect elements not implemented -- you need to add your own code here'
!$$$$$$                  stop
            end select

            ! Assemble into global matrix
            if (.not. banded) then
                do i = 1, 2*nen
                    do j = 1, 2*nen
                        kmat(edof(i), edof(j)) = kmat(edof(i), edof(j)) + ke(i, j)
                    end do
                end do

! Hint: Can you eliminate the loops above by using a different Fortran array syntax?
            else
                !print *, 'ERROR in fea/buildstiff'
                !print *, 'Band form not implemented -- you need to add our own code here'
                !stop
				do i = 1, 2*nen
                    do j = 1, 2*nen
                      if (  (  (edof(i)+1-edof(j)) <= bw) .AND. (  (edof(i)+1-edof(j)) > 0 )) then                        
                        kmat(edof(i)+1 -edof(j), edof(j)) = kmat(edof(i)+1 -edof(j), edof(j)) + ke(i, j)
                      end if
                    end do
            	end do

                !do i = 1,bw
                !  print '(24(f4.2,tr1))', kmat(i,1:neqn)
                !end do
                !pause


                
            end if 
        end do
    end subroutine buildstiff
!
!--------------------------------------------------------------------------------------------------
!
    subroutine enforce

        !! This subroutine enforces the support boundary conditions

        use fedata

        integer :: i, ii, idof
        real(wp) :: penal

        ! Correct for supports
                
        if (.not. banded) then
            if (.not. penalty) then
              !print '(24(f6.2,tr1))', kmat
              
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    p(1:neqn) = p(1:neqn) - kmat(1:neqn, idof) * bound(i, 3)
                    !print *, 'kmat', kmat(10, idof)
                    p(idof) = bound(i, 3)
                    kmat(1:neqn, idof) = 0
                    kmat(idof, 1:neqn) = 0
                    kmat(idof, idof) = 1
					!print *, 'i = ', i
                    
					!print *, 'bound   ', bound(i,3), '   kmat= ', kmat(10, idof)
                    

                end do
            else
                penal = penalty_fac*maxval(kmat)
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    kmat(idof, idof) = kmat(idof, idof) + penal
                    p(idof) = penal * bound(i, 3)  
                end do  
            end if
        else
            if (.not. penalty) then
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    !print *, kmat(9,2)
                    do ii = 1, neqn
                      
                      if (( ii > idof -bw) .AND. (ii < idof) ) then
                    	p( ii ) = p(ii) - kmat(idof-ii + 1 , ii) * bound(i,3) 
                        !print *, 'Upper area - other indice = ',   idof-ii + 1 ,'ii = ', ii                	  
					  elseif ((ii >= idof) .AND. (ii < idof+bw) ) then
                      	p(ii)   = p(ii) - kmat( ii-idof +1   , idof) * bound(i,3)
                        !print *, 'Lower area - other indice = ',   ii - idof + 1,  '  idof = ', idof, 
                      end if
					end do

                    p(idof) = bound(i, 3)
                    idof = int(2*(bound(i,1)-1) + bound(i,2))  
                    kmat(1:bw, idof) = 0
                    do ii=1,bw
                      if (idof-ii+1 >= 1) then
                    	kmat(ii, idof-ii+1) = 0
                      end if
                    end do
                    kmat(1, idof) = 1
                end do
            else
                penal = penalty_fac*maxval(kmat)
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    kmat(1, idof) = kmat(1, idof) + penal
                    p(idof) = penal * bound(i, 3)  
                end do  
            end if
        end if
				
                !print '(f6.2,tr1)', p
                
                
    end subroutine enforce
!
!--------------------------------------------------------------------------------------------------
!
    subroutine recover

        !! This subroutine recovers the element stress, element strain, 
        !! and nodal reaction forces

        use fedata
        use link1
        use plane42rect

        integer :: e, i, nen
        integer :: edof(mdim)
        real(wp), dimension(mdim) :: xe, de
        real(wp), dimension(mdim, mdim) :: ke
        real(wp) :: young, area 
! Hint for continuum elements:
        real(wp):: nu, dens, thk
        real(wp), dimension(3) :: estrain, estress

        ! Reset force vector
        p = 0

        do e = 1, ne
          
            ! Find coordinates etc...
            nen = element(e)%numnode
            do i = 1,nen
                xe(2*i-1) = x(element(e)%ix(i), 1)
                xe(2*i)   = x(element(e)%ix(i), 2)
                edof(2*i-1) = 2 * element(e)%ix(i) - 1
                edof(2*i)   = 2 * element(e)%ix(i)
                de(2*i-1) = d(edof(2*i-1))
                de(2*i)   = d(edof(2*i))
            end do

            ! Find stress and strain
            select case( element(e)%id )
            case( 1 )
                young = mprop(element(e)%mat)%young
                area  = mprop(element(e)%mat)%area
                call link1_ke(xe, young, area, ke)
                p(edof(1:2*nen)) = p(edof(1:2*nen)) + matmul(ke(1:2*nen,1:2*nen), de(1:2*nen))
                call link1_ss(xe, de, young, estress, estrain)
                stress(e, 1:3) = estress
                strain(e, 1:3) = estrain
            case( 2 )
				young = mprop(element(e)%mat)%young
                nu = mprop(element(e)%mat)%nu
                thk  = mprop(element(e)%mat)%thk
                call plane42rect_ke(xe, young, nu, thk, ke)
                
                p(edof(1:2*nen)) = p(edof(1:2*nen)) + matmul(ke(1:2*nen,1:2*nen), de(1:2*nen))
                call plane42rect_ss(xe, de, young, nu, estress, estrain)
                stress(e, 1:3) = estress
                strain(e, 1:3) = estrain

            end select
        end do
    end subroutine recover

end module fea
