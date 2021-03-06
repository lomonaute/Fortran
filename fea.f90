module fea

    !! This module contains procedures that are common to FEM-analysis
    
    implicit none
    save
    private
    public :: displ, initial, buildload, buildstiff, enforce, recover, mmul, eigen, &
    			trans_loading, enforce_p, enforce_mat, build_trans_load, ccd, newmark

contains

!
!--------------------------------------------------------------------------------------------------
!
    subroutine initial
        
        !! This subroutine is mainly used to allocate vectors and matrices
        
        use fedata
        use link1
        use plane42
               
! Hint for continuum elements:nt!eger, parameter :: mdim = 8 
        !integer, dimension(mdim) :: edof
        integer :: e
        ! This subroutine computes the number of global equation,
        ! half bandwidth, etc and allocates global arrays.

        ! Calculate number of equations
        neqn = 2*nn

        if (.not. banded) then
            allocate (kmat(neqn, neqn))
            allocate (mmat(neqn, neqn))
            allocate (left_side (neqn, neqn))
        else
            bw = 0 
			do e=1,ne
            	bw = max( bw, 2*(maxval(element(e)%ix) - minval(element(e)%ix) +1))
                
            end do
            allocate (kmat(bw,neqn))
            allocate (mmat(bw,neqn))
            allocate (left_side (bw, neqn))
            
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
		real(wp) :: von_mises(ne),principal_stress(ne, 3) , bcos
		integer, parameter :: out_unit=20
        
        ! Build load-vector
        call buildload

        ! Build stiffness matrix
        call buildstiff

        ! Remove rigid body modes
        call enforce

        if (.not. banded) then
            ! Factor stiffness matrix
            !print *,'p', p
        	!  
            call factor(kmat)
            !print *,'kmat', kmat    
            ! Solve for displacement vector
            call solve(kmat, p)
        else
          call bfactor(kmat)
          !print *,'p', p
          !print *,'kmat', kmat 
          
          call bsolve(kmat, p)

          
        end if

        ! Transfer results
        d(1:neqn) = p(1:neqn)
        open (unit = 1, file = "displ_static.txt")
        write(1,*) d
        close(1)
        ! Recover stress
        call recover
                
        ! Output results
        call output
		!call stopwatch('stop')
        ! Plot deformed shape
        call plot( undeformed + deformed )

        ! Plot element values
        allocate (plotval(ne))
        print*, 'de', d
        do e = 1, ne
            if (element(e)%id == 1) then
                plotval(e) = stress(e,1)
            else if (element(e)%id == 2) then
                
                plotval(e) = stress(e,1)

            end if
        end do
		
        call plot( elements, eval=plotval, title="Stress", legend=.true. )
        !print*, 'loop'
        
		! Von mises stress
        do e= 1, ne
			von_mises(e) = sqrt(stress(e,1)**2 + stress(e,2)**2 - stress(e,1)*stress(e,2) + 3*stress(e,3)**2)
        end do
   end subroutine displ
!
!--------------------------------------------------------------------------------------------------
!
    subroutine buildload

        !! This subroutine builds the global load vector

        use fedata
        use plane42

        integer :: i, e,j, eface
! Hint for continuum elements:
        integer, dimension(mdim) :: edof
        real(wp), dimension(mdim) :: xe
        real(wp), dimension(mdim) :: re
        real(wp) ::  thk, fe, dens

        ! Build load vector
        p(1:neqn) = 0
        do i = 1, np            
            select case(int(loads(i, 1)))
            case( 1 )
            	! Build nodal load contribution
                p(int(2*loads(i,2)-2+loads(i,3))) = loads(i,4)
                              
            case( 2 )
            	! Build uniformly distributed surface (pressure) load contribution
                e     = loads(i, 2)
                eface = loads(i, 3)
                fe	  = loads(i, 4)
				thk   = mprop(element(e)%mat)%thk
                dens   = mprop(element(e)%mat)%dens
                xe    = 0
                edof  = 0 
            do j = 1, element(i)%numnode
                 xe(2*j-1) = x(element(e)%ix(j),1)
                 xe(2*j  ) = x(element(e)%ix(j),2)
                 edof(2*j-1) = 2 * element(e)%ix(j) - 1
                 edof(2*j)   = 2 * element(e)%ix(j)
            end do
                
            call plane42_re(xe, eface, fe, thk, dens, re)
                
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
        use plane42

        integer :: e, i, j
        integer :: nen 
! Hint for system matrix in band form:
!        integer :: irow, icol
        integer, dimension(mdim) :: edof
        real(wp), dimension(mdim) :: xe
        real(wp), dimension(mdim, mdim) :: ke 
! Hint for modal analysis:
        real(wp), dimension(mdim, mdim) :: me, mHRZ 
        real(wp) :: young, area 
! Hint for modal analysis and continuum elements:
        real(wp) :: nu, dens, thk

        kmat = 0
        mmat = 0

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
                 dens = mprop(element(e)%mat)%dens
            	 call plane42_ke(xe, young, nu, thk, dens, ke, me, mHRZ)

            end select

            ! Assemble into global matrix
            if (.not. banded) then
                do i = 1, 2*nen
                    do j = 1, 2*nen
                      if (lumpedHRZ) then
                      	mmat(edof(i), edof(j)) = mmat(edof(i), edof(j)) + mHRZ(i, j)
                      else
						mmat(edof(i), edof(j)) = mmat(edof(i), edof(j)) + me(i, j)
                      end if
                    kmat(edof(i), edof(j)) = kmat(edof(i), edof(j)) + ke(i, j)
                    end do
                end do

! Hint: Can you eliminate the loops above by using a different Fortran array syntax?
            else
				do i = 1, 2*nen
                    do j = 1, 2*nen
                      if (  (  (edof(i)+1-edof(j)) <= bw) .AND. (  (edof(i)+1-edof(j)) > 0 )) then                        
                        if (lumpedHRZ) then
                    	mmat(edof(i)+1 -edof(j), edof(j)) = mmat(edof(i)+1 -edof(j), edof(j))+mHRZ(i, j)
                      	else
						mmat(edof(i)+1 -edof(j), edof(j)) = mmat(edof(i)+1 -edof(j), edof(j)) + me(i, j)
                      	end if
						kmat(edof(i)+1 -edof(j), edof(j)) = kmat(edof(i)+1 -edof(j), edof(j)) + ke(i, j)
                       end if
                    end do
            	end do
                
            end if 
        end do
        !print*, 'E  =', young
        !print*, 'rho  =', dens
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
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    p(1:neqn) = p(1:neqn) - kmat(1:neqn, idof) * bound(i, 3)
                    p(idof) = bound(i, 3)
                    kmat(1:neqn, idof) = 0.0_wp
                    kmat(idof, 1:neqn) = 0.0_wp
                    kmat(idof, idof) = 1.0_wp
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
                    do ii = 1, neqn
                      
                      if (( ii > idof -bw) .AND. (ii < idof) ) then
                    	p( ii ) = p(ii) - kmat(idof-ii + 1 , ii) * bound(i,3) 
					  elseif ((ii >= idof) .AND. (ii < idof+bw) ) then
                      	p(ii)   = p(ii) - kmat( ii-idof +1   , idof) * bound(i,3)
                      end if
					end do

                    p(idof) = bound(i, 3)
                    idof = int(2*(bound(i,1)-1) + bound(i,2))  
                    kmat(1:bw, idof) = 0.0_wp
                    do ii=1,bw
                      if (idof-ii+1 >= 1) then
                    	kmat(ii, idof-ii+1) = 0.0_wp
                      end if
                    end do
                    kmat(1, idof) = 1.0_wp
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
               
    end subroutine enforce
!
!--------------------------------------------------------------------------------------------------
!
    subroutine recover

        !! This subroutine recovers the element stress, element strain, 
        !! and nodal reaction forces

        use fedata
        use link1
        use plane42

        integer :: e, i, nen
        integer :: edof(mdim)
        real(wp), dimension(mdim) :: xe, de
        real(wp), dimension(mdim, mdim) :: ke
        real(wp), dimension(mdim, mdim) :: me, mHRZ
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
                dens  = mprop(element(e)%mat)%dens
                call plane42_ke(xe, young, nu, thk, dens, ke, me, mHRZ)
                
                p(edof(1:2*nen)) = p(edof(1:2*nen)) + matmul(ke(1:2*nen,1:2*nen), de(1:2*nen))

				call plane42_ss(xe, de, young, nu, estress, estrain)
                stress(e, 1:3) = estress
                strain(e, 1:3) = estrain

            end select
        end do
    end subroutine recover
!
!--------------------------------------------------------------------------------------------------
!
 subroutine mmul(invector, outvector, mtype)
        !! This subroutine builds the mass matrix and the eigenvector

        use fedata
        use link1
        use plane42

		real(wp), dimension(neqn), intent(in) :: invector
        integer, intent(in) :: mtype
		real(wp), dimension(:), intent(out) :: outvector
        
        integer :: e, i, j
        integer :: nen 
! Hint for system matrix in band form:
        integer :: irow, icol
        integer, dimension(mdim) :: edof
        real(wp), dimension(mdim) :: xe
        real(wp), dimension(mdim, mdim) :: ke 
! Hint for modal analysis:
        real(wp), dimension(mdim, mdim) :: me, mHRZ 
        real(wp) :: young, area 
! Hint for modal analysis and continuum elements:
        real(wp) :: nu, dens, thk, mass
        outvector = 0.0_wp
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
            	print*, 'Not implemented for truss structures'
                stop
            case( 2 )
				dens = mprop(element(e)%mat)%dens
                 thk = mprop(element(e)%mat)%thk
                 nu = mprop(element(e)%mat)%nu
                 young = mprop(element(e)%mat)%young                 
            	call plane42_ke(xe, young, nu, thk, dens, ke, me, mHRZ)
            end select
			mass = 0
            do i =1, mdim
              mass = mass +	me(i,i)
            end	do
            do i = 1, mdim              
              do j = 1, mdim
                if (mtype == 1) then
					outvector(edof(i)) = outvector(edof(i)) + ke(i,j) * invector(edof(j))
                	elseif (mtype == 2) then
                	outvector(edof(i)) = outvector(edof(i)) + me(i,j) * invector(edof(j))
                	elseif (mtype == 3) then
                    outvector(edof(i)) = outvector(edof(i)) + mHRZ(i,j) * invector(edof(j))
                	else 
                  		print*, 'Problem in mmul subroutine - wrong value for mtype :', mtype
                    	stop
                end if                        	

			  end do
            end do  
        end do
        
    end subroutine mmul

!
!--------------------------------------------------------------------------------------------------
!
    subroutine eigen

        !! This subroutine calculates displacements

        use fedata
        use numeth
        use processor
        
        logical :: error
		integer :: i, j, p_it, n, idof
		integer, parameter :: p_max=1000, neig = 5, x_val=5
        real(wp), dimension(neqn) :: x_vector, x_vector_p_1, y_vector, y_vector_p_1, d_mode
        real(wp), dimension(neqn, neig) :: z_vector, d_vector
        real(wp), dimension(neig) :: c_var, lambda, omega
        real(wp) :: rp
        real(dp), parameter :: epsilon = 10.0**(-16)
        
        ! Build load-vector
        call buildload

        ! Build stiffness matrix
        call buildstiff

        ! Remove rigid body modes
        call enforce
		
		! Factor stiffness matrix  
        if (.not. banded) then                
            call factor(kmat)
        else
          	call bfactor(kmat)
        end if
		
		! Initialize {D}_0 and {Z}_0
        d_vector = 0.0_wp
        z_vector = 0.0_wp
        
		do i = 1, neig
			!Initial guess for eigenvector {X}_0 and boundary conditions
			x_vector = 0.5*1.0_wp
        	
        	do n = 1, nb
        		idof = int(2 * (bound(n,1)-1) + bound(n,2))
        		x_vector(idof) = bound(n, 3)
        	end do

            call mmul(x_vector, y_vector,2)
            
			! Compute {Z}j = [M]{D}j	
			do j=1, i-1
            	call mmul(d_vector(1:neqn, j), z_vector(1:neqn, j),2)                
            end do

        	! Iteration 
        	do p_it = 1,p_max
            	x_vector_p_1 = x_vector
            	y_vector_p_1 = y_vector
          	! Solve KX = Y
				if (.not. banded) then
            		call solve(kmat, y_vector)
        		else
          			call bsolve(kmat, y_vector)
        		end if
		      
            	! Transfer results
            	x_vector(1:neqn) = y_vector(1:neqn)
			
				do n = 1, nb
        			idof = int(2 * (bound(n,1)-1) + bound(n,2))
        			x_vector(idof) = bound(n, 3)
            	end do
                
            	do j=1, i-1
              		c_var(j) =  dot_product(x_vector, z_vector(1:neqn,j))
              		x_vector = x_vector - c_var(j)*d_vector(1:neqn,j)
            	end do
		              
				! Compute new Y
            	call mmul(x_vector, y_vector,2)
			
            	! Compute rp
            	rp = sqrt(dot_product(x_vector, y_vector)) 
            	y_vector = y_vector / rp
			
            	error = ((sqrt(dot_product(x_vector - x_vector_p_1, x_vector - x_vector_p_1)) &
              		/ sqrt(dot_product(x_vector, x_vector))) < epsilon)
            
            
				if (error) then
              		print *, 'Eigenmode calculation converged'
              		exit
				end if
        	end do       
        
			d_vector(1:neqn,i) = x_vector / rp
        	lambda(i) = dot_product(x_vector, y_vector_p_1) / rp**2 		
        end do

        
        omega = sqrt(lambda)
        !print*, 'omega', omega
        d = reshape(d_vector(1:neqn, 5), (/neqn/))
        !call plot( what = 16, eigenfreq=lambda(2), eigenvector=d, tend =10.d0)

		!open (unit = 1, file = "results.txt")
        !do i=1,5
        !write(1,*) i,  omega(i)
		!write(1,'(2f20.8)') reshape(d_vector(1:neqn, i), (/neqn/))  
		!end do
        !close(1)
   end subroutine eigen

!--------------------------------------------------------------------------------------------------
!
    subroutine trans_loading

        !! Classical Central Difference (CCD) method

        use fedata
        use numeth
        use processor

        integer :: e, i, n_node, method
        real(wp), dimension(:), allocatable :: plotval
!$$$$$$         real(wp), dimension(neqn,neqn) :: cmat
!$$$$$$         real(wp) :: von_mises(ne), principal_stress(ne, 3)
!$$$$$$         real(wp) :: D_n_1(neqn), D_n_2(neqn), D_dot_n_1(neqn), D_2dot_n_1(neqn), d_mode(neqn)
        real(wp) :: D_n_1(neqn), D_dot_n_1(neqn), D_2dot_n_1(neqn), d_mode(neqn)
        real(wp) :: t, omega, transient_contribution
        real(wp), parameter :: delta_t = 0.05
        real(wp), parameter :: alpha = 0.8
        real(wp), parameter :: beta = 0.8
        integer, parameter :: t_steps = 1000
        real(wp) :: storage(6,t_steps)

		method = 3
		n_node = neqn / 2 - 2
        t = 0
        omega = 10*pi
        transient_contribution = 0

		! Build the steady state part of the loadvector
        call buildload
        
        ! Build stiffness matrix
		call buildstiff

        ! Remove rigid body modes
        call enforce

		! transient
        !open (unit = 1, file = "results1.txt")
		
!$$$$$$         call ccd(delta_t, t_steps, alpha, beta, n_node, storage, omega, transient_contribution)
        call newmark(delta_t, t_steps, omega, alpha, beta, method, n_node, transient_contribution, storage)
     	
		open (unit = 1, file = "trans_u_y.txt")	
		write(1,*) storage(2,:)
		close(1)       
             
        ! Recover stress
        call recover
        
        ! call stopwatch('stop')       
        ! Output results
        call output
        
		
        ! Plot deformed shape
        call plot( undeformed + deformed )

        ! Plot element values
        allocate (plotval(ne))
        !print*, 'de', d
        do e = 1, ne
            if (element(e)%id == 1) then
                plotval(e) = stress(e,1)
            else if (element(e)%id == 2) then
                
                plotval(e) = stress(e,1)

            end if
        end do

   end subroutine trans_loading
!
!--------------------------------------------------------------------------------------------------
!
    subroutine enforce_p

        !! This subroutine enforces the support boundary conditions

        use fedata

        integer :: i, ii, idof
        real(wp) :: penal

        ! Correct for supports
                
        if (.not. banded) then
            if (.not. penalty) then
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    p(1:neqn) = p(1:neqn) - kmat(1:neqn, idof) * bound(i, 3)
                    p(idof) = bound(i, 3)
                end do
            else
                penal = penalty_fac*maxval(kmat)
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    p(idof) = penal * bound(i, 3)  
                end do  
            end if
        else
            if (.not. penalty) then
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    do ii = 1, neqn
                      
                      if (( ii > idof -bw) .AND. (ii < idof) ) then
                    	p( ii ) = p(ii) - kmat(idof-ii + 1 , ii) * bound(i,3) 
					  elseif ((ii >= idof) .AND. (ii < idof+bw) ) then
                      	p(ii)   = p(ii) - kmat( ii-idof +1   , idof) * bound(i,3)
                      end if
					end do

                    p(idof) = bound(i, 3)
!$$$$$$                     idof = int(2*(bound(i,1)-1) + bound(i,2)) 
                end do
            else
                penal = penalty_fac*maxval(kmat)
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    p(idof) = penal * bound(i, 3)  
                end do  
            end if
        end if
               
    end subroutine enforce_p
!
!--------------------------------------------------------------------------------------------------
!
    subroutine enforce_mat(mat)

        !! This subroutine enforces the support boundary conditions

        use fedata
		
		real(wp), dimension(:,:), intent(inout) :: mat
        integer :: i, ii, idof
        real(wp) :: penal

        ! Correct for supports
                
        if (.not. banded) then
            if (.not. penalty) then
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    mat(1:neqn, idof) = 0.0_wp
                    mat(idof, 1:neqn) = 0.0_wp
                    mat(idof, idof) = 1.0_wp
                end do
            else
                penal = penalty_fac*maxval(mat)
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    mat(idof, idof) = mat(idof, idof) + penal
                end do  
            end if
        else
            if (.not. penalty) then
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))  
                    mat(1:bw, idof) = 0.0_wp
                    do ii=1,bw
                      if (idof-ii+1 >= 1) then
                    	mat(ii, idof-ii+1) = 0.0_wp
                      end if
                    end do
                    mat(1, idof) = 1.0_wp
                end do
            else
                penal = penalty_fac*maxval(mat)
                do i = 1, nb
                    idof = int(2*(bound(i,1)-1) + bound(i,2))
                    mat(1, idof) = mat(1, idof) + penal
                end do  
            end if
        end if
               
    end subroutine enforce_mat
!
!--------------------------------------------------------------------------------------------------
!
!$$$$$$     subroutine build_trans_load(t, omega, transient_contribution, out_vector)
!$$$$$$ 
!$$$$$$         !! This subroutine builds the global load vector
!$$$$$$ 
!$$$$$$         use fedata
!$$$$$$         use plane42
!$$$$$$ 
!$$$$$$         integer :: i, e, j, eface
!$$$$$$ ! Hint for continuum elements:
!$$$$$$         integer, dimension(mdim) :: edof
!$$$$$$         real(wp), dimension(mdim) :: xe
!$$$$$$         real(wp), dimension(mdim) :: re
!$$$$$$         real(wp) ::  thk, fe, dens
!$$$$$$         real(wp), intent(in) :: t, omega, transient_contribution
!$$$$$$         real(wp), intent(out) :: out_vector(1:neqn)
!$$$$$$ 
!$$$$$$         ! Build load vector
!$$$$$$         out_vector(1:neqn) = 0
!$$$$$$         
!$$$$$$         do i = 1, np            
!$$$$$$             select case(int(loads(i, 1)))
!$$$$$$             case( 1 )
!$$$$$$                 ! Build nodal load contribution
!$$$$$$                 out_vector(int(2*loads(i,2)-2+loads(i,3))) = loads(i,4)* (1 + sin(t * omega) * transient_contribution)
!$$$$$$                               
!$$$$$$             case( 2 )
!$$$$$$                 ! Build uniformly distributed surface (pressure) load contribution
!$$$$$$                 e     = loads(i, 2)
!$$$$$$                 eface = loads(i, 3)
!$$$$$$                 fe    = loads(i, 4) * (1 + sin(t * omega) * transient_contribution)
!$$$$$$                 thk   = mprop(element(e)%mat)%thk
!$$$$$$                 dens   = mprop(element(e)%mat)%dens
!$$$$$$                 xe    = 0
!$$$$$$                 edof  = 0 
!$$$$$$             do j = 1, element(i)%numnode
!$$$$$$                  xe(2*j-1) = x(element(e)%ix(j),1)
!$$$$$$                  xe(2*j  ) = x(element(e)%ix(j),2)
!$$$$$$                  edof(2*j-1) = 2 * element(e)%ix(j) - 1
!$$$$$$                  edof(2*j)   = 2 * element(e)%ix(j)
!$$$$$$             end do
!$$$$$$                 
!$$$$$$             call plane42_tur_blade(xe, eface, fe, thk, dens, re)
!$$$$$$                 
!$$$$$$             do j=1,8 !allocate the loads in {re} onto the {p} vector
!$$$$$$                 out_vector(edof(j)) = out_vector(edof(j)) + re(j);
!$$$$$$             end do
!$$$$$$             
!$$$$$$             case default
!$$$$$$                 print *, 'ERROR in fea/buildload'
!$$$$$$                 print *, 'Load type not known'
!$$$$$$                 stop
!$$$$$$             end select
!$$$$$$         end do     
!$$$$$$     end subroutine build_trans_load
!$$$$$$ !
!$$$$$$ !--------------------------------------------------------------------------------------------------
!$$$$$$ !       
    subroutine build_trans_load(steady_load, t, omega, transient_contribution)
		! This subroutine multiplies the load vector to find the transient load
        use fedata

        real(wp), intent(in) :: steady_load(1:neqn), t, omega, transient_contribution

		p = steady_load * (1 + transient_contribution * sin(t * omega))

    end subroutine build_trans_load
!
!--------------------------------------------------------------------------------------------------
!   
    subroutine ccd(delta_t, t_steps, alpha, beta, n_node, storage, omega, transient_contribution)
    	
		use fedata
        use numeth
        use processor
        
    	!real(wp), intent(in), dimension(:) :: p
        real(wp), intent(in) :: delta_t, alpha, beta, omega, transient_contribution
        integer, intent(in) :: n_node, t_steps
        
		real(wp), dimension(neqn) :: steady_load, D_n_1, D_n_2
        real(wp), dimension(6, t_steps) :: storage
        integer :: i
        
        real(wp) :: D_dot_n_1(neqn), D_2dot_n_1(neqn), t
        real(wp) :: int_vect_1(neqn), int_vect_2(neqn), int_vect_3(neqn), int_vect_4(neqn), int_vect_5(neqn)
        
		D_n_1 = 0
        D_n_2 = 0
        D_dot_n_1 = 0
        D_2dot_n_1 = 0

        storage = 0
        steady_load = p
        
		do i = 1, t_steps
            t = i * delta_t
            p = 0
        	call build_trans_load(steady_load, t, omega, transient_contribution)
		
			! construct D
            call mmul(D_n_1, int_vect_1, 1)
        	if (lumpedHRZ) then
          		call mmul(D_n_1, int_vect_2, 3)
          		call mmul(D_n_2, int_vect_3, 3)
        	else
          		call mmul(D_n_1, int_vect_2, 2)
          		call mmul(D_n_2, int_vect_3, 2)
        	end if
        	call mmul(D_n_2, int_vect_5, 1)
			int_vect_4 = alpha * int_vect_3 + beta * int_vect_5
                        
			p = p - int_vect_1 + 2.0_wp * int_vect_2 / delta_t**2 - int_vect_3 / delta_t**2 &
            	+ (int_vect_4) / (2.0_wp * delta_t)
                
        	call enforce_p
    		D = p
       	
        	if (.not. banded) then
           		left_side = mmat * (1 / delta_t**2) + (alpha * mmat + beta * kmat) / (2 * delta_t)
            	call enforce_mat(left_side)
            	call factor(left_side)
            	call solve(left_side, D)
        	else
          		left_side = mmat * (1 / delta_t**2)+ (alpha * mmat + beta * kmat) / (2 * delta_t)
            	call enforce_mat(left_side)
        		call bfactor(left_side)
            	call bsolve(left_side, D)
        	end if
            
			D_dot_n_1 = 0.5_wp * (D - D_n_2) / delta_t
        	D_2dot_n_1 = (D - 2.0_wp * D_n_1 + D_n_2) / delta_t**2.0_wp

        	storage(1, i) = D(2 * n_node - 1)
        	storage(2, i) = D(2 * n_node)
        	storage(3, i) = D_dot_n_1(2 * n_node - 1)
        	storage(4, i) = D_dot_n_1(2 * n_node)
        	storage(5, i) = D_2dot_n_1(2 * n_node - 1)
        	storage(6, i) = D_2dot_n_1(2 * n_node)

			D_n_2 = D_n_1
        	D_n_1 = D
        end do
        
	end subroutine ccd
!
!--------------------------------------------------------------------------------------------------
!
    subroutine newmark(delta_t, t_steps, omega, alpha, beta, method, &
    				   n_node, transient_contribution, storage)
		
        use fedata
        use numeth
        use processor
        
    	!real(wp), intent(in) :: p(neqn)
        real(wp), intent(in) :: delta_t, alpha, beta, transient_contribution, omega
        integer, intent(in) :: n_node, method, t_steps
        real(wp), intent(out), dimension(6, t_steps) :: storage
        
        real(wp) :: newmark_beta, newmark_gamma, t
        real(wp) :: D_n_1(neqn), D_n_2(neqn), D_dot_n_1(neqn), D_2dot_n_1(neqn)
        real(wp) :: r1(neqn), r2(neqn), r3(neqn)
        real(wp) :: int_vect_1(neqn), int_vect_2(neqn), int_vect_3(neqn)
        real(wp) :: D_dot(neqn), D_2dot(neqn), steady_load(neqn)
        integer :: i
		
		if (method == 0) then  ! average acceleration
        	newmark_beta = 0.25
            newmark_gamma = 0.5
        elseif (method == 1) then  ! linear acceleration
        	newmark_beta = 1.0 / 6.0
            newmark_gamma = 0.5
        elseif (method == 2) then  ! Fox-Goodwin
        	newmark_beta = 1.0 / 12.0
            newmark_gamma = 0.5
        elseif (method == 3) then  ! algorithmically damped
        	newmark_gamma = 0.5
            if (newmark_gamma < 0.5) then
              	newmark_gamma = 0.5
            end if
            newmark_beta = (newmark_gamma + 0.5_wp)**2 * 0.25
        else
            newmark_beta = 0.25
            newmark_gamma = 0.5
            print*, 'Not right method selected in fea/newmark, average acceleration method selected'
		end if
		
		! construct left side
		left_side = mmat * (1.0_wp / (newmark_beta * delta_t**2.0_wp)) &
        			+ (alpha * mmat + beta * kmat) * (newmark_gamma / (newmark_beta * delta_t)) + kmat
		
		if (.not. banded) then
          	call enforce_mat(left_side)
            call factor(left_side)
    	else
        	call enforce_mat(left_side)
      		call bfactor(left_side)
        end if

		D_n_1 = 0
        D_n_2 = 0
        D_dot_n_1 = 0
        
        storage = 0
        steady_load = p

        call mmul(D_n_1, int_vect_1, 1)
        call mmul(D_dot_n_1, int_vect_2, 1)
        
        if (lumpedHRZ) then
          	call mmul(D_dot_n_1, int_vect_3, 3)
        else
          	call mmul(D_dot_n_1, int_vect_3, 2)
        end if
            
		p = p - int_vect_1 - alpha * int_vect_2 - beta * int_vect_3
        call enforce_p

        D_2dot_n_1 = p
        
		call enforce_mat(mmat)
        if (.not. banded) then
            call solve(mmat, D_2dot_n_1)
    	else
            call bsolve(mmat, D_2dot_n_1)
        end if

        storage(1, 1) = D_n_1(2 * n_node - 1)
        storage(2, 1) = D_n_1(2 * n_node)
        storage(3, 1) = D_dot_n_1(2 * n_node - 1)
        storage(4, 1) = D_dot_n_1(2 * n_node)
        storage(5, 1) = D_2dot_n_1(2 * n_node - 1)
        storage(6, 1) = D_2dot_n_1(2 * n_node)

        do i = 2,t_steps
        	t = (i - 1) * delta_t
            
            p = 0
        	call build_trans_load(steady_load, t, omega, transient_contribution)
			
!$$$$$$             print*, ' '
!$$$$$$             print*, p
           
            ! construct right side
            r1 = D_n_1 * (1.0_wp / (newmark_beta * delta_t**2.0_wp)) &
                    + D_dot_n_1 * (1.0_wp / (newmark_beta * delta_t)) &
                    + D_2dot_n_1 * (1.0_wp / (2.0_wp * newmark_beta) - 1.0_wp)
            
            r2 = D_n_1 * (newmark_gamma / (newmark_beta * delta_t)) &
                    + D_dot_n_1 * (newmark_gamma / newmark_beta - 1) &
                    + D_2dot_n_1 * delta_t * (newmark_gamma / (2.0_wp * newmark_beta) - 1.0_wp)
    
            r3 = r2
            
            call mmul(r1, int_vect_1, 1)
            call mmul(r2, int_vect_2, 1)
                    
!$$$$$$             print*, ' '
!$$$$$$             print*, int_vect_1
            
            if (.not. lumpedHRZ) then
              	call mmul(r1, int_vect_1, 2)
            	call mmul(r2, int_vect_2, 2)
    		else
              	call mmul(r1, int_vect_1, 3)
            	call mmul(r2, int_vect_2, 3)
        	end if

            call mmul(r3, int_vect_3, 1)

!$$$$$$             print*, 'int_vect_3'
!$$$$$$             print*, int_vect_3

            p = p + int_vect_1 + alpha * int_vect_2 + beta * int_vect_3

            call enforce_p
    		D = p

!$$$$$$             print*, 'D ='
!$$$$$$             print*, D            

        	if (.not. banded) then
            	call solve(left_side, D)
    		else
            	call bsolve(left_side, D)
        	end if

!$$$$$$             print*, 'D ='
!$$$$$$             print*, D
        
            D_dot = (D - D_n_1) * (newmark_gamma / (newmark_beta * delta_t)) &
                    - D_dot_n_1 * (newmark_gamma / newmark_beta - 1) &
                    - delta_t * D_2dot_n_1 * ((newmark_gamma / (2.0_wp * newmark_beta)) - 1.0_wp)
    
            D_2dot = (D - D_n_1 - delta_t * D_dot_n_1) * (1.0_wp / (newmark_beta * delta_t**2)) &
                    - D_2dot_n_1 * ((1.0_wp / (2.0_wp * newmark_beta)) - 1.0_wp)
            
            storage(1, i) = D(2 * n_node - 1)
            storage(2, i) = D(2 * n_node)
            storage(3, i) = D_dot(2 * n_node - 1)
            storage(4, i) = D_dot(2 * n_node)
            storage(5, i) = D_2dot(2 * n_node - 1)
            storage(6, i) = D_2dot(2 * n_node)
            
            D_n_1 = D
            D_dot_n_1 = D_dot
            D_2dot_n_1 = D_2dot
            
!$$$$$$             print*, ' '
!$$$$$$             print*, D_n_1 - D
            
    	end do  
        
	end subroutine newmark
    
end module fea
