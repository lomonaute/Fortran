module plane42 

    !! This module contains subroutines specific to the plane42 element
    !!        
    !! The plane42 element has 4 nodes. Each node has 2 degrees-of-freedom,
    !! namely, displacement along the \(x\)- and \(y\)-coordinate directions.
    !!
    !! The nodes are numbered counter-clockwise as indicated below. The
    !! figure also shows how the element edges are labelled. For example,
    !! edge 1 is the edge between element node 1 and 2.
    !!
    !!       N4    E3    N3
    !!          o------o 
    !!          |      |
    !!       E4 |      | E2
    !!          |      |
    !!          o------o
    !!       N1    E1    N2
    !!             
    !!
    !! `N1` = element node 1, `N2` = element node 2, etc  
    !! `E1` = element face 1, `E2` = element face 2, etc  
    
    use types
    implicit none
    save
    
    private
    public :: plane42_ke, plane42_re, plane42_vo, plane42_ss, shape, plane42_tur_blade

contains

    subroutine plane42_ke(xe, young, nu, thk, dens, ke, me, mHRZ)

        !! This subroutine constructs the stiffness matrix for
        !! a isoparametric 4-noded quad element.

        real(wp), intent(in) :: young
            !! Young's Modulus for this element
        real(wp), intent(in) :: nu
            !! Poisson's Ratio for this element
        real(wp), intent(in) :: thk
			!! density
        real(wp), intent(in) :: dens
            !! Thickness of this element
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed configuration
            !!
            !! * `xe(1:2)` = \((x,y)\)-coordinates of element node 1
            !! * `xe(3:4)` = \((x,y)\)-coordinates of element node 1
            !! * `xe(5:6)` = \((x,y)\)-coordinates of element node 2
            !! * `xe(7:8)` = \((x,y)\)-coordinates of element node 2
            !!
            !! See also [[plane42rect]]
        real(wp), dimension(:,:), intent(out) :: ke
            !! Stiffness matrix
        real(wp), dimension(:,:), intent(out) :: me
        	!! mass element matrix
        real(wp), dimension(:,:), intent(out) :: mHRZ
        
        real(wp) :: cmat(3,3), fact
        real(wp) :: gaus_co(2), gaus_w
		integer  :: i,j
        real(wp), dimension(3, 8) :: bmat

        real(wp), dimension(2, 2) :: jac
        real(wp), dimension(2, 8) :: n
			!! Jacobian matrix
        real(wp) :: detjac, mass
			!! Determinant of the jacobian matrix

            
        gaus_co = (/ 1.0_wp/sqrt(3.0_wp), -1.0_wp/sqrt(3.0_wp) /)

        gaus_w = 1.0_wp

        ! build constitutive matrix (plane stress)
        cmat = 0.0_wp
        
		fact = young / (1.0_wp - nu**2.0_wp)
        cmat(1,1) = fact
        cmat(1,2) = fact*nu
        cmat(2,1) = fact*nu
        cmat(2,2) = fact
        cmat(3,3) = fact*(1.0_wp-nu)/2.0_wp
		
		ke = 0.0_wp
        me = 0.0_wp
        do i = 1,2
			do j = 1,2

            	call shape(xe, gaus_co(i), gaus_co(j), n, bmat, jac, detjac)
				ke = ke + gaus_w * thk * matmul(transpose(bmat), matmul(cmat, bmat)) * detjac
				me = me + gaus_w * thk * dens * matmul(transpose(n), n) * detjac                
			end do
        end do
		
		
        mHRZ = 0
        mass = 0
        do i =1, int(sqrt(1.0*size(me)))
           mass = mass +	me(i,i)
        end	do
        do i=1,int(sqrt(1.0*size(me)))
          mHRZ(i,i) = me(i,i)*sum(me)/mass
        end do 

    end subroutine plane42_ke
!
!--------------------------------------------------------------------------------------------------
!
    subroutine plane42_re(xe, eface, fe, thk, dens, re)

        !! This subroutine computes the element load vector due
        !! to surface traction (traction is always perpendicular
        !! to element face).

        integer, intent(in) :: eface
            !! Element face where traction (pressure) is applied

        real(wp), intent(in) :: fe
            !! Value of surface traction (pressure)
        real(wp), intent(in) :: thk, dens
            !! Thickness of this element
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed configuration (see also [[plane42rect_ke]])
        real(wp), intent(out) :: re(8)
            !! Element force vector
            !!
            !! * `re(1:2)` = \((f_x^1, f_y^1)\) force at element node 1 in \(x\)- and y-direction
            !! * `re(3:4)` = \((f_x^2, f_y^2)\) force at element node 1 in \(x\)- and y-direction
            !! * etc...
        real(wp) :: eta, xi, gaus_co(2), gaus_w
        real(wp) :: jac(2,2), n(2,8), bmat(3,8), detjac
        
        integer :: i

		gaus_w = 1.0_wp
        gaus_co = (/ 1.0_wp/sqrt(3.0_wp), -1.0_wp/sqrt(3.0_wp) /)
        re = 0
        if (eface == 1) then
          do i=1,2
				eta = -1.0_wp
                xi = gaus_co(i)
				
				call shape(xe, xi, eta, n, bmat, jac, detjac)
				
				re = re + gaus_w * thk * fe * matmul( transpose(n), (/ -jac(1,2), jac(1,1) /) )
            end do
            
        elseif (eface == 2) then
        	do i=1,2
				eta = gaus_co(i)
                xi = 1.0_wp
				
        		call shape(xe, xi, eta, n, bmat, jac, detjac)

				re = re + gaus_w * thk * fe * matmul(transpose(n), (/ -jac(2,2), jac(2,1) /) )
            end do
            
        elseif (eface == 3) then
 			do i=1,2
				eta = 1.0_wp
                xi = gaus_co(i)
				
				call shape(xe, xi, eta, n, bmat, jac, detjac)
				
				re = re + gaus_w * thk * fe * matmul(transpose(n), (/ jac(1,2), -jac(1,1) /) )
            end do 
        elseif (eface == 4) then
			do i=1,2
				eta = gaus_co(i)
                xi = -1.0_wp
                				
        		call shape(xe, xi, eta, n, bmat, jac, detjac)
                
				re = re + gaus_w * thk * fe * matmul(transpose(n), (/ jac(2,2), -jac(2,1) /) )
            end do  
        endif
        
    end subroutine plane42_re
!
!--------------------------------------------------------------------------------------------------
!
    subroutine plane42_vo(xe, acel, dens, thk, re)

        !! This subroutine computes the element load vector due
        !! to volumetric forces

        real(wp), intent(in) :: thk
            !! Thickness of this element
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed configuration (see also [[plane42rect_ke]])
        real(wp), intent(in) :: dens
        	!! Density of the material
        real(wp), dimension(2), intent(in) :: acel
        	!! Volume acceleration
        real(wp), intent(out) :: re(8)
            !! Element force vector
        real(wp), dimension(2) :: gaus_co
        real(wp) :: eta, xi, gaus_w, phi(2)
        real(wp) :: jac(2,2), detjac, n(2,8), bmat(3,8)
        integer :: i, j
		
		gaus_w = 1.0_wp
        gaus_co = (/ -1.0_wp/sqrt(3.0_wp), 1.0_wp/sqrt(3.0_wp) /)
        phi = 4.0_wp * thk * dens * transpose(acel)
        
		re = 0 
        do i = 1,2
          xi = gaus_co(i)
          do j = 1,2
            eta = gaus_co(j)

            call shape(xe, xi, eta, n, bmat, jac, detjac)
!			          
            re = re + gaus_w * thk * matmul(transpose(n), phi) * detjac
          end do
        end do
	
    end subroutine plane42_vo
!
!--------------------------------------------------------------------------------------------------
!
    subroutine plane42_ss(xe, de, young, nu, estress, estrain)

        !! This subrotuine computes the element stress and strain (The location inside the element
        !! where stress and and strain is evaluated, is defined inside the subroutine).

        real(wp), intent(in) :: young
            !! Young's Modulus for this element
        real(wp), intent(in) :: nu 
            !! Poisson's Ratio for this element
        real(wp), dimension(:), intent(in)  :: xe
            !! Nodal coordinates of this element in undeformed configuration (see also [[plane42rect_ke]])
        real(wp), dimension(:), intent(in)  :: de
            !! Displacement field

        real(wp), dimension(:), intent(out) :: estress
            !! Stress at a point inside the element

        real(wp), dimension(:), intent(out) :: estrain

        real(wp), dimension(2,8) :: n
            !! Strain at a point inside the element

        real(wp) :: bmat(3, 8), cmat(3, 3), fact !, princ_stress, direc_stress, von_mises 

        real(wp) :: jac(2,2), detjac
        
        ! Build strain-displacement matrix
		call shape(xe, 0.0_wp, 0.0_wp, n, bmat, jac, detjac)

        ! Compute element strain
        estrain = matmul(bmat, de)
      ! Build constitutive matrix (plane stress)
        cmat = 0
		fact = young / (1.0_wp - nu**2.0_wp)
		cmat(1,1) = fact
        cmat(1,2) = fact*nu
        cmat(2,1) = fact*nu
        cmat(2,2) = fact
        cmat(3,3) = fact*(1.0_wp-nu)/2.0_wp

        ! Compute element stress
        estress = matmul(cmat, estrain)
		! Von mises stress

    end subroutine plane42_ss
!
!--------------------------------------------------------------------------------------------------
!
	subroutine shape(xe, xi, eta, n, bmat, jac, detjac)
		
		real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed configuration (see also [[plane42rect_ke]])
        real(wp), intent(in) :: xi
            !! Xi coordinates in the natural coordinate system
        real(wp), intent(in) :: eta
            !! Eta coordinates in the natural coordinate system
       
		real(wp), dimension(2, 8), intent(out) :: n
		real(wp), dimension(3, 8), intent(out) :: bmat

        real(wp), dimension(2, 2), intent(out) :: jac
			!! Jacobian matrix
        real(wp), intent(out) :: detjac
			!! Determinant of the jacobian matrix
        real(wp), dimension(4, 8) :: n_bar   
        real(wp) :: L(3, 4), gamma(2, 2), gamma_bar(4, 4) !!


        L = reshape( shape=(/ 3, 4 /), order=(/ 2, 1 /), source= & 
        		(/ 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &         
                   0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, &       
                   0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp  /) )
      
        jac(1,1) = 0.25_wp*((eta-1.0_wp)*xe(1)+(1.0_wp-eta)*xe(3)+(1.0_wp+eta)*xe(5)-(1.0_wp+eta)*xe(7))
        jac(1,2) = 0.25_wp*((eta-1.0_wp)*xe(2)+(1.0_wp-eta)*xe(4)+(1.0_wp+eta)*xe(6)-(1.0_wp+eta)*xe(8))
        jac(2,1) = 0.25_wp*((xi-1.0_wp)*xe(1)-(1.0_wp+xi)*xe(3)+(1.0_wp+xi)*xe(5)+(1.0_wp-xi)*xe(7))
		jac(2,2) = 0.25_wp*((xi-1.0_wp)*xe(2)-(1.0_wp+xi)*xe(4)+(1.0_wp+xi)*xe(6)+(1.0_wp-xi)*xe(8))
		
		detjac = jac(1,1)*jac(2,2)-jac(2,1)*jac(1,2)

        
        
        gamma = reshape( shape=(/ 2, 2 /), order=(/ 2, 1 /), source= & 
        		(/ jac(2,2)/detjac, -jac(1,2)/detjac, &         
                   -jac(2,1)/detjac, jac(1,1)/detjac /) ) 

		gamma_bar = reshape( shape=(/ 4, 4 /), order=(/ 2, 1 /), source= & 
        		(/ gamma(1,1), gamma(1,2), 0.0_wp, 0.0_wp, &         
                   gamma(2,1), gamma(2,2), 0.0_wp, 0.0_wp, &  
                   0.0_wp, 0.0_wp, gamma(1,1), gamma(1,2), &       
                   0.0_wp, 0.0_wp, gamma(2,1), gamma(2,2)  /) )
                   
        n_bar = 0.25_wp*reshape( shape=(/ 4, 8 /), order=(/ 2, 1 /), source= & 
        		(/ eta-1.0_wp,  0.0_wp, 1.0_wp-eta,  0.0_wp, 1.0_wp+eta,  0.0_wp, -1.0_wp-eta,   0.0_wp, &         
                   xi -1.0_wp,  0.0_wp, -1.0_wp-xi,  0.0_wp, 1.0_wp+ xi,  0.0_wp,  1.0_wp- xi,   0.0_wp, &
                    0.0_wp, eta-1.0_wp,  0.0_wp, 1.0_wp-eta,  0.0_wp, 1.0_wp+eta,   0.0_wp, -1.0_wp-eta, & 
                    0.0_wp, xi -1.0_wp,  0.0_wp, -1.0_wp-xi,  0.0_wp, 1.0_wp+ xi,   0.0_wp,  1.0_wp- xi /) )
		n(:,:) = 0.0_wp
        n(1, 1) = (1 - xi)*(1 - eta)
        n(2, 2) = n(1, 1)
        n(1, 3) = (1 + xi)*(1 - eta)
        n(2, 4) = n(1, 3)
        n(1, 5) = (1 + xi)*(1 + eta)
        n(2, 6) = n(1, 5)
        n(1, 7) = (1 - xi)*(1 + eta)
        n(2, 8) = n(1,7)
        n(:,:) = 0.25_wp * n(:,:)
        bmat = matmul(L, matmul(gamma_bar, n_bar))

    end subroutine shape
!
!--------------------------------------------------------------------------------------------------
!
    subroutine plane42_tur_blade(xe, eface, fe, thk, dens, re)

        !! This subroutine computes the element load vector due
        !! to surface traction (traction is always perpendicular
        !! to element face).

        integer, intent(in) :: eface
            !! Element face where traction (pressure) is applied

        real(wp), intent(in) :: fe
            !! Value of surface traction (pressure)
        real(wp), intent(in) :: thk, dens
!$$$$$$             !! nodal velocity vector
!$$$$$$         real(wp), intent(in) :: velo(8)
            !! Thickness of this element
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed configuration (see also [[plane42rect_ke]])
        real(wp), intent(out) :: re(8)
            !! Element force vector
            !!
            !! * `re(1:2)` = \((f_x^1, f_y^1)\) force at element node 1 in \(x\)- and y-direction
            !! * `re(3:4)` = \((f_x^2, f_y^2)\) force at element node 1 in \(x\)- and y-direction
            !! * etc...
        real(wp) :: eta, xi, gaus_co(2), gaus_w
        real(wp) :: jac(2,2), n(2,8), bmat(3,8), detjac
        
        integer :: i

		gaus_w = 1.0_wp
        gaus_co = (/ 1.0_wp/sqrt(3.0_wp), -1.0_wp/sqrt(3.0_wp) /)
        re = 0
        if (eface == 1) then
          do i=1,2
				eta = -1.0_wp
                xi = gaus_co(i)
				
				call shape(xe, xi, eta, n, bmat, jac, detjac)
				
				re = re + gaus_w * thk * fe * matmul( transpose(n), (/ -jac(1,2), jac(1,1) /) )
            end do
            
        elseif (eface == 2) then
        	do i=1,2
				eta = gaus_co(i)
                xi = 1.0_wp
				
        		call shape(xe, xi, eta, n, bmat, jac, detjac)

				re = re + gaus_w * thk * fe * matmul(transpose(n), (/ -jac(2,2), jac(2,1) /) )
            end do
            
        elseif (eface == 3) then
 			do i=1,2
				eta = 1.0_wp
                xi = gaus_co(i)
				
				call shape(xe, xi, eta, n, bmat, jac, detjac)
				
				re = re + gaus_w * thk * fe * matmul(transpose(n), (/ jac(1,2), -jac(1,1) /) )
            end do 
        elseif (eface == 4) then
			do i=1,2
				eta = gaus_co(i)
                xi = -1.0_wp
                				
        		call shape(xe, xi, eta, n, bmat, jac, detjac)
                
				re = re + gaus_w * thk * fe * matmul(transpose(n), (/ jac(2,2), -jac(2,1) /) )
            end do  
        endif
        
    end subroutine plane42_tur_blade
    
end module plane42
