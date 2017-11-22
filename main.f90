program main

    !! The "main" routine of the program  
    !!   
    !! This is where execution starts when running the program
    
    use processor
    use fea

    implicit none
	!call stopwatch('star')
    ! Read model data
    call input

    ! Initialize problem
    call initial

    ! Calculate displacements
    if (antype == 'static') then
    	call displ
    elseif (antype == 'modal') then
    	call eigen
    elseif (antype == 'trans') then
    	call trans_loading_ccd
    end if

    ! Close plot window(s)
    call plot( done )
   
end program main
