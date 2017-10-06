program main

    !! The "main" routine of the program  
    !!   
    !! This is where execution starts when running the program
    
    use processor
    use fea

    implicit none

    ! Read model data
    call input

    ! Initialize problem
    call initial

    ! Calculate displacements
    call displ
    
    ! Close plot window(s)
    call plot( done )
   
end program main
