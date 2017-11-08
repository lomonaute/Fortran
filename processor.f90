module processor

    !! This module contains input and output subroutines  
    !!  
    !! ** Note: It is not necessary to alter in this module **

    use fedata
    implicit none

    public

    ! Plot options    
    integer, parameter :: undeformed = 1
    integer, parameter :: deformed = 2
    integer, parameter :: vectors = 4
    integer, parameter :: elements = 8
    integer, parameter :: eigenmode = 16
    integer, parameter :: done = 32
        
    ! Device numbers
    integer, parameter :: psfile = -1
    integer, parameter :: giffile = -2
    integer, parameter :: matlab = -3
    
    ! Counters for output files
    integer :: ips = 0
    integer :: igif = 0
    integer :: imatlab = 0    

    private :: plotmatlabdef, plotmatlabeval, plotmatlabevec, plotmatlabeig
    
contains

    subroutine input

        !!  This subroutine reads in the input in ansys format.

        integer :: k, e
        integer :: nncheck, necheck
        integer :: net
        integer :: ivalue, ivalue2
        integer :: currentet, currentmp, currentr
        character(len=12), dimension(10) :: ename
        character(len=12) :: command, cvalue, cvalue2
        real(wp) :: rvalue, rvalue2
        logical :: fnexist
        integer, dimension(:,:), allocatable :: eix

        ! Reset counters
        nn = 0
        nncheck = 0
        net = 0
        ne = 0
        necheck = 0
        nm = 0
        nb = 0
        np = 0

        ! Note: From Fortran 2003+ command line arguments are accessible, but not in Fortran90
        write (*, *)
        write (*, '("File(s) in current directory: ")')
        call system('dir') ! Windows
!        call system('ls' ) ! Linux
        write (*, *)
        write (*, '("Enter fem input file (<filename> or enter for previous <filename>)")')
        write (*, '("Filename?")')
        read (*, '(a50)') filename
        if (filename == ' ') then
            inquire (file = '.fem_filename', exist = fnexist)
            if (.not. fnexist) then
                write (*, *)
                write (*, '("ERROR: previous filename does not exist")')
                stop
            end if
            open (10, file = '.fem_filename')
            read (10, '(a50)') filename
            close (10)
        end if
        inquire (file = trim(filename), exist = fnexist)
        if (.not. fnexist) then
            write (*, *)
            write (*, '("ERROR: file ", a, " does not exist")') trim(filename)
            stop
        end if

        open (10, file = '.fem_filename')
        write (10, '(a)') trim(filename)
        close (10)

        write (*, *)
        write (*, '("Reading input file ", a)') trim(filename)

        open (10, file = trim(filename))

        ! Pass 1: Scan through file in order to find out how many 
        ! elements, nodes, etc that are defined
        
        ! Read in /PREP7 data
        do
            read (10, '(a)', end=100) command
            if (command == '/PREP7') exit
            if (command == '/' .or. command == ' ' .or. command == '!' ) cycle
        end do
        do
            read (10, '(a1)', end = 200) command
            if (command == '/' .or. command == ' ') cycle

            backspace (10)
            read (10, *) command
            select case (command)
            case ('FINISH', 'finish')
                exit
            case ('N', 'n')
                nn = nn + 1
                backspace (10)
                read (10, *) command, ivalue
                if (ivalue > nncheck) nncheck = ivalue
            case('ET', 'et')
                net = net + 1
            case('EN', 'en')
                ne = ne + 1
                backspace (10)
                read (10, *) command, ivalue
                if (ivalue > necheck) necheck = ivalue
            case('MP', 'mp')
                backspace (10)
                read (10, *) command, cvalue, ivalue
                if (ivalue > nm) nm = ivalue
            case('D', 'd')
                nb = nb + 1
            case('F', 'f')
                np = np + 1
            case('SFE', 'sfe')
                np = np + 1
            end select
        end do

        if (nn == 0) then
            write (*, *) 'ERROR: No nodes defined'
            stop
        else if (nn /= nncheck) then
            write (*, *) 'ERROR: Node number(s) skipped'
            stop
        else if (net == 0) then
            write (*, *) 'ERROR: No element types defined'
            stop
        else if (ne == 0) then
            write (*, *) 'ERROR: No elements defined'
            stop
        else if (ne /= necheck) then
            write (*, *) 'ERROR: Element number(s) skipped'
            stop
        else if (nm == 0) then
            write (*, *) 'ERROR: No material properties defined'
            stop
        else if (nb == 0) then
            write (*, *) 'ERROR: No supports defined'
            stop
        else if (np == 0) then
            write (*, *) 'WARNING: No loads defined'
        end if
        close (10)

        allocate (x(nn, 2))
        allocate (element(ne))
        allocate (eix(ne,4))
        allocate (mprop(nm))
        allocate (bound(nb, 3), loads(np, 4))

        ! Make thickness equal to one as default
        mprop%thk  = 1 
        
        ! And reset all other mprop parameters
        mprop%young = 0
        mprop%nu = 0
        mprop%dens = 0
        mprop%youngy = 0
        mprop%shear = 0  
              
        do e = 1, ne
            element(e)%ix = 0
        end do

        open (10, file = trim(filename))

        nb = 0
        np = 0
        currentet = 0
        currentmp = 0
        currentr = 0
        accel=0.0

        ! Pass 2: Read model definition
        do
            read (10, '(a1)', end=200) command
            if (command == '/' .or. command == ' ' .or. command == '!' ) cycle

            backspace (10)
            read (10, *) command
            select case (command)
            case ('FINISH', 'finish')
                exit
            case ('N', 'n')
                backspace (10)
                read (10, *) command, ivalue, x(ivalue, 1), x(ivalue, 2), rvalue

            case ('ET', 'et')
                backspace (10)
                read (10, *) command, ivalue, cvalue
                select case (cvalue)
                case ('LINK1',  'link1')
                    ename(ivalue) = 'link1'
                case ('PLANE42', 'plane42')
                    ename(ivalue) = 'plane42'
                case ('PLANE42RECT', 'plane42rect')
                    ename(ivalue) = 'plane42'
                case default
                    write (*, *) 'ERROR: Undefined element: ', trim(cvalue)
                    stop
                end select
                currentet = ivalue
            case ('TYPE', 'type')
                backspace (10)
                read (10, *) command, ivalue
                currentet = ivalue
            case ('EN', 'en')
                if (currentet == 0) then
                    write (*, *) 'ERROR: No previous element type pointer defined'
                    stop
                end if
                if (currentmp == 0) then
                    write (*, *) 'ERROR: No previous material property pointer defined'
                    stop
                end if
                backspace (10)
                select case (ename(currentet))
                case ('link1')
                    read (10, *) command, ivalue, (element(ivalue)%ix(k), k = 1, 2)
                    element(ivalue)%mat = currentmp
                    element(ivalue)%id  = 1
                    element(ivalue)%numnode = 2
                case ('plane42')
                    read (10, *) command, ivalue, (element(ivalue)%ix(k), k = 1, 4)
                    element(ivalue)%mat = currentmp
                    element(ivalue)%id  = 2
                    element(ivalue)%numnode = 4
                case default
                    write(*, *) 'ERROR: Unknown element type'
                    stop
                end select
            case ('MP', 'mp')
                backspace (10)
                read (10, *) command, cvalue, ivalue, rvalue
                select case (cvalue)
                case ('EX', 'ex')
                    mprop(ivalue)%young = rvalue
                case ('PRXY', 'prxy')
                    mprop(ivalue)%nu = rvalue
                case ('DENS', 'dens')
                    mprop(ivalue)%dens = rvalue
                case ('EY', 'ey')
                    mprop(ivalue)%youngy = rvalue
                case ('GXY', 'gxy')
                    mprop(ivalue)%shear = rvalue
                case default
                    write (*, *) 'ERROR: Undefined material property: ', trim(cvalue)
                    stop
                end select
                currentmp = ivalue
            case ('MAT', 'mat')
                backspace (10)
                read (10, *) command, ivalue
                currentmp = ivalue
            case ('R', 'r')
                backspace (10)
                read (10, *) command, ivalue, rvalue
                select case (ename(currentet)) 
                case ('link1')
                    ! define truss area
                    mprop(ivalue)%area = rvalue
                case ('plane42')
                    ! define 4-noded quad element thickness
                    mprop(ivalue)%thk = rvalue
                case default
                    write (*, *) 'ERROR: Undefined real constant (r) card'
                    stop
                end select
                currentr = ivalue
            case ('REAL', 'real')
                ! TODO
                ! This command is only partially implemented -- it should not be used
                ! (it has no effect at the moment). The command is supposed to set
                ! the default value of "R", which is to be used for area or thickness
                ! in materials, that do not have this value set explicitly (via a
                ! call of "R" command before calling "MP" on some card).
                backspace (10)
                read (10, *) command, ivalue
                currentr = ivalue
            case ('D', 'd')
                backspace (10)
                nb = nb + 1
                read (10, *) command, ivalue, cvalue, rvalue
                bound(nb, 1) = ivalue
                select case (cvalue)
                case ('UX', 'ux')
                    bound(nb, 2) = 1
                case ('UY', 'uy')
                    bound(nb, 2) = 2
                case default
                    write (*, *) 'ERROR: Nodal dof not "ux" or "uy"'
                    stop
                end select
                bound(nb, 3) = rvalue
            case ('F', 'f')
                backspace (10)
                np = np + 1
                read (10, *) command, ivalue, cvalue, rvalue
                loads(np, 1) = 1
                loads(np, 2) = ivalue
                select case (cvalue)
                case ('FX', 'fx')
                    loads(np, 3) = 1
                case ('FY', 'fy')
                    loads(np, 3) = 2
                case default
                    write (*, *) 'ERROR: Nodal dof not "fx" or "fy"'
                    stop
                end select
                loads(np, 4) = rvalue
            case ('SFE', 'sfe')
                backspace (10)
                np = np + 1
                read (10, *) command, ivalue, ivalue2, cvalue, cvalue2, rvalue
                loads(np, 1) = 2
                loads(np, 2) = ivalue
                loads(np, 3) = ivalue2
                loads(np, 4) = rvalue
            case ('ACEL', 'acel')
                backspace (10)
                read (10, *) command, rvalue, rvalue2
                accel(1) = rvalue
                accel(2) = rvalue2
            case ('ANTYPE','antype')
            	cycle
            case default
                write (*, *) 'ERROR: Unknown keyword: ', trim(command)
                stop
            end select
        end do

        close (10)

        ! Set default analysis type
        antype = 'STATIC' 

        open (10, file = trim(filename))

        ! Read in /SOLU data (if it exists)
        do
          read (10, '(a)', end=400) command
          if (command == '/SOLU' .or. command == '/solu' ) exit
          if (command == '/' .or. command == ' ') cycle
        end do
        do
          read (10, '(a1)', end=300) command
          if (command == '/' .or. command == ' ') cycle

          backspace (10)
          read (10, *) command
          select case (command)
          case ('FINISH', 'finish')
                exit
          case ('ANTYPE', 'antype')
                backspace (10)
                read (10, *) command, cvalue
                select case (cvalue)
                case ('STATIC', 'static')
                    antype = 'static'
                case ('STATIC_NL', 'static_nl')
                    antype = 'static_nl'
                case ('MODAL', 'modal')
                    antype = 'modal'
                case ('ANGLE', 'angle')
                    antype = 'angle'
                case ('TRANS', 'trans')
                    antype = 'trans'
                case default
                    write (*, *) 'ERROR: Only static, static_nl, modal, angle and trans analyses are implemented'
                    stop
                end select
            end select
        end do

        400 continue

        return

        100 write (*, *)
        write (*, '("ERROR: no /PREP7 in input file")')
        stop

        200 write (*, *)
        write (*, '("ERROR: no /PREP7 FINISH in input file")')
        stop

        300 write (*, *)
        write (*, '("ERROR: no /SOLU FINISH in input file")')
        stop
    end subroutine input
!
!--------------------------------------------------------------------------------------------------
!
    subroutine output

        !!  This subroutine writes results into a text file

        integer :: i, e

        write (*, *)
        write (*, '("Writing fem output file to ", a)') trim(filename)//'.txt'
        open (10, file=trim(filename)//'.txt')

        write (10, '(" Node            Displacements                Applied/reaction forces")')
        write (10, '("number     1-direction     2-direction      1-direction     2-direction")')
        do i = 1, nn
            write (10, '(i6, 1x, f15.9, 1x, f15.9, 2x, f15.9, 1x, f15.9)') i, d(2*i-1), d(2*i), p(2*i-1), p(2*i)
        end do
        write (10, '("                                            ___________     ___________")')
        write (10, '("                              Summation ", f15.9, 1x, f15.9)') &
          sum(p(1:neqn-1:2)), sum(p(2:neqn:2))

        write (10, *)
        write (10, '("Element                    Element strain                                   Element Stress")')
        write (10, '("number      1-direction     2-direction     3-direction      1-direction     2-direction     3-direction")')
        do e = 1, ne
            select case( element(e)%id )
            case( 1 )
                write (10, '(1x, i6, 1x, f15.9, 34x, f15.9)') e, strain(e, 1), stress(e, 1)
            case( 2 )
                write (10, '(1x, i6, 3(1x, f15.9), 1x, 3(1x, f15.9))') e, (strain(e, i), i = 1, 3) ,(stress(e, i), i = 1, 3)
            end select
        end do

        close (10)
    end subroutine output
!
!--------------------------------------------------------------------------------------------------
!    
    subroutine stopwatch(oper)

        !! This subroutine computes elapsed wallclock time
        !!
        !! Timing information is written in terminal window

        character(len=4), intent(in) :: oper
            !! Select which "button" to press on your Stopwatch:
            !!
            !! * 'STAR' or 'star' = reset and start the stopwatch
            !! * 'STOP' or 'stop' = print time spent since last 'star' operation (the stop watch is not reset)
        integer time_array_0(8), time_array_1(8)
        real(wp), save :: start_time, finish_time

        select case (oper)
        case ('STAR', 'star')
            call date_and_time(values=time_array_0)
            start_time = time_array_0 (5) * 3600 + time_array_0 (6) * 60 + time_array_0 (7) + 0.001 * time_array_0 (8)
        case ('STOP', 'stop')
            call date_and_time(values=time_array_1)
            finish_time = time_array_1 (5) * 3600 + time_array_1 (6) * 60 + time_array_1 (7) + 0.001 * time_array_1 (8)
            write (6, '(8x, 1a, 1f16.6)') 'elapsed wall clock time:', finish_time - start_time
        case default
            write (*, '("ERROR: in Processor/stopwatch")')
            stop
        end select
    end subroutine stopwatch
!
!--------------------------------------------------------------------------------------------------
!
    subroutine plot( what,                                &
                    ! Output device, etc
                    device, wait, color,                  &
                    ! Titles, legend, etc
                    title, legend, iter,                  &
                    ! For vector plot
                    p1, p2, ang,                          &
                    ! For element plot 
                    eval, mineval, maxeval,               &
                    ! For mode shapes
                    eigenfreq, eigenvector, tend, tdelta, &
                    ! Plotting range
                    xmin, xmax, ymin, ymax )
       
        !! Plotting routine for displaying result on screen, save as graphic file, or
        !! output as a Matlab-script.  
        !!   
        !! This subroutine has many optional arguments, which may or may-not be required
        !! depending on `what` to plot. It is therefore recommended to
        !! call this subroutine using "keyword-arguments", for example like this
        !!
        !!     call plot( deformed + elements, color=.false., eval=my_results )
        !!
        !! In the call above `color` is not required, but included to select a black-and-white plot.
        !! On the other hand `eval` *is* required because `elements` is specified in the argument
        !! `what` (which describes what kind of plot is to be generated).
        
        integer, intent(in) :: what
            !! Select what to plot, options are
            !!
            !! `undeformed`
            !! : Plot the undeformed configuration
            !!
            !! `deformed`
            !! : Plot the deformed configuration
            !!
            !! `vectors`
            !! : Plot principal directions of stress  
            !!   Optional arguments `p1`, `p2`, and `ang` must be defined
            !!
            !! `elements`
            !! : Color elements according to the values in optional argument `eval`,
            !!   which must of course be defined. 
            !!
            !! `eigenmode`
            !! : Plot an eigenmode of the structure using optional arguments `eigenfreq`
            !!   and `eigenvector`, which must be defined. 
            !!
            !! `done`
            !! : Close all plotting windows -- must be called at the end of a program.
            !!
            !! It is possible to combine plot options in several ways, for example
            !!
            !! * `undeformed + deformed`: plot deformed configuration on top of undeformed configuration
            !! * `vectors + deformed`: plot principal directions on deformed geometry
            !! * `vectors + deformed + undeformed`: plot principal directions on deformed geometry on top of
            !!   undeformed configuration
            !! * `elements + deformed`: plot element values on deformed geometry
            !! * `elements + deformed + undeformed`: plot element values on deformed geometry on top of
            !!   undeformed configuration
            !!
            !! If `deformed` is not specified, then undeformed geometry will be used. It is valid to specify
            !! `undeformed` explicitly without `deformed`.
            !!
            !! If `device` is `matlab` then combining several plots will simply generate multiple Matlab-scripts
            !! (i.e. the plots are not combined into one automatically).
        integer, intent(in), optional :: device
            !! Select output device. Available options are
            !!
            !! * Integer in the range from 1 to 8: Open or select window 1-8 and plot there
            !! * `giffile`: save graphics in a GIF-file (pixel based)
            !! * `psfile`: save graphics in a PS-file (vector graphics)
            !! * `matlab`: save graphics as a Matlab-script
            !!
            !! The file name is automatically generated like this <input_file\>\_<0001\>.<suffix\> there
            !! <input_file\> is the name of the input file, <0001\> is a number automatically incremented
            !! by one for every file written of a certain type. The <suffix\> depends on the type of
            !! file output (gif, ps, or m)
        logical, intent(in), optional :: wait
            !! If `wait` is not defined or defined with the value .true. then execution is paused until 
            !! the user presses "enter" on the keyboard (the terminal window must be the active window).  
            !!
            !! Specifying `wait=.false.` can be useful for custom animations.  
            !!  
            !! The wait-mechanism only applies for device = 1 to 8 (i.e. for plotting on screen).
        logical, intent(in), optional :: color
            !! If `color` is not defined or defined with the value .true. then colors will be used for plotting,
            !! otherwise a black-and-white color palette is used.            
        logical, intent(in), optional :: legend
            !! If `legend` is defined with the value .true. then a colorbar is added to the plot. This
            !! is only useful in combination with plot option `elements`.
        character(len=*), intent(in), optional :: title
            !! If `title` is defined, then a title is added to the plot. Optional argument required for plot option `matlab`
        integer, intent(in), optional :: iter
            !! If `iter` is defined, then a subtitle including the iteration number will be added to the plot.            
        real(wp), intent(in), dimension(:), optional :: eval
            !! Values for coloring of elements. Optional argument required for plot option `elements`.
        real(wp), intent(in), optional :: mineval
            !! If `mineval` is defined, then the value will be used instead of min(eval) for the minimum value of the colorbar.
        real(wp), intent(in), optional :: maxeval
            !! If `maxeval` is defined, then the value will be used instead of max(eval) for the maximum value of the colorbar.
        real(wp), intent(in), optional :: eigenfreq
            !! Eigenfrequency of mode shape. Optional argument required for plot option `eigenmode`.
        real(wp), intent(in), dimension(:), optional :: eigenvector
            !! Eigenvector for mode shape.  Optional argument required for plot option `eigenmode`.
        real(wp), intent(in), optional :: tend
            !! Optional end-time for animation of mode shape. If `tend`
            !!
            !! * is defined, then the eigenmode is animated from time t = 0 to t = `tend`
            !! * is *not* defined, then the eigenmode is animated from time t = 0 to t = 5
        real(wp), intent(in), optional :: tdelta
            !! Optional time increment for animation of mode shape. If `tdelta`
            !!
            !! * is defined and positive, then this value will be used as increments of time between
            !!   frames (the animation is running at a speed depending on hardware). 
            !! * is defined and negative, then the eigenmode will be shown in realtime *scaled* by the
            !!   absolute value of `tdelta`. That is, if `tdelta = -0.5` a frequency of 10 Hz will look
            !!   like only 5 Hz. If `tdelta = -3` then a frequency of 10 Hz will look like 30 Hz.
            !! * is *not* defined, then the eigenmode will be shown in realtime.
        real(wp), intent(in), dimension(:), optional :: p1
            !! Value for principal stress along principal direction 1. Optional argument required for plot option `vectors`
        real(wp), intent(in), dimension(:), optional :: p2
            !! Value for principal stress along principal direction 2. Optional argument required for plot option `vectors`
        real(wp), intent(in), dimension(:), optional :: ang
            !! Direction of principal stress 1. Optional argument required for plot option `vectors`
        real(wp), intent(in), optional :: xmin
            !! Optional lower limit of plotting x-range. Overrides automatic limit if defined.
        real(wp), intent(in), optional :: xmax
            !! Optional upper limit of plotting x-range. Overrides automatic limit if defined.
        real(wp), intent(in), optional :: ymin
            !! Optional lower limit of plotting y-range. Overrides automatic limit if defined.
        real(wp), intent(in), optional :: ymax
            !! Optional upper limit of plotting y-range. Overrides automatic limit if defined.
        
        real, dimension(:), pointer :: pgx, pgy
        real, dimension(4) :: xbox, ybox
        real, dimension(5), target :: px, py, pxd, pyd
        real :: clx, cly, mine, maxe, h, w, wsize, aspect, x0, y0, wemax, wxmin, wxmax, wymin, wymax
        real(wp) :: avepx, avepy, ca, sa, sf, time, tstop
        integer(selected_int_kind(15)) :: c0, c1, cr
        integer :: text_thk, vec_thk, edge_thk, truss_thk
        integer :: cid, ncid, ncidmin, ncidmax, dev, cu, cd, cf, ce, i, j, e, nen, p        
        character(len=1) :: ch
        character(len=14) :: strng
        character(len=100) :: varp        
        logical :: col, leg, waituser, init_page
                 
        ! PGPLOT does not provide a module and therefore relies on implicit
        ! interfaces generated by the compiler (the compiler uses type information
        ! from each argument in a given situation). This is OK for subroutines
        ! because they do not return a value (however, the compiler cannot check
        ! for errors in the arguments of a given subroutine -- a mistake will only
        ! show a runtime). For functions we need to tell the compiler what the
        ! return type is. This can be done in 3 ways in Fortran
        ! 
        ! 1)
        ! integer :: pgopen ! The shortest way. Now pgopen can be a function or a variable
        ! 
        ! 2)
        ! integer, external :: pgopen ! A little more detailed. Now pgopen must be a function
        ! 
        ! 3) 
        ! We can define the interface explicitly. Now the compiler also knows what the 
        ! arguments should look like. Ideally, we should provide interfaces for all PGPLOT
        ! functions and subroutines. However, we only need 1 PGPLOT function, so we only
        ! define 1 interface and let the compiler generate implicit interfaces for the
        ! subroutines. Here is what the interface looks like
        interface
            function pgopen( c ) result( k )
                character(len=*), intent(in) :: c
                integer :: k
            end function
        end interface        
        
        ! Save pgplot window identifiers
        integer, save :: panel(8) = 0

	! Set various line thicknesses
	text_thk = 2
        vec_thk = nint(6*scale_thk)
        edge_thk = nint(6*scale_thk)
        truss_thk = nint(12*scale_thk)
        
        ! Is plot request valid ?
        if (iand(what, done+undeformed+deformed+elements+vectors+eigenmode) == 0) then
            write (*, *) 'ERROR: Nothing has been selected for plotting'
            stop
        end if
        
        ! Use default device ?
        if (present(device)) then
            dev = device
        else
            dev = 1
        end if
        
        ! Select plotting device
        init_page = .true.
        select case (dev)        
        case (psfile)
            ips = ips + 1
            write( varp, '(a,"_",i0.4,".ps/cps")' ) trim(filename), ips
            p = pgopen(varp)            
        case (giffile)
            igif = igif + 1
            write( varp, '(a,"_",i0.4,".gif/gif")' ) trim(filename), igif
            p = pgopen(varp)            
        case (matlab)
            if (iand(what, undeformed) /= 0) then
                write (*, *) 'WARNING: Matlab plot of undeformed structure not implemented'
            end if                
            if (iand(what, deformed) /= 0) then
                if (.not.present(title)) then
                    write (*, *) 'ERROR: Cannot plot to Matlab because title is not defined'
                    stop
                end if                
                call plotmatlabdef(title)
            end if
            if (iand(what, elements) /= 0) then
                if (.not.present(title)) then
                    write (*, *) 'ERROR: Cannot plot to Matlab because title is not defined'
                    stop
                end if
                if (.not.present(eval)) then
                    write (*, *) 'ERROR: Cannot plot to Matlab because eval is not defined'
                    stop
                end if                
                call plotmatlabeval(title, eval)
            end if
            if (iand(what, vectors) /= 0) then
                if (.not.present(title)) then
                    write (*, *) 'ERROR: Cannot plot to Matlab because title is not defined'
                    stop
                end if                            
                if (.not.present(p1)) then
                    write (*, *) 'ERROR: Cannot plot to Matlab because p1 is not defined'
                    stop
                end if                
                if (.not.present(p2)) then
                    write (*, *) 'ERROR: Cannot plot to Matlab because p2 is not defined'
                    stop
                end if                
                if (.not.present(ang)) then
                    write (*, *) 'ERROR: Cannot plot to Matlab because ang is not defined'
                    stop
                end if                
                call plotmatlabevec(title, p1, p2, ang)
            end if
            if (iand(what, eigenmode) /= 0) then
                if (.not.present(title)) then
                    write (*, *) 'ERROR: Cannot plot to Matlab because title is not defined'
                    stop
                end if                            
                if (.not.present(eigenfreq)) then
                    write (*, *) 'ERROR: Cannot plot to Matlab because eigenfreq is not defined'
                    stop
                end if                
                if (.not.present(eigenvector)) then
                    write (*, *) 'ERROR: Cannot plot to Matlab because eigenvector is not defined'
                    stop
                end if                
                if (.not.present(tdelta)) then
                    write (*, *) 'ERROR: Cannot plot to Matlab because tdelta is not defined'
                    stop
                end if 
                if (.not.present(tend)) then
                    write (*, *) 'ERROR: Cannot plot to Matlab because tend is not defined'
                    stop
                end if 
                call plotmatlabeig(title, eigenfreq, eigenvector, (/ tdelta, tend /))
            end if
            return
        case ( 1:8 )
            if ( panel(dev) == 0 ) then
                ! Open new pgplot window
                p = pgopen(plotdevice)
                panel(dev) = p                
            else
                ! Select pgplot window
                p = panel(dev)                
                call pgslct(p)
                init_page = .false.
            end if
        case default
            write (*, *) "ERROR: Unknown output device"
            stop
        end select

        ! Do we have a valid pgplot device ?
        if ( .not. p > 0) then
            write (*, *) 'ERROR: pgplot pgopen failed'
            stop            
        end if                

        ! Close window ?
        if( iand(what, done) > 0 .and. dev > 0 ) then
            panel(dev) = 0
            call pgclos
            return
        end if
        
        ! Initialize page for gif, ps, or new screen device ?
        if( init_page ) then
            ! Set size of window ("page" size in inches, aspect ratio)
            wsize = 8.0
            aspect = 0.7                
            call pgpap(wsize, aspect)
            
            ! We don't want PGPLOT asking whether to close a window or not
            call pgask( 0 )     

            ! Set background color to white, text to black
            call pgscr (0, 1.0, 1.0, 1.0)
            call pgscr (1, 0.0, 0.0, 0.0)
            
            ! Start a fresh page
            call pgpage       
        end if        
            
        ! Find size of biggest element
        clx = -huge(1.0)
        cly = -huge(1.0)
        do e = 1, ne
            nen = element(e)%numnode
            do i = 1, nen
                do j = i+1, nen
                   clx = max(real(abs(x(element(e)%ix(i),1)-x(element(e)%ix(j),1))), clx)
                   cly = max(real(abs(x(element(e)%ix(i),2)-x(element(e)%ix(j),2))), cly)
                end do
            end do
        end do
        
        ! Find plotting range ...
        wxmin = huge(1.0)
        wxmax = -huge(1.0)
        wymin = huge(1.0)
        wymax = -huge(1.0)
        do e = 1, ne
            nen = element(e)%numnode
            px(1:nen) = real(x(element(e)%ix(1:nen),1))
            py(1:nen) = real(x(element(e)%ix(1:nen),2))

            ! Get size of undeformed structure if 'undeformed' is requested, or if 'deformed' is not requested
            if (iand(what, undeformed) /= 0 .or. iand(what, deformed) == 0) then
                wxmin = min( minval(px(1:nen)), wxmin )
                wxmax = max( maxval(px(1:nen)), wxmax )
                wymin = min( minval(py(1:nen)), wymin )
                wymax = max( maxval(py(1:nen)), wymax )
            end if
            
            ! Get size of deformed structure if 'deformed' is requested and 'eigenmode' is not requested
            if (iand(what, deformed) /= 0 .and. iand(what, eigenmode) == 0) then
                pxd(1:nen) = real(px(1:nen) + scale_def * d(2*element(e)%ix(1:nen)-1))
                pyd(1:nen) = real(py(1:nen) + scale_def * d(2*element(e)%ix(1:nen)))
                wxmin = min( minval(pxd(1:nen)), wxmin )
                wxmax = max( maxval(pxd(1:nen)), wxmax )
                wymin = min( minval(pyd(1:nen)), wymin )
                wymax = max( maxval(pyd(1:nen)), wymax )
            end if
            
            ! Get size with eigenmode at +/- max amplitude
            if (iand(what, eigenmode) /= 0) then
                pxd(1:nen) = real(px(1:nen) + scale_def * eigenvector(2*element(e)%ix(1:nen)-1))
                pyd(1:nen) = real(py(1:nen) + scale_def * eigenvector(2*element(e)%ix(1:nen)))
                wxmin = min( minval(pxd(1:nen)), wxmin )
                wxmax = max( maxval(pxd(1:nen)), wxmax )
                wymin = min( minval(pyd(1:nen)), wymin )
                wymax = max( maxval(pyd(1:nen)), wymax )
                
                pxd(1:nen) = real(px(1:nen) - scale_def * eigenvector(2*element(e)%ix(1:nen)-1))
                pyd(1:nen) = real(py(1:nen) - scale_def * eigenvector(2*element(e)%ix(1:nen)))
                wxmin = min( minval(pxd(1:nen)), wxmin )
                wxmax = max( maxval(pxd(1:nen)), wxmax )
                wymin = min( minval(pyd(1:nen)), wymin )
                wymax = max( maxval(pyd(1:nen)), wymax )
            end if
        end do
        
        ! Add a little bit of space
        w = wxmax-wxmin
        h = wymax-wymin
        wxmin = wxmin - 0.02*w
        wxmax = wxmax + 0.02*w
        wymin = wymin - 0.02*h
        wymax = wymax + 0.02*h

        ! Apply user-defined plotting range ?
        if (present(xmin)) wxmin = real(xmin)
        if (present(xmax)) wxmax = real(xmax)
        if (present(ymin)) wymin = real(ymin)
        if (present(ymax)) wymax = real(ymax)
        
        ! Save max dimension (x- or y-direction) of elements in the structure
        wemax = max( clx, cly )
        
        ! Handle optional arguments
        if (present(color)) then
            col = color
        else
            col = .true.
        end if        
        
        if (present(legend)) then
            leg = legend
        else
            leg = .false.
        end if        
        
        if (present(wait)) then
            waituser = wait
        else
            waituser = .true.
        end if
        
        if (iand(what, vectors) /= 0) then
            if (.not.present(p1)) then
                write (*, *) 'ERROR: Cannot plot vectors because p1 is not defined'
                stop
            end if
            if (.not.present(p2)) then
                write (*, *) 'ERROR: Cannot plot vectors because p2 is not defined'
                stop
            end if
            if (.not.present(ang)) then
                write (*, *) 'ERROR: Cannot plot vectors because ang is not defined'
                stop
            end if
        end if
        
        if (iand(what, elements) /= 0) then
            if (.not.present(eval)) then
                write (*, *) 'ERROR: Cannot plot elements because eval is not defined'
                stop
            end if
            if (.not.present(mineval)) then
                mine = real(minval(eval))
            else
                mine = real(mineval)
            end if
            if (.not.present(maxeval)) then
                maxe = real(maxval(eval))
            else
                maxe = real(maxeval)
            end if
        end if
        
        if (iand(what, eigenmode) /= 0) then
            if (.not.present(eigenfreq)) then
                write (*, *) 'ERROR: Cannot plot mode shape because eigenfreq is not defined'
                stop
            end if
            if (.not.present(eigenvector)) then
                write (*, *) 'ERROR: Cannot plot mode shape because eigenvector is not defined'
                stop
            end if
        end if
        
        if (leg .and. .not.present(eval)) then
            write (*, *) 'ERROR: Cannot create colorbar because eval is not defined'
            stop
        end if   
                    
        ! Use color or black-and-white ?
        if (col) then
            ! Color palette
            call pgscr(16, 0.0, 0.0, 1.0)
            call pgscr(17, 0.0, 0.2, 0.8)
            call pgscr(18, 0.0, 0.4, 0.6)
            call pgscr(19, 0.0, 0.6, 0.4)
            call pgscr(20, 0.0, 0.8, 0.2)
            call pgscr(21, 0.0, 1.0, 0.0)
            call pgscr(22, 0.0, 1.0, 0.0)
            call pgscr(23, 0.5, 1.0, 0.0)
            call pgscr(24, 1.0, 1.0, 0.0)
            call pgscr(25, 1.0, 0.6, 0.0)
            call pgscr(26, 1.0, 0.3, 0.0)
            call pgscr(27, 1.0, 0.0, 0.0)
            ncidmin = 16
            ncidmax = 27
            if(iand(what, elements ) > 0 ) then                                
                if(iand(what, deformed) > 0 .and. iand(what, undeformed) > 0 ) then
                    cd = 14 ! Dark Gray for deformed structure                
                    cu = 15 ! Light Gray for undeformed structure
                else
                    cd = 0 ! White (background) for deformed structure                
                    cu = 0 ! White (background) for undeformed structure                
                end if
            else                
                cd = 2 ! Red for deformed structure
                cu = 4 ! Blue for undeformed structure
            end if
            cf = 7 ! Yellow for fill
            ce = 12 ! Blue + Magenta for edges
        else
            ! Gray scale palette
            call pgscr (16, 0.9, 0.9, 0.9)
            call pgscr (17, 0.8, 0.8, 0.8)
            call pgscr (18, 0.7, 0.7, 0.7)
            call pgscr (19, 0.6, 0.6, 0.6)
            call pgscr (20, 0.5, 0.5, 0.5)
            call pgscr (21, 0.4, 0.4, 0.4)
            call pgscr (22, 0.3, 0.3, 0.3)
            call pgscr (23, 0.2, 0.2, 0.2)
            call pgscr (24, 0.1, 0.1, 0.1)
            ncidmin = 16
            ncidmax = 24
            if(iand(what, elements ) > 0 ) then
                if(iand(what, deformed) > 0 .and. iand(what, undeformed) > 0 ) then
                    cd = 14 ! Dark Gray for deformed structure
                    cu = 15 ! Light Gray for undeformed structure
                else
                    cd = 0 ! White (background) for deformed structure
                    cu = 0 ! Light Gray for undeformed structure
                end if
            else
                cd = 14 ! Dark Gray for deformed structure
                cu = 15 ! Light Gray for undeformed structure
            end if            
            cf = 15 ! Light Gray for fill
            ce = 1 ! Black (foreground) for edges
        end if
        ncid = ncidmax-ncidmin+1

        ! Beginning buffering
        call pgbbuf

        ! Clear entire window
        call pgeras
        
        ! Plot title ?
        if (present(title)) then
            call pgsvp(0.0, 0.8, 0.9, 1.0)
            call pgswin(0.0, 1.0, 0.0, 1.0)
            call pgsci(1)
            call pgsch(2.0)
            call pgslw(text_thk)
            call pgptxt(0.5, 0.0, 0.0, 0.5, title)
        end if

        ! Plot iteration number ?
        if (present(iter)) then
            call pgsvp(0.0, 0.8, 0.85, 0.9)
            call pgswin(0.0, 1.0, 0.0, 1.0)
            call pgsci(1)
            call pgsch(1.5)
            call pgslw(text_thk)
            write (strng, '(i0)') iter
            call pgptxt(0.5, 0.0, 0.0, 0.5, 'Iteration #'//trim(strng))
        end if
                
        ! Plot legend ?
        if (leg) then
            call pgsvp(0.7, 1.0, 0.1, 0.8)
            call pgswin(0.0, 1.0, 0.0, 1.0)
            call pgsch(1.5)
            call pgslw(text_thk)
            call pgqtxt(0.0, 0.0, 0.0, 1.0, "12345678901234", xbox, ybox)
            h = 1.0/ncid
            w = 1.0/ncid * (0.8-0.1)/(1-0.7) * 0.8
            do i = 1, ncid
                px(1) = 1.0 - w
                py(1) = h*(i-1)
                px(2) = 1.0
                py(2) = py(1)
                px(3) = px(2)
                py(3) = h*i
                px(4) = px(1)
                py(4) = py(3)                
                call pgsci(ncidmin-1+i)
                call pgsfs(1)
                call pgpoly(4,px,py)
                call pgsci(0)
                call pgsfs(2)
                call pgpoly(4,px,py)
                call pgsci(1)
                write (strng, '(2x, es12.3)') mine+(maxe-mine)/ncid*(i-1)
                call pgptxt(px(1)-w/2 - xbox(4), h*(i-1) - 0.5*(ybox(3)+ybox(4)), 0.0, 1.0, strng)                
            end do
            write (strng, '(2x, es12.3)') maxe
            call pgptxt(px(1)-w/2 - xbox(4), h*(i-1) - 0.5*(ybox(3)+ybox(4)), 0.0, 1.0, strng )
        end if
        
        ! Set viewport and plotting range, making sure that aspect ratio 
        ! is 1:1 when plotting the structure
        call pgsvp(0.0, 0.7, 0.05, 0.75)
        call pgwnad( wxmin, wxmax, wymin, wymax )        

        ! Draw a rectangle covering the entire viewport, using background color.
        ! Without this operation, for some unknown reason, bad clipping can
        ! occur making some line look thinner than they should (Windows/ClearWin+)
        call pgsfs(1)
        call pgsci(0)
        call pgrect( wxmin, wxmax, wymin, wxmax )
        
        ! Draw whatever is requested...
        do e = 1, ne
            nen = element(e)%numnode

            ! Get data for undeformed geometry
            px(1:nen) = real(x(element(e)%ix(1:nen),1))
            py(1:nen) = real(x(element(e)%ix(1:nen),2))

            ! Get data for deformed geometry
            pxd(1:nen) = px(1:nen) + real(scale_def * d(2*element(e)%ix(1:nen)-1))
            pyd(1:nen) = py(1:nen) + real(scale_def * d(2*element(e)%ix(1:nen)))

            ! Use deformed geometry for "element values" and "vectors" ?
            if (iand(what, deformed) /= 0) then
                pgx => pxd
                pgy => pyd                    
            else
                pgx => px
                pgy => py
            end if
                    
            ! Plot element values ?
            if (iand(what, elements) /= 0) then
                ! Select color index
                if (abs(maxe-mine) < 1e-10) then
                    cid = ncidmin
                else
                    cid = ncidmin + int(ncid * (eval(e)-mine)/(maxe-mine))                
                    if (cid < ncidmin) then
                        cid = ncidmin
                    else if (cid > ncidmax) then
                        cid = ncidmax
                    end if
                end if
                call pgsci(cid)

                ! Draw element                
                select case (element(e)%id)
                case (1)
                    call pgslw(truss_thk)
                    call pgline(nen,pgx,pgy)
                case (2)
                    call pgslw(edge_thk)
                    call pgsfs(1)
                    call pgpoly(nen,pgx,pgy)
                    call pgsci(0)
                    call pgsfs(2)
                    call pgpoly(nen,pgx,pgy)
                case default
                    write (*, *) 'ERROR: Cannot plot unknown element type'
                    stop
                end select
            end if
            
            ! Plot vectors ?
            if (iand(what, vectors) /= 0) then
                ! Thick lines
                call pgslw(vec_thk)
                
                ! Compute scaling of principal directions
                 sf = sqrt(2.0_wp) / wemax / scale_vec

                ! Find center of element
                avepx = sum(pgx(1:nen))/nen
                avepy = sum(pgy(1:nen))/nen

                ! Plot principal direction 1
                ca = cos(-ang(e))
                sa = sin(-ang(e))
                xbox(1) = real(avepx - 0.5 * abs(p1(e)) * ca/sf)
                ybox(1) = real(avepy - 0.5 * abs(p1(e)) * sa/sf)
                xbox(2) = real(avepx + 0.5 * abs(p1(e)) * ca/sf)
                ybox(2) = real(avepy + 0.5 * abs(p1(e)) * sa/sf)
                if (p1(e) >= 0) then
                    call pgsci(cu)
                else
                    call pgsci(cd)
                end if
                call pgline(2,xbox,ybox)

                ! Plot principal direction 2
                ca = cos(-ang(e)+pi/2)
                sa = sin(-ang(e)+pi/2)
                xbox(1) = real(avepx - 0.5 * abs(p2(e)) * ca/sf)
                ybox(1) = real(avepy - 0.5 * abs(p2(e)) * sa/sf)
                xbox(2) = real(avepx + 0.5 * abs(p2(e)) * ca/sf)
                ybox(2) = real(avepy + 0.5 * abs(p2(e)) * sa/sf)
                if (p2(e) >= 0) then
                    call pgsci(cu)
                else
                    call pgsci(cd)
                end if
                call pgline(2,xbox,ybox)
            end if
            
            ! Plot undeformed structure ?
            if (iand(what, undeformed) /= 0) then
                call pgsci(cu)
                call pgslw(edge_thk)
                select case (element(e)%id)
                case (1)
                    call pgline(nen,px,py)
                case (2)
                    call pgsfs(2)
                    call pgpoly(nen,px,py)
                case default
                    write (*, *) 'ERROR: Cannot plot unknown element type'
                    stop
                end select
            end if
            
            ! Plot deformed structure ?
            if (iand(what, deformed) /= 0) then
                call pgsci(cd)
                call pgslw(edge_thk)
                select case (element(e)%id)
                case (1)
                    if(iand(what, elements) == 0) then
	                call pgline(nen,pxd,pyd)
                    end if
                case (2)
                    call pgsfs(2)
                    call pgpoly(nen,pxd,pyd)
                case default
                    write (*, *) 'ERROR: Cannot plot unknown element type'
                    stop
                end select
            end if
        end do
        
        ! Write user data in plot window
        call pgiden
        
        ! Flush render buffer
        call pgebuf
                
        ! Plot eigenmode ?
        if (iand(what, eigenmode) /= 0) then
            if (present(tend)) then
                if (tend <= 0) then
                    write (*, *) 'ERROR: tend must be positive'
                    stop
                else
                    tstop = tend
                end if
            else
                tstop = 5
            end if
            call system_clock( c0, cr )
            do 
                if( .not. present(tdelta) ) then
                    call system_clock( c1 )
                    time = ( c1-c0 )/real(cr, wp) 
                else if( tdelta > 0 ) then
                    time = time + tdelta
                else if( tdelta < 0 ) then
                    call system_clock( c1 )
                    time = ( c1-c0 )/real(cr, wp)*abs(tdelta)
                else
                    write (*, *) 'ERROR: tdelta must not be zero'
                    stop                    
                end if
                
                ! Clear canvas (using pgeras is slow comparing to just filling a rect)
                call pgbbuf
                call pgsci(0)
                call pgsfs(1)
                call pgrect(wxmin, wxmax, wymin, wymax)
                
                ! Compute scaling factor for eigenvector
                sf = sin(eigenfreq*time)
                
                ! Draw structure (repeating some of the code above...)
                do e = 1, ne
                    nen = element(e)%numnode
                    pxd(1:nen) = real(x(element(e)%ix(1:nen),1) + sf * scale_def * eigenvector(2*element(e)%ix(1:nen)-1))
                    pyd(1:nen) = real(x(element(e)%ix(1:nen),2) + sf * scale_def * eigenvector(2*element(e)%ix(1:nen)))
                    select case (element(e)%id)
                    case (1)
                        call pgslw(truss_thk)
                        call pgsci(ce)
                        call pgline(nen,pxd,pyd)
                    case (2)                        
                        call pgsci(cf)
                        call pgsfs(1)
                        call pgpoly(nen,pxd,pyd)
                        call pgslw(edge_thk)
                        call pgsci(ce)
                        call pgsfs(2)
                        call pgpoly(nen,pxd,pyd)
                    case default
                        write (*, *) 'ERROR: Unknown element type in plotting mode shape'
                        stop
                    end select
                end do
                call pgebuf
                
                if (time > tstop) exit
            end do
        end if
        
        ! Wait for user ?
        select case (dev)
        case (1:8)
            if( waituser ) then
                call pgcurs( x0, y0, ch )
            end if
        case default
            call pgclos 
        end select                
    end subroutine plot
!
!--------------------------------------------------------------------------------------------------
!   
    subroutine plotmatlabdef(title)

        ! Subroutine to plot the un/deformed structure using Matlab

        character(len=*), intent(in) :: title
        character(len=100) :: fdata, fscript
        integer i, j, e

        ! write datafile
        imatlab = imatlab + 1
        write (fdata, '(a,"_plotdeformed_data",i0.4,".m")') trim(filename), imatlab
        open (13, file=trim(fdata))

        ! write nodal coordinates
        write(13,'("X = [")')
        do i = 1,size(x,1)
            write (13,'(3(f15.9,1x))') (x(i,j), j=1,size(x,2))
        end do
        write(13,'("];")')
        write(13,'( )')

        ! write topology matrix
        write(13,'("IX = [")')
        do e = 1, ne
            write (13,'(20(i6,1x))') (element(e)%ix(i) ,i=1,size(element(e)%ix,1))
        end do
        write(13,'("];")')

        ! write deformation vector
        write(13,'("D = [")')
        do i = 1, neqn
            write (13,'(f15.9,1x)') (d(i) )
        end do
        write(13,'("];")')
        close(13)

        ! Create matlab script
        write(fscript, '(a,"_plotdeformed_",i0.4,".m")' ) trim(filename), imatlab
        open (13, file = trim(fscript))
        write(13,*) '% Plotting Un-Deformed and Deformed Structure'
        write(13,*) 'close all'
        write(13,*) 'clear all'
        write(13,*) fdata(1:len_trim(fdata)-2) // ';'
        write(13,*) '% Make plot'
        write(13,*) 'figure'
        write(13,*) 'set(gcf,',"'",'Name',"','", trim(title) ,"'",')'

        ! Element dependent code
        ! NOTE: not possible to mix element types !!!
        if (element(1)%id == 1) then
            write(13,*) 'subplot(2,1,2)'
            write(13,*) 'hold on'
            write(13,*) 'for e = 1:size(IX,1)'
            write(13,*) '   edof = [2*IX(e,1)-1 2*IX(e,1) 2*IX(e,2)-1 2*IX(e,2)];'
            write(13,*) '   xx = X(IX(e,1:2),1) + D(edof(1:2:4));'
            write(13,*) '   yy = X(IX(e,1:2),2) + D(edof(2:2:4));'
            write(13,*) '   plot(xx,yy,',"'",'b',"',","'",'LineWidth',"'",',1.5)'
            write(13,*) 'end'
            write(13,*) 'title(',"'",'Deformed',"'",')'
            write(13,*) 'axis equal'
            write(13,*) 'xaxes = get(gca,',"'",'xlim',"'",');'
            write(13,*) 'yaxes = get(gca,',"'",'ylim',"'",');'
            write(13,*) 'axis off;'
            write(13,*) 'subplot(2,1,1)'
            write(13,*) 'hold on'
            write(13,*) 'for e = 1:size(IX,1)'
            write(13,*) '   xx = X(IX(e,1:2),1);'
            write(13,*) '   yy = X(IX(e,1:2),2);'
            write(13,*) '   plot(xx,yy,',"'",'b',"',","'",'LineWidth',"',",'1.5)'
            write(13,*) 'end'
            write(13,*) 'title(',"'",'Undeformed',"'",')'
            write(13,*) 'axis([min(xaxes(1),min(X(:,1))) max(xaxes(2),max(X(:,1)))...'
            write(13,*) ' min(yaxes(1),min(X(:,2))) max(yaxes(2),max(X(:,2))) ]);'
            write(13,*) 'axis off;'
        else if (element(1)%id == 2) then
            write(13,*) 'subplot(2,1,2)'
            write(13,*) 'hold on'
            write(13,*) 'for e = 1:size(IX,1)'
            write(13,*) '    edof = [2*IX(e,1)-1 2*IX(e,1) 2*IX(e,2)-1 2*IX(e,2)...'
            write(13,*) '        2*IX(e,3)-1 2*IX(e,3) 2*IX(e,4)-1 2*IX(e,4)];'
            write(13,*) '    xx = X(IX(e,1:4),1) + D(edof(1:2:8));'
            write(13,*) '    yy = X(IX(e,1:4),2) + D(edof(2:2:8));'
            write(13,*) '    patch(xx,yy,[1 1 0]);'
            write(13,*) 'end'
            write(13,*) 'title(',"'",'Deformed',"'",')'
            write(13,*) 'axis equal;'
            write(13,*) 'xaxes = get(gca,',"'",'xlim',"'",');'
            write(13,*) 'yaxes = get(gca,',"'",'ylim',"'",');'
            write(13,*) 'axis off;'
            write(13,*) 'subplot(2,1,1)'
            write(13,*) 'hold on'
            write(13,*) 'for e = 1:size(IX,1)'
            write(13,*) '    xx = X(IX(e,1:4),1);'
            write(13,*) '    yy = X(IX(e,1:4),2);'
            write(13,*) '    patch(xx,yy,[1 1 0]);'
            write(13,*) 'end'
            write(13,*) 'title(',"'",'Undeformed',"'",')'
            write(13,*) 'axis([min(xaxes(1),min(X(:,1))) max(xaxes(2),max(X(:,1)))...'
            write(13,*) ' min(yaxes(1),min(X(:,2))) max(yaxes(2),max(X(:,2))) ]);'
            write(13,*) 'axis equal;'
            write(13,*) 'axis off;'
        else
            write (*,'("Unsupported element type in Matlab routine:")')
            write (*,'("plotmatlab -> plotmatlabdef")')
        end if
        
        ! End file and close
        write(13,*) 'hold off'
        write(13,*) 'set(gcf,',"'",'color',"'",',[ 1  1 1]);'
        close(13)

        ! Print to screen:
        print*,'Matlab plotting called: un/deformed'
        print*,'Generated data file: ', trim(fdata)
        print*,'and script file: ', trim(fscript)
        print*,' '
    end subroutine plotmatlabdef
!
!--------------------------------------------------------------------------------------------------
!
    subroutine plotmatlabeval(title,data1)

        ! Subroutine to plot the element values using Matlab

        character(len=*), intent(in) :: title
        real(wp), dimension(:), intent(in) :: data1
        character(len=100) :: fdata, fscript
        integer i, j, e

        ! write datafile
        imatlab = imatlab + 1
        write (fdata, '(a,"_plotelements_data_",i0.4,".m")') trim(filename), imatlab
        open (13, file=trim(fdata))

        ! write nodal coordinates
        write(13,'("X = [")')
        do i = 1, size(x,1)
            write (13,'(3(f15.9,1x))') (x(i,j), j=1,size(x,2))
        end do
        write(13,'("];")')
        write(13,'( )')

        ! write topology matrix
        write(13,'("IX = [")')
        do e = 1, ne
            write (13,'(20(i6,1x))') (element(e)%ix(i) ,i=1,size(element(e)%ix,1))
        end do
        write(13,'("];")')

        ! write data1/element values
        write(13,'("plotval = [")')
        do i = 1, ne
            write (13,'(f32.15)') data1(i)
        end do
        write(13,'("];")')
        close(13)

        ! Create matlab script
        write(fscript, '(a,"_plotelements_",i0.4,".m")' ) trim(filename), imatlab
        open (13, file = trim(fscript))
        write(13,*) '% Plotting Element Values'
        write(13,*) 'close all'
        write(13,*) 'clear all'
        write(13,*) fdata(1:len_trim(fdata)-2) // ';'
        write(13,*) '% Determine colorscale'
        write(13,*) 'cmap = jet;'
        write(13,*) 'cinterp = linspace(min(plotval),max(plotval),size(cmap,1));'
        write(13,*) '% Make plot'
        write(13,*) 'title(',"'", trim(title) ,"'",')'
        write(13,*) 'hold on'
        ! Element dependent code
        ! NOTE: not possible to mix element types !!!
        if (element(1)%id == 1) then
            write(13,*) 'for e = 1:size(IX,1)'
            write(13,*) '   [dummy,arr_pos] = min(abs(cinterp-plotval(e)));'
            write(13,*) '   xx = X(IX(e,1:2),1);'
            write(13,*) '   yy = X(IX(e,1:2),2);'
            write(13,*) '   plot(xx,yy,',"'",'Color',"'",',cmap(arr_pos,:),',"'",'Linewidth',"'",',1.5);'
            write(13,*) 'end'
        else if (element(1)%id == 2) then
            write(13,*) 'for e = 1:size(IX,1)'
            write(13,*) '   [dummy,arr_pos] = min(abs(cinterp-plotval(e)));'
            write(13,*) '   xx = X(IX(e,1:4),1);'
            write(13,*) '   yy = X(IX(e,1:4),2);'
            write(13,*) '   patch(xx,yy,cmap(arr_pos,:));'
            write(13,*) 'end'
        else
            write (*,'("Unsupported element type in Matlab routine:")')
            write (*,'("plot -> plotmatlabelements")')
        end if
        
        ! End file and close
        write(13,*) 'axis equal;'
        write(13,*) 'axis off;'
        write(13,*) 'caxis([min(plotval) max(plotval)]);'
        write(13,*) 'colorbar;'
        write(13,*) 'hold off'
        write(13,*) 'set(gcf,',"'",'color',"'",',[ 1  1 1]);'
        close(13)

        ! Print to screen:
        print*,'Matlab plotting called: elements'
        print*,'Generated data file: ', trim(fdata)
        print*,'and script file: ', trim(fscript)
        print*,' '
    end subroutine plotmatlabeval
!
!--------------------------------------------------------------------------------------------------
!
    subroutine plotmatlabevec(title,ppp1,ppp2,pppang)

        ! Vector field plot using Matlab

        real(wp), dimension(:), intent(in) :: ppp1, ppp2, pppang
        character(len=*), intent(in) :: title
        character(len=100) :: fdata, fscript
        integer i, j, e

        ! write datafile
        imatlab = imatlab + 1
        write (fdata, '(a,"_plotvector_data",i0.4,".m")') trim(filename), imatlab
        open (13, file=trim(fdata))

        ! write nodal coordinates
        write(13,'("X = [")')
        do i = 1, size(x,1)
            write (13,'(3(f15.9,1x))') (x(i,j), j=1,size(x,2))
        end do
        write(13,'("];")')
        write(13,'( )')

        ! write topology matrix
        write(13,'("IX = [")')
        do e = 1, ne
            write (13,'(20(i6,1x))') (element(e)%ix(i) ,i=1,size(element(e)%ix,1))
        end do
        write(13,'("];")')

        ! write data1/element values
        write(13,'("vect = [")')
        do i = 1, ne
            write(13,*) ppp1(i), ppp2(i), pppang(i)
        end do
        write(13,'("];")')
        close(13)

        ! Create matlab script
        write(fscript, '(a,"_plotvector_",i0.4,".m")' ) trim(filename), imatlab
        open (13, file = trim(fscript))
        write(13,*) '% Plotting Vector Field, i.e. principle stresses'
        write(13,*) 'close all'
        write(13,*) 'clear all'
        write(13,*) fdata(1:len_trim(fdata)-2) // ';'
        write(13,*) '% Make plot'
        ! Element dependent code
        ! NOTE: not possible to mix element types !!!
        if (element(1)%id == 1) then
            print*,'LINK1 error: cannot do vector plot of truss structure'
            print*,'Files created are empty !!'
        else if (element(1)%id == 2) then
            write(13,*) '% User scale parameter: See fedata - scale_vec'
            write(13,*) 'scale_vec = 1;'
            write(13,*) '% Define characteristic length'
            write(13,*) 'clx = 0;   cly = 0;'
            write(13,*) 'for e = 1:size(IX,1)'
            write(13,*) '    for i = 1:4'
            write(13,*) '        for j = i+1:4'
            write(13,*) '            clx = max(abs(  X(IX(e,i),1) - X(IX(e,j),1)   ),clx);'
            write(13,*) '            cly = max(abs(  X(IX(e,i),2) - X(IX(e,j),2)   ),cly);'
            write(13,*) '        end'
            write(13,*) '    end'
            write(13,*) 'end'
            write(13,*) 'clmax = max(clx, cly);'
            write(13,*) 'scal = max(max(abs(vect(:,1:2))))*sqrt(10)/clmax / scale_vec;'
            write(13,*) '% Make plot'
            write(13,*) 'figure'
            write(13,*) 'hold on'
            write(13,*) 'for e = 1:size(IX,1)'
            write(13,*) '    xx = X(IX(e,1:4),1);'
            write(13,*) '    yy = X(IX(e,1:4),2);'
            write(13,*) '    patch(xx,yy,[1 1 1]);'
            write(13,*) '    % Find approx. center for isoparametric element'
            write(13,*) '    xc = sum(X(IX(e,:),1))/size(IX,2);'
            write(13,*) '    yc = sum(X(IX(e,:),2))/size(IX,2);'
            write(13,*) '    % Directions for vect(:)'
            write(13,*) '    vec = [cos(-vect(e,3)) sin(-vect(e,3)) ...'
            write(13,*) '        cos(-vect(e,3)+pi/2) sin(-vect(e,3)+pi/2)];'
            write(13,*) '    % Plot magnitude and direction of vect_1'
            write(13,*) '    cc = ',"'", 'b',"'",';'
            write(13,*) '    if vect(e,1) < 0,    cc = ',"'",'r',"'",';     end'
            write(13,*) '    quiver(xc,yc,vec(1),vec(2),abs(vect(e,1))/scal,cc)'
            write(13,*) '    quiver(xc,yc,-vec(1),-vec(2),abs(vect(e,1))/scal,cc)'
            write(13,*) '    % Plot magnitude and direction of vect_2'
            write(13,*) '    cc = ',"'",'b',"'",';'
            write(13,*) '    if vect(e,2) < 0,    cc = ',"'",'r',"'",';     end'
            write(13,*) '    quiver(xc,yc,vec(3),vec(4),abs(vect(e,2))/scal,cc)'
            write(13,*) '    quiver(xc,yc,-vec(3),-vec(4),abs(vect(e,2))/scal,cc)'
            write(13,*) 'end   '
        else
            write (*,'("Unsupported element type in Matlab routine:")')
            write (*,'("plot -> plotmatlabevec")')
        endif
        
        ! End file and close
        write(13,*) 'title( ',"'", trim(title),"'",')'
        write(13,*) 'axis equal;  axis off;  hold off'
        write(13,*) 'set(gcf,',"'",'color',"'",',[ 1  1 1]);'
        close(13)

        ! Print to screen:
        print*,'Matlab plotting called: Vector'
        print*,'Generated data file: ', trim(fdata)
        print*,'and script file: ', trim(fscript)
        print*,' '
    end subroutine plotmatlabevec
!
!--------------------------------------------------------------------------------------------------
!
    subroutine plotmatlabeig(title,freq,evec,times)

        ! Eigenmode plot using Matlab
        ! freq  = eigenvalue
        ! evec  = eigenvector
        ! times = [totaltime, timeinterval]

        real(wp), dimension(:), intent(in) :: evec,times
        real(wp), intent(in) :: freq
        character(len=*), intent(in) :: title
        character(len=100) :: fdata, fscript
        integer i, j, e

        ! write datafile
        imatlab = imatlab + 1
        write (fdata, '(a,"_ploteig_data_",i0.4,".m")') trim(filename), imatlab
        open (13, file=trim(fdata))
        
        ! write nodal coordinates
        write(13,'("X = [")')
        do i = 1, size(x,1)
            write (13,'(3(f15.9,1x))') (x(i,j), j=1,size(x,2))
        end do
        write(13,'("];")')
        write(13,'( )')

        ! write topology matrix
        write(13,'("IX = [")')
        do e = 1, ne
            write (13,'(20(i6,1x))') (element(e)%ix(i) ,i=1,size(element(e)%ix,1))
        end do
        write(13,'("];")')

        ! write deformation vector
        write(13,'("D = [")')
        do i = 1, neqn
            write (13,'(f15.9,1x)') (evec(i) )
        end do
        write(13,'("];")')
        close(13)

        ! Plot script
        ! Create matlab script
        write(fscript, '(a,"_ploteig_",i0.4,".m")' ) trim(filename), imatlab
        open (13, file = trim(fscript))
        write(13,*) '% Plotting Eigenmodes'
        write(13,*) 'close all'
        write(13,*) 'clear all'
        write(13,*) fdata(1:len_trim(fdata)-2) // ';'
        write(13,*) 'freq = ',freq,';'
        write(13,*) 'timeint = ',times(2),';'
        write(13,*) 'timetot = ',times(1),';'
        write(13,*) '% Find max window size.'
        write(13,*) 'lxmin = min(X(:,1));        lxmax = max(X(:,1));'
        write(13,*) 'lymin = min(X(:,2));        lymax = max(X(:,2));'
        write(13,*) 'dxmin = min(D(1:2:end));    dxmax = max(D(1:2:end));'
        write(13,*) 'dymin = min(D(2:2:end));    dymax = max(D(2:2:end));'
        write(13,*) 'lxmin = lxmin - max(abs(dxmin),abs(dxmax))*1.05;'
        write(13,*) 'lxmax = lxmax + max(abs(dxmin),abs(dxmax))*1.05;'
        write(13,*) 'lymin = lymin - max(abs(dymin),abs(dymax))*1.05;'
        write(13,*) 'lymax = lymax + max(abs(dymin),abs(dymax))*1.05;'
        write(13,*) '% Make plot'
        write(13,*) 'figure'
        write(13,*) 'set(gcf,',"'",'color',"'",',[ 1  1 1]);'
        write(13,*) 'times = 0:timeint:timetot;'
        write(13,*) 'for i = 1:length(times)'
        write(13,*) '    tfact = sin(freq*times(i));'
        write(13,*) '    clf;'
        write(13,*) '    hold on'
        write(13,*) '    for e = 1:size(IX,1)'
        ! Element dependent code
        ! NOTE: not possible to mix element types !!!
        if (element(1)%id == 1) then
            write(13,*) '       edof = [2*IX(e,1)-1 2*IX(e,1) 2*IX(e,2)-1 2*IX(e,2)];'
            write(13,*) '       xx = X(IX(e,1:2),1) + tfact*D(edof(1:2:4));'
            write(13,*) '       yy = X(IX(e,1:2),2) + tfact*D(edof(2:2:4));'
            write(13,*) '       plot(xx,yy,',"'",'b',"',","'",'LineWidth',"'",',1.5)'
        elseif (element(1)%id == 2) then
            write(13,*) '       edof = [2*IX(e,1)-1 2*IX(e,1) 2*IX(e,2)-1 2*IX(e,2)...'
            write(13,*) '          2*IX(e,3)-1 2*IX(e,3) 2*IX(e,4)-1 2*IX(e,4)];'
            write(13,*) '       xx = X(IX(e,1:4),1) + tfact*D(edof(1:2:8));'
            write(13,*) '       yy = X(IX(e,1:4),2) + tfact*D(edof(2:2:8));'
            write(13,*) '       patch(xx,yy,[1 1 0]);'
        else
            write (*,'("Unsupported element type in Matlab routine:")')
            write (*,'("plotmatlab -> plotmatlabdef")')
        endif
        
        ! End file and close
        write(13,*) '    end'
        write(13,*) '    axis([lxmin lxmax lymin lymax])'
        write(13,*) '    axis off'
        write(13,*) '    title( ',"'", trim(title),"'",')'
        write(13,*) '    pause(0.01)'
        write(13,*) 'end'
        close(13)

        ! Print to screen:
        print*,'Matlab plotting called: eigenmode'
        print*,'Generated data file: ', trim(fdata)
        print*,'and script file: ', trim(fscript)
        print*,' '
    end subroutine plotmatlabeig

end module processor
