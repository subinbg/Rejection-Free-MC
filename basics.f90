module basics

    use, intrinsic :: iso_fortran_env

    private
    public                                    :: openfiles, closefiles, introduction, seedgen
                                      !rc     = LJ potential cut-off, unit L
                                      !dr_max = maximum deg. of rand. trans., unit L
    real, public                              :: rc, sigma, epsilon, dr_max
    real, public                              :: T, beta, density
    integer, public                           :: N, L
    integer, public                           :: nbin, total_mc_step
    integer, parameter, public                :: single=2


    real,    dimension(:,:), allocatable, public :: positions
    logical, dimension(:)  , allocatable, public :: cluster, moved
    real,    dimension(:)  , allocatable, public :: rbin
    integer, dimension(:)  , allocatable, public :: rhist


    integer, parameter, public                :: mrfile=10, Efile=11, grfile=12


    type, public :: potential_type
        real    :: potential_energy
        logical :: overlap=.false.
        
        contains
            procedure :: add_potential_type
            procedure :: substract_potential_type
            generic   :: operator(+) => add_potential_type
            generic   :: operator(-) => substract_potential_type
    end type potential_type


    type, public :: variable_type
        character(len=30)           :: nam
        real                        :: val              

        !abstract interface 
        !    subroutine func

        !    end subroutine
        !end interface
        
        integer, dimension(:), pointer :: array    => null() ! nrm: normalizer
        !procedure(func), pointer       :: function => null() ! avg: average value
        real                           :: add=0.0            ! msd: average squared value
        logical                        :: scalar = .TRUE.    ! err: estimated error
        real                           :: run_nrm, run_avg, run_msd, run_err
        real                           :: blk_nrm, blk_avg, blk_msd, blk_err
    end type variable_type   


    type, public :: variable_ptr
        type(variable_type), pointer :: p
    end type variable_ptr


    contains

        subroutine introduction
            write(*,*) " "
            write(*,*) " "
            write(*,*) "*********************************************"
            write(*,*) "Monte-Carlo Simulation of Lennard-Jones Fluid"
            write(*,*) "*********************************************"
            write(*,*) " "
        end subroutine introduction


        function seedgen (pid) result (seedg)
            implicit none
            integer(kind=int64) :: seedg
            integer, intent(IN) :: pid
            integer :: s

            call system_clock(s)
            seedg = abs( mod((s*181)*((pid-83)*359), 104729) ) 
        end function seedgen


        function add_potential_type (a,b) result (c)
            implicit none
            type(potential_type)              :: c
            class(potential_type), intent(in) :: a,b

            c%potential_energy = a%potential_energy + b%potential_energy
        end function add_potential_type


        function substract_potential_type (a,b) result (c)
            implicit none
            type(potential_type)              :: c
            class(potential_type), intent(in) :: a,b

            c%potential_energy = a%potential_energy - b%potential_energy
        end function substract_potential_type


        subroutine openfiles
            open(mrfile, file='move_ratio.dat'  , status='new')
            open(Efile , file='total_energy.dat', status='new')
            open(grfile, file='radial.dat'      , status='new')

            write(Efile,  100) "BlockN", "/Energy/", "STD"
            write(mrfile, 100) "BlockN", "/MoveRatio/", "STD"
            write(grfile, 100) "BlockN", "/Distance/", "G(r)"
100         format (a13,a13,a13)
        end subroutine openfiles


        subroutine closefiles
            close(mrfile)
            close(Efile )
            close(grfile)
        end subroutine closefiles


end module basics