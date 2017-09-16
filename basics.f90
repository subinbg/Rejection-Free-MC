module basics

    ! (여기서 '변수'는 단순한 변수 뿐만 아니라 type, subroutine, function 포함)
    ! 모듈 내에서 선언된 변수들은 모듈 내에서 전역변수가 된다.
    ! 모듈 내에서 public으로 선언된 변수들은 모듈과 그 모듈을 사용하는 상위 파일
    ! 모두에서의 전역변수가 된다.
    ! 특별히 public 혹은 private로 선언되지 않은 변수들은 모두 public(default) 이다.
    ! 모듈 내에서 private로 선언된 변수들은 모듈 내에서만 통용되는 전역변수가 된다.
    ! 따라서 private로 선언된 변수를 수정하거나 출력하려면 반드시 module 내에
    ! 그 변수를 수정하거나 출력하게 하는 function/subroutine 등을 만들어야 한다.

    ! 마찬가지로, 본 program에서 선언된 변수들은
    ! 본 program의 contains...에 들어있는 function과 subroutine들에
    ! 영향을 주는 전역변수이다.

    ! 실험 결과, module에서 불러오는 subroutine/function
    ! 본 program에 contain되는 subroutine/function
    ! 모두 call-by-reference 이다. 즉 예컨대 array를 subroutine에 전달한 후
    ! 함수 내부에서 그 array의 내용을 바꿀 경우, 함수 외부에서도 바뀐 내용이 유지된다.
    ! 변수여도 마찬가지이다. 변수 A와 B를 선언한 후
    ! 어떤 함수 swap(A,B)를 부르면 그 함수가 어디서 선언되었건 간에
    ! 실제 변수 A와 B의 값을 변경시킨다.
    ! 따라서 만약에 그러한 상황을 방지하고 싶다면 intent(in) 선언을 하든지
    ! 함수 내의 다른 변수에 복사하여 연산을 수행해야 한다.

    ! 우리에게 program A, module B, module C가 있다고 하자.
    ! 만약 program A와 module C에서 use B를 하고 있다면,
    ! (ifort -c module_B.f90)
    ! (ifort -I module_B.o -c module_C.f90)
    ! (ifort module_C.o program_A.f90)
    ! B에서 선언된 변수가 변경될 경우, 그것이 program A와 module C에 모두 반영된다.
    ! 즉 어디선가 module B의 변수가 변경되었을 경우, program A와 module C에서
    ! 그 변수를 불러오면 변경된 값이 불러져온다.
    ! 따라서 module B는 program A와 module C의 "공통" 전역변수인 것이다.

    ! 한편 module B가 program이 아니라 module C와 module D에서만
    ! use 되고 있다면, module C와 module D가 program에서 use되고 있다고 하더라도
    ! module B가 계속 main program에 올려져 있을 것인가는 알 수 없다.
    ! 이 때 module B의 값을 공통 전역변수마냥 계속 살려두고 싶다면
    ! 변수를 선언할 때 save를 덧붙이거나, module B 또한 program에
    ! use 시켜야만 한다.
    ! (ifort module_C.o module_B.o program_A.f90)

    use, intrinsic :: iso_fortran_env
    ! 여기서의 intrinsic은 fortran 표준 내장함수를 사용한다는 것이다.
    ! complier마다 고유한 내장함수가 있을 수 있으며,
    ! 사용자가 정의한 모듈 중에서 우연히 표준 내장함수와 이름이 같은 것이
    ! 있을 수 있다. 그런 경우를 방지하기 위하여 intrinsic 키워드를 쓴다.
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

        ! interface: procedure(function, subroutine)을 type으로 지정하고 싶을 때 사용
        ! abstract : interface에 func/subrout 이름을 명시하지 않고
        ! 단지 해당 interface의 구조를 갖고 있는 func/subrout와 연결하고 싶을 때 사용
        !abstract interface 
        !    subroutine func ! radial_check에 연결할 생각이었음.

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
            ! type(potential_typ) 으로 선언할 시,
            ! The passed-object dummy argument must be a polymorphic 
            ! dummy data object if the type being defined is extensible 오류 출력.

            c%potential_energy = a%potential_energy + b%potential_energy
        end function add_potential_type


        function substract_potential_type (a,b) result (c)
            implicit none
            type(potential_type)              :: c
            class(potential_type), intent(in) :: a,b
            ! type(potential_typ) 으로 선언할 시,
            ! The passed-object dummy argument must be a polymorphic 
            ! dummy data object if the type being defined is extensible 오류 출력.

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