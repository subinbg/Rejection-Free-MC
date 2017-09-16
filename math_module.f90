module math_module

    use, intrinsic :: iso_fortran_env
    public         :: box_muller, random_translation, metropolis
    real, public   :: pi = 4.0*atan(1.0)

    contains

        subroutine box_muller(ranseed)
            real, intent(inout) :: ranseed
            real                :: U, V, tpi

            tpi=2*pi
            U = 0; V = 0
            do while(.true.)
                call random_number(U)
                call random_number(V)
                if(U>epsilon(V)) exit
            enddo
            ranseed = sqrt(-2*log(U))*cos(tpi*V)
        end subroutine box_muller

        function random_translation (dr_max, vector) result (new_position)
            implicit none
            real, dimension(3)             :: new_position
            real, intent(in)               :: dr_max
            real, dimension(3), intent(in) :: vector
            real, dimension(3)             :: xyz

            call random_number(xyz)
            xyz(:) = 2.0 * xyz(:) - 1.0

            new_position(:) = vector(:) + xyz(:) * dr_max
        end function random_translation

        function metropolis ( delta ) result ( accept )
            implicit none
            logical          :: accept
            real, intent(in) :: delta
            real, parameter  :: exponent_guard = 75.0
            real             :: randcomp

            if ( delta > exponent_guard ) then
                accept = .true.
            else if ( delta < 0.0 ) then
                accept = .false.
            else
                call random_number(randcomp)
                accept = (1 - exp(-delta)) > randcomp
            end if
        end function metropolis

end module math_module