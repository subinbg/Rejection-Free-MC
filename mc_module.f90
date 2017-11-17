module mc_module

    use, intrinsic :: iso_fortran_env
    use            :: basics

    public :: allocate_arrays, deallocate_arrays
    public :: potential_total, potential_atom   , potential_RF
    public :: move           , create           , destroy
    public :: radial_check   

    integer, parameter, private :: lt=-1, gt=1
    logical,            public  :: global1, global2

    contains

        subroutine random_stabilize
            implicit none
            integer :: i
            real    :: j

            do i=1,2000    
                call random_number(j)
            end do
        end subroutine random_stabilize


        subroutine allocate_arrays
            implicit none

            real               :: increment, fraction
            real, dimension(3) :: origin, loc
            integer            :: i,j,k,lattice,number

            if (rc > 0.5) then
                stop 'Error: rc too large'
            end if

            allocate(positions(3,N), cluster(N), moved(N))
            !call random_number(positions)

            cluster(:)      = .false.
            moved(:)        = .false.

            lattice   = int( real(N)**(1.0/3.0) + 1.0 )
            fraction  = L/real(lattice+2)!/2
            origin(:) = fraction
            number = 0

            do i=1,lattice
                do j=1,lattice
                    do k=1,lattice
                        number = number + 1
                        if (number > N) exit

                        loc(:) = fraction * [i-1,j-1,k-1]
                        positions(:,number) =  origin(:) + loc(:)
                        
                    end do
                end do                        
            end do            

            allocate(rbin(nbin+1))
            allocate(rhist(nbin))

            rhist(:)=0
            rbin(1)=0

            increment = L/(2*real(nbin))
            do i=2,nbin+1
                rbin(i) = rbin(i-1) + increment
            end do
        end subroutine allocate_arrays


        subroutine deallocate_arrays
            implicit none
            deallocate(positions)
            deallocate(rbin)
            deallocate(rhist)
            deallocate(cluster)
            deallocate(moved)
        end subroutine deallocate_arrays


        function potential_total () result (totalE)
            implicit none

            type(potential_type) :: totalE
            type(potential_type) :: atom
            integer              :: i, flag=-1

            totalE = potential_type ( potential_energy = 0.0, overlap=.false.)

            do i = 1, N-1
                atom = potential_atom ( positions(:,i), i, gt )
                if ( atom%overlap ) then
                    totalE%overlap = .true.
                    return
                end if
                totalE = totalE + atom
            end do

        end function potential_total


        function potential_atom ( atom_position, i, jrange ) result (atom)
            implicit none

            type(potential_type)           :: atom
            type(potential_type)           :: atom_others
            real, dimension(3), intent(in) :: atom_position
            real, dimension(3)             :: rij
            real                           :: rij_sq, LL
            real, parameter                :: sr2_overlap = 10.0 ! overlap threshold
            integer                        :: j, j1, j2

            integer, intent(in)            :: i
            integer, optional, intent(in)  :: jrange
            
            LL = L**2

            if ( present(jrange) ) then
                SELECT CASE ( jrange )
                    CASE ( lt ) ! j < i
                        j1 = 1
                        j2 = i-1
                    CASE ( gt ) ! j > i
                        j1 = i+1
                        j2 = N
                    case ( single )
                        j1 = i
                        j2 = i
                    CASE default ! should never happen
                        WRITE ( unit = error_unit, fmt='(a,i10)') 'j_range error ', jrange
                        STOP 'Impossible error in potential_atom'
                END SELECT
            else
                j1 = 1
                j2 = N
            end if

            atom = potential_type ( potential_energy = 0.0, overlap=.false.)

            do j = j1,j2
                if (i == j) cycle

                rij(:) = atom_position(:) - positions(:,j)
                rij(:) = rij(:) - anint( rij(:) )
                rij_sq = sum( rij**2 )

                !write(mrfile, *) i,j,positions(:,j)
                !write(mrfile, *) i,j,atom_position(:)

                if (rij_sq < rc**2) then

                    rij_sq = sigma**2/(rij_sq * LL)

                    !if ( rij_sq > sr2_overlap) then
                    !    atom%overlap = .true.
                    !    return
                    !end if

                    atom_others%potential_energy = &
                                        & (rij_sq)**6 - (rij_sq)**3

                    atom = atom + atom_others
                end if
            end do
        end function potential_atom


        function potential_RF ( atom_position, i, global) result (atom)
            implicit none

            type(potential_type)           :: atom
            type(potential_type)           :: atom_others
            real, dimension(3), intent(in) :: atom_position
            real, dimension(3)             :: rij
            real                           :: rij_sq, LL
            real, parameter                :: sr2_overlap = 10.0 ! overlap threshold
            integer                        :: j, j1, j2

            integer, intent(in)            :: i
            logical                        :: global
            
            LL = L**2

            atom = potential_type ( potential_energy = 0.0, overlap=.false.)



            rij(:) = atom_position(:) - positions(:,i)
            rij(:) = rij(:) - anint( rij(:) )
            rij_sq = sum( rij**2 )

            if (rij_sq < rc**2) then
                global = .false.
                rij_sq = sigma**2/(rij_sq * LL)
                atom%potential_energy = (rij_sq)**6 - (rij_sq)**3
            else
                global = .true.
            end if

        end function potential_RF


        subroutine move ( i, ri )
            implicit none
            integer, intent(in) :: i
            real, dimension(3), intent(in) :: ri

            positions(:,i) = ri(:)
        end subroutine move


        subroutine create ( ri )
            implicit none
            real, dimension(3), intent(in) :: ri

            N = N+1
            positions(:,N) = ri(:)
        end subroutine create


        subroutine destroy ( i )
            implicit none
            integer, intent(in) :: i

            positions(:,i) = positions(:,N)
            N = N-1
        end subroutine destroy


        subroutine radial_check
            implicit none
            integer                 :: i,j,ix
            real                    :: rij_n
            real, dimension(3)      :: rij
            real, dimension(nbin+1) :: rij_eval, result
            
            do i=1,N-1
                do j=i+1,N
                    rij = positions(:,i) - positions(:,j)
                    rij(:) = rij(:) - anint( rij(:) )
                    rij_n = sqrt( sum( rij**2 ) )
                    rij_n = rij_n * L

                    if(rij_n < L/2) then
                        rij_eval = rbin(:) - rij_n
                        result = pack([(ix,ix=1,size(rij_eval))],rij_eval.gt.0)
                        rhist(result(1)-1) = rhist(result(1)-1) + 2
                    end if
                end do
            end do
        end subroutine radial_check

end module mc_module