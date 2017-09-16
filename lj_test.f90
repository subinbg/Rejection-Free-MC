program lj_test

    use, intrinsic :: iso_fortran_env!, only: input_unit... only: Fortran2003~
    
    use math_module
    use basics
    use mc_module


    ! MPI Library
    use mpi 
    !!!!!!!!!!!!!

    implicit none

    ! MPI variables!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                 :: myid , ntasks , ierr , islave
    integer , dimension ( MPI_STATUS_SIZE ) :: status
    real                                    :: E_avg, E_err, m_r_avg, m_r_err
    real, allocatable                       :: hist_temp(:)
    integer                                 :: energy_avg_tag=1,energy_err_tag=2
    integer                                 :: moveR_avg_tag =3,moveR_err_tag =4
    integer                                 :: gr_tag=5
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    type(potential_type)                          :: total, atom_old, atom_new
    type(variable_type), target                   :: m_r, E, Cv, gr
    type(variable_ptr), dimension(:), allocatable :: variables

    
    real, dimension(3) :: ri, rij, pivot
    real, allocatable  :: hist_norm(:)
    real               :: m_ratio,delta,rij_sq,rad
    integer            :: nstep, nblock, moves, pindex
    integer            :: i,j,blk,stp,norm_moves
    logical            :: cluster_growth


    ! Initialization of MPI!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call MPI_INIT ( ierr )
    call MPI_COMM_RANK ( MPI_COMM_WORLD , myid , ierr )
    call MPI_COMM_SIZE ( MPI_COMM_WORLD , ntasks , ierr )
    call random_seed(put=[seedgen(myid)])
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    nblock  = 10
    nstep   = 500
    nbin    = 1000
    
    rc      = 0.3
    sigma   = 1
    epsilon = 1
    dr_max  = 0.1
    T       = 1.0
    beta    = 1.0/T
    N       = 1000
    L       = 10
    
    density = real(N)/L**3
    total_mc_step = nblock * nstep * N
    
    call random_stabilize
    call record_preparation

    call allocate_arrays
    allocate(hist_norm(nbin))
    gr%run_nrm = 0.0

    ! Basic Length Unit = L
    positions(:,:) = positions(:,:) / L
    positions(:,:) = positions(:,:) - anint( positions(:,:) )

    total = potential_total()

    ! MPI master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (myid == 0) then
        allocate(hist_temp(nbin))
        call introduction
        call openfiles
        write(*,*) "Number of Parallel Machines =", ntasks
        write(*,*) "Initial Potential Energy =", total%potential_energy
        write(*,*) " "

        total = potential_total()
        if ( total%overlap ) then
            stop 'ERROR in Initialization: Severe Overlap Occurred.'
        end if

        write(*,*) "Initialization Completed."
        write(*,*) "Main Monte Carlo Procedure..."
        write(*,*) " "
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call run_begin
    m_ratio = 0.0

    ! Main Monte-Carlo Procedure
    do blk = 1, nblock

        call blk_begin

        do stp = 1, nstep

            moves = 0
            norm_moves = 0
            
            !do i = 1, N
            call random_number(rad)
            i = floor(rad * real(N))

                call random_number(pivot)
                pivot(:) = 2.0 * pivot(:) - 1.0
                pivot(:) = pivot(:) - anint( pivot(:) )
                
                cluster(:) = .false.
                moved  (:) = .false.
                cluster_growth = .true.
                pindex = i
                
                do while (cluster_growth)
                    
                    cluster(pindex) = .true.

                    ri(:) = 2.0*pivot(:) - positions(:,pindex)
                    ri(:) = ri(:) - anint( ri(:) )

                    do j=1,N
                        if (pindex==j) cycle
                        if (moved(j))  cycle

                        !rij(:) = positions(:,pindex) - positions(:,j)
                        !rij(:) = rij(:) - anint( rij(:) )
                        !rij_sq = sum( rij**2 )

                        atom_old = potential_RF(positions(:,pindex), j, global1)
                        atom_new = potential_RF(ri(:), j, global2)

                        if ( ( global1 ) .and. ( global2 ) ) cycle
                        
                        delta = atom_new%potential_energy - atom_old%potential_energy
                        delta = delta * beta

                        !print *, delta, atom_new%potential_energy, atom_old%potential_energy

                        norm_moves = norm_moves + 1
                        if ( metropolis(delta) .and. (.not.(total%overlap)) ) then
                            cluster(j) = .true.
                            !call move( pindex, ri )
                            moves = moves + 1
                        end if
                    end do

                    moved(pindex) = .true.
                    call move ( pindex, ri )
                    pindex = -1

                    do j=1,N
                        if ( ( cluster(j) ) .and. ( .not. moved(j) ) ) then
                            pindex = j
                            exit
                        end if
                    end do

                    if (pindex < 0) cluster_growth = .false.
                
                end do  

            !end do

            m_ratio = real(moves) / real(norm_moves)

            total = potential_total()
            call calc_variables
            call blk_add


            ! MPI master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (myid == 0 ) then
                write(*,fmt=101) "Completion(%)  BLOCK:", &
                                & real(blk)/real(nblock)*100, "STEP:", real(stp)/real(nstep)*100
101             format(a23, f7.2, a7, f7.2)
            end if
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        end do

        call blk_end

        ! MPI Section!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (myid .ne. 0) then ! Slave node
            call MPI_SEND(  E%blk_avg,1,MPI_REAL,0,energy_avg_tag,MPI_COMM_WORLD,ierr)
            call MPI_SEND(  E%blk_err,1,MPI_REAL,0,energy_err_tag,MPI_COMM_WORLD,ierr)
            call MPI_SEND(m_R%blk_avg,1,MPI_REAL,0, moveR_avg_tag,MPI_COMM_WORLD,ierr)
            call MPI_SEND(m_R%blk_err,1,MPI_REAL,0, moveR_err_tag,MPI_COMM_WORLD,ierr)
        else ! Master node
            do islave = 1, ntasks-1
                call MPI_RECV(  E_avg,1,MPI_REAL,islave,energy_avg_tag,MPI_COMM_WORLD,status,ierr)
                call MPI_RECV(  E_err,1,MPI_REAL,islave,energy_err_tag,MPI_COMM_WORLD,status,ierr)
                call MPI_RECV(m_r_avg,1,MPI_REAL,islave, moveR_avg_tag,MPI_COMM_WORLD,status,ierr)
                call MPI_RECV(m_r_err,1,MPI_REAL,islave, moveR_err_tag,MPI_COMM_WORLD,status,ierr)

                E%blk_avg   = E%blk_avg   + E_avg
                E%blk_err   = E%blk_err   + E_err
                m_R%blk_avg = m_R%blk_avg + m_r_avg
                m_R%blk_err = m_R%blk_err + m_r_err
            end do

            write(Efile ,102) blk, E%blk_avg/real(ntasks)  , E%blk_err/real(ntasks)
            write(mrfile,102) blk, m_r%blk_avg/real(ntasks), m_r%blk_err/real(ntasks)
102         format(I10, f, f)
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (blk > nblock/5)  then 
            call radial_check
            gr%run_nrm = gr%run_nrm + 1.0
            hist_norm(:) = rhist(:) / ( gr%run_nrm*N)  !/(rbin(2)-rbin(1))

            do j = 1,nbin
                hist_norm(j) = hist_norm(j) / ( 4.0/3.0 * pi * density * & 
                                                   & ( rbin(j+1)**3 - rbin(j)**3 ) )
                !write(grfile,102) blk, rbin(j), hist_norm(j)
            end do

            ! MPI Section!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (myid .ne. 0) then ! Slave node
                call MPI_SEND(hist_norm,nbin,MPI_REAL,0,gr_tag,MPI_COMM_WORLD,ierr)
            else ! Master node
                do islave = 1, ntasks-1
                    call MPI_RECV(hist_temp,nbin,MPI_REAL,islave,gr_tag,MPI_COMM_WORLD,status,ierr)
                    hist_norm(:) = hist_norm(:) + hist_temp(:)
                end do

                do j=1,nbin
                    write(grfile,102) blk, rbin(j), hist_norm(j)/real(ntasks)
                end do
            end if
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        end if

    end do

    call run_end

    ! MPI Master!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (myid == 0 ) then
        call closefiles
        deallocate(hist_temp)
        write(*,*) " "
        write(*,*) "All Computation Completed."
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call deallocate_arrays
    deallocate(hist_norm)

    ! MPI Section!!!!!!!!!!!!!!
    call MPI_FINALIZE ( ierr )
    !!!!!!!!!!!!!!!!!!!!!!!!!!!

    contains

        subroutine record_preparation
            implicit none
            allocate(variables(3))
            variables(1)%p => m_r
            variables(2)%p => Cv
            variables(3)%p => E

        end subroutine record_preparation


        subroutine calc_variables
            implicit none

            !type(variable_type), dimension(3) :: variables
            ! 단순하게 이렇게 써 버리면, 이 function과 subroutine을 나가는 순간
            ! variables 라는 변수는 main program에 아무런 영향을 주지 못한다.
            ! variables는 전역변수 m_r, Cv, E로부터 만들어진 함수 지역변수이며
            ! variables는 전역변수들의 주소를 담고 있는 것이 아니기 때문이다.
            ! 참고로 fortran array는 주소 copy를 지원하지 않는다. 즉,
            ! A(:) = 1, B(:) = 2, A = B, B=3 을 하였다고 해서 A 또한 3으로
            ! 이루어진 array가 되지 않는다.

            !type(variable_type), dimension(3), pointer :: variables
            ! 이것도 맞지 않다. 이것은 variable_type의 pointer 세 개를 담는
            ! array가 아니라, variable_type이 세 개 들어가 있는 array를 향한
            ! pointer이기 때문이다.

            ! pointer array를 만들고 싶다면 아래와 같이 하면 된다.
            !type variable_ptr
            !    type(variable_type), pointer :: p
            !end type variable_ptr
            !type(variable_ptr), dimension(3) :: variables
            !variables(1)%p => m_r
            !variables(2)%p => Cv
            !variables(3)%p => E

            m_r%nam = 'Move ratio'
            m_r%val = m_ratio
            
            Cv%nam = 'Cv/N'
            Cv%val = total%potential_energy/(T*sqrt(real(N)))

            E%nam = 'Total energy'
            E%val = 1.5*T + total%potential_energy/real(N) 

            !gr%nam   = 'Radial distribution function'
            !gr%val   = 0.0
            !gr%array => rhist

            !variables = [m_r, E, Cv]
        end subroutine calc_variables


        subroutine run_begin !(variables)
            implicit none
            integer :: i,j
            !type(variable_type), dimension(:) :: variables
            j = size(variables(:),1)
            
            do i=1,j
                !gr%scalar = .false.
                variables(i)%p%run_nrm=0.0
                variables(i)%p%run_avg=0.0
                variables(i)%p%run_msd=0.0
                variables(i)%p%run_err=0.0
            end do
        end subroutine run_begin


        subroutine blk_begin !(variables)
            implicit none
            integer :: i,j
            !type(variable_type), dimension(:) :: variables
            j = size(variables(:),1)

            do i=1,j
                variables(i)%p%blk_nrm=0.0
                variables(i)%p%blk_avg=0.0
                variables(i)%p%blk_msd=0.0
                variables(i)%p%blk_err=0.0
            end do
        end subroutine blk_begin


        subroutine blk_add!( variables )
            implicit none
            !type(variable_type), dimension(:) :: variables
            !type(variable_ptr), dimension(3)  :: variables
            integer :: i,j
            j = size(variables(:),1)

            do i=1,j
                variables(i)%p%blk_avg = variables(i)%p%blk_avg + variables(i)%p%val
                variables(i)%p%blk_msd = variables(i)%p%blk_msd + variables(i)%p%val**2
                variables(i)%p%blk_nrm = variables(i)%p%blk_nrm + 1.0
            end do
        end subroutine blk_add


        subroutine blk_end!( variables )
            implicit none
            !type(variable_type), dimension(:) :: variables
            !type(variable_ptr), dimension(3)  :: variables
            integer :: i,j
            j = size(variables(:),1)

            do i=1,j
                variables(i)%p%blk_avg = variables(i)%p%blk_avg / variables(i)%p%blk_nrm
                variables(i)%p%blk_msd = variables(i)%p%blk_msd / variables(i)%p%blk_nrm
                variables(i)%p%blk_err = variables(i)%p%blk_msd - variables(i)%p%blk_avg**2
                variables(i)%p%blk_err = sqrt( variables(i)%p%blk_err )

                variables(i)%p%run_avg = variables(i)%p%run_avg + variables(i)%p%blk_avg
                variables(i)%p%run_msd = variables(i)%p%run_msd + variables(i)%p%blk_avg**2
                variables(i)%p%run_nrm = variables(i)%p%run_nrm + 1.0
            end do
        end subroutine blk_end


        subroutine run_end!( variables )
            implicit none
            !type(variable_type), dimension(:) :: variables
            !type(variable_ptr), dimension(3)  :: variables
            integer :: i,j
            j = size(variables(:),1)
            
            do i=1,j
                variables(i)%p%run_avg = variables(i)%p%run_avg / variables(i)%p%run_nrm
                variables(i)%p%run_msd = variables(i)%p%run_msd / variables(i)%p%run_nrm
                variables(i)%p%run_err = variables(i)%p%run_msd - variables(i)%p%run_avg**2
                variables(i)%p%run_err = sqrt( variables(i)%p%run_err )
            end do
            !dim = size(variables%scalar, 1)

        end subroutine run_end

end program lj_test
