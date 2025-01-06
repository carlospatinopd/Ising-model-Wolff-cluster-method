program ising_wolff
    implicit none
    
    !------------------------------------------------------------------------!
    !                               parameters                               !
    !------------------------------------------------------------------------!
    integer, parameter :: N = 1e5               ! Monte Carlo steps
    integer, parameter :: L = 50,nx = L, ny = L ! square lattice size
    real(8), parameter :: T = 3                 ! temperature
    real(8), parameter :: p = 1-exp(-2/T)       ! bonds probability
    
    !------------------------------------------------------------------------!
    !                               variables                                !
    !------------------------------------------------------------------------!
    integer :: it, i, j, ip, im, jp, jm
    real(8) :: r
    integer, dimension(nx,nx) :: S

    !------------------------------------------------------------------------!
    !                        Monte Carlo Wolff Method                        !
    !------------------------------------------------------------------------!

    call initial_state(S) ! asign random values to the spins (1 or -1)

    do it = 1,N
        call random_spin(i,j) ! randomly choose a spin
        call pbc(ip,im,jp,jm) ! neighbors with periodic boundary conditions

    end do

    !------------------------------------------------------------------------!
    !                        functions and subroutines                       !
    !------------------------------------------------------------------------!

    contains

    ! random initial state
    subroutine initial_state(S)
        integer, intent(out) :: S(nx,ny)

        do i = 1,L
            do j = 1,L
                call random_number(r)
                S(i,j) = 2*int(r+0.5d0)-1
            end do
        end do
    end subroutine initial_state

    ! random initial state
    subroutine random_spin(i,j)
        integer, intent(out) :: i, j
        call random_number(r)
        i = 1+int(r*nx)
        call random_number(r)
        j = 1+int(r*ny)
    end subroutine random_spin

    ! neighbors and periodic boudary conditions
    subroutine pbc(ip,im,jp,jm)
        integer, intent(out) :: ip, im, jp, jm

        if (i == nx) then
            ip = 1
        else
            ip = i+1
        end if
        if (i == 1) then
            im = nx
        else
            im = i-1
        end if
        if (j == ny) then
            jp = 1
        else
            jp = j+1
        end if
        if (j == 1) then
            jm = ny
        else
            jm = j-1
        end if
    end subroutine pbc


end program ising_wolff