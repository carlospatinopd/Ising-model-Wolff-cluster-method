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
    integer :: it, i, j, Si, n_add, ic
    real(8) :: r
    integer, dimension(nx,nx) :: S, cluster
    integer, dimension(2,4) :: s_add

    !------------------------------------------------------------------------!
    !                        Monte Carlo Wolff Method                        !
    !------------------------------------------------------------------------!

    call initial_state(S) ! asign random values to the spins (1 or -1)

    do it = 1,N
        cluster = 1
        call random_spin(i,j) ! randomly choose a spin
        cluster(i,j) = -1
        Si = S(i,j) ! value of the choosen spin
        call cluster_formation(i,j,n_add,s_add)
        do while (n_add > 0)
            do ic = 1,n_add
                call cluster_formation(s_add(1,ic),s_add(1,ic),n_add,s_add)
                cluster(s_add(1,ic),s_add(2,ic)) = -1
            end do
        end do
        do i = 1,nx
            do j = 1,ny
                S(i,j) = S(i,j)*cluster(i,j)
            end do 
        end do 
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

    ! neighbors with periodic boudary conditions for cluster formation
    subroutine cluster_formation(i,j,n_add,s_add)
        integer, intent(in) :: i, j
        integer :: ip, im, jp, jm
        integer, intent(out) :: s_add(2,4), n_add

        if (i == nx) then; ip = 1; else; ip = i+1; end if
        if (i == 1) then; im = nx; else; im = i-1; end if
        if (j == ny) then; jp = 1; else; jp = j+1; end if
        if (j == 1) then; jm = ny; else;  jm = j-1; end if

        n_add = 0
        call random_number(r)
        if (S(ip,j) == Si .AND. r >= p) then
            n_add = n_add+1
            s_add(1,n_add) = ip
            s_add(2,n_add) = j
        end if
        call random_number(r)
        if (S(im,j) == Si .AND. r >= p) then
            n_add = n_add+1
            s_add(1,n_add) = im
            s_add(2,n_add) = j
        end if
        call random_number(r)
        if (S(i,jp) == Si .AND. r >= p) then
            n_add = n_add+1
            s_add(1,n_add) = i
            s_add(2,n_add) = jp
        end if
        call random_number(r)
        if (S(i,jm) == Si .AND. r >= p) then
            n_add = n_add+1
            s_add(1,n_add) = i
            s_add(2,n_add) = jm
        end if

        Si = S(i,j)
    end subroutine cluster_formation

end program ising_wolff