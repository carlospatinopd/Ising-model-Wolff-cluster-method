program ising_wolff
    implicit none
    
    !------------------------------------------------------------------------!
    !                               parameters                               !
    !------------------------------------------------------------------------!
    integer, parameter :: N = 1e5                ! Monte Carlo steps
    integer, parameter :: L = 50,nx = L, ny = L   ! square lattice size
    real(8), parameter :: T = 2.2d0               ! temperature
    real(8), parameter :: p = 1.0d0-exp(-2.0d0/T) ! bonds probability
    
    !------------------------------------------------------------------------!
    !                               variables                                !
    !------------------------------------------------------------------------!
    integer :: it, i, j, Si, n_add, ic
    real(8) :: r
    integer, dimension(nx,nx) :: S ! array of spins
    logical, dimension(nx,ny) :: C ! array of clustered spins
    integer, dimension(2,4) :: s_add

    !------------------------------------------------------------------------!
    !                        Monte Carlo Wolff Method                        !
    !------------------------------------------------------------------------!

    call initial_state(S) ! asign random values to the spins (1 or -1)

    do it = 1,N
        C = .FALSE.
        call random_spin(i,j) ! randomly choose a spin
        C(i,j) = .TRUE.
        Si = S(i,j) ! value of the choosen spin
        call cluster_formation(i,j,n_add,s_add,C)
        do while (n_add > 0)
            do ic = 1,n_add
                Si=S(s_add(1,ic),s_add(2,ic))
                call cluster_formation(s_add(1,ic),s_add(2,ic),n_add,s_add,C)
            end do
        end do
        do i = 1,nx
            do j = 1,ny
                if (C(i,j)) S(i,j) = -S(i,j)
            end do 
        end do 
    end do

    open (unit = 100, file = "final_state.txt", status = "replace")
        do i = 1,L
            do j = 1,L
                write(100,*) S(i,j)
            end do
        end do
    close(100)

    !------------------------------------------------------------------------!
    !                        functions and subroutines                       !
    !------------------------------------------------------------------------!

    contains

    ! random initial state
    subroutine initial_state(S)
        integer, intent(out) :: S(nx,ny)
        open (unit = 100, file = "initial_state.txt", status = "replace")
            do i = 1,L
                do j = 1,L
                    call random_number(r)
                    S(i,j) = 2*int(r+0.5d0)-1
                    write(100,*) S(i,j)
                end do
            end do
        close(100)
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
    subroutine cluster_formation(i,j,n_add,s_add,C)
        integer, intent(in) :: i, j
        integer, intent(inout) :: s_add(2,4), n_add
        logical, intent(inout) :: C(nx,ny)
        integer :: ip, im, jp, jm

        ip = mod(i,nx)+1
        im = mod(i-2,nx)+1
        jp = mod(j,ny)+1
        jm = mod(j-2,ny)+1

        n_add = 0
        call random_number(r)
        if (S(ip,j) == Si .AND. r < p .AND. .NOT. C(ip,j)) then
            n_add = n_add+1
            s_add(1,n_add) = ip
            s_add(2,n_add) = j
            C(ip,j) = .TRUE.
        end if
        call random_number(r)
        if (S(im,j) == Si .AND. r < p .AND. .NOT. C(im,j)) then
            n_add = n_add+1
            s_add(1,n_add) = im
            s_add(2,n_add) = j
            C(im,j) = .TRUE.
        end if
        call random_number(r)
        if (S(i,jp) == Si .AND. r < p .AND. .NOT. C(i,jp)) then
            n_add = n_add+1
            s_add(1,n_add) = i
            s_add(2,n_add) = jp
            C(i,jp) = .TRUE.
        end if
        call random_number(r)
        if (S(i,jm) == Si .AND. r < p .AND. .NOT. C(i,jm)) then
            n_add = n_add+1
            s_add(1,n_add) = i
            s_add(2,n_add) = jm
            C(i,jm) = .TRUE.
        end if
    end subroutine cluster_formation

end program ising_wolff