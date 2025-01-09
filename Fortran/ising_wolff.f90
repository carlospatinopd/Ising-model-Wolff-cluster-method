program ising_wolff
    implicit none
    
    !------------------------------------------------------------------------!
    !                               parameters                               !
    !------------------------------------------------------------------------!

    integer, parameter :: N = 1e3                 ! Monte Carlo steps
    integer, parameter :: L = 50,nx = L, ny = L   ! square lattice size
    real(8), parameter :: T = 2.0d0              ! temperature
    real(8), parameter :: p = 1.0d0-exp(-2.0d0/T) ! bonds probability
    integer, parameter :: f = 1e0                 ! saving data frecuency
    
    !------------------------------------------------------------------------!
    !                               variables                                !
    !------------------------------------------------------------------------!

    integer :: it, i, j, Si, n_add, ic, E, M, ip, jp
    real(8) :: r
    integer, dimension(nx,nx) :: S ! array of spins
    logical, dimension(nx,ny) :: C ! array of clustered spins
    integer, dimension(2,nx*ny) :: s_add

    !------------------------------------------------------------------------!
    !                        Monte Carlo Wolff Method                        !
    !------------------------------------------------------------------------!

    call initial_state(S) ! asign random values to the spins (1 or -1)

    open (unit = 100, file = "state_evolution.txt", status = "replace")
    open (unit = 101, file = "magnetization.txt", status = "replace")
    open (unit = 102, file = "energy.txt", status = "replace")
    E = 0
    M = 0
    do it = 1,N
        C = .FALSE.
        call random_spin(i,j) ! randomly choose a spin
        C(i,j) = .TRUE.
        Si = S(i,j) ! value of the choosen spin
        n_add = 1
        s_add(1,1) = i
        s_add(2,1) = j
        
        do ic = 1,n_add
            call cluster_formation(Si,s_add(1,ic),s_add(2,ic),n_add,s_add,C)
        end do

        ! cluster flip
        do i = 1,nx
            do j = 1,ny
                if (C(i,j)) S(i,j) = -S(i,j)
            end do 
        end do 

        ! save the data
        if (mod(it,f) == 0) then
            M = 0
            E = 0
            do i = 1,nx
                ip = mod(i,nx)+1
                do j = 1,ny
                    jp = mod(j,ny)+1
                    M = M+S(i,j)
                    E = E-(S(i,j)*(S(ip,j)+S(i,jp)))
                    write(100,*) S(i,j)
                end do 
            end do
            write(101,*) M
            write(102,*) E/2
        end if
    end do
    close(100)
    close(101)
    close(102)

    !------------------------------------------------------------------------!
    !                        functions and subroutines                       !
    !------------------------------------------------------------------------!

    contains

    ! random initial state
    subroutine initial_state(S)
        integer, intent(out) :: S(nx,ny)

        open (unit = 100, file = "parameters.txt", status = "replace")
        open (unit = 101, file = "initial_state.txt", status = "replace")
            write(100,*) L, N, f
            do i = 1,L
                do j = 1,L
                    call random_number(r)
                    S(i,j) = 2*int(r+0.5d0)-1
                    write(101,*) S(i,j)
                end do
            end do
        close(100)
        close(101)

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
    subroutine cluster_formation(Si,i,j,n_add,s_add,C)
        integer, intent(in) :: Si, i, j
        integer, intent(inout) :: s_add(2,nx*ny), n_add
        logical, intent(inout) :: C(nx,ny)
        integer :: ip, im, jp, jm

        ip = mod(i,nx)+1
        im = mod(i-2,nx)+1
        jp = mod(j,ny)+1
        jm = mod(j-2,ny)+1

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