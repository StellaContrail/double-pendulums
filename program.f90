! Physics constant values
module constant
    implicit none
    double precision,parameter :: l = 0.3d0, m = 0.5d0, g = 9.8d0, pi = acos(-1.0d0)
end module

! Main module
module double_pendulum
    use constant
    implicit none
contains
    ! Runge-kutta method of four order
    subroutine rungekutta(x, v, dt)
        ! x -> θ, v -> dθ/dt
        ! abbreviated because of Fortran cannot handle more than 132 words in line
        double precision,intent(inout) :: x(:), v(:)
        double precision,intent(in) :: dt
        ! ex. k1(1) -> k1 of x1, k1(2) -> k1 of x2
        double precision :: k1(4), k2(4), k3(4), k4(4)
        k1(1) = dt * f_x(v(1))
        k1(2) = dt * f_x(v(2))
        k1(3) = dt * f_v1(x(1), x(2), v(1), v(2))
        k1(4) = dt * f_v2(x(1), x(2), v(1), v(2))

        k2(1) = dt * f_x(v(1) + 0.5d0*k1(3))
        k2(2) = dt * f_x(v(2) + 0.5d0*k1(4))
        k2(3) = dt * f_v1(x(1) + 0.5d0*k1(1), x(2) + 0.5d0*k1(2), v(1) + 0.5d0*k1(3), v(2) + 0.5d0*k1(4))
        k2(4) = dt * f_v2(x(1) + 0.5d0*k1(1), x(2) + 0.5d0*k1(2), v(1) + 0.5d0*k1(3), v(2) + 0.5d0*k1(4))

        k3(1) = dt * f_x(v(1) + 0.5d0*k2(3))
        k3(2) = dt * f_x(v(2) + 0.5d0*k2(4))
        k3(3) = dt * f_v1(x(1) + 0.5d0*k2(1), x(2) + 0.5d0*k2(2), v(1) + 0.5d0*k2(3), v(2) + 0.5d0*k2(4))
        k3(4) = dt * f_v2(x(1) + 0.5d0*k2(1), x(2) + 0.5d0*k2(2), v(1) + 0.5d0*k2(3), v(2) + 0.5d0*k2(4))

        k4(1) = dt * f_x(v(1) + k3(3))
        k4(2) = dt * f_x(v(2) + k3(4))
        k4(3) = dt * f_v1(x(1) + k3(1), x(2) + k3(2), v(1) + k3(3), v(2) + k3(4))
        k4(4) = dt * f_v2(x(1) + k3(1), x(2) + k3(2), v(1) + k3(3), v(2) + k3(4))
        
        x(1) = x(1) + (k1(1) + 2d0*(k2(1) + k3(1)) + k4(1)) / 6d0
        x(2) = x(2) + (k1(2) + 2d0*(k2(2) + k3(2)) + k4(2)) / 6d0
        v(1) = v(1) + (k1(3) + 2d0*(k2(3) + k3(3)) + k4(3)) / 6d0
        v(2) = v(2) + (k1(4) + 2d0*(k2(4) + k3(4)) + k4(4)) / 6d0
    end subroutine

    ! Calculate positions of balls from θ
    subroutine get_position(theta, q)
        double precision,intent(in) :: theta(:)
        double precision,intent(out) :: q(:,:)
        double precision sin1, sin2, cos1, cos2
        sin1 = sin(theta(1))
        cos1 = cos(theta(1))
        sin2 = sin(theta(2))
        cos2 = cos(theta(2))

        q(1,1) = l * sin1
        q(1,2) = -l * cos1
        q(2,1) = l * (sin1 + sin2)
        q(2,2) = -l * (cos1 + cos2)
    end subroutine

    ! For θ (x)
    double precision function f_x(v)
        double precision,intent(in) :: v
        f_x = v
    end function
    ! For dθ1/dt (v1)
    double precision function f_v1(x1, x2, v1, v2)
        double precision,intent(in) :: x1, x2, v1, v2
        double precision sin1, sin2, s, c
        double precision numerator, denominator
        sin1 = sin(x1)
        sin2 = sin(x2)
        s = sin(x1-x2);
        c = cos(x1-x2);

        numerator = g*(c*sin2-2d0*sin1)-s*l*(v2**2d0+c*v1**2d0)
        denominator = l*(2d0-c**2d0)
        f_v1 = numerator/denominator
    end function
    ! for dθ2/dt (v2)
    double precision function f_v2(x1, x2, v1, v2)
        double precision,intent(in) :: x1, x2, v1, v2
        double precision sin1, sin2, s, c
        double precision numerator, denominator
        sin1 = sin(x1)
        sin2 = sin(x2)
        s = sin(x1-x2)
        c = cos(x1-x2)

        numerator = 2d0*g*(c*sin1-sin2) + s*l*(2d0*v1**2d0+c*v2**2d0)
        denominator = l * (2d0-c**2d0)
        f_v2 = numerator/denominator
    end function
end module

program main
    use constant
    use double_pendulum
    implicit none
    integer,parameter :: fo = 10
    integer i, n ! n = 10000
    double precision :: theta(2) = (/ 0.5d0*pi, 0.5d0*pi /), theta_d(2) = 0.0d0
    ! q is generalized coordinates
    ! q(x:1 or v:2, 1 or 2)
    double precision :: q(2, 2) = 0d0, dt = 0.0d0, t = 0.0d0, observeTime = 100d0
    ! [caution] if observeTime is too big (or dt is), θ (or dθ/dt) could be extremely a large number, and result would be bunch of NaN
    open(fo, file='output.txt')
    write (*, '(a)', advance='no') "Observe Time(sec) :"
    read (*, *) observeTime
    write (*, '(a)', advance='no') "Division :"
    read (*, *) n
    dt = observeTime/dble(n)
    write (*, *)
    write (*, *) "x1 y1 x2 y2 t"
    ! Output what's initially inside
    call get_position(theta, q)
    write (fo, '(7(f15.6, x))') q(1,1), q(1,2), q(2,1), q(2,2), t, theta_d(1), theta_d(2)

    do i = 1, n
        call rungekutta(theta, theta_d, dt)
        call get_position(theta, q)
        t = t + dt
        write (fo, '(7(f15.6, x))') q(1,1), q(1,2), q(2,1), q(2,2), t, theta_d(1), theta_d(2)
        if (i < 13) then
            write (*, '(5(f10.6, x))') q(1,1), q(1,2), q(2,1), q(2,2), t
        else if (i == 13) then
            write (*, *) "( mid progress was skipped )"
        else if (i > n - 13) then
            write (*, '(5(f10.6, x))') q(1,1), q(1,2), q(2,1), q(2,2), t
        end if
    end do

    close(fo)
end program