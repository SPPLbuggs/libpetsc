    module eqn_lib
    use props
    implicit none
    
    real(8), allocatable :: f(:,:,:), f_mi(:,:,:)
    
    contains
    
    subroutine feval(g, i, j, n, m, dof, f, b_temp)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j, n, m, dof
    real(8), intent(in) :: f(n,m,dof)
    real(8), intent(out) :: b_temp(dof)
    
    call feqn1(g, i, j, f(:,:,1), f(:,:,2), b_temp(1))
    
    if (g%dof > 1) call feqn2(g, i, j, f(:,:,1), f(:,:,2), f_mi(:,:,2), b_temp(2))
    
    end subroutine
    
    subroutine feqn1(g, i, j, f1, f2, b)
    type(grid), intent(in) :: g
    integer, intent(in)  :: i, j
    real(8), intent(in)  :: f1(:,:), f2(:,:)
    real(8), intent(out) :: b
    real(8) :: dfdx = 0, dfdy = 0
    
    if (g%nx > 1) then
        dfdx = ((f1(i+1,j) - f1(i,j)) / g%dx(i) &
               -(f1(i,j) - f1(i-1,j)) / g%dx(i-1)) &
               / g%dlx(i-1)
    end if
    
    if (g%ny > 1) then
        dfdy = ((f1(i,j+1) - f1(i,j)) / g%dy(j) &
               -(f1(i,j) - f1(i,j-1)) / g%dy(j-1)) &
               / g%dly(j-1)
    end if
    
    b = dfdx + dfdy + (f2(1,1) - f2(i,j)) * 100.
    
    end subroutine
    
    subroutine feqn2(g, i, j, f1, f2, f2_mi, b)
    type(grid), intent(in) :: g
    integer, intent(in)  :: i, j
    real(8), intent(in)  :: f1(:,:), f2(:,:), f2_mi(:,:)
    real(8), intent(out) :: b
    real(8) :: dfdx = 0, dfdy = 0, vx(2) = 0, vy(2) = 0, mu, Te
    
    mu = 0.1
    Te  = 0.5
    
    if (g%nx > 1) then
        if (g%type_x(i-1,j-1) == -1) then
            vx(1) = 0
            vx(2) = 0.5 * mu * (f2(i+1,j) + f2(i,j)) * (f1(i+1,j) - f1(i,j) &
                    - Te * log(f2(i+1,j) / f2(i,j))) / g%dx(i)
        else if (g%type_x(i-1,j-1) == 1) then
            vx(1) = 0.5 * mu * (f2(i,j) + f2(i-1,j)) * (f1(i,j) - f1(i-1,j) &
                    - Te * log(f2(i,j) / f2(i-1,j))) / g%dx(i-1)
            vx(2) = 0
        else
            vx(1) = 0.5 * mu * (f2(i,j) + f2(i-1,j)) * (f1(i,j) - f1(i-1,j) &
                    - Te * log(f2(i,j) / f2(i-1,j))) / g%dx(i-1)
            vx(2) = 0.5 * mu * (f2(i+1,j) + f2(i,j)) * (f1(i+1,j) - f1(i,j) &
                    - Te * log(f2(i+1,j) / f2(i,j))) / g%dx(i)
        end if
        
        dfdx = (vx(2) - vx(1)) / g%dlx(i-1)
    end if
    
    if (g%ny > 1) then
        if (g%type_y(i-1,j-1) == -1) then
            vy(1) = 0
            vy(2) = 0.5 * mu * (f2(i,j+1) + f2(i,j)) * (f1(i,j+1) - f1(i,j) &
                    - Te * log(f2(i,j+1) / f2(i,j))) / g%dy(j)
        else if (g%type_y(i-1,j-1) == 1) then
            vy(1) = 0.5 * mu * (f2(i,j) + f2(i,j-1)) * (f1(i,j) - f1(i,j-1) &
                    - Te * log(f2(i,j) / f2(i,j-1))) / g%dy(j-1)
            vy(2) = 0
        else
            vy(1) = 0.5 * mu * (f2(i,j) + f2(i,j-1)) * (f1(i,j) - f1(i,j-1) &
                    - Te * log(f2(i,j) / f2(i,j-1))) / g%dy(j-1)
            vy(2) = 0.5 * mu * (f2(i,j+1) + f2(i,j)) * (f1(i,j+1) - f1(i,j) &
                    - Te * log(f2(i,j+1) / f2(i,j))) / g%dy(j)
        end if
        
        dfdy = (vy(2) - vy(1)) / g%dly(i-1)
    end if
    
    b = f2(i,j) - f2_mi(i,j) + g%dt * (dfdx + dfdy)
    
    end subroutine
    end module
    
    
    
    
    
