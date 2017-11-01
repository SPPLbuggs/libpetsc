    module eqn_lib
    use props
    implicit none
    
    real(8), allocatable :: f(:,:,:)
    
    contains
    
    subroutine feval(g, i, j, n, m, dof, f, b_temp)
    type(grid), intent(in) :: g
    integer, intent(in) :: i, j, n, m, dof
    real(8), intent(in) :: f(n,m,dof)
    real(8), intent(out) :: b_temp(dof)
    
    call feqn1(g, i, j, f(:,:,1), b_temp(1))
    
    end subroutine
    
    subroutine feqn1(g, i, j, f, b)
    type(grid), intent(in) :: g
    integer, intent(in)  :: i, j
    real(8), intent(in)  :: f(:,:)
    real(8), intent(out) :: b
    real(8) :: V, dfdx = 0, dfdy = 0
    
    !v = sin(2 * 3.14159 * g%t)
    V = 1
    
    if (g%nx > 1) then
        if (g%type_x(i-1,j-1) == -1) then
            dfdx = f(i+1,j) - 2 * f(i,j) + V
        else if (g%type_x(i-1,j-1) == 1) then
            dfdx = f(i-1,j) - 2 * f(i,j)
        else
            dfdx = f(i+1,j) - 2 * f(i,j) + f(i-1,j)
        end if
    end if
    
    if (g%ny > 1) then
        if (g%type_y(i-1,j-1) == -1) then
            dfdy = f(i,j+1) - f(i,j)
        else if (g%type_y(i-1,j-1) == 1) then
            dfdy = f(i,j-1) - f(i,j)
        else
            dfdy = f(i,j+1) - 2 * f(i,j) + f(i,j-1)
        end if
    end if
    
    b = dfdx + dfdy
    
    end subroutine
    end module
