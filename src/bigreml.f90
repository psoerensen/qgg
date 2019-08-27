!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module global
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    module kinds

    implicit none
    
    integer, parameter :: real64 = selected_real_kind(15, 307)
    integer, parameter :: int32 = selected_int_kind(9)

    end module kinds


    module global

    use kinds
    use iso_c_binding
    use f2cio

    implicit none
    public
    real(real64), allocatable :: V(:,:),P(:,:),G(:,:)
    integer(int32), allocatable :: indxg(:)
    integer(int32) :: ng
    integer(kind=c_int64_t) :: pos14, nbytes14, offset14, i14
    integer(c_int):: cfres, nbytes
    character(len=1000, kind=c_char) :: filename
    character(len=20, kind=c_char) :: mode

    
    end module global

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module funcs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    module bigfuncs

    use kinds
    use global
    use iso_c_binding
    use f2cio
 
    implicit none
    
    contains


    function crossprod(a,b) result(c)
    implicit none
    external dgemm
    real(real64), dimension(:,:), intent(in)  :: a,b
    real(real64) :: c(size(a,1),size(b,2))
    call dgemm('n', 'n', size(a,1), size(b,2), size(a,2), 1.0D0, a, size(a,1), b, size(a,2), 0.0D0, c, size(a,1))
    end function crossprod

    function matvec(a,b) result(c)
    implicit none
    external dgemm
    real(real64), dimension(:,:), intent(in)  :: a
    real(real64), dimension(:), intent(in)  :: b
    real(real64) :: c(size(a,1))
    call dgemm('n', 'n', size(a,1), 1, size(a,2), 1.0D0, a, size(a,1), b, size(a,2), 0.0D0, c, size(a,1))
    end function matvec

    function diag(a) result(b)
    implicit none
    real(real64), dimension(:,:), intent(in)  :: a
    real(real64) :: b(size(a,1))
    integer(int32) :: i
    do i=1,size(a,1)
      b(i) = a(i,i)
    enddo
    end function diag

    function inverse(a) result(b)
    implicit none
    external dpotrf, dpotri
    real(real64), dimension(:,:), intent(in)  :: a
    real(real64) :: b(size(a,1),size(a,2))
    integer(int32) :: info,i,j
    info=0
    call dpotrf('U',size(a,1),a,size(a,1),info)             ! cholesky decompostion of a
    call dpotri('U',size(a,1),a,size(a,1),info)             ! inverse of a
    ! copy to lower
    b=a
    do i=1,size(a,1)
      do j=i,size(a,1)
        b(j,i) = b(i,j)
      enddo
    enddo
    end function inverse

 
    function readG(row,fp) result(gr)
    implicit none
    type(c_ptr), value :: fp
    integer(int32), intent(in) :: row
    integer(int32) :: i
    real(real64) :: gr(size(V,1)),grw(ng)
    i14=row
    pos14 = (i14-1)*nbytes14 
    cfres=cseek(fp,pos14,0)            
    cfres=fread_real(grw,8,ng,fp)
    do i=1,size(V,1)
      gr(i) = grw(indxg(i))
    enddo
    end function readG
  

    
    end module bigfuncs
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module subs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    module bigsubs

    use kinds
    use global
    use bigfuncs
    use iso_c_binding
    use f2cio
    
    implicit none
    
    contains
 
    subroutine chol(a)
    implicit none
    external dpotrf
    real(real64), dimension(:,:), intent(inout)  :: a
    integer(int32) :: info,i,j
    info=0
    call dpotrf('U',size(a,1),a,size(a,1),info)             ! cholesky decompostion of a
    ! copy to lower
    do i=1,size(a,1)
      do j=i,size(a,1)
        a(j,i) = a(i,j)
      enddo
    enddo
    end subroutine chol

    subroutine chol2inv(a)
    implicit none
    external dpotri
    real(real64), dimension(:,:), intent(inout)  :: a
    integer(int32) :: info,i,j
    info=0
    call dpotri('U',size(a,1),a,size(a,1),info)             ! inverse of a
    ! copy to lower
    do i=1,size(a,1)
      do j=i,size(a,1)
        a(j,i) = a(i,j)
      enddo
    enddo
    end subroutine chol2inv


    subroutine computeV(weights,fnames)
    implicit none
    integer(int32) :: i,j,r
    real(real64), dimension(:),intent(in) :: weights
    real(real64) :: gr(size(V,1)), grw(ng)
    character(len=*, kind=c_char), dimension(:),intent(in) :: fnames
 
    type(c_ptr):: fp

    do r=1,size(weights)

    if (r<size(weights)) then

    filename = trim(adjustl(fnames(r))) // C_NULL_CHAR
    mode =  'rb' // C_NULL_CHAR
    fp = fopen(filename, mode)

    do i=1,size(V,1)
      gr = readG(indxg(i),fp)
      !print*, i, indxg(i), gr(1:3)
      do j=1,size(V,1)
        !V(i,j)=V(i,j)+grw(indxg(j))*weights(r)
        V(i,j)=V(i,j)+gr(j)*weights(r)
      enddo
    enddo

    cfres=fclose(fp)
    
    else
    
    do j=1,size(V,1)
        V(j,j)=V(j,j)+weights(r)  ! residual
    enddo
    
    endif
    enddo
  
    end subroutine computeV

    
    end module bigsubs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine reml(n,nf,nr,tol,maxit,ncores,fnr,ngr,indx,y,X,theta,ai,b,varb,u,Vy,Py,llik,trPG,trVG)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use global
    use bigsubs
    use bigfuncs
    use iso_c_binding
    use f2cio

    implicit none

    ! input and output variables
    integer(int32) :: n,nf,nr,maxit,ngr,indx(n),ncores,nchar
    real(real64) :: tol
    real(real64)  :: y(n),X(n,nf),theta(nr)
    character(len=1000) :: rfnames(nr-1)
    character(len=14) :: fnr
    
    ! local variables
    integer(int32) :: i,j,it
    real(real64) :: theta0(nr),ai(nr,nr),s(nr),trPG(nr),trVG(nr),delta(nr),b(nf),u(n,nr)
    real(real64) :: VX(n,nf),XVX(nf,nf),VXXVX(n,nf),Vy(n),Py(n),Pu(n,nr),gr(n,nr-1),varb(nf,nf) 
    real(real64) :: llik, ldV, ldXVX, yPy

    type(c_ptr):: fileunit(nr-1), fp
    
    logical :: exst
    
    ! allocate variables
    allocate(indxg(n))
 
    indxg=indx   
    ng=ngr

    nbytes14 = 8.0d0*dble(ng) 

    filename = 'param.qgg' // C_NULL_CHAR
    mode =  'r' // C_NULL_CHAR
    fp = fopen(filename, mode)
    do i=1,nr-1
      cfres=fgets_char(filename,100,fp)
      nchar=index(filename, '.grm')
      rfnames(i) = filename(2:(nchar+3)) 
    enddo
    cfres=fclose(fp)

    call omp_set_num_threads(ncores)

    do it = 1, maxit
    
    ! compute V and save to file
    allocate(V(n,n))
    V=0.0D0

    call computeV(theta,rfnames)

    ! compute inverse of V (store in V) and log determinant of V
    ldV=0.0D0
    call chol(V)                ! cholesky decomposition of V
    ldV = sum(log(diag(V)**2))  ! log determinant of V
    call chol2inv(V)            ! inverse V using cholesky

    ! compute P and save to file
    VX = crossprod(V,X)                   ! n*nf = n*n n*nf 
    XVX = crossprod(transpose(X),VX)
    call chol(XVX)                        ! cholesky decomposition of XVX
    ldXVX = sum(log(diag(XVX)**2))        ! log determinant of XVX
    call chol2inv(XVX)                    ! inverse XVX using cholesky
    VXXVX = crossprod(VX,XVX)             ! n*nf = n*nf nf*nf
    
    b=matmul(transpose(VXXVX),y)
    
    Vy = matvec(V,y)

    ! compute P (stored in V), trVG and trPG
    trVG = 0.0D0
    trPG = 0.0D0
    trVG(nr) = sum(diag(V))  ! residual effects

    do i=1,nr-1
      filename = trim(adjustl(rfnames(i))) // C_NULL_CHAR
      mode =  'rb' // C_NULL_CHAR
      fileunit(i) = fopen(filename, mode)
    enddo

    do j=1,n
      do i=1,nr-1            ! random effects excluding residual
        gr(1:n,i) = readG(j,fileunit(i))
        trVG(i) = trVG(i) + sum(gr(1:n,i)*V(1:n,j)) 
      enddo
      V(1:n,j) = V(1:n,j) - matvec(VXXVX(:,1:nf),VX(j,1:nf))  ! update the j'th column of V to be the j'th column of P 
      do i=1,nr-1
        trPG(i) = trPG(i) + sum(gr(1:n,i)*V(1:n,j))
      enddo 
    enddo
    trPG(nr) = sum(diag(V))   ! residual effects
 
    ! compute Py and yPy (P stored in V)
    Py = matvec(V,y)
    yPy = sum(y*Py)
    
    ! compute u (unscaled)
    do i=1,n
      do j=1,nr-1           ! random effects excluding residual
        gr(1:n,j) = readG(i,fileunit(j))
        u(i,j) = sum(gr(1:n,j)*Py)
      enddo
    enddo
    u(:,nr) = Py            ! random residual effects

    ! compute Pu  (P stored in V)
    do i=1,nr                     ! random effects including residual
      Pu(:,i) = matvec(V,u(:,i))
    enddo

    do i=1,nr-1
      cfres=fclose(fileunit(i))
    enddo

    deallocate(V) 
    
    ! compute average information and first derivatives
    ai=0.0D0
    s=0.0D0
    do i=1,nr
      do j=i,nr
         ai(i,j) = 0.5D0*sum(u(:,i)*Pu(:,j))
         ai(j,i) = ai(i,j)
      enddo
      if (i<nr) then 
        s(i) = -0.5*( trPG(i) - sum(u(:,i)*Py) )
      else
        s(i) = -0.5*( trPG(i) - sum(Py*Py) )
      endif
    enddo

    ! compute u (scaled)
    do i=1,nr-1             ! random effects excluding residual
      u(:,i) = u(:,i)*theta(i)
    enddo
    
    ! compute theta and asd
    ai = inverse(ai)
    theta0 = theta + matmul(ai,s)
    where (theta0<tol) theta0 = tol
    !where (abs(theta0)<tol) theta0 = tol   ! check this if theta include covariances 
    delta = theta - theta0
    theta = theta0

    ! compute restricted log likelihood
    llik = -0.5D0*( ldV + ldXVX + yPy )
    
    !print *, theta

    if (it.eq.maxit) exit
    if ( maxval(abs(delta))<tol ) exit
    
    enddo

    ! scale residuals
    u(1:n,nr) = u(1:n,nr)*theta(nr)
    varb = XVX
 
    deallocate(indxg) 
 
    end subroutine reml
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
