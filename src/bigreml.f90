!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module global
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    module global

    implicit none
    public
    real*8, allocatable :: V(:,:),P(:,:),G(:,:)
    integer*4, allocatable :: indxg(:)
    integer*4 :: ng
    
    end module global

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module funcs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    module bigfuncs

    use global
    
    implicit none
    
    contains

    function crossprod(a,b) result(c)
    implicit none
    external dgemm
    real*8, dimension(:,:), intent(in)  :: a,b
    real*8 :: c(size(a,1),size(b,2))
    call dgemm('n', 'n', size(a,1), size(b,2), size(a,2), 1.0D0, a, size(a,1), b, size(a,2), 0.0D0, c, size(a,1))
    end function crossprod

    function matvec(a,b) result(c)
    implicit none
    external dgemm
    real*8, dimension(:,:), intent(in)  :: a
    real*8, dimension(:), intent(in)  :: b
    real*8 :: c(size(a,1))
    call dgemm('n', 'n', size(a,1), 1, size(a,2), 1.0D0, a, size(a,1), b, size(a,2), 0.0D0, c, size(a,1))
    end function matvec

    function diag(a) result(b)
    implicit none
    real*8, dimension(:,:), intent(in)  :: a
    real*8 :: b(size(a,1))
    integer*4 :: i
    do i=1,size(a,1)
      b(i) = a(i,i)
    enddo
    end function diag

    function inverse(a) result(b)
    implicit none
    external dpotrf, dpotri
    real*8, dimension(:,:), intent(in)  :: a
    real*8 :: b(size(a,1),size(a,2))
    integer*4 :: info,i,j
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

    function readG(row,funit) result(gr)
    !function readG(row,fname) result(gr)
    implicit none
    !character(len=*), intent(in) :: fname
    integer*4, intent(in) :: row,funit
    integer*4 :: i
    real*8 :: gr(size(V,1)),grw(ng)
    !open (unit=12,file=trim(adjustl(fname)) , status='old', form='unformatted', access='direct', recl=8*ng)
    !read (unit=12, rec=indxg(row)) grw
    read (unit=funit, rec=indxg(row)) grw
    !close (unit=12)
    do i=1,size(V,1)
      gr(i) = grw(indxg(i))
    enddo

    end function readG


  
    
    end module bigfuncs
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module subs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    module bigsubs

    use global
    !use mymkl
    
    implicit none
    
    contains
 
    subroutine chol(a)
    implicit none
    external dpotrf
    real*8, dimension(:,:), intent(inout)  :: a
    integer*4 :: info,i,j
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
    real*8, dimension(:,:), intent(inout)  :: a
    integer*4 :: info,i,j
    info=0
    call dpotri('U',size(a,1),a,size(a,1),info)             ! inverse of a
    ! copy to lower
    do i=1,size(a,1)
      do j=i,size(a,1)
        a(j,i) = a(i,j)
      enddo
    enddo
    end subroutine chol2inv

    subroutine loadG(fname)
    implicit none
    integer*4 :: i,j,k
    real*8 :: grw(ng)
    character(len=*), intent(in) :: fname
    logical :: exst

    inquire(file=trim(adjustl(fname)), exist=exst)
    if(.not.(exst)) then
       print *, 'Trying to open file:'
       print*, fname
       print *, 'file does not exist'
       stop
    endif
    
    open (unit=12,file=trim(adjustl(fname)), status='old', form='unformatted', access='direct', recl=8*ng)
    do i=1,size(G,1)
      read (unit=12, rec=indxg(i)) grw
      do j=1,size(G,1)
        G(i,j) = grw(indxg(j))
      enddo
    enddo
    close (unit=12)
    end subroutine loadG

 
    subroutine computeV(weights,fnames)
    implicit none
    integer*4 :: i,j,k,r
    real*8, dimension(:),intent(in) :: weights
    real*8 :: grw(ng)
    character(len=*), dimension(:),intent(in) :: fnames
    logical :: exst

    do r=1,size(weights)
    if (r<size(weights)) then
    inquire(file=trim(adjustl(fnames(r))), exist=exst)
    if(.not.(exst)) then
       print *, 'Trying to open file:'
       print*, fnames(r)
       print *, 'file does not exist'
       stop
    endif

    open (unit=12,file=trim(adjustl(fnames(r))), status='old', form='unformatted', access='direct', recl=8*ng)
    do i=1,size(V,1)
      read (unit=12,rec=indxg(i)) grw
      do j=1,size(V,1)
        V(i,j)=V(i,j)+grw(indxg(j))*weights(r)
      enddo
    enddo
    close (unit=12)

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
    !subroutine reml(n,nf,nr,tol,maxit,ncores,ngr,indx,y,X,theta,ai,b,varb,u,Vy,Py,llik,trPG,trVG)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use global
    !use mymkl
    use bigsubs
    use bigfuncs

    implicit none

    ! input and output variables
    integer*4 :: n,nf,nr,maxit,ngr,indx(n),ncores
    real*8 :: tol
    real*8  :: y(n),X(n,nf),theta(nr)
    character(len=1000) :: rfnames(nr-1)
    character(len=14) :: fnr
    
    ! local variables
    integer*4 :: i,j,k,it,row,fileunit(nr)
    real*8 :: theta0(nr),ai(nr,nr),s(nr),trPG(nr),trVG(nr),delta(nr),b(nf),u(n,nr),asd(nr,nr)
    real*8 :: VX(n,nf),XVX(nf,nf),VXXVX(n,nf),Vy(n),Py(n),Pu(n,nr),gr(n,nr-1),varb(nf,nf) 
    real*8 :: llik, ldV, ldXVX, yPy, ymean
    
    logical :: exst
    
    ! allocate variables
    allocate(indxg(n))
 
    indxg=indx   
    ng=ngr

    open (unit=10, file=trim(adjustl(fnr)), status='old')

    !read G filenames and check they exist
    do i=1,nr-1
    read(unit=10,fmt=*) rfnames(i)
    inquire(file=trim(adjustl(rfnames(i))), exist=exst)
    if(.not.(exst)) then
       print *, 'Trying to open file:'
       print*, rfnames(i)
       print *, 'file does not exist'
       stop
    endif
    enddo

    close(unit=10, status='delete')

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
    fileunit(i)=10+i
    open (unit=fileunit(i),file=trim(adjustl(rfnames(i))) , status='old', form='unformatted', access='direct', recl=8*ng)
    enddo

    do j=1,n
      do i=1,nr-1            ! random effects excluding residual
        gr(1:n,i) = readG(j,fileunit(i))
        !gr(1:n,i) = readG(j,rfnames(i))
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
        !gr(1:n,j) = readG(i,rfnames(j))
        u(i,j) = sum(gr(1:n,j)*Py)
      enddo 
      !u(i,1:(nr-1)) = matvec(transpose(gr),Py)
    enddo
    u(:,nr) = Py            ! random residual effects

    ! compute Pu  (P stored in V)
    do i=1,nr                     ! random effects including residual
      Pu(:,i) = matvec(V,u(:,i))
    enddo

    do i=1,nr-1
    close(fileunit(i))
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
    
    print *, theta

    if (it.eq.maxit) exit
    if ( maxval(abs(delta))<tol ) exit
    
    enddo

    ! scale residuals
    u(1:n,nr) = u(1:n,nr)*theta(nr)
    varb = XVX
 
    deallocate(indxg) 
 
    end subroutine reml
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
