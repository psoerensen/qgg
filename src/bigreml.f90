!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module global
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    module kinds

    implicit none
    
    integer, parameter :: real64 = selected_real_kind(15, 307)
    integer, parameter :: int32 = selected_int_kind(9)

    end module kinds


  module f2ciobig
  
  use iso_c_binding

  implicit none
  private
  public :: fopen, fclose, fcread, fcwrite, fwrite_real, fread_real, fgets_char, cseek 

     
  interface

     function fopen(filename, mode) bind(C,name='fopen')
       !filename: file name to associate the file stream to 
       !mode: null-terminated character string determining file access mode 
       import
       implicit none
       type(c_ptr) fopen
       character(kind=c_char), intent(in) :: filename(*)
       character(kind=c_char), intent(in) :: mode(*)
     end function fopen

     function fclose(fp) bind(C,name='fclose')
       !fp: the file stream to close 
       import
       implicit none
       integer(c_int) fclose
       type(c_ptr), value :: fp
     end function fclose
     
     function fcread(buffer,size,nbytes,fp) bind(C,name='fread')
       ! buffer: pointer to the array where the read objects are stored 
       ! size: size of each object in bytes 
       ! count: the number of the objects to be read 
       ! fp: the stream to read 
       import
       implicit none
       integer(c_int) fcread
       integer(kind=c_int), value :: size
       integer(kind=c_int), value :: nbytes
       integer(kind=c_int8_t), dimension(nbytes) :: buffer 
       type(c_ptr), value :: fp
     end function fcread
     
     function cseek(fp,offset,origin) bind(C,name='fseek')
       !fp: file stream to modify 
       !offset: number of characters to shift the position relative to origin 
       !origin: position to which offset is added (SEEK_SET, SEEK_CUR, SEEK_END) 
       import
       implicit none
       integer(c_int) cseek
       type(c_ptr), value :: fp
       integer(kind=c_int64_t), value :: offset
       integer(kind=c_int), value :: origin
     end function cseek
     
     function fcwrite(buffer,size,nbytes,fp) bind(C,name='fwrite')
       ! buffer: pointer to the array where the write objects are stored 
       ! size: size of each object in bytes 
       ! count: the number of the objects to be written 
       ! fp: the stream to write 
       import
       implicit none
       integer(c_int) fcwrite
       integer(kind=c_int), value :: size
       integer(kind=c_int), value :: nbytes
       integer(kind=c_int8_t), dimension(nbytes) :: buffer 
       type(c_ptr), value :: fp
     end function fcwrite

     function fwrite_real(buffer,size,nbytes,fp) bind(C,name='fwrite')
       ! buffer: pointer to the array where the write objects are stored 
       ! size: size of each object in bytes 
       ! count: the number of the objects to be written 
       ! fp: the stream to write 
       import
       implicit none
       integer(c_int) fwrite_real
       integer(kind=c_int), value :: size
       integer(kind=c_int), value :: nbytes
       real(c_double), dimension(nbytes) :: buffer 
       type(c_ptr), value :: fp
     end function fwrite_real

     function fread_real(buffer,size,nbytes,fp) bind(C,name='fread')
       ! buffer: pointer to the array where the read objects are stored 
       ! size: size of each object in bytes 
       ! count: the number of the objects to be read 
       ! fp: the stream to read 
       import
       implicit none
       integer(c_int) fread_real
       integer(kind=c_int), value :: size
       integer(kind=c_int), value :: nbytes
       real(kind=c_double), dimension(nbytes) :: buffer 
       type(c_ptr), value :: fp
     end function fread_real

     function fgets_char(buffer,nbytes,fp) bind(C,name='fgets')
       ! buffer: pointer to the array where the read objects are stored 
       ! count: the number of the objects to be read 
       ! fp: the stream to read 
       import
       implicit none
       integer(c_int) fgets_char
       integer(kind=c_int), value :: nbytes
       character(kind=c_char) :: buffer(1000) 
       type(c_ptr), value :: fp
     end function fgets_char



  end interface
     
  end module f2ciobig


    module global

    use kinds
    use iso_c_binding
    use f2ciobig

    implicit none
    public
    real(c_double), allocatable :: V(:,:),P(:,:),G(:,:)
    integer(c_int), allocatable :: indxg(:)
    integer(c_int) :: ng
    integer(kind=c_int64_t) :: pos14, nbytes14, offset14, i14
    integer(c_int):: cfres, nbytes

    
    end module global

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module funcs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    module bigfuncs

    use kinds
    use global
    use iso_c_binding
    use f2ciobig
 
    implicit none
    
    contains


    function crossprod(a,b) result(c)
    implicit none
    external dgemm
    real(c_double), dimension(:,:), intent(in)  :: a,b
    real(c_double) :: c(size(a,1),size(b,2))
    call dgemm('n', 'n', size(a,1), size(b,2), size(a,2), 1.0D0, a, size(a,1), b, size(a,2), 0.0D0, c, size(a,1))
    end function crossprod

    function matvec(a,b) result(c)
    implicit none
    external dgemm
    real(c_double), dimension(:,:), intent(in)  :: a
    real(c_double), dimension(:), intent(in)  :: b
    real(c_double) :: c(size(a,1))
    call dgemm('n', 'n', size(a,1), 1, size(a,2), 1.0D0, a, size(a,1), b, size(a,2), 0.0D0, c, size(a,1))
    end function matvec

    function diag(a) result(b)
    implicit none
    real(c_double), dimension(:,:), intent(in)  :: a
    real(c_double) :: b(size(a,1))
    integer(c_int) :: i
    do i=1,size(a,1)
      b(i) = a(i,i)
    enddo
    end function diag

    function inverse(a) result(b)
    implicit none
    external dpotrf, dpotri
    real(c_double), dimension(:,:), intent(in)  :: a
    real(c_double) :: b(size(a,1),size(a,2))
    integer(c_int) :: info,i,j
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
    integer(c_int), intent(in) :: row
    integer(c_int) :: i
    real(c_double) :: gr(size(V,1)),grw(ng)
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
    use f2ciobig
    
    implicit none
    
    contains
 
    subroutine chol(a)
    implicit none
    external dpotrf
    real(c_double), dimension(:,:), intent(inout)  :: a
    integer(c_int) :: info,i,j
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
    real(c_double), dimension(:,:), intent(inout)  :: a
    integer(c_int) :: info,i,j
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
    integer(c_int) :: i,j,r
    real(c_double), dimension(:),intent(in) :: weights
    real(c_double) :: gr(size(V,1))
    character(len=*, kind=c_char), dimension(:),intent(in) :: fnames
    character(len=1000, kind=c_char) :: filename4
    character(len=20, kind=c_char) :: mode

    type(c_ptr):: fp

    do r=1,size(weights)

    if (r<size(weights)) then

    filename4 = trim(adjustl(fnames(r))) // C_NULL_CHAR
    mode =  'rb' // C_NULL_CHAR
    fp = fopen(filename4, mode)

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
!    subroutine reml(n,nf,nr,tol,maxit,ncores,fnr,ngr,indx,y,X,theta,ai,b,varb,u,Vy,Py,llik,trPG,trVG)
    subroutine reml(n,nf,nr,tol,maxit,ncores,ngr,indx,y,X,theta,ai,b,varb,u,Vy,Py,llik,trPG,trVG,ncharsg,fnGCHAR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use global
    use bigsubs
    use bigfuncs
    use iso_c_binding
    use f2ciobig

    implicit none

    ! input and output variables
    integer(c_int) :: n,nf,nr,maxit,ngr,indx(n),ncores,nchar,ncharsg(nr-1),fnGCHAR(nr-1,1000)
    real(c_double) :: tol
    real(c_double)  :: y(n),X(n,nf),theta(nr)
    character(len=1000, kind=c_char)::  rfnames(nr-1)
    character(len=1000, kind=c_char) :: filename2,filename3
    character(len=20, kind=c_char) :: mode

    
    ! local variables
    integer(c_int) :: i,j,it
    real(c_double) :: theta0(nr),ai(nr,nr),s(nr),trPG(nr),trVG(nr),delta(nr),b(nf),u(n,nr)
    real(c_double) :: VX(n,nf),XVX(nf,nf),VXXVX(n,nf),Vy(n),Py(n),Pu(n,nr),gr(n,nr-1),varb(nf,nf) 
    real(c_double) :: llik, ldV, ldXVX, yPy

    type(c_ptr):: fileunit(nr-1)
    
    !if (c_double /= kind(1.0d0))  &
    !error stop 'Default REAL isn''t interoperable with FLOAT!!'

    ! allocate variables
    allocate(indxg(n))
 
    indxg=indx   
    ng=ngr

    nbytes14 = int(8.0d0*dble(ng),kind=c_int64_t) 

    !filename1 = 'param.qgg' // C_NULL_CHAR
    !mode =  'r' // C_NULL_CHAR
    !fp = fopen(filename1, mode)
    !do i=1,nr-1
    !  cfres=fgets_char(filename2,1000,fp)
    !  nchar=index(filename2, '.grm')
    !  rfnames(i) = filename2(1:(nchar+3)) 
    !enddo
    !cfres=fclose(fp)

    do i=1,nr-1
      nchar = ncharsg(i)
      do j=1,nchar
        filename2(j:j) = char(fnGCHAR(i,j))
      enddo
      rfnames(i) = filename2(1:nchar)
    enddo

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
      filename3 = trim(adjustl(rfnames(i))) // C_NULL_CHAR
      mode =  'rb' // C_NULL_CHAR
      fileunit(i) = fopen(filename3, mode)
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
