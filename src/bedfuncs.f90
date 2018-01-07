
    module bedfuncs

    implicit none
    
    contains

    function raw2int(n,nbytes,raw) result(g)
    implicit none

    integer, intent(in) :: nbytes,n
    integer, parameter :: byte = selected_int_kind(1) 
    integer(byte), intent(in) :: raw(nbytes)
    integer*4 :: i,j,k,rawbits,g(n) 
    integer, dimension(4) :: rawcodes
 
    rawcodes = (/ 0, 3, 1, 2 /)
    
    g=0
    k=0
    do i=1,nbytes 
      do j=0,6,2
        k = k + 1
        rawbits = ibits(raw(i), j, 2)
        g(k) = rawcodes(rawbits+1)
        if (k==n) exit 
      enddo
      if (k==n) exit
    enddo

    end function raw2int

    function raw2real(n,nbytes,raw) result(w)
    implicit none

    integer, intent(in) :: nbytes,n
    integer, parameter :: byte = selected_int_kind(1) 
    integer(byte), intent(in) :: raw(nbytes)
    integer*4 :: i,j,k,rawbits
    real*8 :: w(n)
    real*8, dimension(4) :: rawcodes
 
    rawcodes = (/ 0.0D0, 3.0D0, 1.0D0, 2.0D0 /)
    
    w=0.0D0
    k=0
    do i=1,nbytes 
      do j=0,6,2
        k = k + 1
        rawbits = ibits(raw(i), j, 2)
        w(k) = rawcodes(rawbits+1)
        if (k==n) exit 
      enddo
      if (k==n) exit
    enddo

    end function raw2real


    function scale(nr,g) result(w)
    implicit none

    integer, intent(in) :: nr
    real*8, intent(in) :: g(nr)
    real*8 :: mean,sd,tol,nsize,w(nr)

    tol=0.0001D0
    w=g
    nsize=real(count(w<3.0D0))
    mean=sum(w, mask=w<3.0D0)/nsize
    where(w<3.0D0) 
    w=w-mean
    elsewhere
    w=0.0D0
    end where
    sd=sqrt(sum(w**2)/(nsize-1))
    if(sd>tol) w=w/sd
    if(sd<tol) w=0.0D0

    end function scale



    end module bedfuncs


  !==============================================================================================================
  subroutine readbed(n,nr,rws,nc,cls,scaled,W,nbytes,fnRAW)	
  !==============================================================================================================
  ! use call fseek in combination with stream access
  ! this should allow reading bedfiles using recl and direct access 
  use bedfuncs 
  
  implicit none
  
  integer*4 :: n,nr,nc,rws(nr),cls(nc),scaled,nbytes 
  real*8 :: W(nr,nc),gsc(nr)
  character(len=1000) :: fnRAW

  integer*4, parameter :: byte = selected_int_kind(1) 
  integer(byte) :: raw(nbytes)
  integer*4 :: i,stat,g(n)

  integer, parameter :: k14 = selected_int_kind(14) 
  integer (kind=k14) :: pos

  open(unit=13, file=fnRAW, status='old', access='direct', form='unformatted', recl=nbytes)
  !open(unit=13, file=fnRAW, status='old', access='stream', form='unformatted', action='read')

  W=0.0D0  
  do i=1,nc 
    read(13, iostat=stat, rec=cls(i)) raw
    !pos = 1 + (cls(i)-1)*nbytes
    !read(13, pos=pos) raw
    if (stat /= 0) exit
    g = raw2int(n,nbytes,raw)
    if (scaled==0) then
      where(g==3) g=0
      W(1:nr,i) = dble(g(rws))
    endif
    if (scaled==1) then
      gsc=dble(g(rws))
      W(1:nr,i)=scale(nr,gsc)
    endif
  enddo 

  close(unit=13)

  end subroutine readbed




  !==============================================================================================================
  subroutine qcbed(n,nr,rws,nc,cls,af,nmiss,n0,n1,n2,nbytes,fnRAW,ncores)	
  !==============================================================================================================

  use bedfuncs 
  
  implicit none
  
  integer*4 :: n,nr,nc,rws(nr),cls(nc),nbytes,g(n),grws(nr),ncores 
  real*8 :: n0(nc),n1(nc),n2(nc),ntotal,af(nc),nmiss(nc)
  character(len=1000) :: fnRAW

  integer, parameter :: byte = selected_int_kind(1) 
  integer(byte) :: raw(nbytes)
  integer :: i, stat

  integer, parameter :: k14 = selected_int_kind(14) 
  integer (kind=k14) :: pos

  call omp_set_num_threads(ncores)

  af=0.0D0
  nmiss=0.0D0
  ntotal=real(nr)  
  
  open(unit=13, file=fnRAW, status='old', access='direct', form='unformatted', recl=nbytes)
  !open(unit=13, file=fnRAW, status='old', access='stream', form='unformatted', action='read')

  !$omp parallel do
  do i=1,nc 
    read(13, iostat=stat, rec=cls(i)) raw
    if (stat /= 0) exit
    !pos = 1 + (cls(i)-1)*nbytes
    !read(13, pos=pos) raw
    g = raw2int(n,nbytes,raw)
    grws = g(rws)
    nmiss(i)=dble(count(grws==3))
    n0(i)=dble(count(grws==0))
    n1(i)=dble(count(grws==1)) 
    n2(i)=dble(count(grws==2))
    af(i)=(n1(i)+2.0D0*n2(i))/(2.0D0*(ntotal-nmiss(i)))
  enddo 
  close(unit=13)

  end subroutine qcbed

  !==============================================================================================================
  subroutine mafbed(n,nr,rws,nc,cls,af,nbytes,fnRAW)	
  !==============================================================================================================

  use bedfuncs 
  
  implicit none
  
  integer*4 :: n,nr,nc,rws(nr),cls(nc),nbytes,g(n),grws(nr) 
  real*8 :: ntotal,af(nc),nmiss(nc)
  character(len=1000) :: fnRAW

  integer, parameter :: byte = selected_int_kind(1) 
  integer(byte) :: raw(nbytes)
  integer :: i, stat

  af=0.0D0
  nmiss=0.0D0
  ntotal=real(nr)  
  
  open(unit=13, file=fnRAW, status='old', access='direct', form='unformatted', recl=nbytes)

  do i=1,nc 
    read(13, iostat=stat, rec=cls(i)) raw
    if (stat /= 0) exit
    g = raw2int(n,nbytes,raw)
    grws = g(rws)
    nmiss(i)=dble(count(grws==3))
    af(i)=sum(grws, mask=grws<3)/(2.0D0*(ntotal-nmiss(i)))
  enddo 

  close(unit=13)

  end subroutine mafbed


  !==============================================================================================================
  subroutine prsbed(n,nr,rws,nc,cls,scaled,nbytes,fnRAW,msize,ncores,nprs,s,prs)	
  !==============================================================================================================
  ! C = A*B (dimenions: mxn = mxk kxn) 
  ! call dgemm("n","n",m,n,k,1.0d0,a,m,b,k,0.0d0,c,m)
  ! C = A*B + C
  ! call dgemm("n","n",m,n,k,1.0d0,a,m,b,k,1.0d0,c,m)

  use bedfuncs 
  
  implicit none
  
  integer*4 :: i,j,k,n,nr,nc,rws(nr),cls(nc),scaled,nbytes,ncores,msize,nprs 
  real*8 :: W(nr,msize)
  character(len=1000) :: fnRAW
  real*8 :: prs(nr,nprs),s(nc,nprs)

  call omp_set_num_threads(ncores)

  prs = 0.0D0
  W = 0.0D0 

  do i=1,nc,msize
  
  if((i+msize-1)<nc) then 
    call readbed(n,nr,rws,msize,cls(i:(i+msize-1)),scaled,W,nbytes,fnRAW)
    k = msize
    call dgemm("n","n",nr,nprs,k,1.0d0,W,nr,s(i:(i+msize-1),:),k,1.0d0,prs,nr)
  endif
  if((i+msize-1)>=nc) then
    call readbed(n,nr,rws,size(cls(i:nc)),cls(i:nc),scaled,W,nbytes,fnRAW)
    k = size(cls(i:nc))
    call dgemm("n","n",nr,nprs,k,1.0d0,W(:,1:size(cls(i:nc))),nr,s(i:nc,:),k,1.0d0,prs,nr)
  endif  
  
  print*,'Finished block',i

  enddo
  end subroutine prsbed



!==============================================================================================================
  subroutine solvebed(n,nr,rws,nc,cls,scaled,nbytes,fnRAW,ncores,nit,lambda,tol,y,g,e,s,mean,sd)
!==============================================================================================================

  use bedfuncs 

  !implicit none
  
  integer*4 :: i,j,k,n,nr,nc,rws(nr),cls(nc),scaled,nbytes,nit,it,ncores,incx,incy
  real*8 :: y(n),e(n),raww(n),w(n),g(n)
  real*8 :: dww(nc),s(nc),os(nc),lambda(nc),mean(nc),sd(nc)
  real*8 :: lhs,rhs,snew,tol,sigma
  character(len=1000) :: fnRAW

  integer, parameter :: byte = selected_int_kind(1) 
  integer(byte) :: raw(nbytes)
  integer :: stat

  call omp_set_num_threads(ncores)

  open(unit=13, file=fnRAW, status='old', access='direct', form='unformatted', recl=nbytes)

  inx = 1
  incy = 1

  ! genotypes coded 0,1,2,3=missing => where 0,1,2 means 0,1,2 copies of alternative allele 
  do i=1,nc
  read(13, iostat=stat, rec=cls(i)) raw
  if (stat /= 0) exit
  raww = raw2real(n,nbytes,raw)
  where (raww<3.0D0)
  w = (raww-mean(i))/sd(i)
  elsewhere
  w = 0.0D0
  end where
  dww(i)=0.0D0
  !$omp parallel do
  do j=1,nr
    dww(i)=dww(i)+w(rws(j))*w(rws(j))
  enddo  
  !$omp end parallel do
  !dww(i)=dot_product(w(rws),w(rws))
  if(s(i).eq.0.0D0) then
    s(i)=0.0D0  
    !$omp parallel do
    do j=1,nr
      s(i)=s(i)+w(rws(j))*y(rws(j))
    enddo  
    !$omp end parallel do
    s(i)=s(i)/dww(i))/nc
    !s(i)=(dot_product(w(rws),y(rws))/dww(i))/nc
    !s(i)=(ddot(nr,w(rws),incx,y(rws),incy)/dww(i))/nc
  endif     
  enddo
  close (unit=13)
  os=s

  e=0.0D0
  e(rws)=y(rws)

  do it=1,nit
  g=0.0D0
  open (unit=13,file=fnRAW, status='old', form='unformatted', access='direct', recl=nbytes)
  do i=1,nc
  read(13, iostat=stat, rec=cls(i)) raw
  if (stat /= 0) exit
  raww = raw2real(n,nbytes,raw)
  where (raww<3.0D0)
  w = (raww-mean(i))/sd(i)
  elsewhere
  w = 0.0D0
  end where
  lhs=dww(i)+lambda(i)
  !rhs=ddot(nr,w(rws),incx,e(rws),incy) + dww(i)*s(i)
  rhs=dww(i)*s(i)
  !$omp parallel do
  do j=1,nr
    rhs=rhs+w(rws(j))*e(rws(j))
  enddo  
  !$omp end parallel do
  !rhs=dot_product(w(rws),e(rws)) + dww(i)*s(i)
  snew=rhs/lhs
  
  !$omp parallel do
  do j=1,nr
    e(rws(j))=e(rws(j))-w(rws(j))*(snew-s(i))
  enddo  
  !$omp end parallel do
  !e(rws)=e(rws) - w(rws)*(snew-s(i))
  s(i)=snew
  !g=g+w*s(i)
  enddo
  close (unit=13)
  print*,(sum((s-os)**2)/sqrt(dble(nc)))
  if( (sum((s-os)**2)/sqrt(dble(nc)))<tol) exit  
  os=s  
  enddo
  
  nit=it
  tol=sum((s-os)**2)
  
  end subroutine solvebed



      subroutine bedgemm(imax,ncores)
      implicit none
      !integer,parameter::imax=1024*4
      integer :: imax,ncores
      real*8::flop
      real*8, dimension(:,:), allocatable::a,b,c
      !integer::i,j,m,k,n
      real*8::time0,start,finish
      external dgemm      

      allocate(a(imax,imax),b(imax,imax),c(imax,imax)) 

      call omp_set_num_threads(ncores)

      flop=real(imax)*real(imax)*real(imax)*2.0D0
      c(:,:) = 0.0d0
      a(:,:) = 1.0d0
      b(:,:) = 2.0d0
 
      start = time()
      call dgemm("n","n",imax,imax,imax,1.0d0,a,imax,b,imax,0.0d0,c,imax)
      finish = time()
      time0 = finish - start
      print*, "omp time:",time0, flop/time0/1000000000.0D0,"Gflops"
      deallocate(a,b,c)

      end subroutine bedgemm


  !==============================================================================================================
  subroutine grmbed(n,nr,rws,nc,cls,scaled,nbytes,fnRAW,msize,ncores,G)	
  !==============================================================================================================
  use bedfuncs 
  
  implicit none
  
  integer*4 :: i,j,n,nr,nc,rws(nr),cls(nc),scaled,nbytes,ncores,msize 
  real*8 :: G(nr,nr), W(nr,msize), traceG
  character(len=1000) :: fnRAW

  call omp_set_num_threads(ncores)

  G = 0.0D0
  W = 0.0D0 

  do i=1,nc,msize
  
  if((i+msize-1)<nc) then 
    call readbed(n,nr,rws,msize,cls(i:(i+msize-1)),scaled,W,nbytes,fnRAW)
    call dsyrk('u', 'n', nr, msize, 1.0D0, W, nr, 1.0D0, G, nr)
  endif
  if((i+msize-1)>=nc) then
    call readbed(n,nr,rws,size(cls(i:nc)),cls(i:nc),scaled,W,nbytes,fnRAW)
    call dsyrk('u', 'n', nr, size(cls(i:nc)), 1.0D0, W(:,1:size(cls(i:nc))), nr, 1.0D0, G, nr)
  endif  
  
  print*,'Finished block',i

  enddo
 
  traceG = 0.0D0
  do i=1,size(G,1)
      traceG = traceG + G(i,i)
  enddo
  traceG = traceG/real(nr) 
 
  do i=1,size(G,1)
    do j=i,size(G,1)
      G(i,j) = G(i,j)/traceG
      G(j,i) = G(i,j)
    enddo
  enddo

  end subroutine grmbed

