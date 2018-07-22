!==============================================================================================================
! functions and subroutines used for bedfiles 
!==============================================================================================================
!
!	https://www.cog-genomics.org/plink2/formats
!	00 Homozygous for first allele in .bim file
!	01 Missing genotype
!	10 Heterozygous
!	11 Homozygous for second allele in .bim file
    

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
    ! 00 01 10 11

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
    ! 00 01 10 11
    
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

    tol=0.00001D0
    w=g
    !nsize=real(count(w<3.0D0))
    nsize=dble(count(w<3.0D0))
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
  subroutine bed2raw(n,m,cls,nbytes,append,fnBED,fnRAW)	
  !==============================================================================================================

  use bedfuncs 
  
  implicit none
  
  integer*4 :: n,m,cls(m),nbytes,nchar,append  
  character(len=1000) :: fnRAW
  character(len=1000) :: fnBED

  integer*4, parameter :: byte = selected_int_kind(1) 
  integer(byte) :: raw(nbytes), magic(3)
  integer*4 :: stat,i,offset

  integer, parameter :: k14 = selected_int_kind(14) 
  integer (kind=k14) :: pos, nbytes14, offset14,i14

  offset=3
  nchar=index(fnBED, '.bed')
  open(unit=13, file=fnBED(1:(nchar+3)), status='old', access='stream', form='unformatted', action='read')

  nchar=index(fnRAW, '.raw')
  if (append==0) open(unit=14, file=fnRAW(1:(nchar+3)), status='new', access='stream', action='write')
  if (append==1) open(unit=14, file=fnRAW(1:(nchar+3)), status='old', access='stream', action='write',position='append')
  
  nbytes14 = nbytes
  offset14 = offset
  
  read(13) magic
  do i=1,m 
    if(cls(i)==1) then
    i14=i
    pos = 1 + offset14 + (i14-1)*nbytes14
    read(13, pos=pos) raw
      write(14) raw
      print*, 'writing record', i, 'to file' 
    endif
  enddo 

  close(unit=13)
  close(unit=14)

  end subroutine bed2raw



  !==============================================================================================================
  subroutine readbed(n,nr,rws,nc,cls,scaled,W,nbytes,fnRAW)	
  !==============================================================================================================

  use bedfuncs 
  
  implicit none
  
  integer*4 :: n,nr,nc,rws(nr),cls(nc),scaled,nbytes,nchar,offset  
  real*8 :: W(nr,nc),gsc(nr),gr(n),n0,n1,n2,nmiss,af,ntotal
  character(len=1000) :: fnRAW

  integer*4, parameter :: byte = selected_int_kind(1) 
  integer(byte) :: raw(nbytes)
  integer*4 :: i, stat

  integer, parameter :: k14 = selected_int_kind(14) 
  integer (kind=k14) :: pos,nbytes14,offset14,i14

  offset=0
  nchar=index(fnRAW, '.bed')
  if(nchar>0) offset=3
  if(nchar==0) nchar=index(fnRAW, '.raw')

  open(unit=13, file=fnRAW(1:(nchar+3)), status='old', access='direct', form='unformatted', recl=nbytes)

  ntotal=dble(nr)  

  W=0.0D0  

  do i=1,nc
    read(13, iostat=stat, rec=cls(i)) raw
    gr = raw2real(n,nbytes,raw)
    if (scaled==0) then
      where(gr==3.0D0) gr=0.0D0
      W(1:nr,i) = gr(rws)
    endif
    if (scaled==1) then
      gsc=gr(rws)
      W(1:nr,i)=scale(nr,gsc)
    endif
    if (scaled==2) then
      nmiss=dble(count(gr==3.0D0))
      n0=dble(count(gr==0.0D0))
      n1=dble(count(gr==1.0D0)) 
      n2=dble(count(gr==2.0D0))
      af=(n1+2.0D0*n2)/(2.0D0*(ntotal-nmiss))
      W(1:nr,i) = gr(rws)
      where(W(1:nr,i)==0.0D0) W(1:nr,i)=-2.0D0*(af)*(1.0D0-af)
      where(W(1:nr,i)==1.0D0) W(1:nr,i)=1.0D0 - 2.0D0*(af)*(1.0D0-af)
      where(W(1:nr,i)==2.0D0) W(1:nr,i)=-2.0D0*(af)*(1.0D0-af)
      where(W(1:nr,i)==3.0D0) W(1:nr,i)=0.0D0
    endif
  enddo 

  close(unit=13)

  end subroutine readbed


  !==============================================================================================================
  subroutine readbedstream(n,nr,rws,nc,cls,scaled,W,nbytes,fnRAW,ncores)	
  !==============================================================================================================

  use bedfuncs 
  
  implicit none
  
  integer*4 :: n,nr,nc,rws(nr),cls(nc),scaled,nbytes,ncores,nchar,offset  
  real*8 :: W(nr,nc),gsc(nr,ncores),gr(n,ncores),n0,n1,n2,nmiss,af,ntotal
  character(len=1000) :: fnRAW

  integer*4, parameter :: byte = selected_int_kind(1) 
  integer(byte) :: raw(nbytes,nc)
  integer*4 :: i,j
  integer, external :: omp_get_thread_num

  integer, parameter :: k14 = selected_int_kind(14) 
  integer (kind=k14) :: pos(nc),nbytes14,offset14,i14

  call omp_set_num_threads(ncores)

  offset=0
  nchar=index(fnRAW, '.bed')
  if(nchar>0) offset=3
  if(nchar==0) nchar=index(fnRAW, '.raw')
  open(unit=13, file=fnRAW(1:(nchar+3)), status='old', access='stream', form='unformatted', action='read')

  ntotal=dble(nr)  

  nbytes14 = nbytes
  offset14 = offset
  do i=1,nc 
    i14=cls(i)
    pos(i) = 1 + offset14 + (i14-1)*nbytes14
    read(13, pos=pos(i)) raw(1:nbytes,i)
  enddo

  W=0.0D0
  !$omp parallel do private(i,j,nmiss,n0,n1,n2,af)
  do i=1,nc
    j=omp_get_thread_num()+1 
    gr(1:n,j) = raw2real(n,nbytes,raw(1:n,i))

    if (scaled==0) then
      where(gr(1:n,j)==3.0D0) gr(:,j)=0.0D0
      W(1:nr,i) = gr(rws,j)
    endif
    if (scaled==1) then
      gsc(1:nr,j)=gr(rws,j)
      W(1:nr,i)=scale(nr,gsc(1:nr,j))
    endif
    if (scaled==2) then
      nmiss=dble(count(gr(rws,j)==3.0D0))
      n0=dble(count(gr(rws,j)==0.0D0))
      n1=dble(count(gr(rws,j)==1.0D0)) 
      n2=dble(count(gr(rws,j)==2.0D0))
      af=(n1+2.0D0*n2)/(2.0D0*(ntotal-nmiss))
      W(1:nr,i) = gr(rws,j)
      where(W(1:nr,i)==0.0D0) W(1:nr,i)=-2.0D0*(af)*(1.0D0-af)
      where(W(1:nr,i)==1.0D0) W(1:nr,i)=1.0D0 - 2.0D0*(af)*(1.0D0-af)
      where(W(1:nr,i)==2.0D0) W(1:nr,i)=-2.0D0*(af)*(1.0D0-af)
      where(W(1:nr,i)==3.0D0) W(1:nr,i)=0.0D0
    endif


  enddo 
  !$omp end parallel do

  close(unit=13)

  end subroutine readbedstream

  !==============================================================================================================
  subroutine qcbed(n,nr,rws,nc,cls,af,nmiss,n0,n1,n2,nbytes,fnRAW,ncores)	
  !==============================================================================================================

  use bedfuncs 
  
  implicit none
  
  integer*4 :: n,nr,nc,rws(nr),cls(nc),nbytes,g(n),grws(nr),ncores,nchar,offset 
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
  !ntotal=real(nr)  
  ntotal=dble(nr)  
  
  offset=0
  nchar=index(fnRAW, '.bed')
  
  if(nchar>0) offset=3
  if(nchar==0) nchar=index(fnRAW, '.raw')
  !print*,'nchar',fnRAW(1:(nchar+3)),nchar,offset

  open(unit=13, file=fnRAW(1:(nchar+3)), status='old', access='direct', form='unformatted', recl=nbytes)
  !open(unit=13, file=fnRAW(1:(nchar+3)), status='old', access='stream', form='unformatted', action='read')

  do i=1,nc 
    read(13, iostat=stat, rec=cls(i)) raw
    !if (stat /= 0) exit
    !pos = 1 + offset + (cls(i)-1)*nbytes
    !read(13, pos=pos) raw
    g = raw2int(n,nbytes,raw)
    grws = g(rws)
    nmiss(i)=dble(count(grws==3))
    n0(i)=dble(count(grws==0))
    n1(i)=dble(count(grws==1)) 
    n2(i)=dble(count(grws==2))
    if ( nmiss(i)<ntotal ) af(i)=(n1(i)+2.0D0*n2(i))/(2.0D0*(ntotal-nmiss(i)))
    !af(i)=(n1(i)+2.0D0*n2(i))/(2.0D0*(ntotal-nmiss(i)))
  enddo 

  close(unit=13)

  end subroutine qcbed

  !==============================================================================================================
  subroutine mafbed(n,nr,rws,nc,cls,af,nmiss,n0,n1,n2,nbytes,fnRAW,ncores)	
  !==============================================================================================================

  use bedfuncs 
  
  implicit none
  
  integer*4 :: n,nr,nc,rws(nr),cls(nc),nbytes,g(n,ncores),grws(nr,ncores),ncores,nchar,offset 
  real*8 :: n0(nc),n1(nc),n2(nc),ntotal,af(nc),nmiss(nc)
  character(len=1000) :: fnRAW
  integer, external :: omp_get_thread_num

  integer, parameter :: byte = selected_int_kind(1) 
  integer(byte) :: raw(nbytes,nc)
  integer :: i,j

  integer, parameter :: k14 = selected_int_kind(14) 
  integer (kind=k14) :: pos(nc),nbytes14,offset14,i14

  call omp_set_num_threads(ncores)

  offset=0
  nchar=index(fnRAW, '.bed')
  if(nchar>0) offset=3
  if(nchar==0) nchar=index(fnRAW, '.raw')

  open(unit=13, file=fnRAW(1:(nchar+3)), status='old', access='stream', form='unformatted', action='read')

  nbytes14 = nbytes
  offset14 = offset
  do i=1,nc 
    i14=cls(i)
    pos(i) = 1 + offset14 + (i14-1)*nbytes14
    read(13, pos=pos(i)) raw(1:nbytes,i)
  enddo

  af=0.0D0
  nmiss=0.0D0
  ntotal=dble(nr)  

  !$omp parallel do private(i,j)
  do i=1,nc
    j=omp_get_thread_num()+1 
    !read(13, pos=pos(i)) raw(1:nbytes,j)
    !g(1:n,j) = raw2int(n,nbytes,raw(1:n,j))
    g(1:n,j) = raw2int(n,nbytes,raw(1:n,i))
    grws(1:nr,j) = g(rws,j)
    nmiss(i)=dble(count(grws(1:nr,j)==3))
    n0(i)=dble(count(grws(1:nr,j)==0))
    n1(i)=dble(count(grws(1:nr,j)==1)) 
    n2(i)=dble(count(grws(1:nr,j)==2))
    if ( nmiss(i)<ntotal ) af(i)=(n1(i)+2.0D0*n2(i))/(2.0D0*(ntotal-nmiss(i)))
  enddo 
  !$omp end parallel do

  close(unit=13)

  end subroutine mafbed


  !==============================================================================================================
  subroutine mmbed(n,nr,rws,nc,cls,scaled,nbytes,fnRAW,msize,ncores,nprs,s,prs)	
  !==============================================================================================================
  ! C = A*B (dimenions: mxn = mxk kxn) 
  ! call dgemm("n","n",m,n,k,1.0d0,a,m,b,k,0.0d0,c,m)
  ! C = A*B + C
  ! call dgemm("n","n",m,n,k,1.0d0,a,m,b,k,1.0d0,c,m)

  use bedfuncs 
  
  implicit none
  
  integer*4 :: i,j,k,n,nr,nc,rws(nr),cls(nc),scaled,nbytes,ncores,msize,nprs,ncw 
  real*8 :: W(nr,msize)
  character(len=1000) :: fnRAW
  real*8 :: prs(nr,nprs),s(nc,nprs)

  call omp_set_num_threads(ncores)

  prs = 0.0D0
  W = 0.0D0 

  do i=1,nc,msize

    if((i+msize-1)<nc) ncw = size(cls(i:(i+msize-1)))
    if((i+msize-1)>=nc) ncw = size(cls(i:nc))          

    call readbed(n,nr,rws,ncw,cls(i:(i+ncw-1)),scaled,W(:,1:ncw),nbytes,fnRAW)
    call dgemm("n","n",nr,nprs,ncw,1.0d0,W(:,1:ncw),nr,s(i:(i+ncw-1),:),ncw,1.0d0,prs,nr)
 
    print*,'Finished block',i

  enddo
  end subroutine mmbed


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
  subroutine grmbed(n,nr,rws,nc,cls1,cls2,scaled,nbytes,fnRAW,msize,ncores,fnG,gmodel)	
  !==============================================================================================================
  use bedfuncs 
  
  implicit none
  
  integer*4 :: i,j,n,nr,nc,rws(nr),cls1(nc),cls2(nc),scaled,nbytes,ncores,msize,nchar,ncw,gmodel 
  real*8 :: G(nr,nr), W1(nr,msize), traceG
  character(len=1000) :: fnRAW,fnG
  real*8, allocatable :: W2(:,:)

  call omp_set_num_threads(ncores)

  G = 0.0D0
  W1 = 0.0D0
  if(gmodel==3) then 
    allocate(W2(nr,msize))
    W2 = 0.0D0
  endif

  do i=1,nc,msize

    if((i+msize-1)<nc) ncw = size(cls1(i:(i+msize-1)))
    if((i+msize-1)>=nc) ncw = size(cls1(i:nc))          
  
    select case (gmodel)

      case (1) ! additive 
      call readbed(n,nr,rws,ncw,cls1(i:(i+ncw-1)),scaled,W1(:,1:ncw),nbytes,fnRAW)
      call dsyrk('u', 'n', nr, ncw, 1.0D0, W1(:,1:ncw), nr, 1.0D0, G, nr)

      case (2) ! dominance
      scaled=2
      call readbed(n,nr,rws,ncw,cls1(i:(i+ncw-1)),scaled,W1(:,1:ncw),nbytes,fnRAW)
      call dsyrk('u', 'n', nr, ncw, 1.0D0, W1(:,1:ncw), nr, 1.0D0, G, nr)

      case (3) ! epistasis
      call readbed(n,nr,rws,ncw,cls1(i:(i+ncw-1)),scaled,W1(:,1:ncw),nbytes,fnRAW)
      call readbed(n,nr,rws,ncw,cls2(i:(i+ncw-1)),scaled,W2(:,1:ncw),nbytes,fnRAW)
      do j=1,ncw
        W1(:,j) = W1(:,j)*W2(:,j)
      enddo
      call dsyrk('u', 'n', nr, ncw, 1.0D0, W1(:,1:ncw), nr, 1.0D0, G, nr)

    end select
  
    print*,'Finished block',i

  enddo
 
  traceG = 0.0D0
  do i=1,size(G,1)
      traceG = traceG + G(i,i)
  enddo
  !traceG = traceG/real(nr) 
  traceG = traceG/dble(nr) 
 
  do i=1,size(G,1)
    do j=i,size(G,1)
      G(i,j) = G(i,j)/traceG
      G(j,i) = G(i,j)
    enddo
  enddo

  nchar=index(fnG, '.grm')
  open(unit=10, file=fnG(1:(nchar+3)), status='unknown', access='stream', form='unformatted', action='write')
  do j=1,size(G,1)
    write(unit=10) G(1:size(G,1),j)
  enddo
  close(10)

  end subroutine grmbed


!==============================================================================================================
  subroutine solvebed(n,nr,rws,nc,cls,scaled,nbytes,fnRAW,ncores,nit,lambda,tol,y,g,e,s,mean,sd)
!==============================================================================================================

  use bedfuncs 

  !implicit none
  
  integer*4 :: i,j,k,n,nr,nc,rws(nr),cls(nc),scaled,nbytes,nit,it,ncores,nchar,offset
  real*8 :: y(n),e(n),raww(n),w(n),g(n)
  real*8 :: dww(nc),s(nc),os(nc),lambda(nc),mean(nc),sd(nc)
  real*8 :: lhs,rhs,snew,tol,sigma,dots
  character(len=1000) :: fnRAW
  real*8, external  :: ddot

  integer, parameter :: k14 = selected_int_kind(14) 
  integer (kind=k14) :: pos

  integer, parameter :: byte = selected_int_kind(1) 
  integer(byte) :: raw(nbytes)
  integer :: stat

  call omp_set_num_threads(ncores)

  offset=0
  nchar=index(fnRAW, '.bed')
  if(nchar>0) offset=3
  if(nchar==0) nchar=index(fnRAW, '.raw')
  open(unit=13, file=fnRAW(1:(nchar+3)), status='old', access='direct', form='unformatted', recl=nbytes)
  !open(unit=13, file=fnRAW(1:(nchar+3)), status='old', access='stream', form='unformatted', action='read')

  ! genotypes coded 0,1,2,3=missing => where 0,1,2 means 0,1,2 copies of alternative allele 
  do i=1,nc
    read(13, iostat=stat, rec=cls(i)) raw
    !if (stat /= 0) exit
    !pos = 1 + offset + (cls(i)-1)*nbytes
    !read(13, pos=pos) raw

    raww = raw2real(n,nbytes,raw)
    where (raww<3.0D0)
      w = (raww-mean(i))/sd(i)
    elsewhere
      w = 0.0D0
    end where
    dww(i)=0.0D0
    dww(i)=dot_product(w(rws),w(rws))
    if(s(i).eq.0.0D0) then
      s(i)=0.0D0  
      !s(i)=(dot_product(w(rws),y(rws))/dww(i))/nc
      s(i)=(ddot(nr,w(rws),1,y(rws),1)/dww(i))/nc
    endif     
  enddo
  close (unit=13)
  os=s

  e=0.0D0
  e(rws)=y(rws)

  do it=1,nit
    g=0.0D0
    open(unit=13,file=fnRAW(1:(nchar+3)), status='old', form='unformatted', access='direct', recl=nbytes)
    !open(unit=13, file=fnRAW(1:(nchar+3)), status='old', access='stream', form='unformatted', action='read')
    do i=1,nc
      read(13, iostat=stat, rec=cls(i)) raw
      !if (stat /= 0) exit
      !pos = 1 + offset + (cls(i)-1)*nbytes
      !read(13, pos=pos) raw
      raww = raw2real(n,nbytes,raw)
        where (raww<3.0D0)
        w = (raww-mean(i))/sd(i)
      elsewhere
        w = 0.0D0
      end where
      lhs=dww(i)+lambda(i)
      !rhs=ddot(nr,w(rws),1,e(rws),1) 
      !rhs=rhs + dww(i)*s(i)
      !rhs=dot_product(w(rws),e(rws)) 
      !rhs=rhs + dww(i)*s(i)
      dots = 0.0D0
      !$omp parallel &
      !$omp   shared ( w,e ) &
      !$omp   private ( j )
      !$omp do reduction ( + : dots )
      do j=1,nr
        dots = dots + w(rws(j))*e(rws(j))
      end do
      !$omp end do
      !$omp end parallel
      rhs=dww(i)*s(i)+dots

      snew=rhs/lhs
  
      !e(rws)=e(rws) - w(rws)*(snew-s(i))
      !do j=1,nr
      !call daxpy(nr, (snew-s(i)), w(rws(j)), 1, e(rws(j)), 1)
      !enddo
      !$omp parallel do
      do j=1,nr
         e(rws(j))=e(rws(j))-w(rws(j))*(snew-s(i))
      enddo  
      !$omp end parallel do

      s(i)=snew
      g=g+w*s(i)
    enddo
    close (unit=13)
    print*,(sum((s-os)**2)/sqrt(dble(nc)))
    if( (sum((s-os)**2)/sqrt(dble(nc)))<tol) exit  
    os=s  
  enddo
  
  nit=it
  tol=sum((s-os)**2)
  
  end subroutine solvebed





  
!==============================================================================================================
  subroutine mtsolvebed(n,nr,rws,nc,cls,scaled,nbytes,fnRAW,ncores,nit,lambda,tol,nt,y,g,e,s,mean,sd)
!==============================================================================================================

  use bedfuncs 

  !implicit none
  
  integer*4 :: i,j,k,n,nr,nc,nt,t,rws(nr),cls(nc),scaled,nbytes,nit,it,ncores,nchar,offset
  real*8 :: y(n,nt),e(n,nt),g(n,nt),crit(nt)
  real*8 :: raww(n),w(n)
  real*8 :: dww(nc),s(nc,nt),os(nc,nt),lambda(nt),mean(nc),sd(nc)
  real*8 :: lhs(nt),rhs(nt),snew(nt),dots(nt),tol,sigma
  character(len=1000) :: fnRAW
  real*8, external  :: ddot
  integer*4, external :: omp_get_thread_num

  integer, parameter :: k14 = selected_int_kind(14) 
  integer (kind=k14) :: pos(nc),nbytes14,offset14,i14

  integer, parameter :: byte = selected_int_kind(1) 
  integer(byte) :: raw(nbytes)
  integer :: stat

  call omp_set_num_threads(ncores)

  offset=0
  nchar=index(fnRAW, '.bed')
  if(nchar>0) offset=3
  if(nchar==0) nchar=index(fnRAW, '.raw')

  open(unit=13, file=fnRAW(1:(nchar+3)), status='old', access='stream', form='unformatted', action='read')

  nbytes14 = nbytes
  offset14 = offset

  do i=1,nc
    i14=cls(i)
    pos(i) = 1 + offset14 + (i14-1)*nbytes14
  enddo

  ! genotypes coded 0,1,2,3=missing => where 0,1,2 means 0,1,2 copies of alternative allele 
  !!$omp parallel do private(t,i)
  do i=1,nc
    read(13, pos=pos(i)) raw
    !call flush(13)
    raww = raw2real(n,nbytes,raw)
    where (raww<3.0D0)
      w = (raww-mean(i))/sd(i)
    elsewhere
      w = 0.0D0
    end where
    dww(i)=0.0D0
    dww(i)=dot_product(w(rws),w(rws))
    !$omp parallel do private(t)
    do t=1,nt
      if(s(i,t).eq.0.0D0) then
        s(i,t)=(ddot(nr,w(rws),1,y(rws,t),1)/dww(i))/nc
      endif
    enddo     
    !$omp end parallel do
  enddo
  !!$omp end parallel do

  close (unit=13)
  os=s

  e=0.0D0
  e(rws,1:nt)=y(rws,1:nt)

  do it=1,nit
    g=0.0D0
    open(unit=13, file=fnRAW(1:(nchar+3)), status='old', access='stream', form='unformatted', action='read')

    do i=1,nc

     read(13, pos=pos(i)) raw
      raww = raw2real(n,nbytes,raw)
        where (raww<3.0D0)
        w = (raww-mean(i))/sd(i)
      elsewhere
        w = 0.0D0
      end where

      !$omp parallel do private(t,j)
      do t=1,nt
        lhs(t)=dww(i)+lambda(t)
        !dots(t) = dot_product(w(rws),e(rws,t))
        dots(t)=0.0D0
        do j=1,nr
          dots(t) = dots(t) + w(rws(j))*e(rws(j),t)
        end do
        rhs(t)=dww(i)*s(i,t)+dots(t)
        snew(t)=rhs(t)/lhs(t)
        e(rws,t)=e(rws,t)-w(rws)*(snew(t)-s(i,t))
        s(i,t)=snew(t)
        !g(1:n,t)=g(1:n,t)+w*s(i,t)
      enddo
      !$omp end parallel do

    enddo
    close (unit=13)
    do t=1,nt
      crit(t) = sum((s(1:nc,t)-os(1:nc,t))**2)/sqrt(dble(nc))
      print*,crit(t)
    enddo
    os=s  
    if( maxval(crit)<tol ) exit

  enddo
  
  nit=it
  !tol=sum((s-os)**2)
  
  end subroutine mtsolvebed




  !==============================================================================================================
  subroutine ldbed(n,nr,rws,nc,cls,msize,nbytes,fnRAW,ncores,fnLD)	
  !==============================================================================================================

  use bedfuncs 
  
  implicit none
  
  integer*4 :: n,nr,nc,rws(nr),cls(nc),nbytes,ncores,nchar,offset,msize,nrw
  real*8 :: ld(nc,2*msize+1),w1(nr*4,ncores),w2(nr*4,ncores),dots(msize,ncores)
  character(len=1000) :: fnRAW,fnLD
  integer, external :: omp_get_thread_num

  integer, parameter :: byte = selected_int_kind(1) 
  integer(byte) :: raw(nbytes),raww(nr,nc)
  integer :: i,j,k,thread

  integer, parameter :: k14 = selected_int_kind(14) 
  integer (kind=k14) :: pos(nc),nbytes14,offset14,i14

  call omp_set_num_threads(ncores)

  nrw =nr*4

  offset=0
  nchar=index(fnRAW, '.bed')
  if(nchar>0) offset=3
  if(nchar==0) nchar=index(fnRAW, '.raw')

  open(unit=13, file=fnRAW(1:(nchar+3)), status='old', access='stream', form='unformatted', action='read')

  nbytes14 = nbytes
  offset14 = offset
  do i=1,nc 
    i14=cls(i)
    pos(i) = 1 + offset14 + (i14-1)*nbytes14
    read(13, pos=pos(i)) raw(rws)
    raww(1:nr,i)=raw(rws)
  enddo

  ld=0.0D0
  ld(1:nc,msize+1) = 1.0D0
  !$omp parallel do private(i,j,k,thread)
  do i=1,nc
    thread=omp_get_thread_num()+1 
    w1(1:nrw,thread) = raw2real(nrw,nr,raww(rws,i))
    w1(1:nrw,thread)=scale(nrw,w1(1:nrw,thread))
    dots=0.0D0
    do j=1,msize
      k = i+j
      if(k<(nc+1)) then 
        w2(1:nrw,thread) = raw2real(nrw,nr,raww(rws,k))
        w2(1:nrw,thread)=scale(nrw,w2(1:nrw,thread))
        dots(j,thread) = dot_product(w1(1:nrw,thread),w2(1:nrw,thread))/dble(nrw)
        ld(i,msize+1+j)=dots(j,thread)
      endif
    enddo
    dots=0.0D0
    do j=1,msize
      k = i-j
      if(k>1) then 
        w2(1:nrw,thread) = raw2real(nrw,nr,raww(rws,k))
        w2(1:nrw,thread)=scale(nrw,w2(1:nrw,thread))
        dots(j,thread) = dot_product(w1(1:nrw,thread),w2(1:nrw,thread))/dble(nrw)
        ld(i,msize+1-j)=dots(j,thread)
      endif
    enddo
  enddo 
  !$omp end parallel do

  close(unit=13)

  nchar=index(fnLD, '.ld')
  open(unit=14, file=fnLD(1:(nchar+2)), status='unknown', access='stream', form='unformatted', action='write')
  do i=1,nc
    write(unit=14) ld(i,1:(2*msize+1))
  enddo
  close(14)

  end subroutine ldbed




  !==============================================================================================================
  subroutine readbinold(n,nr,rws,w,nbytes,fnBIN)	
  !==============================================================================================================

  implicit none
  
  integer*4 :: n,nr,rws(nr),nbytes,nchar,i  
  real*8 :: w(nr),raw(n)
  character(len=1000) :: fnBIN


  nchar=index(fnBIN, '.bin')
  !open(unit=13, file=fnBIN(1:(nchar+3)), status='old', access='direct', form='unformatted', action='read')

  do i =1,1000 
  open(unit=13, file=fnBIN(1:(nchar+3)), status='old', access='direct', form='unformatted', recl=n*8)
  read(13, rec=1) raw
  close(unit=13)
  w = raw(rws)
  print*,i 
  enddo

  end subroutine readbinold


  !==============================================================================================================
  subroutine readbin(n,nr,rws,nc,cls,W,nbytes,fnBIN)	
  !==============================================================================================================

  implicit none
  
  integer*4 :: n,nr,rws(nr),nc,cls(nc),nbytes,nchar,i,j,k  
  real*8 :: W(nr,nc),raw(n)
  character(len=1000) :: fnBIN

  integer, parameter :: k14 = selected_int_kind(14) 
  integer (kind=k14) :: pos(nc),nbytes14,offset14,i14

  nchar=index(fnBIN, '.bin')
  open(unit=13, file=fnBIN(1:(nchar+3)), status='old', access='stream', form='unformatted', action='read')
  nbytes14 = nbytes*n
  offset14 = 0
  do i=1,nc 
    i14=cls(i)
    pos(i) = 1 + offset14 + (i14-1)*nbytes14
    read(13, pos=pos(i)) raw(1:n)
    W(1:nr,i) = raw(rws)
    !call flush(13)
  enddo

  !nchar=index(fnBIN, '.bin')
  !open(unit=13, file=fnBIN(1:(nchar+3)), status='old', access='direct', form='unformatted', recl=n*8)
  !do i=1,nc 
  !  read(13, rec=cls(i)) raw
  !  W(1:nr,i) = raw(rws)
  !  call flush(13)
  !enddo


  close(unit=13)

  end subroutine readbin



  !==============================================================================================================
  subroutine psets(m,stat,nsets,setstat,msets,p,np,ncores)	
  !==============================================================================================================

  implicit none
  
  integer*4 :: m,nsets,msets(nsets),p(nsets),np,ncores   
  integer*4 :: i,j,k,k1,k2,seed(ncores),maxm,thread,multicore   
  real*8 :: stat(m),setstat(nsets),u,pstat
  integer, external :: omp_get_thread_num

  p=0

  multicore=0
  if (ncores>1) multicore=1
 
  maxm = m - maxval(msets) - 1

  select case (multicore)

  case (0) ! no multicore 

  do i=1,nsets
    do j=1,np
      call random_number(u)
      k1 = 1 + floor(maxm*u)  ! sample: k = n + floor((m+1-n)*u) n, n+1, ..., m-1, m
      k2 = k1+msets(i)-1
      pstat = sum(stat(k1:k2))
      if (pstat < setstat(i)) p(i) = p(i) + 1
    enddo
  enddo   

  case (1) ! multicore 

  call omp_set_num_threads(ncores)

  !$omp parallel do private(i,j,k1,k2,u,pstat)
  do i=1,nsets
    thread=omp_get_thread_num()+1 
    do j=1,np
      call random_number(u)
      k1 = 1 + floor(maxm*u)  ! sample: k = n + floor((m+1-n)*u) n, n+1, ..., m-1, m
      k2 = k1+msets(i)-1
      pstat = sum(stat(k1:k2))
      if (pstat < setstat(i)) p(i) = p(i) + 1
    enddo
  enddo   
  !$omp end parallel do

  end select


  end subroutine psets


!==============================================================================================================
! functions for memory mapping of files 
!==============================================================================================================
! https://www.pgroup.com/userforum/viewtopic.php?p=8139&sid=900dee3e8bacb79da27dc14d5e644cf2


   module mmapfuncs



    interface
    subroutine memcpy(dest, src, n) bind(C,name='memcpy')
    use iso_c_binding
    integer(c_intptr_t), value:: dest
    integer(c_intptr_t), value:: src
    integer(c_size_t), value :: n
    end subroutine memcpy
    end interface

    !interface
    !integer(c_intptr_t) function mmap(addr,len,prot,flags,fildes,off) result(result) bind(c,name='mmap') 
    !use iso_c_binding 
    !integer(c_intptr_t), value :: addr 
    !integer(c_size_t), value :: len
    !integer(c_int), value :: prot 
    !integer(c_int), value :: flags 
    !integer(c_int), value :: fildes 
    !integer(c_size_t), value :: off 
    !end function mmap 
    !end interface

    interface
    integer(c_int) function munmap(addr, len) bind(c,name='munmap')
    use iso_c_binding 
    integer(c_intptr_t), value :: addr 
    integer(c_size_t), value :: len
    end function munmap
    end interface

    interface 
    type(c_ptr) function mmap(addr,len,prot,flags,fildes,off) bind(c,name='mmap') 
    use iso_c_binding 
    integer(c_int), value :: addr 
    integer(c_size_t), value :: len 
    integer(c_int), value :: prot 
    integer(c_int), value :: flags 
    integer(c_int), value :: fildes 
    integer(c_size_t), value :: off 
    end function mmap 
    end interface 

    end module


   
    !==============================================================================================================
    subroutine fmmap(n,m,nr,rws,nc,cls,W,nbytes,fnBIN)	
   !==============================================================================================================

    use mmapfuncs 
    use iso_c_binding 

    implicit none


    type(c_ptr) :: cptr
    integer(c_intptr_t) :: adr
    integer(c_size_t) :: len, off, n_size_t, m_size_t, nbytes_size_t 
    integer,parameter :: prot_read=1	  
    integer,parameter :: map_private=2    
    integer,parameter :: map_share=1    

    integer :: fd,nchar,i,nbytes,null 
    integer :: n,m,nr,rws(nr),nc,cls(nc) 
    real*8, pointer :: x(:) 
    real*8 :: mapx(n) 
    real*8 :: W(nr,nc) 

    integer*4 :: fildes, getfd

    character(len=1000) :: fnBIN

    integer, parameter :: k14 = selected_int_kind(14) 
    integer (kind=k14) :: pos(nc),nbytes14,offset14,i14,k,k1,k2,n14,m14,nc14
    integer (kind=c_long) :: nbytes_c_long,i_c_long, off_c_long

    nchar=index(fnBIN, '.bin')

    open(unit=13, file=fnBIN(1:(nchar+3)), status='old', access='stream', form='unformatted', action='read')

    fd = fnum( unit=13 )

    nbytes14=nbytes
    n14=n
    m14=m
    nc14=nc
    !len = n14*m14*nbytes14
    len = n14*nc14*nbytes14

    off=0
    offset14=minval(cls) 
    off=(offset14-1)*nbytes14

    !len1 = n*nbytes
    !null=0
    !do i = 1,nc
    !off = 0 + (i-1)*n*nbytes
    !adr = mmap(loc(null),len,prot_read,map_private,fd,off)
    !print*,'was here',off,len,len1
    !print*, loc(mapx)
    !call memcpy(loc(mapx), adr, len1)
    !print*,'was here'
    !W(1:nr,i)=mapx(1:n) 
    !print*,'was here'
    !k = munmap(adr, len1)
    !print*,k
    !enddo

    
    off = 0 !+ (cls(i)-1)*n*nbytes
    cptr = mmap(0,len,prot_read,map_share,fd,off) 
    !cptr = mmap(0,len,prot_read,map_private,fd,off) 
    call c_f_pointer(cptr,x,[len]) 
    do i = 1,nc
      !i14 =cls(i)
      i14 =i
      k1=(i14-1)*n14+1
      k2=i14*n14
      mapx(1:n)=x(k1:k2) 
      W(1:nr,i)=mapx(rws) 
    enddo
    
    !nbytes_c_long = nbytes 
    !do i = 1,nc
    !  i_c_long = cls(i) 
    !  off_c_long = 0 + (i_c_long-1)*nbytes_c_long
    !  cptr = mmap(0,len,prot_read,map_private,fd,off_c_long) 
    !  !k1=(i-1)*nr+1
    !  !k2=i*nr
    !  call c_f_pointer(cptr,x,[len]) 
    !  W(1:nr,i)=x(rws) 
    !  !k = munmap(cptr, len)  
    !enddo



    end subroutine fmmap	
