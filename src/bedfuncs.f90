
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

  use bedfuncs 
  
  implicit none
  
  integer*4 :: n,nr,nc,rws(nr),cls(nc),scaled,nbytes 
  real*8 :: W(nr,nc),gsc(nr)
  character(len=1000) :: fnRAW

  integer*4, parameter :: byte = selected_int_kind(1) 
  integer(byte) :: raw(nbytes)
  integer*4 :: i,stat,g(n)
  
  open(unit=13, file=fnRAW, status='old', access='direct', form='unformatted', recl=nbytes)

  W=0.0D0  
  do i=1,nc 
    read(13, iostat=stat, rec=cls(i)) raw
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
  subroutine qcbed(n,nr,rws,nc,cls,af,nmiss,n0,n1,n2,nbytes,fnRAW)	
  !==============================================================================================================

  use bedfuncs 
  
  implicit none
  
  integer*4 :: n,nr,nc,rws(nr),cls(nc),nbytes,g(n),grws(nr) 
  real*8 :: n0(nc),n1(nc),n2(nc),ntotal,af(nc),nmiss(nc)
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
  subroutine prsbed(n,nr,rws,nc,cls,scaled,nprs,s,prs,nbytes,fnRAW)	
  !==============================================================================================================
  
  use bedfuncs 

  implicit none
  
  integer*4 :: i,j,k,n,nr,nc,rws(nr),cls(nc),scaled,nprs,nbytes,g(n)
  real*8 :: prs(nr,nprs),gsc(nr),s(nc,nprs)
  character(len=1000) :: fnRAW
  integer, parameter :: byte = selected_int_kind(1) 
  integer(byte) :: raw(nbytes)
  integer :: stat

  prs=0.0D0  
  open(unit=13, file=fnRAW, status='old', access='direct', form='unformatted', recl=nbytes)

  do i=1,nc
    read(13, iostat=stat, rec=cls(i)) raw
    if (stat /= 0) exit
    g = raw2int(n,nbytes,raw)
    if (scaled==0) then
      where(g==3) g=0
      do j=1,nprs
        prs(1:nr,j)=prs(1:nr,j) + dble(g(rws))*s(i,j)
      enddo
    endif
    if (scaled==1) then
      gsc=scale(nr,dble(g(rws)))
      do j=1,nprs
        prs(1:nr,j)=prs(1:nr,j)+gsc*s(i,j)
      enddo
    endif
  enddo
  close (unit=13)

  end subroutine prsbed



!==============================================================================================================
  subroutine bedsolve(n,nr,rws,nc,cls,scaled,nbytes,fnRAW,nit,lambda,tol,y,g,e,s,mean,sd)
!==============================================================================================================

  use bedfuncs 

  implicit none
  
  integer*4 :: i,j,k,n,nr,nc,rws(nr),cls(nc),scaled,nbytes,nit,it
  real*8 :: y(n),e(n),raww(n),w(n),g(n)
  real*8 :: dww(nc),s(nc),os(nc),lambda(nc),mean(nc),sd(nc)
  real*8 :: lhs,rhs,snew,tol,sigma
  character(len=1000) :: fnRAW

  integer, parameter :: byte = selected_int_kind(1) 
  integer(byte) :: raw(nbytes)
  integer :: stat

  open(unit=13, file=fnRAW, status='old', access='direct', form='unformatted', recl=nbytes)

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
  dww(i)=dot_product(w(rws),w(rws)) 
  if(s(i).eq.0.0D0) then 
    s(i)=(dot_product(w(rws),y(rws))/dww(i))/nc
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
  rhs=dot_product(w(rws),e(rws)) + dww(i)*s(i)
  snew=rhs/lhs
  
  e(rws)=e(rws) - w(rws)*(snew-s(i))
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
  
  end subroutine bedsolve
