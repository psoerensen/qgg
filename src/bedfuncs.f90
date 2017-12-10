
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
  
  open(13, file=fnRAW, status='old', access='direct', form='unformatted', recl=nbytes)

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

  close(13)

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
  integer :: i, stat, rawbits

  af=0.0D0
  nmiss=0.0D0
  ntotal=real(nr)  
  
  open(13, file=fnRAW, status='old', access='direct', form='unformatted', recl=nbytes)

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

  close(13)

  end subroutine qcbed

  !==============================================================================================================
  subroutine prsbed(n,nr,rws,nc,cls,scaled,nprs,s,prs,nbytes,fnRAW)	
  !==============================================================================================================
  
  use bedfuncs 

  implicit none
  
  integer*4 :: i,j,k,n,nr,nc,rws(nr),cls(nc),scaled,nprs,nbytes,g(n)
  real*8 :: prs(nr,nprs),grws(nr),gsc(nr),s(nc,nprs)
  character(len=1000) :: fnRAW
  integer, parameter :: byte = selected_int_kind(1) 
  integer(byte) :: raw(nbytes)
  integer :: stat, rawbits

  prs=0.0D0  
  open (unit=13,file=fnRAW, status='old', form='unformatted', access='direct', recl=n)
  do i=1,nc
    read(13, iostat=stat, rec=cls(i)) raw
    if (stat /= 0) exit
    g = raw2int(n,nbytes,raw)

    if (scaled==0) then
      where(g==3) g=0
      grws=dble(g(rws))
      do j=1,nprs
        do k=1,nr
          prs(k,j)=prs(k,j)+grws(k)*s(i,j)
        enddo
      enddo
    endif

    if (scaled==1) then
      grws=dble(g(rws))
      gsc=scale(nr,grws)
      do j=1,nprs
        do k=1,nr
          prs(k,j)=prs(k,j)+gsc(k)*s(i,j)
        enddo
      enddo
    endif

  enddo
  close (unit=13)

  end subroutine prsbed
