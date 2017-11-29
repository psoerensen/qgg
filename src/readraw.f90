  !==============================================================================================================
  subroutine readraw(n,nc,cls,scaled,W,fnRAW)	
  !==============================================================================================================
  ! http://www.tek-tips.com/viewthread.cfm?qid=890231
  
  implicit none
  
  integer*4 :: i,n,nc,cls(nc),scaled
  real*8 :: W(n,nc),raww(n),wsc(n),nsize,n0,mean,sd
  real*8 :: tol
  character*1 :: raw(n)
  character(len=1000) :: fnRAW

  tol=0.0001D0

  W=0.0D0  
  open (unit=13,file=fnRAW, status='old', form='unformatted', access='direct', recl=n)
  if (scaled==0) then
   do i=1,nc
    read (unit=13,rec=cls(i)) raw
    raww=ichar(raw)
    where (raww.ne.0.0D0)
      wsc=raww-1.0D0
    elsewhere
      wsc = 0.0D0
    end where
    W(1:n,i)=wsc
   enddo
  endif
  if (scaled==1) then
   do i=1,nc
    read (unit=13,rec=cls(i)) raw
    raww=ichar(raw)
    nsize=count(raww/=0.0D0)
    mean=sum(raww)/nsize
    where (raww.ne.0.0D0)
      wsc=(raww-mean)
    elsewhere
      wsc = 0.0D0
    end where
    sd=sqrt(sum(wsc**2)/nsize)
    where (wsc.ne.0.0D0)
      wsc=wsc/sd
    elsewhere
      wsc = 0.0D0
    end where
    if(sd<tol) wsc=0.0D0
    W(1:n,i)=wsc
   enddo
  endif
  close (unit=13)

  end subroutine readraw
