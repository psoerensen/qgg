  !==============================================================================================================
  subroutine prsraw(n,nr,rws,nc,cls,scaled,nprs,s,prs,fnRAW)	
  !==============================================================================================================
  ! http://www.tek-tips.com/viewthread.cfm?qid=890231
  
  implicit none
  
  integer*4 :: i,j,n,nr,nc,rws(nr),cls(nc),scaled,nprs
  real*8 :: prs(nr,nprs),raww(n),wsc(n),s(nc,nprs),nsize,n0,mean,sd
  real*8 :: tol
  character*1 :: raw(n)
  character(len=1000) :: fnRAW

  tol=0.0001D0

  prs=0.0D0  
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
    do j=1,nprs
      prs(1:nr,j)=prs(1:nr,j)+wsc(rws)*s(i,j)
    enddo
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
    do j=1,nprs
      prs(1:nr,j)=prs(1:nr,j)+wsc(rws)*s(i,j)
    enddo
   enddo
  endif
  close (unit=13)

  end subroutine prsraw
