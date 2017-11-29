  !==============================================================================================================
  subroutine mafraw(n,nr,rws,nc,cls,af,nmiss,fnRAW)	
  !==============================================================================================================
  
  implicit none
  
  integer*4 :: i,n,nc,nr,rws(nr),cls(nc),raww(nr)
  real*8 :: n0,n1,n2,ntotal,af(nc),nmiss(nc)
  character*1 :: raw(n)
  character(len=1000) :: fnRAW

  af=0.0D0
  nmiss=0
  ntotal=real(n)  
  open (unit=13,file=fnRAW, status='old', form='unformatted', access='direct', recl=n)
  do i=1,nc
   read (unit=13,rec=cls(i)) raw
   raww=ichar(raw(rws))
   nmiss(i)=count(raww==0)
   n0=count(raww==1)
   n1=count(raww==2)
   n2=count(raww==3)
   af(i)=(n1+2*n2)/(2*(ntotal-nmiss(i)))
  enddo 

  end subroutine mafraw
