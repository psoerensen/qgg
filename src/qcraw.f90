  !==============================================================================================================
  subroutine qcraw(n,nr,rws,nc,cls,af,nmiss,n0,n1,n2,fnRAW)	
  !==============================================================================================================
  
  implicit none
  
  integer*4 :: i,n,nc,nr,rws(nr),cls(nc),raww(nr)
  real*8 :: n0(nc),n1(nc),n2(nc),ntotal,af(nc),nmiss(nc)
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
   n0(i)=count(raww==1)
   n1(i)=count(raww==2) 
   n2(i)=count(raww==3)
   af(i)=(n1(i)+2*n2(i))/(2*(ntotal-nmiss(i)))
  enddo 

  end subroutine qcraw
