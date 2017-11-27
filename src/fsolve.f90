!==============================================================================================================
  subroutine fsolve(n,m,cls,nr,rws,fnRAW,nit,lambda,tol,y,g,e,s,mean,sd)
!==============================================================================================================

  implicit none
  
  integer*4 :: i,j,k,m,n,nr,it,nit,cls(m),rws(nr)
  real*8 :: y(n),e(n),raww(n),w(n),g(n)
  real*8 :: dww(m),s(m),os(m),lambda(m),mean(m),sd(m)
  real*8 :: lhs,rhs,snew,tol,sigma
  character*1 :: raw(n)
  character(len=1000) :: fnRAW
  
  open (unit=13,file=fnRAW, status='old', form='unformatted', access='direct', recl=n)
  ! genotypes coded 0=missing,1,2,3 => where 1,2,3 means 0,1,2 copies of alternative allele 
  do i=1,m
  read (unit=13,rec=cls(i)) raw
  raww = ichar(raw)
  where (raww.ne.0.0D0)
  w = (raww-1.0D0-mean(i))/sd(i)
  elsewhere
  w = 0.0D0
  end where
  dww(i)=dot_product(w(rws),w(rws)) 
  if(s(i).eq.0.0D0) then 
    s(i)=(dot_product(w(rws),y(rws))/dww(i))/m
  endif     
  enddo
  close (unit=13)
  os=s

  e=0.0D0
  e(rws)=y(rws)

  do it=1,nit
  g=0.0D0
  open (unit=13,file=fnRAW, status='old', form='unformatted', access='direct', recl=n)
  do i=1,m
  read (unit=13,rec=cls(i)) raw
  raww = ichar(raw)		
  where (raww.ne.0.0D0)
  w = (raww-1.0D0-mean(i))/sd(i)
  elsewhere
  w = 0.0D0
  end where
  lhs=dww(i)+lambda(i)
  rhs=dot_product(w,e) + dww(i)*s(i)
  snew=rhs/lhs
  
  !e=e - w*(snew-s(i))
  e(rws)=e(rws) - w(rws)*(snew-s(i))
  s(i)=snew
  g=g+w*s(i)
  enddo
  close (unit=13)
  print*,(sum((s-os)**2)/sqrt(real(m)))
  if( (sum((s-os)**2)/sqrt(real(m)))<tol) exit  
  os=s  
  enddo
  
  nit=it
  tol=sum((s-os)**2)
  
  end subroutine fsolve
