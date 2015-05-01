program Concat_Arrays
implicit none
 
  real, dimension(10,10) :: a 
  real, dimension(10) :: b = [ 4, 5, 6, 4, 5, 2, 7, 4, 9, 1]
  real, dimension(:,:), allocatable :: c
 
  
  a = 1.0d0
  !a(1,1:3) = [1, 2, 3]
  !a(2,:) = [1, 2, 3]
  !a(3,:) = [1, 2, 3]
  !a(1,1) = 1
  !a(2,1) = 2
  !a(3,1) = 3
  !a(1,2) = 1
  !a(2,2) = 2
  !a(3,2) = 3
  !a(1,3) = 1
  !a(2,3) = 2
  !a(3,3) = 3
  
  allocate(c(size(a,1)+1, size(a,2)+1))
  c = 0.0d0 
  c(1:size(a,1),1:size(a,2)) = a
  c(size(a,1)+1, 1:size(a,2)) = b
  c(1:size(a,1), size(a,2)+1) = b
  
  !c(1:size(a)) = a
  !c(size(a)+1:size(a)+size(b)) = b
  
  write(*,*) size(a,1)
  write(*,*) a
 
end program Concat_Arrays
