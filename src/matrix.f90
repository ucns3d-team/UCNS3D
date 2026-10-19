module lapck
   use mpiinfo
   use declaration
   implicit none

contains


   subroutine transposematrix(qff, qtff, ideg)
      implicit none
      integer, intent(in) :: ideg
      real, allocatable, dimension(:, :), intent(in) :: qff
      real, allocatable, dimension(:, :), intent(inout) :: qtff
      integer i, j
      do i = 1, ideg
         do j = 1, ideg
            qtff(i, j) = qff(j, i)
         end do
      end do
   end subroutine


   subroutine qrdecomposition(lscqm, qff, rff, ideg, vqr)
      implicit none
      integer, intent(in) :: ideg
      real, allocatable, intent(in) :: lscqm(:, :)
      real, allocatable, intent(inout) :: qff(:, :), rff(:, :)
      real, intent(inout) :: vqr(:)
      real :: xnorm, vnorm, projection
      integer :: i, j, k, l

      qff(1:ideg, 1:ideg) = 0.0
      rff(1:ideg, 1:ideg) = lscqm(1:ideg, 1:ideg)
      do i = 1, ideg
         qff(i, i) = 1.0
      end do
      do l = 1, ideg
         ! Same Householder sign and normalization as the original routine.
         xnorm = 0.0
         do k = l, ideg
            vqr(k) = rff(k, l)
            xnorm = xnorm + vqr(k)**2
         end do
         vqr(l) = vqr(l) + sign(1.0, vqr(l)) * sqrt(xnorm)
         vnorm = 0.0
         do k = l, ideg
            vnorm = vnorm + vqr(k)**2
         end do
         vnorm = sqrt(vnorm)
         vqr(l:ideg) = vqr(l:ideg) / vnorm

         ! R <- (I - 2 v v^T) R: no dense reflection matrix or matrix product.
         do j = 1, ideg
            projection = 0.0
            do k = l, ideg
               projection = projection + vqr(k) * rff(k, j)
            end do
            do i = l, ideg
               rff(i, j) = rff(i, j) - 2.0 * vqr(i) * projection
            end do
         end do
         ! Q <- Q (I - 2 v v^T).
         do i = 1, ideg
            projection = 0.0
            do k = l, ideg
               projection = projection + qff(i, k) * vqr(k)
            end do
            do j = l, ideg
               qff(i, j) = qff(i, j) - 2.0 * projection * vqr(j)
            end do
         end do
      end do
   end subroutine qrdecomposition

end module
