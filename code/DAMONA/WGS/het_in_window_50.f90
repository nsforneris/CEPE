!ifort -O3 -mcmodel=medium het_in_window.f90 -o het_in_wd

program heterozygosity_ratio
    implicit none
    integer, parameter   :: n_indiv = 132
    integer, parameter   :: max_snp = 600000

    integer              :: i, j, k 
    character(len=20)    :: idval
    integer              :: nsnp, chrom, posval, left, right, total_snps, total_length, miss, het_count
    character(len=1)     :: ref, alt

    integer, allocatable :: pos(:)
    integer, allocatable :: gen(:, :)
    character(len=20), allocatable :: id(:)    

    ! Window size: 50 KB (±25 KB around the SNP)
    integer, parameter :: window_size = 50000
    integer, parameter :: half_window = window_size / 2  ! 25 KB

    ! Allocate memory
    allocate(id(max_snp), pos(max_snp))
    allocate(gen(max_snp, n_indiv))
    
    open(9, file="genfile.txt", status="old")

    ! Read SNP data
    nsnp = 0
    do
       read(9, *, end=10) chrom, idval, posval, ref, alt, (gen(nsnp + 1, j), j = 1, n_indiv)
       nsnp = nsnp + 1
       id(nsnp) = idval
       pos(nsnp) = posval
    end do
10 close(9)
   print*, 'Finished reading VCF data. Total SNPs:', nsnp

   
   ! Open files for output
   open(11, file="out_het_snp_info.txt", status="replace")   ! File for SNP IDs
   open(12, file="out_het_counts.txt", status="replace")     ! File for heterozygosity counts

   ! Write SNP IDs to snp_info.txt

   ! Process each SNP
   do i = 1, nsnp
      ! Define fixed window around SNP (±25 KB)
      left = i
      right = i

      ! Find leftmost SNP within 25 KB
      do while (left > 1 .and. pos(i) - pos(left) <= half_window)
          left = left - 1
      end do

      ! Find rightmost SNP within 25 KB
      do while (right < nsnp .and. pos(right) - pos(i) <= half_window)
          right = right + 1
      end do

      ! Compute the total number of SNPs in this window
      total_snps = right - left + 1

      ! Compute the actual length this window
      total_length = (pos(right) - pos(left))

      ! Ignore SNP if less than 100 SNPs in window
      if (total_snps < 100) cycle

      ! Print SNP info in output file
      write(11, '(i0,1x,a,1x,i0,1x,i0)') i, id(i), total_snps, total_length  

      ! Compute heterozygosity ratio for each individual
      do j = 1, n_indiv
          het_count = 0
          miss = 0
          do k = left, right
              if (gen(k, j) .eq. 1) het_count = het_count + 1
              if (gen(k, j) .ne. 9) miss = miss + 1
          end do
          if (miss .gt. 50) then
               write(12, '(i0,1x,i0,1x,f10.8)') i, j, (het_count*1.D0) / (miss*1.D0)
          else
               cycle
          end if
      end do

      ! Print message after every 10,000 SNPs processed
      if (mod(i, 10000) == 0) then
          print*, 'Processed', i, 'SNPs.'
      end if

   end do

   ! Close output file
   close(11)
   close(12)

end program heterozygosity_ratio
