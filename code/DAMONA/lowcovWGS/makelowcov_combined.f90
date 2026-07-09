program VCFtoLoFo

 ! assumes input are snps (there is no warning for polymorphic sites)
 ! Generates a low-coverage WGS dataset from 'nam' animals, where 'subset_size' = 'nam'/3.
 ! The animals are divided into three groups of equal size:
 ! the first group is assigned 2× coverage (pcov1 = 0.1), the second group a 5x coverage 
 ! (pcov2 = 0.25), and the third group 10x coverage (pcov3 = 0.5).
 
 ! compile : ifort -O3 -mcmodel=medium makelowcov.f90 -o makelowcov
 !./makelowcov

 implicit none

 ! variable declaration
 integer, parameter            :: ik8   = selected_int_kind(8)
 integer, parameter            :: nsnp  = 8417679  !max number of snp in the sequence data
 integer, parameter            :: nam   = 132      !number of animals
 integer, parameter            :: subset_size = 44 !number of animals with one particular coverage
 integer, parameter            :: seed  = 451983
 real, parameter               :: pcov1 = 0.10     !prop1 of initial counts
 real, parameter               :: pcov2 = 0.25     !prop2 of initial counts
 real, parameter               :: pcov3 = 0.50     !prop3 of initial counts

 integer                       :: seed_array(1)   
 integer                       :: dosage1(nam),dosage2(nam),tmp1(nsnp),tmp2(nsnp)
 integer                       :: i,j,k,pp
 integer(ik8)                  :: nreads, suma
 real                          :: rand_val,rtmp
 real                          :: psamp(nam)      !Array to hold coverage values
 character(len=6)              :: cc
 character(len=1)              :: rr,aa
 character(len=100)            :: fmt0
 logical                       :: cond

 type onevcf
       character(len=6)        :: chr
       character(len=1)        :: ref, alt, stat
       integer                 :: pos
       integer                 :: ad1(nam), ad2(nam)
 end type onevcf
 type(onevcf)                  :: vcf(nsnp)

 ! initialize random seed 
 seed_array(1) = seed
 call random_seed(put=seed_array)

 ! initialize psamp array
 psamp(1:subset_size)                 = pcov1  !first 44 individuals 0.10
 psamp(subset_size+1:2*subset_size)   = pcov2  !next  44 individuals 0.25
 psamp(2*subset_size+1:3*subset_size) = pcov3  !last  44 individuals 0.50

 ! Fisher-Yates shuffle (to randomly permute psamp)
 do i = nam, 2, -1
    call random_number(rand_val)
    j = int(rand_val * i) + 1         !pick a random index in range [1, i]
    ! swap elements in psamp array
    rtmp = psamp(i)
    psamp(i) = psamp(j)
    psamp(j) = rtmp
 end do

 ! open file
 open(10, file = 'Damona132_snp_AD', status = 'old')
 open(11, file = 'out_vcf_lcWGS_2x5x10x', status = 'replace')
 open(12, file = 'out_reads_per_animal_pre', status = 'replace')
 open(13, file = 'out_reads_per_animal_pos_2x5x10x', status = 'replace')

 ! reading info from file
 do i = 1, nsnp
   read(10,*) cc, pp, rr, aa, ( (dosage1(j), dosage2(j)), j = 1, nam)
   vcf(i)%chr = cc
   vcf(i)%pos = pp
   vcf(i)%ref = rr
   vcf(i)%alt = aa
   vcf(i)%ad1 = dosage1
   vcf(i)%ad2 = dosage2
 end do
 close(10)
 print*, 'finished reading'
 !
 ! Min number of reads reported on the real sequence datafile
 dosage1 = 0
 dosage2 = 0
 do j = 1, nam
    dosage1(j) = minval(vcf%ad1(j), mask = vcf%ad1(j) .gt. 0)
    dosage2(j) = minval(vcf%ad2(j), mask = vcf%ad2(j) .gt. 0)
 end do
 print*, "min number of reads reported for the REF allele: ", minval(dosage1)
 print*, "min number of reads reported for the ALT allele: ", minval(dosage2)
 !
 ! for each animal replace with new reads
 do j = 1, nam
     nreads = sum(vcf%ad1(j)) + sum(vcf%ad2(j))
     write(12,'(i0,1x,i0,1x,i0)') nreads, dosage1(j), dosage2(j)
     tmp1 = vcf%ad1(j)
     tmp2 = vcf%ad2(j)
     call subsamplead(psamp(j), nsnp, nreads, tmp1, tmp2)
     vcf%ad1(j) = tmp1
     vcf%ad2(j) = tmp2
 end do
 close(12)
 print*, 'finished subsampling'

 ! write new vcf with reduced reads
 write(fmt0,'(a20,i0,a15)') '(a6,1x,i0,2(1x,a1),',nam,'(1x,i0,1x,i0))'
 vcf%stat = '0'
 k = 0
 do i = 1, nsnp
   suma = sum( vcf(i)%ad2 )                        !reads supporting the alt allele
   cond = all( (vcf(i)%ad1 + vcf(i)%ad2) .eq. 0 )  !no call status
   if( cond .or. (suma .le. 2)  ) cycle
   vcf(i)%stat = '1'
   k = k + 1
   write(11,fmt0) vcf(i)%chr,vcf(i)%pos,vcf(i)%ref,vcf(i)%alt,((vcf(i)%ad1(j),vcf(i)%ad2(j)),j=1,nam)
 end do
 close(11)
 print*, k, ' variants remained'

 ! for each animal print number of reads before and after reducing coverage
 do j = 1, nam
     nreads = sum( (vcf%ad1(j) + vcf%ad2(j)), mask = vcf%stat .eq. '1' )
     write(13,'(i0,1x,f4.2)') nreads, psamp(j)
 end do
 close(13)

 end program VCFtoLoFo

 subroutine subsamplead(prop, snps, reads, ad1, ad2)
 implicit none
 real, intent(in)          :: prop
 integer, intent(in)       :: snps, reads
 integer, intent(inout)    :: ad1(*), ad2(*)
 integer                   :: b(reads,2)
 integer                   :: bord(reads)
 integer, allocatable      :: subord(:)
 integer                   :: suma, totr, tota, subreads
 integer                   :: i, ii, j

 do i = 1, reads
    bord(i) = i
 end do
 ! creates array b with all the sampled alleles for an animal
 ! accross all snp
 b = 0
 suma = 0
 do i = 1, snps
     totr = ad1(i)
     tota = ad2(i)
     do ii = 1, totr
           suma = suma + 1
           b(suma, 1) =  i    !pointer to snp info
           b(suma, 2) =  0    !ref
     end do
     do ii = 1, tota
           suma = suma + 1
           b(suma, 1) =  i    !pointer to snp info
           b(suma, 2) =  1    !alt
     end do
 end do

 ! random sample n=subreads dosages
 subreads = int(prop*reads)
 allocate( subord(subreads) )
 subord = 0
 call ransam(bord, subord, reads, subreads)

 ! replace old dosages with new dosages
 do i = 1, snps
    ad1(i) = 0
    ad2(i) = 0
 end do
 do i = 1, subreads
    j = subord(i)
    ad1(  b(j,1)  ) = ad1( b(j,1) ) + (-1)*b(j,2) + 1
    ad2(  b(j,1)  ) = ad2( b(j,1) ) + b(j,2)
 end do

 deallocate( subord )
 return
 end subroutine subsamplead

 subroutine ransam(x, a, n, k)
 implicit none
 integer, intent(in)               :: n, k
 integer, intent(in)               :: x(*)
 integer, intent(inout)            :: a(*)
 integer                           :: j,l,m
 real                              :: rand_val

 m = 0
 do j = 1, n
    call random_number(rand_val)
    l = int(real(n-j+1, kind=8) * rand_val) + 1
    if (l .gt. (k-m) ) cycle
    m = m + 1
    a(m) = x(j)
    if (m .ge. k) exit
 end do
 return
 end subroutine ransam
