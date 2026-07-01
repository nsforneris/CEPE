program VCFtoLoFo

 ! assumes input are biallelic snps (there is no warning for polymorphic sites)
 ! arguments used for each coverage: 
 ! 0.1x: seed = 541986; pcov = 0.005
 ! 0.2x: seed = 541987; pcov = 0.01
 ! 0.5x: seed = 541985; pcov = 0.025
 ! 1x:   seed = 541984; pcov = 0.05
 ! 2x:   seed = 541983; pcov = 0.10
 ! 5x:   seed = 541988; pcov = 0.25

 ! compile : ifort -O3 -mcmodel=medium makelowcov.f90 -o makelowcov
 !./makelowcov

 implicit none

 ! variable declaration
 integer, parameter            :: ik8 = selected_int_kind(8)
 integer, parameter            :: nsnp = 8417679  !max number of snp in the sequence data
 integer, parameter            :: nam = 132       !number of animals
 real, parameter               :: pcov = 0.025    !prop of initial counts
 integer, parameter            :: seed = 541985

 character(len=6)              :: cc
 character(len=1)              :: rr, aa
 character(len=100)            :: fmt0
 logical                       :: cond
 integer                       :: dosage1(nam),dosage2(nam),tmp1(nsnp),tmp2(nsnp)
 integer(ik8)                  :: nreads, suma
 integer                       :: i,j,k,pp

 type onevcf
       character(len=6)        :: chr
       character(len=1)        :: ref, alt, stat
       integer                 :: pos
       integer                 :: ad1(nam), ad2(nam)
 end type onevcf
 type(onevcf)                  :: vcf(nsnp)

 call srand(seed)

 ! --------------------------------------------------------------------------
 ! 1) Reading of the read counts (AD) from the full sequence (WGS) for all
 !    animals
 ! --------------------------------------------------------------------------
 open(10, file = 'Damona132_snp_AD', status = 'old')
 
 open(11, file = 'out_vcf_lcWGS_0.5x', status = 'replace')
 open(12, file = 'out_reads_per_animal_pre', status = 'replace')
 open(13, file = 'out_reads_per_animal_pos_0.5x', status = 'replace')

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
 
 ! --------------------------------------------------------------------------
 ! 2) Min number of reads reported on the real sequence datafile
 ! --------------------------------------------------------------------------
 dosage1 = 0
 dosage2 = 0
 do j = 1, nam
   dosage1(j) = minval(vcf%ad1(j), mask = vcf%ad1(j) .gt. 0)
   dosage2(j) = minval(vcf%ad2(j), mask = vcf%ad2(j) .gt. 0)
 end do
 print*, "min number of reads reported for the REF allele: ", minval(dosage1)
 print*, "min number of reads reported for the ALT allele: ", minval(dosage2)

 ! --------------------------------------------------------------------
 ! 3) Actual downsampling: for each animal, the total number of reads
 !    is taken and a proportion "pcov" is sub-sampled (defined above,
 !    e.g. 0.10 = 10% of the original coverage), preserving the counts
 !    heterogeneity among individuals(subsamplead subrutine)
 ! --------------------------------------------------------------------
 do j = 1, nam
     nreads = sum(vcf%ad1(j)) + sum(vcf%ad2(j))
     write(12,'(i0,1x,i0,1x,i0)') nreads, dosage1(j), dosage2(j)
     tmp1 = vcf%ad1(j)
     tmp2 = vcf%ad2(j)
     call subsamplead(pcov, nsnp, nreads, tmp1, tmp2)
     vcf%ad1(j) = tmp1
     vcf%ad2(j) = tmp2
 end do
 close(12)
 print*, 'finished subsampling'

 ! --------------------------------------------------------------------
 ! 4) Writing of the VCF with low coverage. Variants with no reads in
 !    any animal, or with 2 or fewer reads supporting the ALT allel, are
 !    discarded (they are not considered identifiable in variant calling)
 ! --------------------------------------------------------------------
 write(fmt0,'(a20,i0,a15)') '(a6,1x,i0,2(1x,a1),',nam,'(1x,i0,1x,i0))'
 vcf%stat = '0'
 k = 0
 do i = 1, nsnp
   suma = sum( vcf(i)%ad2 )                        !reads supporting the alt allele at variant i
   cond = all( (vcf(i)%ad1 + vcf(i)%ad2) .eq. 0 )  !all samples with zero reads
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
   write(13,'(i0)') nreads
 end do
 close(13)

end program VCFtoLoFo

! --------------------------------------------------------------------
! Subrutine subsamplead: performs random sub-sampling of reads
! (without replacement) for one animal, accross all SNPS in the panel.
! It builds a vector "b" with one entry per individual read
! (indicating which SNP it belongs to and whether it is the ref or alt
! allele), randomly selects "subreads = prop*reads" of those reads,
! and reconstructs the ad1/ad2 counts per SNP from the subsample.
! --------------------------------------------------------------------
subroutine subsamplead(prop, snps, reads, ad1, ad2)
 real, intent(in)          :: prop
 integer, intent(in)       :: snps, reads
 integer, intent(inout)    :: ad1(*), ad2(*)
 integer                   :: b(reads,2)
 integer                   :: bord(reads)
 integer, allocatable      :: subord(:)
 integer                   :: suma, totr, tota, subreads

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
           !print "(i0,1x,i0)", b(suma, 1), b(suma, 2)
     end do
     do ii = 1, tota
           suma = suma + 1
           b(suma, 1) =  i    !pointer to snp info
           b(suma, 2) =  1    !alt
           !print "(i0,1x,i0)", b(suma, 1), b(suma, 2)
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

! --------------------------------------------------------------------
! Subrutine ransam: random sub-sampling of k elements without
! replacement from a vector of n elemnts, used by subsamplead to choose
! which reads survive the downsampling.
! --------------------------------------------------------------------
subroutine ransam(x, a, n, k)
 integer, intent(in)               :: n, k
 integer, intent(in)               :: x(*)
 integer, intent(inout)            :: a(*)

 m = 0
 do 50 j = 1, n
   l = int( (float(n-j+1)) * rand(0) ) + 1
   if (l .gt. (k-m) ) go to 50
   m = m + 1
   a(m) = x(j)
   if (m .ge. k) go to 99
 50 continue
 99 return
end subroutine ransam
