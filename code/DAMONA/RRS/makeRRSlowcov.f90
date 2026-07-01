program VCFtoLoFo

  ! usage: after running RRS_panels_generation_lowcov.R
  ! assumes input are biallelic snps (there is no warning for polymorphic sites)
  ! replace file '9' ('...positions...') if running for RRS30K instead of RRS15K
  ! replace pcov = 0.25 if running for RRS15K@5x or RRS30K@5x instead of 2x
  ! adjust output file names
  ! compile : ifort -O3 -mcmodel=medium makeRRSlowcov.f90 -o makeRRSlowcov
  !./makeRRSlowcov

  implicit none

  ! variable declaration
  integer, parameter            :: ik8 = selected_int_kind(8)
  integer, parameter            :: max_nsnp = 8417679  !max number of snp in the sequence data
  integer, parameter            :: max_gbs = 60000     !max number of snp in the RRS data (not exact value)
  integer, parameter            :: nam = 132           !number of animals
  real, parameter               :: pcov = 0.10         !prop of initial counts
  integer, parameter            :: seed = 12345677

  character(len=6)              :: cc
  character(len=1)              :: rr, aa
  character(len=100)            :: fmt0
  character(len=20)             :: chr_pos
  logical                       :: cond
  integer                       :: dosage1(nam),dosage2(nam)
  integer, allocatable          :: tmp1(:),tmp2(:)
  integer(ik8)                  :: nsnp, nsnp_filt, suma, nreads
  integer                       :: i, j, k, pp, pos_underscore, correction_val
  integer                       :: filter_chr(max_gbs), filter_pos(max_gbs)

  type onevcf
    character(len=6)        :: chr
    character(len=1)        :: ref, alt, stat
    integer                 :: pos, chr_num
    integer                 :: ad1(nam), ad2(nam)
  end type onevcf

  type(onevcf), allocatable     :: vcf(:), vcf_filt(:)

  allocate(vcf(max_nsnp))
  call srand(seed)

  ! --------------------------------------------------------------------------
  ! 1) Reading of the read counts (AD) from the full sequence (WGS) for all
  !    animals, and of the list of positions for the RRS panel to be kept
  !    (e.g. RRS15K). The dropout correction file is read further below,
  !    after filtering the VCF down to the RRS panel.
  ! --------------------------------------------------------------------------
  open(8,  file = 'reads_to_correct.txt', status = 'old')
  open(9,  file = 'RRS15K_positions_to_downsample.txt', status = 'old')
  open(10, file = 'Damona132_snp_AD', status = 'old')
  open(11, file = 'out_vcf_RRS15K_2x', status = 'replace')
  open(12, file = 'out_reads_per_animal_pre_RRS15K', status = 'replace')
  open(13, file = 'out_reads_per_animal_pos_RRS15K_2x', status = 'replace')

  nsnp = 0
  do
    read(10, *, end=10) cc, pp, rr, aa, ((dosage1(j), dosage2(j)), j = 1, nam)
    nsnp = nsnp + 1
    vcf(nsnp)%chr = cc
    read(vcf(nsnp)%chr(4:), *) vcf(nsnp)%chr_num  ! Convert to integer
    vcf(nsnp)%pos = pp
    vcf(nsnp)%ref = rr
    vcf(nsnp)%alt = aa
    vcf(nsnp)%ad1 = dosage1
    vcf(nsnp)%ad2 = dosage2
  end do
  10 close(10)
  print*, 'Finished reading VCF data. Total SNPs:', nsnp

  suma = 0
  do
    read(9, *, end=200) cc, pp, chr_pos
    suma = suma + 1
    read(cc(4:), *) filter_chr(suma)  ! Extract chr_num
    filter_pos(suma) = pp
  end do
  200 close(9)
  print*, 'Finished reading RRS filtering file. Found', suma, 'positions to keep.'

  ! --------------------------------------------------------------------
  ! 2) Filtering of the full VCF: only the SNPs belonging to the RRS
  !    panel (e.g. RRS15K_positions_to_downsample.txt) are kept
  ! --------------------------------------------------------------------
  allocate(vcf_filt(suma))
  nsnp_filt = 0
  do i = 1, nsnp
    do j = 1, suma
      if (vcf(i)%chr_num == filter_chr(j) .and. vcf(i)%pos == filter_pos(j)) then
        nsnp_filt = nsnp_filt + 1
        vcf_filt(nsnp_filt) = vcf(i)  ! Keep only matching SNPs
        exit  ! No need to check further once a match is found
      end if
    end do
  end do
  nsnp = nsnp_filt  ! Update number of SNPs after filtering
  print*, 'Filtered VCF data. Remaining SNPs:', nsnp

  ! Free old VCF array and replace it with the filtered one
  deallocate(vcf)
  allocate(vcf(nsnp))
  vcf = vcf_filt(1:nsnp)  ! Copy filtered data back

  ! Free temporary storage
  deallocate(vcf_filt)

  ! --------------------------------------------------------------------
  ! 3) Application of the allelic dropout correction (calculated in script
  !    RRS_panels_generation_lowcov.R, reads_to_correct.txt file): for
  !    each affected SNP and sample, the read count is set to 0 for the
  !    ALT allele (correction_val=0, heterozygous at the cut site) or for
  !    both alleles (correction_val=1, ALT/ALT homozygous at the cut site
  !    -> fragment set to missing)
  ! --------------------------------------------------------------------
  do
    read(8, *, end=100) chr_pos, j, correction_val
    ! Extract chromosome number and position from chr_pos
    pos_underscore = index(chr_pos, "_")
    read(chr_pos(:pos_underscore-1), *) k   ! Read chromosome number
    read(chr_pos(pos_underscore+1:), *) pp  ! Read position
    ! Find corresponding SNP in VCF data
    do i = 1, nsnp
      if (vcf(i)%chr_num == k .and. vcf(i)%pos == pp) then
        if (correction_val == 0) then
          vcf(i)%ad2(j) = 0
        else if (correction_val == 1) then
          vcf(i)%ad1(j) = 0
          vcf(i)%ad2(j) = 0
        end if
        exit ! No need to check further once a match is found
      end if
    end do
  end do
  100 close(8)
  print*, 'Finished applying corrections.'

  ! Min number of reads reported on the real sequence datafile
  dosage1 = 0
  dosage2 = 0
  do j = 1, nam
    dosage1(j) = minval(vcf%ad1(j), mask = vcf%ad1(j) .gt. 0)
    dosage2(j) = minval(vcf%ad2(j), mask = vcf%ad2(j) .gt. 0)
  end do
  print*, "min number of reads reported for the REF allele: ", minval(dosage1)
  print*, "min number of reads reported for the ALT allele: ", minval(dosage2)

  ! --------------------------------------------------------------------
  ! 4) Actual downsampling: for each animal, the total number of reads
  !    (already corrected for dropout) is taken and a proportion "pcov"
  !    is sub-sampled (defined above, e.g. 0.10 = 10% of the original
  !    coverage), preserving the counts heterogeneity among individuals
  !    (subsamplead subrutine)
  ! --------------------------------------------------------------------
  ! for each animal replace with new reads
  allocate(tmp1(nsnp))
  allocate(tmp2(nsnp))

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
  ! 5) Writing of the VCF with low coverage. Variants with no reads in
  !    any animal, or with 2 or fewer reads supporting the ALT allel, are
  !    discarded (they are not considered identifiable in variant calling)
  ! --------------------------------------------------------------------
  write(fmt0,'(a20,i0,a15)') '(a6,1x,i0,2(1x,a1),',nam,'(1x,i0,1x,i0))'
  vcf%stat = '0'
  k = 0
  do i = 1, nsnp
    suma = sum( vcf(i)%ad2 )       !reads supporting the ALT allele at variant i
    cond = all( (vcf(i)%ad1 + vcf(i)%ad2) .eq. 0 )  !all samples with no call status
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
! replacement from a vector of n elemnts (sequential reservoir 
! sampling algorithm), used by subsamplead to choose which reads 
! survive the downsampling.
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

