module assignment_prob
  implicit none
  private count_pdf
contains

  subroutine count_pdf(FilNamPDF, ncount)
    integer :: i, ios
    integer :: UnitPDF = 1022
    integer           , intent(out) :: ncount
    character(len=*)  , intent(in)  :: FilNamPDF
    !character(len=120), intent(in)  :: FilNamPDF
    real(8)  :: energy, pdf

    open(UnitPDF,file = FilNamPDF )
    ncount = 0
    ios = 0
    do
      read(UnitPDF,*,iostat=ios) energy, pdf
      if (ios /= 0)  exit
      ncount = ncount + 1
    enddo
    close(UnitPDF)
  end subroutine
  
  subroutine read_pdf(FilNamPDF,energy,pdf)
  !subroutine read_pdf(FilNamPDF,ncount,energy,pdf)
    integer :: i, j, k
    integer :: ncount
    integer :: UnitPDF =1001, UnitOut = 123
    integer :: ios
    character(len=*)    , intent(in)  :: FilNamPDF
    real(8), allocatable, intent(out) :: energy(:), pdf(:)
    real(8) :: MaxPDF
    real(8), parameter :: LowLimit = 0.0000000001

    print*," *Start to read :", FilNamPDF

    !***count data points in order to allocate energy and pdf array
    call count_pdf(FilNamPDF, ncount)
    print*, ncount
    allocate(energy(ncount), pdf(ncount))
    energy(:) = 0.0d0
    pdf(:)    = 0.0d0

    open(UnitPDF, file = FilNamPDF, status="old")
    do i = 1, ncount
      read(UnitPDF,*) energy(i), pdf(i) 
!      write(1111, *) energy(i), pdf(i)
    enddo

    !If input data of prob distirib is normalized (the max val is 1), the following line must not be required. 
    !If the input is not normalized, this line must be active.
    MaxPDF = maxval(pdf)
    print*, "Max value of pdf = ", MaxPDF

!***Normalization of PDF based on max value of PDF
    do i = 1, ncount
      pdf(i) = exp(pdf(i) - MaxPDF)
      if (pdf(i) <= LowLimit) pdf(i) = 0.0d0
      !write(1111, *) energy(i), pdf(i)
    enddo

    close(UnitPDF)

  end subroutine


  subroutine assign_pdf(enePDF, pdf, potential, prob)
    integer              :: i ,j ,k
    real(8), intent(in)  :: enePDF(:), pdf(:)
    real(4), intent(in)  :: potential
    real(8), intent(out) :: prob
    real(4)              :: MultProb
    integer              :: NdataPDF

    NdataPDF = size(pdf)

    do i = 1, NdataPDF - 1
      if(potential >= enePDF(i) .and. potential < enePDF(i+1)) then 
        MultProb = pdf(i) + pdf(i+1)
        MultProb = MultProb * 0.5d0
        prob     = MultProb
        write(*,'("Potential                          :", f15.3)') potential
        write(*,'("The range where a conf. is detected:", 2f15.3)') enePDF(i), enePDF(i+1)
        write(*,'("Neighbor Prob (i and i+1)          :", 2f12.8)') pdf(i), pdf(i+1) 
        write(*,'("Probability (Non-log scale)        :", f12.8)') prob
        exit
      endif
    enddo
  end subroutine
end module
