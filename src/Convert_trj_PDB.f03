module Convert_trj_PDB 
  use variables 
  implicit none

  private count_atomtype
  public  read_pdb
  public  analyze_trj 
  public  outputPDB

contains
  
  subroutine count_atomtype(PDB, atomtype)
    type(var_PDB), intent(in)    :: PDB
    character(len=*), intent(in) :: atomtype
    integer :: iatom, icount

    icount = 0
    do iatom = 1, PDB%n_atoms
      if (trim(adjustl(PDB%AtomName(iatom))) == atomtype) then 
        !print*, len_trim(adjustl(PDB%AtomName(iatom))), trim(PDB%AtomName(iatom))
        print*, PDB%ResNum(iatom), PDB%AtomName(iatom)
        icount = icount + 1
      endif
    enddo
    print*, "No of CA=", icount

  end subroutine

  subroutine read_pdb(FilName, PDB)
    type(var_PDB), intent(inout) :: PDB 
    character(len=*), intent(in) :: FilName
    character(len=6)   :: atom
    integer, parameter :: FilUnit = 11
    integer, parameter :: OutPDBUnit  = 12
    integer :: i
    integer :: NotAtomLine
  
    open(FilUnit, file=FilName, status="old")
    print*,"  File Name: ",FilName
    allocate(PDB%AtomNum(PDB%n_atoms) , & 
             PDB%AtomName(PDB%n_atoms), & 
             PDB%ResName(PDB%n_atoms) , &
             PDB%ChainId(PDB%n_atoms) , &
             PDB%ResNum(PDB%n_atoms)  , & 
             PDB%x(PDB%n_atoms)       , & 
             PDB%y(PDB%n_atoms)       , & 
             PDB%z(PDB%n_atoms) )
    PDB%AtomNum(:) = 0
    PDB%ResNum(:)  = 0
    PDB%x(:) = 0; PDB%y(:) = 0; PDB%z(:) = 0
  
  !***Detect first atom line
    do
      read(FilUnit,*) atom
      if(atom /= "ATOM" .and. atom /= "HETATM") then
        NotAtomLine = NotAtomLine + 1  
        cycle
      elseif(atom == "ATOM" .or. atom == "HETATM") then
        exit
      endif
    enddo
    backspace(FilUnit)
  !***start reading
    do i = 1, PDB%n_atoms 
      read(FilUnit,"(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3)") &
           atom       , &
           PDB%AtomNum(i) , &
           PDB%AtomName(i), &
           PDB%ResName(i) , &
           PDB%ChainId(i) , &
           PDB%ResNum(i)  , &
           PDB%x(i),        &
           PDB%y(i),        &
           PDB%z(i)  
           if (atom == "TER") stop "TER lines are included, remove these."
    enddo
    close(FilUnit)
    !call count_atomtype(PDB,"CA")

  end subroutine

!------------------------------------------------------
  subroutine analyze_trj(filename,fnameOut, PDB, A_snapshot, isOutputTrj)
      use apply_PBC 
      implicit none
      type(var_PDB)     , intent(inout) :: PDB 
      type(var_snapshot), intent(inout) :: A_snapshot
      logical           , optional      :: isOutputTrj
      character(len=*)                  :: filename, fnameOut
      integer                           :: unit = 11, unit_outPDB = 10
      integer                           :: ios, iatom, n_confs

      integer(4)                        :: istp,iyn15v,iyn15h
      real(4)                           :: sitime, sec, et, kinetic, temperature, rmsf, rmsd 

      !***for coordinated in a restart file
      double precision, allocatable :: cord(:,:)


      !***default of isOutputTrj  is .false.
      if (.not. present(isOutputTrj)) isOutputTrj = .false.
  
      open(unit, file = filename, form="unformatted", status="old")
      allocate(A_snapshot%x(PDB%n_atoms),A_snapshot%y(PDB%n_atoms),A_snapshot%z(PDB%n_atoms))
      
      !***Initialization
      A_snapshot%iconf = 0; A_snapshot%x(:) = 0; A_snapshot%y(:) = 0; A_snapshot%z(:) = 0
      
!      call prepare_apply_PBC(PDB%n_atoms,PDB%ResNum,PDB%n_residues)

      !***Analyzing trajectory
      do 
          read(unit, iostat=ios) istp,sitime,sec,et,kinetic,temperature,&
                                 A_snapshot%potential,rmsf,iyn15v,iyn15h,rmsd
          if (ios /= 0) exit
          read(unit) (A_snapshot%x(iatom), A_snapshot%y(iatom), A_snapshot%z(iatom), iatom = 1, PDB%n_atoms)

          !@@@@ for restart file
!          allocate(cord(1:3,PDB%n_atoms))
!          read(unit)
!          read(unit)
!          read(unit)
!          read(unit) (cord(1:3,iatom),iatom=1,PDB%n_atoms)
!          print*, "double precision : ", cord(1,1)
!          A_snapshot%x(:) = real(cord(1,:))
!          A_snapshot%y(:) = real(cord(2,:))
!          A_snapshot%z(:) = real(cord(3,:))
!          print*, "signle precision : ", A_snapshot%x(1)
          !@@@@

          print*, "step no:", istp, "Potential=", A_snapshot%potential
 
          A_snapshot%iconf = A_snapshot%iconf + 1

!          call ReturnAtom(PDB%n_atoms,PDB%n_chains,PDB%n_residues,A_snapshot%x, A_snapshot%y, A_snapshot%z, & 
!                          A_snapshot%cellsize,PDB%ResNum)

          if (isOutputTrj) call outputPDB(unit_outPDB,fnameOut, PDB, A_snapshot)

          print*,"Conf NO:", A_snapshot%iconf 
      enddo

      n_confs = A_snapshot%iconf
      write(*,'("#Number of conformation ",i8)') n_confs 

  end subroutine

  subroutine DetectStrangeCoord(coords, threshold)
    real(4), intent(in)  :: coords(3)
    real(4), optional :: threshold

    if ( .not. present(threshold)) then 
      threshold = 100.0
      print*, "Threshold=", threshold
    endif

    if (coords(1) >= threshold .or. coords(2) >= threshold .or. coords(3) >= threshold) then
      print*, "STOP because coordinates are too big to realize PDB format."
      print*, "Coordinates=", coords, "   Threshold=", threshold
      stop
    endif
  end subroutine 

  subroutine outputPDB(unit, fnameOut, PDB, a_snapshot)
      integer      , intent(in) :: unit 
      character(len=*), intent(in) :: fnameOut !@
      integer :: i
      type(var_PDB), intent(in) :: PDB 
      type(var_snapshot), optional, intent(in) :: a_snapshot 
      logical :: IsCorrect
      if (present(a_snapshot)) then 
          open(unit, file = fnameOut) ! "frames.pdb")
          write(unit,"('MODEL', i7)") a_snapshot%iconf
          write(unit,"('#Potential, kcal/mol ', f13.3)") a_snapshot%potential 
          do i = 1, PDB%n_atoms 
            call DetectStrangeCoord( (/a_snapshot%x(i), a_snapshot%y(i), a_snapshot%z(i)/),threshold=200.0)
            write(unit,"(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3,a26)") &
              "ATOM  ", &
              PDB%AtomNum(i), &
              PDB%AtomName(i),&
              PDB%ResName(i),&
              PDB%ChainId(i), &
              PDB%ResNum(i), &
              a_snapshot%x(i), a_snapshot%y(i), a_snapshot%z(i),&
              PDB%blnk
          enddo
          write(unit,"(a6)") "ENDMDL"
      else
          open(unit, file =fnameOut, status="replace")
          do i = 1, PDB%n_atoms 
            write(unit,"(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3,a26)") &
              "ATOM  ", &
              PDB%AtomNum(i), &
              PDB%AtomName(i),&
              PDB%ResName(i),&
              PDB%ChainId(i), &
              PDB%ResNum(i), &
              PDB%x(i),        &
              PDB%y(i),        &
              PDB%z(i),        & 
              PDB%blnk
          enddo
      endif
      !close(unit)
  end subroutine

  
end module
