module apply_PBC 
  implicit none
  logical, private :: IfUsed !judge if prepare_apply_PBC has used before ReturnAtom.
  integer, private :: i  
  integer, private, allocatable :: first_atom(:)
  integer, private, allocatable :: last_atom(:)
contains

  subroutine prepare_apply_PBC(n_atoms, residue_numbers, n_residues) 
    integer, intent(in)  :: residue_numbers(n_atoms)
    integer, intent(in)  :: n_residues
    integer, intent(in)  :: n_atoms

    allocate(first_atom(n_residues), last_atom(n_residues))
    do i = 1, n_atoms
      if( i == 1 ) then
        first_atom(1) = 1
      elseif( i /= 1 .and. residue_numbers(i) /= residue_numbers(i-1)) then
        first_atom(residue_numbers(i)) = i
        last_atom(residue_numbers(i-1)) = i - 1
      endif
    enddo
    last_atom(residue_numbers(i-1)) = i - 1
    IfUsed = .true.
  end subroutine


  subroutine ReturnAtom(n_atoms, nchain,n_residues,codx,cody,codz,cellsize)
    logical :: judge
    integer :: iresidue, j, jj, k, kk, n_residues
    integer :: imove(3)
    integer :: NonProtorudedResidueNo
    integer :: iatmst, iatmen
    integer,intent(in) :: n_atoms
    integer :: nchain
    integer :: ilandmark(nchain), iChainStart(nchain), iChainEnd(nchain) 
    real(4), intent(inout) :: codx(n_atoms), cody(n_atoms), codz(n_atoms)
    real(4) :: rdiff
    double precision, intent(in) :: cellsize(3) 
    if (IfUsed .neqv. .true.) then 
      stop "subroutine, prepare_apply_PBC, has not used yet before subroutine, ReturnAtom."
    endif
  !*****************************************************
  !Move atoms at each residues if atoms go beyond a box.
  !*****************************************************
    imove(:) = 0
    do iresidue = 1, n_residues
      iatmst = first_atom(iresidue)                    !First atom number of a residues 
      iatmen = last_atom(iresidue)                     !Last atom number of a residues
      do jj = iatmst + 1, iatmen              !From next atom to last atom in a residue.
        !***x-direction 
        rdiff = codx(iatmst) - codx(jj)       !calc. diff. betw. coordinate of first atom and jjth atom.
        imove(1) = idnint(rdiff / cellsize(1))       !imove is 0 or 1. if 1, move back atom. if 0, remain position of atom. 
        codx(jj) = codx(jj) + (imove(1) * cellsize(1)) !Update a coordinate of atom
  
        !***y-direction
        rdiff = cody(iatmst) - cody(jj) !calc. diff. betw. coordinate of first atom and jjth atom.
        imove(2) = idnint(rdiff / cellsize(2))
        cody(jj) = cody(jj) + (imove(2) * cellsize(2)) 

        !***z-direction
        rdiff = codz(iatmst) - codz(jj) !calc. diff. betw. coordinate of first atom and jjth atom.
        imove(3) = idnint(rdiff / cellsize(3))
        codz(jj) = codz(jj) + (imove(3) * cellsize(3)) 

      enddo

      !***Detection of non-protoruded residue Number
      judge = all(imove == 0 ) !If all of imove components are 0, return true. 
      if (judge .eqv. .true.) then 
        NonProtorudedResidueNo = iresidue 
        !print*,"Protoruded residue is",  NonProtorudedResidueNo
      endif

    enddo

  
    !*****************************************************
    !Move residue if a entire residue goes beyond a box.
    !*****************************************************
    !Landmark residues must be specified 
    !on the basis of the residues of previous conformation
    !which did not go beyond a box.

    !***peptide system only 
    ilandmark(1)   = NonProtorudedResidueNo
    iChainStart(1) = 1 
    iChainEnd(1)   = n_residues 

    !**For ET1
    !ilandmark(1) = 5
    !ilandmark(2) = NonProtorudedResidueNo
    !iChainStart(1) = 1 
    !iChainStart(2) = 21 
    !iChainEnd(1) = 20
    !iChainEnd(2) = 40
    
    !***new
    imove(:) = 0

    do i = 1, nchain
      do iresidue = ilandmark(i)-1, iChainStart(i), -1
        iatmst = first_atom(iresidue)
        iatmen = last_atom(iresidue)
  
        !*************************************************************************************************
        !All atoms in a residue will be moved, 
        !if the distance btw 1st atom of residue(i+1) and residue(i) is over half of the length of x-cell.
        !*************************************************************************************************
        rdiff  = codx(first_atom(iresidue+1)) - codx(first_atom(iresidue)) !If neiboring res. is far to iresi+1 res., rdiff is large. 
        imove(1)  = idnint(rdiff / cellsize(1))
        !imove  = idnint(rdiff / cellsize(1))
        do jj = iatmst, iatmen
          codx(jj) = codx(jj) + (imove(1) * cellsize(1))
          !codx(jj) = codx(jj) + (imove * cellsize(1))
        enddo
  
        rdiff  = cody(first_atom(iresidue + 1)) - cody(first_atom(iresidue))
        imove(2)  = idnint(rdiff / cellsize(2))
        !imove  = idnint(rdiff / cellsize(2))
        do jj = iatmst, iatmen
          cody(jj) = cody(jj) + (imove(2) * cellsize(2))
          !cody(jj) = cody(jj) + (imove * cellsize(2))
        enddo
        
        rdiff  = codz(first_atom(iresidue + 1)) - codz(first_atom(iresidue))
        imove(3)  = idnint(rdiff / cellsize(3))
        !imove  = idnint(rdiff / cellsize(3))
        do jj = iatmst, iatmen
          codz(jj) = codz(jj) + (imove(3) * cellsize(3))
          !codz(jj) = codz(jj) + (imove * cellsize(3))
        enddo

      enddo
    enddo


    imove(:) = 0
    do i = 1, nchain
      do iresidue = ilandmark(i)+1, iChainEnd(i)
        iatmst = first_atom(iresidue)
        iatmen = last_atom(iresidue)
  
        rdiff = codx(first_atom(iresidue-1)) - codx(first_atom(iresidue))
        imove(1) = idnint(rdiff / cellsize(1))
        !imove = idnint(rdiff / cellsize(1))
        do jj = iatmst, iatmen
          codx(jj) = codx(jj) + (imove(1) * cellsize(1))
          !codx(jj) = codx(jj) + (imove * cellsize(1))
        enddo
  
        rdiff = cody(first_atom(iresidue-1)) - cody(first_atom(iresidue))
        imove(2) = idnint(rdiff / cellsize(2))
        !imove = idnint(rdiff / cellsize(2))
        do jj = iatmst, iatmen
          cody(jj) = cody(jj) + (imove(2) * cellsize(2))
          !cody(jj) = cody(jj) + (imove * cellsize(2))
        enddo
  
        rdiff = codz(first_atom(iresidue-1)) - codz(first_atom(iresidue))
        imove(3) = idnint(rdiff / cellsize(3))
        !imove = idnint(rdiff / cellsize(3))
        do jj = iatmst, iatmen
          codz(jj) = codz(jj) + (imove(3) * cellsize(3))
          !codz(jj) = codz(jj) + (imove * cellsize(3))
        enddo
      enddo

    enddo
    
    !write(*,'(a)') "Note) Atoms have been returned (module periodic was uesd)."
  end subroutine


end module
