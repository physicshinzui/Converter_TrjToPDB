module apply_PBC 
  implicit none
  logical, private :: IfUsed !judge if prepare_apply_PBC has used before ReturnAtom.
  integer, private :: i  
  integer, private, allocatable :: first_atom(:)
  integer, private, allocatable :: last_atom(:)
  !TODO: put variables which are shared with some subroutines.
contains

  !Bad subroutine name. this should be "save_firstLastAtomNo_eachResidue".
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

!------
  subroutine get_residue_no_for_landmark(n_atoms,nchains,codx,cody,codz,atom_ranges,ilandmark,res_nums,loc)
    !This routine searches for two landmark residues that are defined as the minimum disatnce of the pair atoms 
    !consisting of chain1's and chain2's. 

    integer, intent(in) :: n_atoms
    integer, intent(in) :: nchains
    real(4), intent(in) :: codx(n_atoms),cody(n_atoms),codz(n_atoms)
    integer, intent(in) :: atom_ranges(:,:)
    integer, intent(in) :: res_nums(:)
    integer, intent(out):: ilandmark(nchains)     
    integer, intent(out):: loc(2)

    integer :: i, j, icou, no_of_pairs, jcou

    real(4), allocatable :: dr(:,:) 
    real(4) :: dx, dy, dz, dr_min

    !TODO: 
    !chain%atomNo_begin, chain%atomNo_end <= looks more readable than atom_range(:,:)? 
    !
    !type chain 
    !  integer, intent(in) :: atomNo_begin, atomNo_end
    !end type
    !type(chain), allocatable :: chain(:)

    print*, "Range1(mol1):",atom_ranges(1,1:2)
    print*, "Range2(mol2):",atom_ranges(2,1:2)
    allocate(dr(atom_ranges(1,2),atom_ranges(2,2)-(atom_ranges(2,1)-1)) )

    print*,"size of dr (x,y) = ", size(dr,1), size(dr,2)
    dr(:,:) = 0.0d0

    select case(nchains)
      case(1)
        print*, "'get_residue_no_for_landmark' routine was not used."

      case(2)
        print*, "This routine is using."

        !***Calculate distances between atoms of mol1 and mol2
        icou = 0
        do i = atom_ranges(1,1), atom_ranges(1,2)
          jcou = 0
          do j = atom_ranges(2,1), atom_ranges(2,2)
            jcou = jcou + 1
            dx = abs(codx(i) - codx(j))
            dy = abs(cody(i) - cody(j))
            dz = abs(codz(i) - codz(j))
            dr(i,jcou) = sqrt(dx**2 + dy**2 + dz**2)
            icou = icou + 1
          enddo
        enddo

        !***Caution: Too specialized for s100b system
        loc(1:2) = minloc(dr)
        loc(2)   = loc(2) + atom_ranges(1,2)
        print*, "Min dist = ", minval(dr), "Min loc = ",loc
        print*, "no of pairs = ", icou
        ilandmark(1) = res_nums(loc(1))
        ilandmark(2) = res_nums(loc(2))
        print*,"Landmark residue no =", ilandmark
        print*,""
        !***

    end select

  end subroutine

  !This name looks bad, "PbcCorrection" is a better name??
  subroutine ReturnAtom(n_atoms , nchains, n_residues, & 
                        codx    , cody   , codz      , & 
                        cellsize, res_nums)

    integer :: iresidue, jj, n_residues
    integer :: nchains
    integer :: imove(3)
    integer :: atom_no_begin, atom_no_end
    integer :: loc(2)
    real(4) :: rdiff
    integer :: atom_ranges(1:2,1:2) !Specialized only for two chain systems

    integer         , intent(in)    :: n_atoms
    integer         , intent(in)    :: res_nums(:) 
    double precision, intent(in)    :: cellsize(1:3) 
    real(4)         , intent(inout) :: codx(n_atoms), cody(n_atoms), codz(n_atoms)

    integer :: ilandmark(nchains), iChainBegin(nchains), iChainEnd(nchains) 


    if (IfUsed .neqv. .true.) then 
      stop "subroutine, prepare_apply_PBC, has not used yet before subroutine, ReturnAtom."
    endif

    ilandmark(:)=0; iChainBegin(:)=0; iChainEnd(:)=0

    !Caution: specialized only for S100B-CTD 
    iChainBegin(1) = 1 
    iChainEnd(1)   = 94 
    iChainBegin(2) = 95 
    iChainEnd(2)   = 108
    atom_ranges(:,:)=0 !@ 
    atom_ranges(1, 1:2)=(/ 1, 1466 /)   !atom range of chain1
    atom_ranges(2, 1:2)=(/ 1467, 1692/) !           of chain2

    !Pick the landmard atom and residue by calculating distances between atoms
    !of mol1 and mol2 if a system consists of two chains.  
    call get_residue_no_for_landmark(n_atoms,nchains,codx,cody,codz,atom_ranges,ilandmark,res_nums,loc)

    !*****************************************************
    !Move atoms at each residues if atoms go beyond a box.
    !*****************************************************
    imove(:) = 0

    !If an atom in a residue is over the box, make it back.  
    do iresidue = 1, n_residues
      atom_no_begin = first_atom(iresidue)                    !First atom number of a residues 
      atom_no_end   = last_atom(iresidue)                     !Last atom number of a residues
      
      !Caution: this if-cause must be used with the following "Temporary code".
      if (iresidue == ilandmark(1)) then 
        atom_no_begin = loc(1)
      elseif (iresidue == ilandmark(2)) then
        atom_no_begin = loc(2)
      endif

      !Caution: Maybe Here are problems. 2017-2-10
      do jj = atom_no_begin + 1, atom_no_end              !From next atom to last atom in a residue.

        !***x-direction 
        rdiff = codx(atom_no_begin) - codx(jj)        !calc. diff. betw. coordinate of first atom and jjth atom.
        imove(1) = idnint(rdiff / cellsize(1)) !imove is 0 or 1. if 1, move back atom. if 0, remain position of atom. 
        codx(jj) = codx(jj) + (imove(1) * cellsize(1)) !Update a coordinate of atom
  
        !***y-direction
        rdiff = cody(atom_no_begin) - cody(jj) !calc. diff. betw. coordinate of first atom and jjth atom.
        imove(2) = idnint(rdiff / cellsize(2))
        cody(jj) = cody(jj) + (imove(2) * cellsize(2)) 

        !***z-direction
        rdiff = codz(atom_no_begin) - codz(jj) !calc. diff. betw. coordinate of first atom and jjth atom.
        imove(3) = idnint(rdiff / cellsize(3))
        codz(jj) = codz(jj) + (imove(3) * cellsize(3)) 
      enddo
     
      !Temporary code: This is too complicated. => But, this works! 
      !This lines should be merged into a function.
      if (iresidue == ilandmark(1) .or. iresidue == ilandmark(2)) then

        print*, "Residue No corresponds to ilandmark(1) or ilandmark(2): ", iresidue

        do jj = atom_no_begin -1, first_atom(iresidue), -1
          !***x-direction 
          rdiff = codx(atom_no_begin) - codx(jj)        !calc. diff. betw. coordinate of first atom and jjth atom.
          imove(1) = idnint(rdiff / cellsize(1)) !imove is 0 or 1. if 1, move back atom. if 0, remain position of atom. 
          codx(jj) = codx(jj) + (imove(1) * cellsize(1)) !Update a coordinate of atom
  
          !***y-direction
          rdiff = cody(atom_no_begin) - cody(jj) !calc. diff. betw. coordinate of first atom and jjth atom.
          imove(2) = idnint(rdiff / cellsize(2))
          cody(jj) = cody(jj) + (imove(2) * cellsize(2)) 

          !***z-direction
          rdiff = codz(atom_no_begin) - codz(jj) !calc. diff. betw. coordinate of first atom and jjth atom.
          imove(3) = idnint(rdiff / cellsize(3))
          codz(jj) = codz(jj) + (imove(3) * cellsize(3)) 
        enddo
      endif

    enddo

    !*****************************************************
    !Move residue if a entire residue goes beyond a box.
    !*****************************************************
    imove(:) = 0
    do i = 1, nchains
      do iresidue = ilandmark(i)-1, iChainBegin(i), -1
        atom_no_begin = first_atom(iresidue)
        atom_no_end = last_atom(iresidue)
  
        !*************************************************************************************************
        !All atoms in a residue will be moved, 
        !if the distance btw 1st atom of residue(i+1) and residue(i) is over half of the length of x-cell.
        !*************************************************************************************************
        
        !If neiboring res. is far to iresi+1 res., rdiff is large. 
        rdiff  = codx(first_atom(iresidue+1)) - codx(first_atom(iresidue)) 
        imove(1)  = idnint(rdiff / cellsize(1))
        do jj = atom_no_begin, atom_no_end
          codx(jj) = codx(jj) + (imove(1) * cellsize(1))
        enddo
  
        rdiff  = cody(first_atom(iresidue+1)) - cody(first_atom(iresidue))
        imove(2)  = idnint(rdiff / cellsize(2))
        do jj = atom_no_begin, atom_no_end
          cody(jj) = cody(jj) + (imove(2) * cellsize(2))
        enddo
        
        rdiff  = codz(first_atom(iresidue+1)) - codz(first_atom(iresidue))
        imove(3)  = idnint(rdiff / cellsize(3))
        do jj = atom_no_begin, atom_no_end
          codz(jj) = codz(jj) + (imove(3) * cellsize(3))
        enddo

      enddo
    enddo

    imove(:) = 0
    do i = 1, nchains
      do iresidue = ilandmark(i)+1, iChainEnd(i)
        atom_no_begin = first_atom(iresidue)
        atom_no_end = last_atom(iresidue)
  
        rdiff = codx(first_atom(iresidue-1)) - codx(first_atom(iresidue))
        imove(1) = idnint(rdiff / cellsize(1))
        do jj = atom_no_begin, atom_no_end
          codx(jj) = codx(jj) + (imove(1) * cellsize(1))
        enddo
  
        rdiff = cody(first_atom(iresidue-1)) - cody(first_atom(iresidue))
        imove(2) = idnint(rdiff / cellsize(2))
        do jj = atom_no_begin, atom_no_end
          cody(jj) = cody(jj) + (imove(2) * cellsize(2))
        enddo
  
        rdiff = codz(first_atom(iresidue-1)) - codz(first_atom(iresidue))
        imove(3) = idnint(rdiff / cellsize(3))
        do jj = atom_no_begin, atom_no_end
          codz(jj) = codz(jj) + (imove(3) * cellsize(3))
        enddo
      enddo
    enddo
    
    !write(*,'(a)') "Note) Atoms have been returned (module periodic was uesd)."
    
    !2017-2-21
    !After this routine, we should translate com into the origin for easiness of
    !viewing structures in pymol. 


  end subroutine



end module
