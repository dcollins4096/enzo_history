module grafic_io

  use grafic_types
#ifdef ENZO
  use hdf_io
  use enzo_io
#endif  
  implicit none
!  integer, parameter :: nblock = 10
  integer, parameter :: RECORD_FLAG=4 !bytes 
  integer, private, save :: myid=0, nprocs=1
  
  interface grafic_write
     module procedure grafic_write_single, grafic_write_double
  end interface

  interface grafic_read
     module procedure grafic_read_single, grafic_read_double
  end interface

contains

  subroutine init_grafic_io
    integer ierr
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  end subroutine init_grafic_io

  ! This routine write a fftw_mpi slice of a cube 
  ! This slice has interior dimension nz, exterior nx.
  ! Beware that interior dimension is padded: nz <- 2(nz/2+1)
  ! Padded region will not be written.
  subroutine grafic_write_double(buffer, local_nz, local_z_start, nz, ny, nx, &
       filename, dim, ndims, padding_in, white_in, mask, level)
    
    ! Arguments
    implicit none
    include 'mpif.h'
    real(dp), allocatable, dimension(:), intent(inout) :: buffer
    integer, intent(in) :: local_nz, local_z_start, nz, ny, nx
    integer, intent(in) :: dim, ndims
    character(len=128), intent(in) :: filename
    integer, optional :: level
    logical, optional, dimension(:,:,:), allocatable, intent(inout) :: mask
    !integer, intent(in) :: filenamelen
    logical, optional :: padding_in
    logical, optional :: white_in

    ! Local variables
    real(sp), allocatable, dimension(:) :: tampon
    integer :: record_size
    integer(RECORD_FLAG) :: taille_tampon
    integer :: n2x
    integer(i8b) :: index, index2, offset, length, toto
    integer :: idummy=0
    integer :: i, j, k, ierr
    integer :: this_level
    logical :: padding, white

    ! Variables for Enzo I/O
    character(128) :: enzoname, dset_name
    integer :: local_dims(3), dimsf(3), nxx, lnz, mask_dims(3)
    real(dp), allocatable, dimension(:) :: work
    logical, dimension(:,:,:), allocatable :: mmask

    padding=.false.
    white=.false.
    this_level=0
    if (present(padding_in)) padding = padding_in
    if (present(padding_in) .and. myid.eq.0) write(*,*) 'padding is on'
    if (present(white_in)) white = white_in
    if (present(level)) this_level = level
    ! 2*4 for leading and trailing record sizes 
    n2x = 2*(nx/2+1)

#ifdef ENZO
    if (padding) then 
       nxx = n2x
       dset_name = filename
    else
       nxx = nx
       dset_name = filename !get_dsetname(filename)
    endif
    dimsf = (/nxx, ny, nz/)

    ! If there's no mask, create one with all .true.
    if (.not.present(mask)) then
       allocate(mmask(n2x,ny,local_nz))
       local_dims = (/n2x, ny, local_nz/)
       mmask = .true.
    else
       mask_dims = shape(mask)
       allocate(mmask(mask_dims(1), mask_dims(2), mask_dims(3)))
       local_dims = mask_dims
       mmask = mask
    endif

    CALL crop_buffer(buffer, work, local_dims, mmask)
    CALL init_phdf5(filename)
    if (dim .eq. 0) then
       CALL create_hdf5_dset(dset_name, H5T_NATIVE_DOUBLE, 3, dimsf, &
            dim, ndims, level=this_level)
    else
       CALL open_hdf5_dset(dset_name)
    endif
    
    CALL write_phdf5_data(work, local_z_start, local_nz, local_dims(3), &
         local_dims(2), local_dims(1), dim, ndims, masked_in=.true.)
    CALL close_phdf5_file
    deallocate(work)
    deallocate(mmask)
#else

    if (padding) then
       write(*,*) 'padding is true in grafic_write'
       allocate(tampon(ny*n2x))
       taille_tampon = ny*n2x*sp ! in bytes
    else
       allocate(tampon(ny*nx))
       taille_tampon = ny*nx*sp ! in bytes
    endif

    if (white) then
       offset = 4*4 + 2*RECORD_FLAG
    else
       offset = 3*4+8*sp + 2*RECORD_FLAG ! First record
    endif

    if (padding) then
       offset = offset + local_z_start*int(ny*n2x*sp + 2*RECORD_FLAG,8)
    else
       offset = offset + local_z_start*int(ny*nx*sp + 2*RECORD_FLAG,8)
    endif

    do k=1,local_nz
       
       index2=1
       ! This loop to (eventually) get rid of fftw padding zone
       if (padding) then
          do j=1,ny
             do i=1,n2x
                index = ((k-1)*ny+j-1)*n2x+i
                tampon(index2) = real(buffer(index),kind=sp)
                index2 = index2 + 1
             enddo
          enddo
       else
          do j=1,ny
             do i=1,nx
                index = ((k-1)*ny+j-1)*n2x+i
                tampon(index2) = real(buffer(index),kind=sp)
                index2 = index2 + 1
             enddo
          enddo
       endif
       ! First write f... Fortran record length
       length=RECORD_FLAG
       call f77_parallel_write(trim(filename),len_trim(filename), length &
            & ,offset,taille_tampon)
       offset = offset+RECORD_FLAG
       if (padding) then
          length = ny*n2x*sp
       else
          length = ny*nx*sp
       endif
       call f77_parallel_write(trim(filename),len_trim(filename), length, &
            & offset, tampon)
       offset = offset+length
       ! Again ..
       length=RECORD_FLAG
       call f77_parallel_write(trim(filename),len_trim(filename), length &
            & ,offset,taille_tampon)
       offset = offset+RECORD_FLAG
    enddo
    deallocate(tampon)
#endif
    if (present(mask)) then
       if (allocated(mask)) deallocate(mask)
    endif

  end subroutine grafic_write_double
       
  subroutine grafic_write_single(buffer, local_nz, local_z_start, nz, ny, nx, &
       filename, dim, ndims, padding_in, white_in, mask, level)
    
    ! Arguments
    implicit none
    include 'mpif.h'
    real(sp), dimension(:), allocatable, intent(inout) :: buffer
    integer, intent(in) :: local_nz, local_z_start, nz, ny, nx
    integer, intent(in) :: dim, ndims
    character(len=128), intent(in) :: filename
    integer, optional :: level
    logical, optional, dimension(:,:,:), allocatable, intent(inout) :: mask
    !integer, intent(in) :: filenamelen
    logical, optional :: padding_in
    logical, optional :: white_in

    ! Local variables
    real(sp), allocatable, dimension(:) :: tampon
    integer :: record_size
    integer(RECORD_FLAG) :: taille_tampon
    integer :: n2x
    integer(i8b) :: index, index2, offset, length, toto
    integer :: idummy=0
    integer :: i,j,k
    integer :: this_level
    logical :: padding, white

    ! Variables for Enzo I/O
    character(128) :: enzoname, dset_name
    integer :: local_dims(3), dimsf(3), nxx, lnz, mask_dims(3)
    real(sp), allocatable, dimension(:) :: work
    logical, dimension(:,:,:), allocatable :: mmask

    ! If there's no mask, create one with all .true.
    allocate(mmask(nx,ny,local_nz))
    if (.not.present(mask)) then
       mmask = .true.
    else
       mmask = mask
    endif

    padding=.false.
    white=.false.
    this_level=0
    if (present(padding_in)) padding = padding_in
    if (present(padding_in)) write(*,*) 'padding is on'
    if (present(white_in)) white = white_in
    if (present(level)) this_level = level
    ! 2*4 for leading and trailing record sizes 
    n2x = 2*(nx/2+1)

#ifdef ENZO
    if (padding) then 
       nxx = n2x
       dset_name = filename
    else
       nxx = nx
       dset_name = filename !get_dsetname(filename)
    endif
    dimsf = (/nxx, ny, nz/)

    ! If there's no mask, create one with all .true.
    if (.not.present(mask)) then
       allocate(mmask(n2x,ny,local_nz))
       local_dims = (/n2x, ny, local_nz/)
       mmask = .true.
    else
       mask_dims = shape(mask)
       allocate(mmask(mask_dims(1), mask_dims(2), mask_dims(3)))
       local_dims = mask_dims
       mmask = mask
    endif

    CALL crop_buffer(buffer, work, local_dims, mmask)
    CALL init_phdf5(filename)
    if (dim .eq. 0) then
       CALL create_hdf5_dset(dset_name, H5T_NATIVE_REAL, 3, dimsf, &
            dim, ndims, level=this_level)
    else
       CALL open_hdf5_dset(dset_name)
    endif

    CALL write_phdf5_data(work, local_z_start, local_nz, local_dims(3), &
         local_dims(2), local_dims(1), dim, ndims, masked_in=.true.)
    CALL close_phdf5_file
    deallocate(work)
    deallocate(mmask)
#else

    if (padding) then
       write(*,*) 'padding is true in grafic_write'
       allocate(tampon(ny*n2x))
       taille_tampon = ny*n2x*sp ! in bytes
    else
       allocate(tampon(ny*nx))
       taille_tampon = ny*nx*sp ! in bytes
    endif

    if (white) then
       offset = 4*4 + 2*RECORD_FLAG
    else
       offset = 3*4+8*sp + 2*RECORD_FLAG ! First record
    endif

    if (padding) then
       offset = offset + local_z_start*int(ny*n2x*sp + 2*RECORD_FLAG,8)
    else
       offset = offset + local_z_start*int(ny*nx*sp + 2*RECORD_FLAG,8)
    endif

    do k=1,local_nz
       
       index2=1
       ! This loop to (eventually) get rid of fftw padding zone
       if (padding) then
       do j=1,ny
          do i=1,n2x
             index = ((k-1)*ny+j-1)*n2x+i
             tampon(index2) = real(buffer(index),kind=sp)
             index2 = index2 + 1
          enddo
       enddo
       else
          do j=1,ny
             do i=1,nx
                index = ((k-1)*ny+j-1)*n2x+i
                tampon(index2) = real(buffer(index),kind=sp)
                index2 = index2 + 1
             enddo
          enddo
       endif

       ! First write f... Fortran record length
       length=RECORD_FLAG
       call f77_parallel_write(trim(filename),len_trim(filename), length &
            & ,offset,taille_tampon)
       offset = offset+RECORD_FLAG
       if (padding) then
          length = ny*n2x*sp
       else
          length = ny*nx*sp
       endif
       call f77_parallel_write(trim(filename),len_trim(filename), length, &
            & offset, tampon)
       offset = offset+length
       ! Again ..
       length=RECORD_FLAG
       call f77_parallel_write(trim(filename),len_trim(filename), length &
            & ,offset,taille_tampon)
       offset = offset+RECORD_FLAG
    enddo
    deallocate(tampon)
#endif
    deallocate(mmask)
    if (present(mask)) then
       if (allocated(mask)) deallocate(mask)
    endif

  end subroutine grafic_write_single

  subroutine grafic_read_double(buffer, local_nz, local_z_start, nz, ny, nx, &
       filename, dim, padding_in, white_in, serial_in)

    ! Arguments
    implicit none
    include 'mpif.h'
    integer, intent(in) :: local_nz, local_z_start, nz, ny, nx
    integer, intent(in) :: dim
    real(dp), dimension(local_nz*ny*2*(nx/2+1)), intent(out) :: buffer
    character(len=128), intent(in) :: filename
    logical, optional :: padding_in ! Read padded zone or not
    logical, optional :: white_in ! Different header for white noise file
    logical, optional :: serial_in  ! Serial HDF5 

    ! Local variables
    real(sp), allocatable, dimension(:) :: tampon
    integer :: record_size
    integer(RECORD_FLAG) :: taille_tampon
    integer :: n2x
    integer(i8b) :: index, index2, offset, length
    integer :: idummy=0
    integer :: i,j,k
    logical :: padding, white, serial

    ! Variables for Enzo I/O
    character(128) :: enzoname, dset_name
    integer :: nxx

    padding = .false.
    white = .false.
    serial = .false.
    if (present(padding_in)) padding = padding_in
    if (present(white_in)) white = white_in
    if (present(serial_in)) serial = serial_in
    ! 2*4 for leading and trailing record sizes 
    n2x = 2*(nx/2+1)

#ifdef ENZO
    if (padding) then
       nxx = n2x
       dset_name = filename
    else
       nxx = nx
       dset_name = get_dsetname(filename)
    endif
    CALL init_phdf5(filename, serial_in=serial)
    CALL open_hdf5_dset(dset_name)  ! dataset name = enzoname = filename
    CALL read_phdf5_data(buffer, local_z_start, local_nz, ny, nxx, dim, &
         padding_in=padding)
    CALL close_phdf5_file
#else

    if (padding) then
       allocate(tampon(ny*n2x))
    else 
       allocate(tampon(ny*nx))
    endif
    taille_tampon = 0

    if (white) then
       offset = 4*4 + 2*RECORD_FLAG
    else
       offset = 3*4+8*sp + 2*RECORD_FLAG ! First record
    endif
    if (padding) then
       offset = offset + local_z_start*int(ny*n2x*sp + 2*RECORD_FLAG,8)
    else 
       offset = offset + local_z_start*int(ny*nx*sp+2*RECORD_FLAG,8)
    endif

    do k=1,local_nz

       length = RECORD_FLAG
       ! First read f... Fortran record length
       call f77_parallel_read(trim(filename),len_trim(filename), length &
            & ,offset,taille_tampon)
       if (padding) then
           if (taille_tampon /= ny*n2x*sp) then
             print*,'Record size ',taille_tampon,' is different from expected', &
                  & ny*n2x*sp
             stop
          endif
       else
           if (taille_tampon /= ny*nx*sp) then
             print*,'Record size ',taille_tampon,' is different from expected', &
                  & ny*nx*sp
             stop
          endif
       endif
       offset = offset+RECORD_FLAG

       if (padding) then
          length = ny*n2x*sp
       else 
          length = ny*nx*sp
       endif
       call f77_parallel_read(trim(filename),len_trim(filename), length, &
            & offset, tampon)
       offset = offset+length
       ! Again ..
       length=RECORD_FLAG
       call f77_parallel_read(trim(filename),len_trim(filename), length &
            & ,offset,taille_tampon)
       offset = offset+RECORD_FLAG

       ! This loop to (eventually) get rid of fftw padding zone
       index2=1
       if (padding) then
          do j=1,ny
             do i=1,n2x ! Read padding (Nyquist) : reading in k space ...
                index=((k-1)*ny+j-1)*n2x+i
                buffer(index) = real(tampon(index2),kind=dp)
                index2 = index2 + 1
             enddo
          enddo
       else
          do j=1,ny
             do i=1,nx ! Do not read padding : reading in real space ...
                index = ((k-1)*ny+j-1)*n2x+i
                buffer(index) = real(tampon(index2),kind=dp)
                index2 = index2 + 1
             enddo
          enddo
       endif

    enddo
    deallocate(tampon)
#endif

  end subroutine grafic_read_double

  subroutine grafic_read_single(buffer, local_nz, local_z_start, nz, ny, nx, &
       filename, dim, padding_in, white_in)

    ! Arguments
    implicit none
    include 'mpif.h'
    integer, intent(in) :: local_nz, local_z_start, nz, ny, nx, dim
    real(sp), dimension(local_nz*ny*2*(nx/2+1)), intent(out) :: buffer
    character(len=128), intent(in) :: filename
    logical, optional :: padding_in ! Read padded zone or not
    logical, optional :: white_in ! Different header for white noise file
    ! Local variables
    real(sp), allocatable, dimension(:) :: tampon
    integer :: record_size
    integer(RECORD_FLAG) :: taille_tampon
    integer :: n2x
    integer(i8b) :: index, index2, offset, length
    integer :: idummy=0
    integer :: i,j,k
    logical :: padding
    logical :: white

    ! Variables for Enzo I/O
    character(128) :: enzoname, dset_name
    integer nxx

    padding = .false.
    white = .false.
    if (present(padding_in)) padding = padding_in
    if (present(white_in)) white = white_in
    ! 2*4 for leading and trailing record sizes 
    n2x = 2*(nx/2+1)

#ifdef ENZO
    nxx = nx
    if (padding) then
       nxx = n2x
       dset_name = filename
    else
       dset_name = get_dsetname(filename)
    endif
    CALL init_phdf5(filename)
    CALL open_hdf5_dset(dset_name)
    CALL read_phdf5_data(buffer, local_z_start, local_nz, ny, nxx, dim, &
         padding_in=padding)
    CALL close_phdf5_file
#else

    if (padding) then
       allocate(tampon(ny*n2x))
    else 
       allocate(tampon(ny*nx))
    endif
    taille_tampon = 0

    if (white) then
       offset = 4*4 + 2*RECORD_FLAG
    else
       offset = 3*4+8*sp + 2*RECORD_FLAG ! First record
    endif
    if (padding) then
       offset = offset + local_z_start*int(ny*n2x*sp + 2*RECORD_FLAG,8)
    else 
       offset = offset + local_z_start*int(ny*nx*sp+2*RECORD_FLAG,8)
    endif

    do k=1,local_nz
       
       length = RECORD_FLAG
       ! First read f... Fortran record length
       call f77_parallel_read(trim(filename),len_trim(filename), length &
            & ,offset,taille_tampon)
       if (padding) then
           if (taille_tampon /= ny*n2x*sp) then
             print*,'Record size ',taille_tampon,' is different from expected', &
                  & ny*n2x*sp
             stop
          endif
       else
           if (taille_tampon /= ny*nx*sp) then
             print*,'Record size ',taille_tampon,' is different from expected', &
                  & ny*nx*sp
             stop
          endif
       endif
       offset = offset+RECORD_FLAG

       if (padding) then
          length = ny*n2x*sp
       else 
          length = ny*nx*sp
       endif
       call f77_parallel_read(trim(filename),len_trim(filename), length, &
            & offset, tampon)
       offset = offset+length
       ! Again ..
       length=RECORD_FLAG
       call f77_parallel_read(trim(filename),len_trim(filename), length &
            & ,offset,taille_tampon)
       offset = offset+RECORD_FLAG

       ! This loop to (eventually) get rid of fftw padding zone
       index2=1
       if (padding) then
          do j=1,ny
             do i=1,n2x ! Read padding (Nyquist) : reading in k space ...
                index=((k-1)*ny+j-1)*n2x+i
                buffer(index) = real(tampon(index2),kind=sp)
                index2 = index2 + 1
             enddo
          enddo
       else
          do j=1,ny
             do i=1,nx ! Do not read padding : reading in real space ...
                index = ((k-1)*ny+j-1)*n2x+i
                buffer(index) = real(tampon(index2),kind=sp)
                index2 = index2 + 1
             enddo
          enddo
       endif

    enddo
    deallocate(tampon)
#endif

  end subroutine grafic_read_single


  subroutine grafic_read_header(filename, head_taille, head_cosmo)
    
    type(taille), intent(out) :: head_taille
    type(cosmo), intent(out) :: head_cosmo
    character(len=128), intent(in) :: filename
    logical :: ok
    character(128) :: enzoname
    integer :: filepart, n_filepart

    inquire(file=filename,exist=ok)
    if (.not. ok) then
       print*,'File ',trim(filename),' does not exist, aborting.'
       stop
    endif

#ifdef ENZO
    CALL grafic_hdf5_read_header(filename, head_taille, head_cosmo)
#else
    open(2,file=filename,form='unformatted',status='old')
    read(2) head_taille%nx, head_taille%ny, head_taille%nz, head_taille%dx, &
         head_taille%lx, head_taille%ly, head_taille%lz, head_cosmo%astart, &
         head_cosmo%omegam, head_cosmo%omegav, head_cosmo%h0
    close(2)
#endif
    head_taille%n0 = head_taille%nx
  end subroutine grafic_read_header

  subroutine grafic_read_header_white(filename, nx, ny, nz, iseed)
    
    integer, intent(out) :: nx,ny,nz,iseed
    character(len=128), intent(in) :: filename
    logical :: ok

    inquire(file=filename,exist=ok)
    if (.not. ok) then
       print*,'File ',trim(filename),' does not exist, aborting.'
       stop
    endif

#ifdef ENZO
    CALL grafic_hdf5_read_header_white(filename, nx, ny, nz, iseed)
#else
    open(2,file=filename,form='unformatted',status='old')
    read(2) nx, ny, nz, iseed
    close(2)
#endif

  end subroutine grafic_read_header_white
  
  subroutine grafic_write_header(filename, head_taille, head_cosmo)
    
    type(taille), intent(in) :: head_taille
    type(cosmo), intent(in) :: head_cosmo
    character(len=128), intent(in) :: filename
    character(128) :: enzoname
    integer :: filepart, n_filepart

#ifdef ENZO
    CALL grafic_hdf5_write_header(filename, head_taille, head_cosmo)
#else
    open(2,file=filename,form='unformatted',status='unknown')
    write(2) head_taille%nx, head_taille%ny, head_taille%nz, head_taille%dx, &
         head_taille%lx, head_taille%ly, head_taille%lz, head_cosmo%astart, &
         head_cosmo%omegam, head_cosmo%omegav, head_cosmo%h0
    close(2)
#endif

  end subroutine grafic_write_header

  subroutine grafic_write_header_white(filename, nx, ny, nz, iseed)
    
    integer, intent(in) :: nx,ny,nz,iseed
    character(len=128), intent(in) :: filename

#ifdef ENZO
    CALL grafic_hdf5_write_header_white(filename, nx, ny, nz, iseed)
#else    
    open(2,file=filename,form='unformatted',status='unknown')
    write(2) nx, ny, nz, iseed
    close(2)
#endif

  end subroutine grafic_write_header_white

       
end module grafic_io
