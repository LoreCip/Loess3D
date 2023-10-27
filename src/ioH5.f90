module ioH5
    
  use hdf5
  use iso_fortran_env, only: RK => real64, RKS => real32
  use iso_c_binding, only: c_loc

  implicit none

  private
  integer, parameter :: compression_level = 9

  public :: write_to_hdf5,    &
            read_from_hdf5,   &
            open_hdf5file,    &
            create_hdf5file,  &
            close_hdf5file,   &
            open_hdf5group,   &
            create_hdf5group, &
            close_hdf5group,  &
            hid_t

  interface write_to_hdf5
      module procedure write_real_kind4,              & ! real(4)
                       write_real_kind8,              & ! real(8)
                       write_real_1d_array_kind4,     & ! 1d real(4) arrays
                       write_real_1d_array_kind8,     & ! 1d real(8) arrays
                       write_real_2d_array_kind4,     & ! 2d real(4) arrays
                       write_real_2d_array_kind8,     & ! 2d real(8) arrays
                       write_real_3d_array_kind4,     & ! 3d real(4) arrays
                       write_real_3d_array_kind8,     & ! 3d real(8) arrays
                       write_integer_kind4,           & ! integer(4)
                       write_integer_1d_array_kind1,  & ! 1d integer(1) arrays
                       write_integer_1d_array_kind4,  & ! 1d integer(4) arrays
                       write_string                     ! strings
  end interface write_to_hdf5
      
  interface read_from_hdf5
      module procedure read_real_kind4,              & ! real(4)
                       read_real_kind8,              & ! real(8)
                       read_real_1d_array_kind4,     & ! 1d real(4) arrays
                       read_real_1d_array_kind8,     & ! 1d real(8) arrays
                       read_real_2d_array_kind4,     & ! 2d real(4) arrays
                       read_real_2d_array_kind8,     & ! 2d real(8) arrays
                       read_real_3d_array_kind4,     & ! 3d real(4) arrays
                       read_real_3d_array_kind8,     & ! 3d real(8) arrays
                       read_integer_kind4,           & ! integer(4)
                       read_integer_1d_array_kind1,  & ! 1d integer(1) arrays
                       read_integer_1d_array_kind4,  & ! 1d integer(4) arrays
                       read_string                     ! strings
  end interface read_from_hdf5
      
contains
      
  subroutine create_hdf5group(file_id, groupname, group_id, error)
      character(len=*), intent(in)  :: groupname
      integer(hid_t),   intent(in)  :: file_id
      integer(hid_t),   intent(out) :: group_id
      integer,          intent(out) :: error
  
      call h5gcreate_f(file_id, groupname, group_id, error)
  
  end subroutine create_hdf5group

  subroutine open_hdf5group(file_id, groupname, group_id, error)
      character(len=*), intent(in)  :: groupname
      integer(hid_t),   intent(in)  :: file_id
      integer(hid_t),   intent(out) :: group_id
      integer,          intent(out) :: error
      call h5gopen_f(file_id, groupname, group_id, error)
  end subroutine open_hdf5group

  subroutine close_hdf5group(group_id, error)
      integer(hid_t),   intent(in)  :: group_id
      integer,          intent(out) :: error
      call h5gclose_f(group_id, error)
  end subroutine close_hdf5group

  subroutine check(error, error_msg)
  
      integer, intent(in) :: error
      character(len=*), intent(in) :: error_msg
  
      if (error /= 0) then
          stop error_msg
      end if
    
  end subroutine check

  subroutine create_hdf5file(filename, file_id, error)
      character(len=*), intent(in)  :: filename
      integer(hid_t),   intent(out) :: file_id
      integer,          intent(out) :: error
      integer :: filter_info
      integer :: filter_info_both
      logical :: avail
      
      ! initialise hdf5
      call h5open_f(error)
      call check(error, "Cannot initialise hdf5")
      
      ! check if gzip compression is available.
      call h5zfilter_avail_f(h5z_filter_deflate_f, avail, error)
      if (.not.avail) then
        stop "gzip filter not available"
      endif
    
      call h5zget_filter_info_f(h5z_filter_deflate_f, filter_info, error)
      filter_info_both=ior(h5z_filter_encode_enabled_f, h5z_filter_decode_enabled_f)
      if (filter_info /= filter_info_both) then
        stop "gzip filter not available for encoding and decoding"
      endif
    
      ! create file
      call h5fcreate_f(filename, h5f_acc_trunc_f, file_id, error)
      call check(error, "Cannot create hdf5 file")
        
  end subroutine create_hdf5file

  subroutine open_hdf5file(filename, file_id, error)
       character(len=*), intent(in)  :: filename
       integer(hid_t),   intent(out) :: file_id
       integer,          intent(out) :: error
       ! initialise hdf5
       call h5open_f(error)
       call check(error, "Cannot initialise hdf5")

       ! open file
       call h5fopen_f(filename, h5f_acc_rdwr_f, file_id, error)
       call check(error, "Cannot opne hdf5 file")

  end subroutine open_hdf5file

  subroutine close_hdf5file(file_id, error)
      integer(hid_t), intent(in)  :: file_id
      integer,        intent(out) :: error
      ! close file
      call h5fclose_f(file_id, error)
      call check(error, "Cannot close hdf5 file")
      
      ! close hdf5
      call h5close_f(error)
      call check(error, "Cannot close hdf5")
      
  end subroutine close_hdf5file

  subroutine write_real_kind4(x, name, id, error)
      real(RKS),   intent(in)  :: x
      character(*),   intent(in)  :: name
      integer(hid_t), intent(in)  :: id
      integer,        intent(out) :: error
      integer(hsize_t), parameter  :: xshape(0) = 0
      integer(hid_t) :: dspace_id
      integer(hid_t) :: dset_id
      integer(hid_t) :: dtype_id

      dtype_id = h5t_native_real
      ! create dataspace
      call h5screate_f(h5s_scalar_f, dspace_id, error)
      call check(error, "Cannot create hdf5 dataspace")
      
      ! create dataset in file
      call h5dcreate_f(id, name, dtype_id, dspace_id, dset_id, error)
      call check(error, "Cannot create hdf5 dataset")
      
      ! write to file
      call h5dwrite_f(dset_id, dtype_id, x, xshape, error)
      call check(error, "Cannot write to hdf5 file")
      
      ! close dataset
      call h5dclose_f(dset_id, error)
      call check(error, "Cannot close hdf5 dataset")
      
      ! close dataspace
      call h5sclose_f(dspace_id, error)
      call check(error, "Cannot close hdf5 dataspace")

  end subroutine write_real_kind4

  subroutine write_real_kind8(x, name, id, error)
       real(RK),   intent(in)  :: x
       character(*),   intent(in)  :: name
       integer(hid_t), intent(in)  :: id
       integer,        intent(out) :: error
       integer(hsize_t), parameter  :: xshape(0) = 0
       integer(hid_t) :: dspace_id
       integer(hid_t) :: dset_id
       integer(hid_t) :: dtype_id

       dtype_id = h5t_native_double
       ! create dataspace
       call h5screate_f(h5s_scalar_f, dspace_id, error)
       call check(error, "Cannot create hdf5 dataspace")

       ! create dataset in file
       call h5dcreate_f(id, name, dtype_id, dspace_id, dset_id, error)
       call check(error, "Cannot create hdf5 dataset")

       ! write to file
       call h5dwrite_f(dset_id, dtype_id, x, xshape, error)
       call check(error, "Cannot write to hdf5 file")

       ! close dataset
       call h5dclose_f(dset_id, error)
       call check(error, "Cannot close hdf5 dataset")

       ! close dataspace
       call h5sclose_f(dspace_id, error)
       call check(error, "Cannot close hdf5 dataspace")

  end subroutine write_real_kind8

  subroutine write_real_1d_array_kind4(x, name, id, error)
      real(RKS),   intent(in)  :: x(:)
      character(*),   intent(in)  :: name
      integer(hid_t), intent(in)  :: id
      integer,        intent(out) :: error
      integer, parameter :: ndims = 1
      integer(hsize_t)   :: xshape(ndims)
      integer(hsize_t)   :: chunk(ndims)
      integer(hid_t)     :: dspace_id
      integer(hid_t)     :: dset_id
      integer(hid_t)     :: prop_id
      integer(hid_t)     :: dtype_id

      xshape = shape(x)
      chunk = shape(x)
      dtype_id = h5t_native_real

      ! create dataspace
      call h5screate_simple_f(ndims, xshape, dspace_id, error)
      call check(error, "Cannot create hdf5 dataspace")

      ! create the dataset creation property list, add the gzip
      ! compression filter and set the chunk size.
      call h5pcreate_f(h5p_dataset_create_f, prop_id, error)
      call h5pset_deflate_f(prop_id, compression_level, error)
      call h5pset_chunk_f(prop_id, ndims, chunk, error)
      call check(error, "Cannot create hdf5 property list")
      
      ! create dataset in file
      call h5dcreate_f(id, name, dtype_id, dspace_id, dset_id, error, prop_id)
      call check(error, "Cannot create hdf5 dataset")

      ! write to file
      call h5dwrite_f(dset_id, dtype_id, x, xshape, error)
      call check(error, "Cannot write to hdf5 file")
      
      ! close property list
      call h5pclose_f(prop_id, error)
      call check(error, "Cannot close hdf5 property list")

      ! close dataset
      call h5dclose_f(dset_id, error)
      call check(error, "Cannot close hdf5 dataset")
      
      ! close dataspace
      call h5sclose_f(dspace_id, error)
      call check(error, "Cannot close hdf5 dataspace")

  end subroutine write_real_1d_array_kind4

  subroutine write_real_1d_array_kind8(x, name, id, error)
      real(RK),   intent(in)  :: x(:)
      character(*),   intent(in)  :: name
      integer(hid_t), intent(in)  :: id
      integer,        intent(out) :: error
      integer, parameter :: ndims = 1
      integer(hsize_t)   :: xshape(ndims)
      integer(hsize_t)   :: chunk(ndims)
      integer(hid_t)     :: dspace_id
      integer(hid_t)     :: dset_id
      integer(hid_t)     :: prop_id
      integer(hid_t)     :: dtype_id

      xshape = shape(x)
      chunk = shape(x)
      dtype_id = h5t_native_double

      ! create dataspace
      call h5screate_simple_f(ndims, xshape, dspace_id, error)
      call check(error, "Cannot create hdf5 dataspace")
      
      ! create the dataset creation property list, add the gzip
      ! compression filter and set the chunk size.
      call h5pcreate_f(h5p_dataset_create_f, prop_id, error)
      call h5pset_deflate_f(prop_id, compression_level, error)
      call h5pset_chunk_f(prop_id, ndims, chunk, error)
      call check(error, "Cannot create hdf5 property list")
      
      ! create dataset in file
      call h5dcreate_f(id, name, dtype_id, dspace_id, dset_id, error, prop_id)
      call check(error, "Cannot create hdf5 dataset")
      
      ! write to file
      call h5dwrite_f(dset_id, dtype_id, x, xshape, error)
      call check(error, "Cannot write to hdf5 file")
      
      ! close property list
      call h5pclose_f(prop_id, error)
      call check(error, "Cannot close hdf5 property list")
      
      ! close dataset
      call h5dclose_f(dset_id, error)
      call check(error, "Cannot close hdf5 dataset")
      
      ! close dataspace
      call h5sclose_f(dspace_id, error)
      call check(error, "Cannot close hdf5 dataspace")
      
  end subroutine write_real_1d_array_kind8

  subroutine write_real_2d_array_kind4(x, name, id, error)
      real(RKS),   intent(in)  :: x(:,:)
      character(*),   intent(in)  :: name
      integer(hid_t), intent(in)  :: id
      integer,        intent(out) :: error
      integer, parameter :: ndims = 2
      integer(hsize_t)   :: xshape(ndims)
      integer(hsize_t)   :: chunk(ndims)
      integer(hid_t)     :: dspace_id
      integer(hid_t)     :: dset_id
      integer(hid_t)     :: prop_id
      integer(hid_t)     :: dtype_id

      xshape = shape(x)
      chunk = shape(x)
      dtype_id = h5t_native_real

      ! create dataspace
      call h5screate_simple_f(ndims, xshape, dspace_id, error)
      call check(error, "Cannot create hdf5 dataspace")

      ! create the dataset creation property list, add the gzip
      ! compression filter and set the chunk size.
      call h5pcreate_f(h5p_dataset_create_f, prop_id, error)
      call h5pset_deflate_f(prop_id, compression_level, error)
      call h5pset_chunk_f(prop_id, ndims, chunk, error)
      call check(error, "Cannot create hdf5 property list")

      ! create dataset in file
      call h5dcreate_f(id, name, dtype_id, dspace_id, dset_id, error, prop_id)
      call check(error, "Cannot create hdf5 dataset")

      ! write to file
      call h5dwrite_f(dset_id, dtype_id, x, xshape, error)
      call check(error, "Cannot write to hdf5 file")

      ! close property list
      call h5pclose_f(prop_id, error)
      call check(error, "Cannot close hdf5 property list")

      ! close dataset
      call h5dclose_f(dset_id, error)
      call check(error, "Cannot close hdf5 dataset")

      ! close dataspace
      call h5sclose_f(dspace_id, error)
      call check(error, "Cannot close hdf5 dataspace")

  end subroutine write_real_2d_array_kind4

  subroutine write_real_2d_array_kind8(x, name, id, error)
      real(RK),   intent(in)  :: x(:,:)
      character(*),   intent(in)  :: name
      integer(hid_t), intent(in)  :: id
      integer,        intent(out) :: error
      integer, parameter :: ndims = 2
      integer(hsize_t)   :: xshape(ndims)
      integer(hsize_t)   :: chunk(ndims)
      integer(hid_t)     :: dspace_id
      integer(hid_t)     :: dset_id
      integer(hid_t)     :: prop_id
      integer(hid_t)     :: dtype_id

      xshape = shape(x)
      chunk = shape(x)
      dtype_id = h5t_native_double

      ! create dataspace
      call h5screate_simple_f(ndims, xshape, dspace_id, error)
      call check(error, "Cannot create hdf5 dataspace")

      ! create the dataset creation property list, add the gzip
      ! compression filter and set the chunk size.
      call h5pcreate_f(h5p_dataset_create_f, prop_id, error)
      call h5pset_deflate_f(prop_id, compression_level, error)
      call h5pset_chunk_f(prop_id, ndims, chunk, error)
      call check(error, "Cannot create hdf5 property list")

      ! create dataset in file
      call h5dcreate_f(id, name, dtype_id, dspace_id, dset_id, error, prop_id)
      call check(error, "Cannot create hdf5 dataset")

      ! write to file
      call h5dwrite_f(dset_id, dtype_id, x, xshape, error)
      call check(error, "Cannot write to hdf5 file")

      ! close property list
      call h5pclose_f(prop_id, error)
      call check(error, "Cannot close hdf5 property list")

      ! close dataset
      call h5dclose_f(dset_id, error)
      call check(error, "Cannot close hdf5 dataset")

      ! close dataspace
      call h5sclose_f(dspace_id, error)
      call check(error, "Cannot close hdf5 dataspace")

  end subroutine write_real_2d_array_kind8

  subroutine write_real_3d_array_kind4(x, name, id, error)
      real(RKS),   intent(in)  :: x(:,:,:)
      character(*),   intent(in)  :: name
      integer(hid_t), intent(in)  :: id
      integer,        intent(out) :: error
      integer, parameter :: ndims = 3
      integer(hsize_t)   :: xshape(ndims)
      integer(hsize_t)   :: chunk(ndims)
      integer(hid_t)     :: dspace_id
      integer(hid_t)     :: dset_id
      integer(hid_t)     :: prop_id
      integer(hid_t)     :: dtype_id

      xshape = shape(x)
      chunk = shape(x)
      dtype_id = h5t_native_real

      ! create dataspace
      call h5screate_simple_f(ndims, xshape, dspace_id, error)
      call check(error, "Cannot create hdf5 dataspace")

      ! create the dataset creation property list, add the gzip
      ! compression filter and set the chunk size.
      call h5pcreate_f(h5p_dataset_create_f, prop_id, error)
      call h5pset_deflate_f(prop_id, compression_level, error)
      call h5pset_chunk_f(prop_id, ndims, chunk, error)
      call check(error, "Cannot create hdf5 property list")

      ! create dataset in file
      call h5dcreate_f(id, name, dtype_id, dspace_id, dset_id, error, prop_id)
      call check(error, "Cannot create hdf5 dataset")

      ! write to file
      call h5dwrite_f(dset_id, dtype_id, x, xshape, error)
      call check(error, "Cannot write to hdf5 file")

      ! close property list
      call h5pclose_f(prop_id, error)
      call check(error, "Cannot close hdf5 property list")

      ! close dataset
      call h5dclose_f(dset_id, error)
      call check(error, "Cannot close hdf5 dataset")

      ! close dataspace
      call h5sclose_f(dspace_id, error)
      call check(error, "Cannot close hdf5 dataspace")

  end subroutine write_real_3d_array_kind4

  subroutine write_real_3d_array_kind8(x, name, id, error)
      real(RK),   intent(in)  :: x(:,:,:)
      character(*),   intent(in)  :: name
      integer(hid_t), intent(in)  :: id
      integer,        intent(out) :: error
      integer, parameter :: ndims = 3
      integer(hsize_t)   :: xshape(ndims)
      integer(hsize_t)   :: chunk(ndims)
      integer(hid_t)     :: dspace_id
      integer(hid_t)     :: dset_id
      integer(hid_t)     :: prop_id
      integer(hid_t)     :: dtype_id
      xshape = shape(x)
      chunk = shape(x)
      dtype_id = h5t_native_double
      ! create dataspace
      call h5screate_simple_f(ndims, xshape, dspace_id, error)
      call check(error, "Cannot create hdf5 dataspace")

      ! create the dataset creation property list, add the gzip
      ! compression filter and set the chunk size.
      call h5pcreate_f(h5p_dataset_create_f, prop_id, error)
      call h5pset_deflate_f(prop_id, compression_level, error)
      call h5pset_chunk_f(prop_id, ndims, chunk, error)
      call check(error, "Cannot create hdf5 property list")

      ! create dataset in file
      call h5dcreate_f(id, name, dtype_id, dspace_id, dset_id, error, prop_id)
      call check(error, "Cannot create hdf5 dataset")

      ! write to file
      call h5dwrite_f(dset_id, dtype_id, x, xshape, error)
      call check(error, "Cannot write to hdf5 file")

      ! close property list
      call h5pclose_f(prop_id, error)
      call check(error, "Cannot close hdf5 property list")

      ! close dataset
      call h5dclose_f(dset_id, error)
      call check(error, "Cannot close hdf5 dataset")

      ! close dataspace
      call h5sclose_f(dspace_id, error)
      call check(error, "Cannot close hdf5 dataspace")

  end subroutine write_real_3d_array_kind8

  subroutine write_integer_kind4(x, name, id, error)
      integer(RKS), intent(in)  :: x
      character(*),    intent(in)  :: name
      integer(hid_t),  intent(in)  :: id
      integer,         intent(out) :: error
      integer(hsize_t), parameter  :: xshape(0) = 0
      integer(hid_t) :: dspace_id
      integer(hid_t) :: dset_id
      integer(hid_t) :: dtype_id

      dtype_id = h5t_native_integer
      
      ! create dataspace
      call h5screate_f(h5s_scalar_f, dspace_id, error)
      call check(error, "Cannot create hdf5 dataspace")

      ! create dataset in file
      call h5dcreate_f(id, name, dtype_id, dspace_id, dset_id, error)
      call check(error, "Cannot create hdf5 dataset")

      ! write to file
      call h5dwrite_f(dset_id, dtype_id, x, xshape, error)
      call check(error, "Cannot write to hdf5 file")

      ! close dataset
      call h5dclose_f(dset_id, error)
      call check(error, "Cannot close hdf5 dataset")

      ! close dataspace
      call h5sclose_f(dspace_id, error)
      call check(error, "Cannot close hdf5 dataspace")
  end subroutine write_integer_kind4

  subroutine write_integer_1d_array_kind4(x, name, id, error)
      integer(RKS), intent(in)  :: x(:)
      character(*),    intent(in)  :: name
      integer(hid_t),  intent(in)  :: id
      integer,         intent(out) :: error

      integer, parameter :: ndims = 1
      integer(hsize_t)   :: xshape(ndims)
      integer(hsize_t)   :: chunk(ndims)
      integer(hid_t)     :: dspace_id
      integer(hid_t)     :: dset_id
      integer(hid_t)     :: prop_id
      integer(hid_t)     :: dtype_id

      xshape = shape(x)
      chunk = shape(x)
      dtype_id = h5t_native_integer

      ! create dataspace
      call h5screate_simple_f(ndims, xshape, dspace_id, error)
      call check(error, "Cannot create hdf5 dataspace")

      ! create the dataset creation property list, add the gzip
      ! compression filter and set the chunk size.
      call h5pcreate_f(h5p_dataset_create_f, prop_id, error)
      call h5pset_deflate_f(prop_id, compression_level, error)
      call h5pset_chunk_f(prop_id, ndims, chunk, error)
      call check(error, "Cannot create hdf5 property list")

      ! create dataset in file
      call h5dcreate_f(id, name, dtype_id, dspace_id, dset_id, error, prop_id)
      call check(error, "Cannot create hdf5 dataset")

      ! write to file
      call h5dwrite_f(dset_id, dtype_id, x, xshape, error)
      call check(error, "Cannot write to hdf5 file")

      ! close property list
      call h5pclose_f(prop_id, error)
      call check(error, "Cannot close hdf5 property list")

      ! close dataset
      call h5dclose_f(dset_id, error)
      call check(error, "Cannot close hdf5 dataset")

      ! close dataspace
      call h5sclose_f(dspace_id, error)
      call check(error, "Cannot close hdf5 dataspace")

  end subroutine write_integer_1d_array_kind4

  subroutine write_integer_1d_array_kind1(x, name, id, error)
      integer(kind=1), intent(in)  :: x(:)
      character(*),    intent(in)  :: name
      integer(hid_t),  intent(in)  :: id
      integer,         intent(out) :: error
      integer, parameter :: ndims = 1
      integer(hsize_t)   :: xshape(ndims)
      integer(hsize_t)   :: chunk(ndims)
      integer(hid_t)     :: dspace_id
      integer(hid_t)     :: dset_id
      integer(hid_t)     :: prop_id
      integer(hid_t)     :: dtype_id

      xshape = shape(x)
      chunk = shape(x)
      dtype_id = h5t_std_i8le

      ! create dataspace
      call h5screate_simple_f(ndims, xshape, dspace_id, error)
      call check(error, "Cannot create hdf5 dataspace")

      ! create the dataset creation property list, add the gzip
      ! compression filter and set the chunk size.
      call h5pcreate_f(h5p_dataset_create_f, prop_id, error)
      call h5pset_deflate_f(prop_id, compression_level, error)
      call h5pset_chunk_f(prop_id, ndims, chunk, error)
      call check(error, "Cannot create hdf5 property list")

      ! create dataset in file
      call h5dcreate_f(id, name, dtype_id, dspace_id, dset_id, error, prop_id)
      call check(error, "Cannot create hdf5 dataset")

      ! write to file
      call h5dwrite_f(dset_id, dtype_id, x, xshape, error)
      call check(error, "Cannot write to hdf5 file")

      ! close property list
      call h5pclose_f(prop_id, error)
      call check(error, "Cannot close hdf5 property list")

      ! close dataset
      call h5dclose_f(dset_id, error)
      call check(error, "Cannot close hdf5 dataset")

      ! close dataspace
      call h5sclose_f(dspace_id, error)
      call check(error, "Cannot close hdf5 dataspace")

  end subroutine write_integer_1d_array_kind1

  subroutine write_string(str, name, id, error)
      character(*),    intent(in), target :: str
      character(*),    intent(in)  :: name
      integer(hid_t),  intent(in)  :: id
      integer,         intent(out) :: error
      integer, parameter :: ndims = 0
      integer(hsize_t)   :: sshape(ndims)
      integer(hid_t)     :: dspace_id
      integer(hid_t)     :: dset_id
      integer(size_t)    :: slength
      integer(hid_t)     :: filetype
      type(c_ptr)        :: cpointer

      slength = len(str)
      sshape  = shape(str)

      ! create file datatypes. save the string as fortran string
      call h5tcopy_f(h5t_fortran_s1,filetype, error)
      call h5tset_size_f(filetype, slength, error)
      call check(error, "Cannot create hdf5 datatype")

      ! create dataspace
      call h5screate_simple_f(ndims, sshape, dspace_id, error)
      call check(error, "Cannot create hdf5 dataspace")

      ! create the dataset in file
      call h5dcreate_f(id, name, filetype, dspace_id, dset_id, error)
      call check(error, "Cannot create hdf5 dataset")

      ! find c pointer
      cpointer = c_loc(str(1:1))
      ! write to file
      call h5dwrite_f(dset_id, filetype, cpointer, error)
      call check(error, "Cannot  write to hdf5 file")

      ! close dataset
      call h5dclose_f(dset_id, error)
      call check(error, "Cannot close hdf5 dataset")

      ! close dataspace
      call h5sclose_f(dspace_id, error)
      call check(error, "Cannot close hdf5 dataspace")

      ! close datatype
      call h5tclose_f(filetype, error)
      call check(error, "Cannot close hdf5 datatype")

  end subroutine write_string

  subroutine read_real_kind4(x, name, id, got, error)
      real(RKS),   intent(out) :: x
      character(*),   intent(in)  :: name
      integer(hid_t), intent(in)  :: id
      logical,        intent(out) :: got
      integer,        intent(out) :: error
      integer(hsize_t), parameter  :: xshape(0) = 0
      integer(hid_t) :: dset_id
      integer(hid_t) :: dtype_id

      dtype_id = h5t_native_real

      ! check if dataset exists
      call h5lexists_f(id, name, got, error)
      if (.not.got) return

      ! open dataset
      call h5dopen_f(id, name, dset_id, error)
      call check(error, "Cannot open hdf5 dataset")

      ! read dataset
      call h5dread_f(dset_id, dtype_id, x, xshape, error)
      call check(error, "Cannot read hdf5 dataset")

      ! close dataset
      call h5dclose_f(dset_id, error)
      call check(error, "Cannot close hdf5 dataset")

      if (error /= 0) got = .false.
  end subroutine read_real_kind4

  subroutine read_real_kind8(x, name, id, got, error)
      real(RK),   intent(out) :: x
      character(*),   intent(in)  :: name
      integer(hid_t), intent(in)  :: id
      logical,        intent(out) :: got
      integer,        intent(out) :: error
      integer(hsize_t), parameter  :: xshape(0) = 0
      integer(hid_t) :: dset_id
      integer(hid_t) :: dtype_id

      dtype_id = h5t_native_double

      ! check if dataset exists
      call h5lexists_f(id, name, got, error)
      if (.not.got) return
      
      ! open dataset
      call h5dopen_f(id, name, dset_id, error)
      call check(error, "Cannot open hdf5 dataset")

      ! read dataset
      call h5dread_f(dset_id, dtype_id, x, xshape, error)
      call check(error, "Cannot read hdf5 dataset")

      ! close dataset
      call h5dclose_f(dset_id, error)
      call check(error, "Cannot close hdf5 dataset")

      if (error /= 0) got = .false.
  end subroutine read_real_kind8

  subroutine read_real_1d_array_kind4(x, name, id, got, error)
      real(RKS),   intent(out) :: x(:)
      character(*),   intent(in)  :: name
      integer(hid_t), intent(in)  :: id
      logical,        intent(out) :: got
      integer,        intent(out) :: error
      integer, parameter :: ndims = 1
      integer(hsize_t)   :: xshape(ndims)
      integer(hid_t)     :: dset_id
      integer(hid_t)     :: dtype_id

      xshape = shape(x)
      dtype_id = h5t_native_real

      ! check if dataset exists
      call h5lexists_f(id, name, got, error)
      if (.not.got) return

      ! open dataset
      call h5dopen_f(id, name, dset_id, error)
      call check(error, "Cannot open hdf5 dataset")

      ! read dataset
      call h5dread_f(dset_id, dtype_id, x, xshape, error)
      call check(error, "Cannot read hdf5 dataset")

      ! close dataset
      call h5dclose_f(dset_id, error)
      call check(error, "Cannot close hdf5 dataset")

      if (error /= 0) got = .false.
  end subroutine read_real_1d_array_kind4

  subroutine read_real_1d_array_kind8(x, name, id, got, error)
      real(RK),   intent(out) :: x(:)
      character(*),   intent(in)  :: name
      integer(hid_t), intent(in)  :: id
      logical,        intent(out) :: got
      integer,        intent(out) :: error
      integer, parameter :: ndims = 1
      integer(hsize_t)   :: xshape(ndims)
      integer(hid_t)     :: dset_id
      integer(hid_t)     :: dtype_id

      xshape = shape(x)
      dtype_id = h5t_native_double

      ! check if dataset exists
      call h5lexists_f(id, name, got, error)
      if (.not.got) return
      
      ! open dataset
      call h5dopen_f(id, name, dset_id, error)
      call check(error, "Cannot open hdf5 dataset")

      ! read dataset
      call h5dread_f(dset_id, dtype_id, x, xshape, error)
      call check(error, "Cannot read hdf5 dataset")

      ! close dataset
      call h5dclose_f(dset_id, error)
      call check(error, "Cannot close hdf5 dataset")

      if (error /= 0) got = .false.
  end subroutine read_real_1d_array_kind8

  subroutine read_real_2d_array_kind4(x, name, id, got, error)
      real(RKS),   intent(out) :: x(:,:)
      character(*),   intent(in)  :: name
      integer(hid_t), intent(in)  :: id
      logical,        intent(out) :: got
      integer,        intent(out) :: error
      integer, parameter :: ndims = 2
      integer(hsize_t)   :: xshape(ndims)
      integer(hid_t)     :: dset_id
      integer(hid_t)     :: dtype_id

      xshape = shape(x)
      dtype_id = h5t_native_real

      ! check if dataset exists
      call h5lexists_f(id, name, got, error)
      if (.not.got) return

      ! open dataset
      call h5dopen_f(id, name, dset_id, error)
      call check(error, "Cannot open hdf5 dataset")

      ! read dataset
      call h5dread_f(dset_id, dtype_id, x, xshape, error)
      call check(error, "Cannot read hdf5 dataset")

      ! close dataset
      call h5dclose_f(dset_id, error)
      call check(error, "Cannot close hdf5 dataset")

      if (error /= 0) got = .false.
  end subroutine read_real_2d_array_kind4

  subroutine read_real_2d_array_kind8(x, name, id, got, error)
      real(RK),   intent(out) :: x(:,:)
      character(*),   intent(in)  :: name
      integer(hid_t), intent(in)  :: id
      logical,        intent(out) :: got
      integer,        intent(out) :: error
      integer, parameter :: ndims = 2
      integer(hsize_t)   :: xshape(ndims)
      integer(hid_t)     :: dset_id
      integer(hid_t)     :: dtype_id

      xshape = shape(x)
      dtype_id = h5t_native_double

      ! check if dataset exists
      call h5lexists_f(id, name, got, error)
      if (.not.got) return

      ! open dataset
      call h5dopen_f(id, name, dset_id, error)
      call check(error, "Cannot open hdf5 dataset")

      ! read dataset
      call h5dread_f(dset_id, dtype_id, x, xshape, error)
      call check(error, "Cannot read hdf5 dataset")

      ! close dataset
      call h5dclose_f(dset_id, error)
      call check(error, "Cannot close hdf5 dataset")

      if (error /= 0) got = .false.
  end subroutine read_real_2d_array_kind8

  subroutine read_real_3d_array_kind4(x, name, id, got, error)
      real(RKS),   intent(out) :: x(:,:,:)
      character(*),   intent(in)  :: name
      integer(hid_t), intent(in)  :: id
      logical,        intent(out) :: got
      integer,        intent(out) :: error
      integer, parameter :: ndims = 3
      integer(hsize_t)   :: xshape(ndims)
      integer(hid_t)     :: dset_id
      integer(hid_t)     :: dtype_id

      xshape = shape(x)
      dtype_id = h5t_native_real

      ! check if dataset exists
      call h5lexists_f(id, name, got, error)
      if (.not.got) return

      ! open dataset
      call h5dopen_f(id, name, dset_id, error)
      call check(error, "Cannot open hdf5 dataset")

      ! read dataset
      call h5dread_f(dset_id, dtype_id, x, xshape, error)
      call check(error, "Cannot read hdf5 dataset")

      ! close dataset
      call h5dclose_f(dset_id, error)
      call check(error, "Cannot close hdf5 dataset")

      if (error /= 0) got = .false.
  end subroutine read_real_3d_array_kind4

  subroutine read_real_3d_array_kind8(x, name, id, got, error)
      real(RK),   intent(out) :: x(:,:,:)
      character(*),   intent(in)  :: name
      integer(hid_t), intent(in)  :: id
      logical,        intent(out) :: got
      integer,        intent(out) :: error
      integer, parameter :: ndims = 3
      integer(hsize_t)   :: xshape(ndims)
      integer(hid_t)     :: dset_id
      integer(hid_t)     :: dtype_id

      xshape = shape(x)
      dtype_id = h5t_native_double

      ! check if dataset exists
      call h5lexists_f(id, name, got, error)
      if (.not.got) return

      ! open dataset
      call h5dopen_f(id, name, dset_id, error)
      call check(error, "Cannot open hdf5 dataset")

      ! read dataset
      call h5dread_f(dset_id, dtype_id, x, xshape, error)
      call check(error, "Cannot read hdf5 dataset")

      ! close dataset
      call h5dclose_f(dset_id, error)
      call check(error, "Cannot close hdf5 dataset")

      if (error /= 0) got = .false.
  end subroutine read_real_3d_array_kind8

  subroutine read_integer_kind4(x, name, id, got, error)
      integer(RKS), intent(out) :: x
      character(*),    intent(in)  :: name
      integer(hid_t),  intent(in)  :: id
      logical,         intent(out) :: got
      integer,         intent(out) :: error
      integer(hsize_t), parameter  :: xshape(0) = 0
      integer(hid_t) :: dset_id
      integer(hid_t) :: dtype_id

      dtype_id = h5t_native_integer

      ! check if dataset exists
      call h5lexists_f(id, name, got, error)
      if (.not.got) return
      
      ! open dataset
      call h5dopen_f(id, name, dset_id, error)
      call check(error, "Cannot open hdf5 dataset")

      ! read dataset
      call h5dread_f(dset_id, dtype_id, x, xshape, error)
      call check(error, "Cannot read hdf5 dataset")

      ! close dataset
      call h5dclose_f(dset_id, error)
      call check(error, "Cannot close hdf5 dataset")

      if (error /= 0) got = .false.
  end subroutine read_integer_kind4

  subroutine read_integer_1d_array_kind1(x, name, id, got, error)
      integer(kind=1), intent(out) :: x(:)
      character(*),    intent(in)  :: name
      integer(hid_t),  intent(in)  :: id
      logical,         intent(out) :: got
      integer,         intent(out) :: error
      integer, parameter :: ndims = 1
      integer(hsize_t)   :: xshape(ndims)
      integer(hid_t)     :: dset_id
      integer(hid_t)     :: dtype_id

      xshape = shape(x)
      dtype_id = h5t_std_i8le

      ! check if dataset exists
      call h5lexists_f(id, name, got, error)
      if (.not.got) return
      
      ! open dataset
      call h5dopen_f(id, name, dset_id, error)
      call check(error, "Cannot open hdf5 dataset")

      ! read dataset
      call h5dread_f(dset_id, dtype_id, x, xshape, error)
      call check(error, "Cannot read hdf5 dataset")

      ! close dataset
      call h5dclose_f(dset_id, error)
      call check(error, "Cannot close hdf5 dataset")

      if (error /= 0) got = .false.
  end subroutine read_integer_1d_array_kind1

  subroutine read_integer_1d_array_kind4(x, name, id, got, error)
      integer(RKS), intent(out) :: x(:)
      character(*),    intent(in)  :: name
      integer(hid_t),  intent(in)  :: id
      logical,         intent(out) :: got
      integer,         intent(out) :: error
      integer, parameter :: ndims = 1
      integer(hsize_t)   :: xshape(ndims)
      integer(hid_t)     :: dset_id
      integer(hid_t)     :: dtype_id

      xshape = shape(x)
      dtype_id = h5t_native_integer

      ! check if dataset exists
      call h5lexists_f(id, name, got, error)
      if (.not.got) return

      ! open dataset
      call h5dopen_f(id, name, dset_id, error)
      call check(error, "Cannot open hdf5 dataset")

      ! read dataset
      call h5dread_f(dset_id, dtype_id, x, xshape, error)
      call check(error, "Cannot read hdf5 dataset")

      ! close dataset
      call h5dclose_f(dset_id, error)
      call check(error, "Cannot close hdf5 dataset")

      if (error /= 0) got = .false.
  end subroutine read_integer_1d_array_kind4

  subroutine read_string(str, name, id, got, error)
      character(*),   intent(out) :: str
      character(*),   intent(in)  :: name
      integer(hid_t), intent(in)  :: id
      logical,        intent(out) :: got
      integer,        intent(out) :: error
      integer,         parameter :: dim0 = 1
      integer(size_t), parameter :: sdim = 100
      integer(hsize_t) :: dims(1) = (/dim0/)
      integer(hsize_t) :: maxdims(1)
      integer(hid_t) :: filetype, memtype, space, dset
      character(len=sdim), allocatable, target :: rdata(:)
      integer(size_t) :: size
      type(c_ptr) :: f_ptr

      ! check if dataset exists
      call h5lexists_f(id, name, got, error)
      if (.not.got) return

      call h5dopen_f(id, name, dset, error)
      call check(error, "Cannot open hdf5 dataset")

      ! get the datatype and its size.
      call h5dget_type_f(dset, filetype, error)
      call h5tget_size_f(filetype, size, error)
      call check(error, "Cannot get hdf5 datatype or size")

      ! make sure the declared length is large enough,
      ! the c string contains the null character.
      if (size > sdim+1) then
        print*,'error: character len is too small'
        stop
      end if

      ! get dataspace.
      call h5dget_space_f(dset, space, error)
      call h5sget_simple_extent_dims_f(space, dims, maxdims, error)
      call check(error, "Cannot get hdf5 dataspace")

      allocate(rdata(1:dims(1)))

      ! create the memory datatype.
      call h5tcopy_f(h5t_fortran_s1,memtype, error)
      call h5tset_size_f(memtype, sdim, error)
      call check(error, "Cannot get hdf5 memory datatype")
      
      ! read the data.
      f_ptr = c_loc(rdata(1)(1:1))

      call h5dread_f(dset, memtype, f_ptr, error,space)
      call check(error, "Cannot read hdf5 dataset")

      ! close and release resources.
      call h5dclose_f(dset, error)
      call h5sclose_f(space, error)
      call h5tclose_f(filetype, error)
      call h5tclose_f(memtype, error)
      call check(error, "Cannot close hdf5 dataset")

      str = rdata(1)

      deallocate(rdata)
      
      if (error /= 0) got = .false.
  end subroutine read_string

end module ioH5
