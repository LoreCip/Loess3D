module ioH5

    use hdf5
    use iso_fortran_env, only: RK => real64

    implicit none
    
    private

    integer, parameter :: compression_level = 9

contains

    subroutine create_hdf5_file(file_name, n, m, l, totL, X, Y, Z, Oin, file_id)

        integer,                   intent(in) :: n, m, l, totL
        real(RK), dimension(totL), intent(in) :: X, Y, Z, Oin
      
        character(len=*), intent(in)  :: file_name
        integer(hid_t)  , intent(out) :: file_id
        
        integer :: status

        ! Open the HDF5 library
        call h5open_f(status)

        ! Create a new HDF5 file
        call h5fcreate_f(trim(file_name), H5F_ACC_TRUNC_F, file_id, status)

        ! Check for errors
        if (status /= 0) then
          write(*,*) 'Error creating HDF5 file:', trim(file_name)
          return
        end if

        ! SAVE THE GRID

        call write_real_1d_array(X, "Yq",      file_id, status)
        call write_real_1d_array(Y, "logtemp", file_id, status)
        call write_real_1d_array(Z, "lognb",   file_id, status)

        call write_real_1d_array(Oin, "f_in",  file_id, status)

        ! SAVE SOME PARAMETERS

        call write_integer(n, "n", file_id, status)
        call write_integer(m, "m", file_id, status)
        call write_integer(l, "l", file_id, status)
        
        ! Check for errors
        if (status /= 0) then
            write(*,*) 'Error saving to HDF5 file:', trim(file_name)
            return
        end if
    end subroutine create_hdf5_file

    subroutine close_hdf5_file(file_id)

      integer(hid_t), intent(in) :: file_id
      integer :: status

      ! Close the HDF5 file
      call h5fclose_f(file_id, status)

      ! Check for errors
      if (status /= 0) then
        write(*,*) 'Error closing HDF5 file'
      end if

      ! Close the HDF5 library
      call h5close_f(status)

  end subroutine close_hdf5_file

    subroutine write_integer(x, name, id, error)
        
        integer, intent(in)  :: x
        CHARACTER(*),    intent(in)  :: name
        integer(HID_T),  intent(in)  :: id
        integer,         intent(out) :: error       
        integer(HSIZE_T), parameter  :: xshape(0) = 0
        integer(HID_T) :: dspace_id
        integer(HID_T) :: dset_id
        integer(HID_T) :: dtype_id

        dtype_id = H5T_NATIVE_INTEGER  

        ! Create dataspace
        call h5screate_f(H5S_SCALAR_F, dspace_id, error)
        if (error /= 0) then
          write(*,'("cannot create HDF5 dataspace",/)')
          stop
        end if    

        ! Create dataset in file
        call h5dcreate_f(id, name, dtype_id, dspace_id, dset_id, error)
        if (error /= 0) then
          write(*,'("cannot create HDF5 dataset",/)')
          stop
        end if       

        ! Write to file
        call h5dwrite_f(dset_id, dtype_id, x, xshape, error)
        if (error /= 0) then
          write(*,'("cannot write to HDF5 file",/)')
          stop
        end if       

        ! Close dataset
        call h5dclose_f(dset_id, error)
        if (error /= 0) then
          write(*,'("cannot close HDF5 dataset",/)')
          stop
        end if       
        
        ! Close dataspace
        call h5sclose_f(dspace_id, error)
        if (error /= 0) then
          write(*,'("cannot close HDF5 dataspace",/)')
          stop
        end if

    end subroutine write_integer

    subroutine write_real_1d_array(x, name, id, error)
        real(RK),   intent(in)  :: x(:)
        character(*),   intent(in)  :: name
        integer(HID_T), intent(in)  :: id
        integer,        intent(out) :: error
      
        integer, parameter :: ndims = 1
        integer(HSIZE_T)   :: xshape(ndims)
        integer(HSIZE_T)   :: chunk(ndims)
        integer(HID_T)     :: dspace_id
        integer(HID_T)     :: dset_id
        integer(HID_T)     :: prop_id
        integer(HID_T)     :: dtype_id
      
        xshape = shape(x)
        chunk = shape(x)
        dtype_id = H5T_NATIVE_DOUBLE
      
        ! Create dataspace
        call h5screate_simple_f(ndims, xshape, dspace_id, error)
        if (error /= 0) then
          write(*,'("cannot create HDF5 dataspace",/)')
          RETURN
        end if
      
        ! Create the dataset creation property list, add the gzip
        ! compression filter and set the chunk size.
        call h5pcreate_f(H5P_DATASET_CREATE_F, prop_id, error)
        call h5pset_deflate_f(prop_id, compression_level, error)
        call h5pset_chunk_f(prop_id, ndims, chunk, error)
        if (error /= 0) then
          write(*,'("cannot create HDF5 property list",/)')
          stop
        end if
      
        ! Create dataset in file
        call h5dcreate_f(id, name, dtype_id, dspace_id, dset_id, error, prop_id)
        if (error /= 0) then
          write(*,'("cannot create HDF5 dataset",/)')
          stop
        end if
      
        ! Write to file
        call h5dwrite_f(dset_id, dtype_id, x, xshape, error)
        if (error /= 0) then
          write(*,'("cannot write to HDF5 file",/)')
          stop
        end if
      
        ! Close property list
        call h5pclose_f(prop_id, error)
        if (error /= 0) then
          write(*,'("cannot close HDF5 property list",/)')
          stop
        end if
      
        ! Close dataset
        call h5dclose_f(dset_id, error)
        if (error /= 0) then
          write(*,'("cannot close HDF5 dataset",/)')
          stop
        end if
      
        ! Close dataspace
        call h5sclose_f(dspace_id, error)
        if (error /= 0) then
          write(*,'("cannot close HDF5 dataspace",/)')
          stop
        end if
      
      end subroutine write_real_1d_array
    
end module