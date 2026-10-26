program test_div_grad_fft
  use mpi

  use m_allocator
  use m_base_backend
  use m_common, only: pi, get_argument, DIR_C, VERT, CELL, &
                      RDR_C2Z, RDR_C2X, RDR_Z2X
  use m_config, only: domain_config_t, solver_config_t
  use m_field, only: field_t
  use m_mesh
  use m_solver, only: solver_t

#ifdef CUDA
  use m_cuda_allocator
  use m_cuda_backend
  use m_cuda_common, only: SZ
#else
  use m_omp_backend
  use m_omp_common, only: SZ
#endif

  implicit none

  class(base_backend_t), pointer :: backend
  class(allocator_t), pointer :: allocator
  type(allocator_t), pointer :: host_allocator
  type(mesh_t), target :: mesh
  type(solver_t) :: solver
  class(field_t), pointer :: f, f_ii, f_final, f_check, dfdx, dfdy, dfdz, temp, f_final_check
  class(field_t), pointer :: f_host
  integer :: i, j, k, dims(3)
  real(dp) :: xloc(3), x, y, z, error_norm

#ifdef CUDA
  type(cuda_backend_t), target :: cuda_backend
  type(cuda_allocator_t), target :: cuda_allocator
  integer :: ndevs, devnum
#else
  type(omp_backend_t), target :: omp_backend
#endif

  type(allocator_t), target :: omp_allocator

  type(domain_config_t) :: domain_cfg
  type(solver_config_t) :: solver_cfg
  character(32) :: backend_name
  integer :: nrank, nproc, ierr
  logical :: use_2decomp

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'

#ifdef CUDA
  ierr = cudaGetDeviceCount(ndevs)
  ierr = cudaSetDevice(mod(nrank, ndevs)) ! round-robin
  ierr = cudaGetDevice(devnum)
  backend_name = "CUDA"
#else
  backend_name = "OMP"
#endif

  call domain_cfg%read(nml_file=get_argument(1))
  call solver_cfg%read(nml_file=get_argument(1))

  if (product(domain_cfg%nproc_dir) /= nproc) then
    if (nrank == 0) print *, 'nproc_dir specified in the input file does &
                              &not match the total number of ranks, falling &
                              &back to a 1D decomposition along Z-dir instead.'
    domain_cfg%nproc_dir = [1, 1, nproc]
  end if

  ! Decide whether 2decomp is used or not
  print *, "-- debug --"
  print *, solver_cfg%poisson_solver_type
  print *, trim(backend_name)
  use_2decomp = solver_cfg%poisson_solver_type == 'FFT' &
                .and. trim(backend_name) == 'OMP'

  mesh = mesh_t(domain_cfg%dims_global, domain_cfg%nproc_dir, &
                domain_cfg%L_global, domain_cfg%BC_x, domain_cfg%BC_y, &
                domain_cfg%BC_z, domain_cfg%stretching, domain_cfg%beta, &
                use_2decomp=use_2decomp)

#ifdef CUDA
  cuda_allocator = cuda_allocator_t(mesh%get_dims(VERT), SZ)
  allocator => cuda_allocator
  if (nrank == 0) print *, 'CUDA allocator instantiated'

  omp_allocator = allocator_t(mesh%get_dims(VERT), SZ)
  host_allocator => omp_allocator

  cuda_backend = cuda_backend_t(mesh, allocator)
  backend => cuda_backend
  if (nrank == 0) print *, 'CUDA backend instantiated'
#else
  omp_allocator = allocator_t(mesh%get_dims(VERT), SZ)
  allocator => omp_allocator
  host_allocator => omp_allocator
  if (nrank == 0) print *, 'OpenMP allocator instantiated'

  omp_backend = omp_backend_t(mesh, allocator)
  backend => omp_backend
  if (nrank == 0) print *, 'OpenMP backend instantiated'
#endif

  print *, domain_cfg%BC_x, "BC_x"
  print *, domain_cfg%BC_y, "BC_y"
  print *, domain_cfg%BC_z, "BC_Z"
  solver = solver_t(backend, mesh, host_allocator)

  if (nrank == 0) print *, 'solver instantiated'

  f => allocator%get_block(DIR_C, CELL)
  f_check => allocator%get_block(DIR_X)

  f_host => host_allocator%get_block(DIR_C)
  dims = solver%mesh%get_dims(CELL)
  do k = 1, dims(3)
    do j = 1, dims(2)
      do i = 1, dims(1)
        xloc = solver%mesh%get_coordinates(i, j, k, CELL)
        !x = xloc(1) + solver%mesh%geo%d(1)/2._dp
        !y = xloc(2) + solver%mesh%geo%d(2)/2._dp
        !z = xloc(3) + solver%mesh%geo%d(3)/2._dp
        ! x = solver%mesh%geo%midp_coords(i, 1)
        ! y = solver%mesh%geo%midp_coords(j, 2)
        ! z = solver%mesh%geo%midp_coords(k, 3)

        !print*, i, j, k, x, y, z
        !f_host%data(i, j, k) = cos(x/solver%mesh%geo%L(2)*2*pi)
        !f_host%data(i, j, k) = cos(y*63*pi) !shoulndt see mode 63 on this mesh
        f_host%data(i, j, k) = cos(2._dp*xloc(1)*pi)
      end do
    end do
  end do
print*, 'cos\n'
do i = 1, dims(1)
  print*, solver%mesh%geo%midp_coords(i, 1), (f_host%data(i, 4, 8)-f_host%data(1, 4, 8))/(-4*pi*pi)
end do
! do j = 1, dims(2)
!   print*, solver%mesh%geo%midp_coords(j, 2), (f_host%data(4, j, 8)-f_host%data(4, 1, 8))/(-4*pi*pi)
! end do
  call backend%set_field_data(f, f_host%data, DIR_C)
  call f%set_data_loc(CELL)
  call host_allocator%release_block(f_host)

  call backend%reorder(f_check, f, RDR_C2X)

  temp => backend%allocator%get_block(DIR_C)
  ! solve poisson equation with FFT based approach
  call backend%poisson_fft%solve_poisson(f, temp)

  call backend%allocator%release_block(temp)


f_host => host_allocator%get_block(DIR_C)

call backend%get_field_data(f_host%data, f)
print*, 'last\n'
do i = 1, dims(1)
  print*, solver%mesh%geo%midp_coords(i, 1), f_host%data(i, 4, 8)-f_host%data(1, 4, 8)
end do
! do j = 1, dims(2)
!   print*, solver%mesh%geo%midp_coords(j, 2), f_host%data(4, j, 8)-f_host%data(4, 1, 8)
! end do
call host_allocator%release_block(f_host)
!return

  ! reorder back to our specialist data structure from 3D Cartesian
  f_ii => backend%allocator%get_block(DIR_Z)
  call backend%reorder(f_ii, f, RDR_C2Z)
  call backend%allocator%release_block(f)

  dfdx => backend%allocator%get_block(DIR_X)
  dfdy => backend%allocator%get_block(DIR_X)
  dfdz => backend%allocator%get_block(DIR_X)
  call solver%gradient_p2v(dfdx, dfdy, dfdz, f_ii)

  call backend%allocator%release_block(f_ii)

  f_final => backend%allocator%get_block(DIR_Z)

  call solver%divergence_v2p(f_final, dfdx, dfdy, dfdz)

  call backend%allocator%release_block(dfdx)
  call backend%allocator%release_block(dfdy)
  call backend%allocator%release_block(dfdz)

  f_final_check => backend%allocator%get_block(DIR_X)
  call backend%reorder(f_final_check, f_final, RDR_Z2X)
  call backend%allocator%release_block(f_final)

  call backend%vecadd(-1._dp, f_check, 1._dp, f_final_check)

  f_host => host_allocator%get_block(DIR_C)

  call backend%get_field_data(f_host%data, f_final_check)
  call backend%allocator%release_block(f_final_check)
  dims = solver%mesh%get_dims(CELL)
  error_norm = norm2(f_host%data(1:dims(1), 1:dims(2), 1:dims(3)))/product(dims)
  print*, 'norm', error_norm

!print*, 'last\n'
!do j = 1, dims(2)
!  print*, solver%mesh%geo%midp_coords(j, 2), f_host%data(4, j, 8)
!end do

  call backend%allocator%release_block(f_check)
  call host_allocator%release_block(f_host)

  if (error_norm .gt. 1d-9) then
    error stop 'div(grad(fft())) test failed!'
  else
    if (mesh%par%is_root()) print *, "TEST PASS"
  end if

  call MPI_Finalize(ierr)

end program test_div_grad_fft
