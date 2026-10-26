! program test_poisson_nyquist_sweep
!   !! Poisson Solver Nyquist / Near-Nyquist Sweep Test
!   !!
!   !! Single swept direction, configurable by one parameter:
!   !!
!   !!   SWEEP_DIR = SWEEP_DIR_Y   (default)
!   !!     Grid: 128 x 65 x 32,  y=dirichlet/dirichlet
!   !!
!   !!   SWEEP_DIR = SWEEP_DIR_X
!   !!     Grid:  65 x 64 x 32,  x=dirichlet/dirichlet
!   !!
!   !! Both use 64 intervals in the swept direction (domain size = 1),
!   !! so results are directly comparable.
!   !!
!   !! Wavenumbers tested: n = 1, 2, 3, 32, 33, 64
!   !!
!   !! Expected failures (XFAIL):
!   !!   n=32  compact6 modified-wavenumber error at k=N/2  -> XFAIL
!   !!   n=33  one mode above Nyquist, aliases onto n=31    -> XFAIL
!   !!   n=64  exact DST/DCT basis mode, solver handles it  -> PASS (not XFAIL)
!   !!
!   !! Each 1D test:
!   !!   Check 1: Poisson solution vs analytical p = -cos(n*pi*s)/(n*pi)^2
!   !!   Check 2: div(grad(p)) recovers original RHS f  (round-trip L2)
!   !!
!   !! 3D test  --  cos(32*pi*x) * cos(32*pi*y) * cos/sin(32*pi*z)
!   !! -----------------------------------------------------------------------
!   !! Grid: 65 x 65 x 64  (64 intervals in every direction, domain = [0,1]^3)
!   !! BCs:  x=dirichlet, y=dirichlet, z=periodic
!   !!
!   !! Two RHS functions, both with analytical solution p = -f / (3*(32*pi)^2):
!   !!   TEST_COS_XYZ:  f = cos(32*pi*x) * cos(32*pi*y) * cos(32*pi*z)
!   !!   TEST_SIN_XYZ:  f = cos(32*pi*x) * cos(32*pi*y) * sin(32*pi*z)
!   !!
!   !! z is periodic so cos(32*pi*z) is the Nyquist FFT mode (k=32 out of 32).
!   !! Both variants are expected to fail for the same reason as the 1D n=32
!   !! case: compact6 modified-wavenumber error accumulates in all three
!   !! directions simultaneously.

!   use iso_fortran_env, only: stderr => error_unit
!   use mpi

!   use m_allocator,        only: allocator_t, field_t
!   use m_base_backend,     only: base_backend_t
!   use m_common,           only: dp, pi, DIR_C, DIR_X, DIR_Y, DIR_Z, &
!                                 VERT, CELL, RDR_C2Z, RDR_C2X, RDR_Z2X
!   use m_mesh,             only: mesh_t
!   use m_solver,           only: allocate_tdsops
!   use m_tdsops,           only: dirps_t
!   use m_vector_calculus,  only: vector_calculus_t

! #ifdef CUDA
!   use cudafor
!   use m_cuda_allocator,   only: cuda_allocator_t
!   use m_cuda_backend,     only: cuda_backend_t
!   use m_cuda_common,      only: SZ
! #else
!   use m_omp_backend,      only: omp_backend_t
!   use m_omp_common,       only: SZ
! #endif

!   implicit none

!   ! ==================================================================
!   ! *** CHANGE THIS ONE PARAMETER TO SWITCH DIRECTION ***
!   !
!   !   SWEEP_DIR_Y  ->  y-sweep, grid 128 x 65 x 32, y-Dirichlet
!   !   SWEEP_DIR_X  ->  x-sweep, grid  65 x 64 x 32, x-Dirichlet
!   ! ==================================================================
!   integer, parameter :: SWEEP_DIR_X = 1
!   integer, parameter :: SWEEP_DIR_Y = 2

!   integer, parameter :: SWEEP_DIR = SWEEP_DIR_Y   ! <--- change here

!   ! ==================================================================
!   ! Everything below derives from SWEEP_DIR -- do not edit
!   ! ==================================================================

!   real(dp), parameter :: ERROR_TOLERANCE = 1.0e-11_dp

!   integer, parameter :: N_SWEEP = 6
!   integer, parameter :: SWEEP_NS(N_SWEEP) = [1, 2, 3, 32, 33, 64]

!   ! n=32 and n=33 fail due to compact6 modified-wavenumber error / aliasing
!   ! n=64 passes (exact DST/DCT basis mode)
!   logical, parameter :: SWEEP_XFAIL(N_SWEEP) = &
!     [.false., .false., .false., .true., .true., .false.]

!   ! Grid: 64 intervals in swept direction, modest size in others
!   integer, parameter :: DIMS_Y(3) = [128, 65, 32]   ! y-sweep
!   integer, parameter :: DIMS_X(3) = [ 65, 64, 32]   ! x-sweep

!   ! ------------------------------------------------------------------
!   ! 3D test: cos(32*pi*x)*cos(32*pi*y)*cos/sin(32*pi*z)
!   ! Grid 65 x 65 x 64: 64 intervals in all directions, domain [0,1]^3
!   ! BCs: x=dirichlet, y=dirichlet, z=periodic
!   ! ------------------------------------------------------------------
!   integer, parameter :: N_3D = 2
!   integer, parameter :: TEST_COS_XYZ = 1   ! z-term = cos(32*pi*z)
!   integer, parameter :: TEST_SIN_XYZ = 2   ! z-term = sin(32*pi*z)
!   ! Both expected to fail: compact6 modified-wavenumber error in all 3 dirs
!   logical, parameter :: XFAIL_3D(N_3D) = [.true., .true.]

!   ! ------------------------------------------------------------------
!   ! MPI
!   ! ------------------------------------------------------------------
!   integer :: nrank, nproc, ierr

!   ! ------------------------------------------------------------------
!   ! Results
!   ! ------------------------------------------------------------------
!   logical  :: results(N_SWEEP)
!   real(dp) :: poisson_errs(N_SWEEP), divgrad_errs(N_SWEEP)

!   logical  :: results_3d(N_3D)
!   real(dp) :: poisson_errs_3d(N_3D), divgrad_errs_3d(N_3D)

!   integer :: i
!   logical :: allpass
!   character(len=1) :: dir_label

!   ! ------------------------------------------------------------------
!   ! Backend / mesh (declared here, built in run_config)
!   ! ------------------------------------------------------------------
!   class(base_backend_t), pointer :: backend
!   class(allocator_t),    pointer :: allocator
!   type(allocator_t),     pointer :: host_allocator

!   type(mesh_t),          target  :: mesh
!   type(dirps_t),         pointer :: xdirps, ydirps, zdirps
!   type(vector_calculus_t)        :: vector_calculus

! #ifdef CUDA
!   type(cuda_backend_t),   target :: cuda_backend
!   type(cuda_allocator_t), target :: cuda_allocator
! #else
!   type(omp_backend_t),    target :: omp_backend
! #endif
!   type(allocator_t),      target :: omp_allocator

!   integer :: nproc_dir(3), dims_global(3)
!   logical :: use_2decomp
!   character(len=20) :: bc_x(2), bc_y(2), bc_z(2)

!   ! ==================================================================
!   ! MPI init
!   ! ==================================================================
!   call MPI_Init(ierr)
!   call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
!   call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

!   if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'

! #ifdef CUDA
!   block
!     integer :: ndevs, devnum
!     ierr = cudaGetDeviceCount(ndevs)
!     ierr = cudaSetDevice(mod(nrank, ndevs))
!     ierr = cudaGetDevice(devnum)
!   end block
! #endif

!   ! ==================================================================
!   ! Derive grid and BCs from SWEEP_DIR
!   ! ==================================================================
!   if (SWEEP_DIR == SWEEP_DIR_Y) then
!     dims_global = DIMS_Y
!     bc_x = ['periodic ', 'periodic ']
!     bc_y = ['dirichlet', 'dirichlet']
!     bc_z = ['periodic ', 'periodic ']
!     dir_label = 'y'
!   else
!     dims_global = DIMS_X
!     bc_x = ['dirichlet', 'dirichlet']
!     bc_y = ['periodic ', 'periodic ']
!     bc_z = ['periodic ', 'periodic ']
!     dir_label = 'x'
!   end if

!   ! ==================================================================
!   ! Header
!   ! ==================================================================
!   if (nrank == 0) then
!     write(stderr,'(A)') ''
!     write(stderr,'(A)') &
!       '  ======================================================='
!     write(stderr,'(2X,A,A,A,I0,A,I0,A,I0)') &
!       'Nyquist sweep  dir=', dir_label, &
!       '  grid ', dims_global(1), ' x ', dims_global(2), &
!       ' x ', dims_global(3)
!     write(stderr,'(A)') &
!       '  ======================================================='
!     write(stderr,'(4X,A,A,A,A,A,A,A,A,A,A)') &
!       'BC: x=[', trim(bc_x(1)), ',', trim(bc_x(2)), &
!       ']  y=[', trim(bc_y(1)), ',', trim(bc_y(2)), &
!       ']  z=[', trim(bc_z(1)), ','//trim(bc_z(2))//']'
!     write(stderr,'(4X,A,ES8.1)') 'Tolerance: ', ERROR_TOLERANCE
!     write(stderr,'(A)') ''
!   end if

!   ! ==================================================================
!   ! Build mesh
!   ! ==================================================================
!   nproc_dir = [1, 1, nproc]

! #ifdef CUDA
!   use_2decomp = .false.
! #else
!   use_2decomp = .true.
! #endif

!   mesh = mesh_t(dims_global, nproc_dir, [1.0_dp, 1.0_dp, 1.0_dp], &
!                 bc_x, bc_y, bc_z, use_2decomp=use_2decomp)

!   ! ==================================================================
!   ! Build backend & allocators
!   ! ==================================================================
! #ifdef CUDA
!   cuda_allocator = cuda_allocator_t(mesh%get_dims(VERT), SZ)
!   allocator      => cuda_allocator
!   omp_allocator  = allocator_t(mesh%get_dims(VERT), SZ)
!   host_allocator => omp_allocator
!   cuda_backend   = cuda_backend_t(mesh, allocator)
!   backend        => cuda_backend
! #else
!   omp_allocator  = allocator_t(mesh%get_dims(VERT), SZ)
!   allocator      => omp_allocator
!   host_allocator => omp_allocator
!   omp_backend    = omp_backend_t(mesh, allocator)
!   backend        => omp_backend
! #endif

!   ! ==================================================================
!   ! Build tdsops and Poisson FFT
!   ! ==================================================================
!   allocate(xdirps, ydirps, zdirps)
!   xdirps%dir = DIR_X
!   ydirps%dir = DIR_Y
!   zdirps%dir = DIR_Z
!   call allocate_tdsops(xdirps, backend, mesh, 'compact6', 'compact6', &
!                        'classic', 'compact6')
!   call allocate_tdsops(ydirps, backend, mesh, 'compact6', 'compact6', &
!                        'classic', 'compact6')
!   call allocate_tdsops(zdirps, backend, mesh, 'compact6', 'compact6', &
!                        'classic', 'compact6')

!   vector_calculus = vector_calculus_t(backend)
!   call backend%init_poisson_fft(mesh, xdirps, ydirps, zdirps)

!   ! ==================================================================
!   ! Run sweep
!   ! ==================================================================
!   do i = 1, N_SWEEP
!     call run_cos_dir_test(backend, host_allocator, mesh, &
!                           xdirps, ydirps, zdirps, vector_calculus, &
!                           SWEEP_DIR, SWEEP_NS(i), &
!                           results(i), poisson_errs(i), divgrad_errs(i))
!   end do

!   ! ==================================================================
!   ! Summary table
!   ! ==================================================================
!   if (nrank == 0) then
!     write(stderr,'(A)') ''
!     write(stderr,'(A)') &
!       '  -------------------------------------------------------'
!     write(stderr,'(A)') '  SUMMARY'
!     write(stderr,'(A)') &
!       '  -------------------------------------------------------'
!     write(stderr,'(2X,A4,2X,A18,A14,A14,2X,A6,2X,A8,2X,A)') &
!       '  n ', 'Function          ', &
!       '  Poisson L2  ', '  DivGrad L2  ', &
!       'Result', 'Expected', 'Verdict'
!     write(stderr,'(A)') ''

!     do i = 1, N_SWEEP
!       call print_row(SWEEP_NS(i), dir_label, &
!                      poisson_errs(i), divgrad_errs(i), &
!                      results(i), SWEEP_XFAIL(i))
!     end do

!     write(stderr,'(A)') ''
!     write(stderr,'(A)') '  Notes:'
!     write(stderr,'(4X,A)') &
!       'n=32  compact6 modified-wavenumber error at k=N/2 (XFAIL)'
!     write(stderr,'(4X,A)') &
!       'n=33  above Nyquist, aliases onto n=31             (XFAIL)'
!     write(stderr,'(4X,A)') &
!       'n=64  exact DST/DCT basis mode, passes cleanly     (PASS)'
!     write(stderr,'(4X,A)') &
!       'Verdict OK = result matches expectation'
!     write(stderr,'(4X,A)') &
!       'Verdict !! = unexpected outcome'
!     write(stderr,'(A)') ''
!   end if

!   ! ==================================================================
!   ! 3D test: cos(32*pi*x)*cos(32*pi*y)*cos/sin(32*pi*z)
!   ! ==================================================================
!   if (nrank == 0) then
!     write(stderr,'(A)') ''
!     write(stderr,'(A)') &
!       '  ======================================================='
!     write(stderr,'(A)') &
!       '  3D test  |  65 x 65 x 64  |  x,y=Dirichlet  z=periodic'
!     write(stderr,'(A)') &
!       '  ======================================================='
!     write(stderr,'(4X,A)') &
!       'f = cos(32*pi*x) * cos(32*pi*y) * cos/sin(32*pi*z)'
!     write(stderr,'(4X,A)') &
!       'Analytical: p = -f / (3*(32*pi)^2)'
!     write(stderr,'(4X,A,ES8.1)') 'Tolerance: ', ERROR_TOLERANCE
!     write(stderr,'(A)') ''
!   end if

!   call run_3d_test(nproc, nrank, N_3D, XFAIL_3D, &
!                    results_3d, poisson_errs_3d, divgrad_errs_3d)

!   if (nrank == 0) then
!     write(stderr,'(A)') ''
!     write(stderr,'(A)') '  3D SUMMARY'
!     write(stderr,'(A)') &
!       '  -------------------------------------------------------'
!     write(stderr,'(2X,A24,A14,A14,2X,A6,2X,A8,2X,A)') &
!       'Function                ', &
!       '  Poisson L2  ', '  DivGrad L2  ', &
!       'Result', 'Expected', 'Verdict'
!     write(stderr,'(A)') ''
!     call print_row_3d('cos(32px)*cos(32py)*cos(32pz)', &
!                       poisson_errs_3d(1), divgrad_errs_3d(1), &
!                       results_3d(1), XFAIL_3D(1))
!     call print_row_3d('cos(32px)*cos(32py)*sin(32pz)', &
!                       poisson_errs_3d(2), divgrad_errs_3d(2), &
!                       results_3d(2), XFAIL_3D(2))
!     write(stderr,'(A)') ''
!     write(stderr,'(A)') '  Notes:'
!     write(stderr,'(4X,A)') &
!       'z=periodic: cos/sin(32*pi*z) is the Nyquist FFT mode'
!     write(stderr,'(4X,A)') &
!       'All 3 directions at modified-wavenumber error regime simultaneously'
!     write(stderr,'(4X,A)') &
!       'cos vs sin in z isolates whether z-parity affects the error'
!     write(stderr,'(A)') ''
!   end if

!   ! ==================================================================
!   ! Final verdict
!   ! ==================================================================
!   allpass = .true.
!   do i = 1, N_SWEEP
!     if (results(i) .neqv. (.not. SWEEP_XFAIL(i))) allpass = .false.
!   end do
!   do i = 1, N_3D
!     if (results_3d(i) .neqv. (.not. XFAIL_3D(i))) allpass = .false.
!   end do

!   if (allpass) then
!     if (nrank == 0) &
!       write(stderr,'(A)') 'ALL TESTS PASSED (accounting for XFAIL).'
!   else
!     error stop 'SOME TESTS FAILED UNEXPECTEDLY.'
!   end if

!   call MPI_Finalize(ierr)

! contains

!   ! ================================================================
!   ! Run one cos(n*pi*s) Poisson test, s = x or y per sweep_dir
!   ! ================================================================
!   subroutine run_cos_dir_test(backend, host_allocator, mesh, &
!                                xdirps, ydirps, zdirps, vector_calculus, &
!                                sweep_dir, n_wave, &
!                                test_passed, poisson_err_out, divgrad_err_out)
!     class(base_backend_t), pointer, intent(in)  :: backend
!     type(allocator_t),     pointer, intent(in)  :: host_allocator
!     type(mesh_t),                   intent(in)  :: mesh
!     type(dirps_t),         pointer, intent(in)  :: xdirps, ydirps, zdirps
!     type(vector_calculus_t),        intent(in)  :: vector_calculus
!     integer,                        intent(in)  :: sweep_dir, n_wave
!     logical,                        intent(out) :: test_passed
!     real(dp),                       intent(out) :: poisson_err_out
!     real(dp),                       intent(out) :: divgrad_err_out

!     class(field_t), pointer :: f_device, f_reference, f_result
!     class(field_t), pointer :: host_field, host_analytical, temp
!     class(field_t), pointer :: dpdx, dpdy, dpdz, gradient_input

!     integer  :: dims(3), ii, jj, kk
!     real(dp) :: coords(3), n_pi, n_pi_sq, coord_val
!     real(dp) :: poisson_error_norm, div_grad_error_norm
!     logical  :: poisson_passed, div_grad_passed

!     dims    = mesh%get_dims(CELL)
!     n_pi    = real(n_wave, dp) * pi
!     n_pi_sq = n_pi * n_pi

!     if (mesh%par%is_root()) then
!       write(stderr,'(4X,A,I0,A,A,A)') &
!         'cos(', n_wave, '*pi*', &
!         merge('x', 'y', sweep_dir == SWEEP_DIR_X), ') ---'
!     end if

!     ! ------------------------------------------------------------------
!     ! Fill RHS: f = cos(n*pi*coord[sweep_dir])
!     ! ------------------------------------------------------------------
!     f_device    => backend%allocator%get_block(DIR_C, CELL)
!     f_reference => backend%allocator%get_block(DIR_X)
!     host_field  => host_allocator%get_block(DIR_C)

!     do kk = 1, dims(3)
!       do jj = 1, dims(2)
!         do ii = 1, dims(1)
!           coords    = mesh%get_coordinates(ii, jj, kk, CELL)
!           coord_val = merge(coords(1), coords(2), sweep_dir == SWEEP_DIR_X)
!           host_field%data(ii, jj, kk) = cos(n_pi * coord_val)
!         end do
!       end do
!     end do

!     call backend%set_field_data(f_device, host_field%data, DIR_C)
!     call f_device%set_data_loc(CELL)
!     call host_allocator%release_block(host_field)

!     ! Keep original RHS in DIR_X layout for the div-grad check
!     call backend%reorder(f_reference, f_device, RDR_C2X)

!     ! ------------------------------------------------------------------
!     ! Solve Poisson: nabla^2 p = f
!     ! ------------------------------------------------------------------
!     temp => backend%allocator%get_block(DIR_C)
!     call backend%poisson_fft%solve_poisson(f_device, temp)
!     call backend%allocator%release_block(temp)

!     ! ------------------------------------------------------------------
!     ! Check 1: p_numerical vs p_analytical = -cos(n*pi*s) / (n*pi)^2
!     ! ------------------------------------------------------------------
!     host_field => host_allocator%get_block(DIR_C)
!     call backend%get_field_data(host_field%data, f_device)

!     ! Remove arbitrary additive constant
!     host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) = &
!       host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) &
!       - host_field%data(1, 1, 1)

!     host_analytical => host_allocator%get_block(DIR_C)
!     do kk = 1, dims(3)
!       do jj = 1, dims(2)
!         do ii = 1, dims(1)
!           coords    = mesh%get_coordinates(ii, jj, kk, CELL)
!           coord_val = merge(coords(1), coords(2), sweep_dir == SWEEP_DIR_X)
!           host_analytical%data(ii, jj, kk) = -cos(n_pi * coord_val) / n_pi_sq
!         end do
!       end do
!     end do

!     host_analytical%data(1:dims(1), 1:dims(2), 1:dims(3)) = &
!       host_analytical%data(1:dims(1), 1:dims(2), 1:dims(3)) &
!       - host_analytical%data(1, 1, 1)

!     host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) = &
!       host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) &
!       - host_analytical%data(1:dims(1), 1:dims(2), 1:dims(3))

!     poisson_error_norm = norm2( &
!       host_field%data(1:dims(1), 1:dims(2), 1:dims(3))) &
!       / product(dims)

!     call host_allocator%release_block(host_analytical)
!     call host_allocator%release_block(host_field)

!     poisson_passed = (poisson_error_norm <= ERROR_TOLERANCE)

!     ! ------------------------------------------------------------------
!     ! Check 2: div(grad(p)) recovers original RHS f
!     ! ------------------------------------------------------------------
!     gradient_input => backend%allocator%get_block(DIR_Z)
!     call backend%reorder(gradient_input, f_device, RDR_C2Z)
!     call backend%allocator%release_block(f_device)

!     dpdx => backend%allocator%get_block(DIR_X)
!     dpdy => backend%allocator%get_block(DIR_X)
!     dpdz => backend%allocator%get_block(DIR_X)

!     call vector_calculus%gradient_c2v( &
!       dpdx, dpdy, dpdz, gradient_input, &
!       xdirps%stagder_p2v, xdirps%interpl_p2v, &
!       ydirps%stagder_p2v, ydirps%interpl_p2v, &
!       zdirps%stagder_p2v, zdirps%interpl_p2v)

!     call backend%allocator%release_block(gradient_input)

!     f_result => backend%allocator%get_block(DIR_Z)

!     call vector_calculus%divergence_v2c( &
!       f_result, dpdx, dpdy, dpdz, &
!       xdirps%stagder_v2p, xdirps%interpl_v2p, &
!       ydirps%stagder_v2p, ydirps%interpl_v2p, &
!       zdirps%stagder_v2p, zdirps%interpl_v2p)

!     call backend%allocator%release_block(dpdx)
!     call backend%allocator%release_block(dpdy)
!     call backend%allocator%release_block(dpdz)

!     f_device => backend%allocator%get_block(DIR_X)
!     call backend%reorder(f_device, f_result, RDR_Z2X)
!     call backend%allocator%release_block(f_result)

!     call backend%vecadd(-1.0_dp, f_reference, 1.0_dp, f_device)

!     host_field => host_allocator%get_block(DIR_C)
!     call backend%get_field_data(host_field%data, f_device)

!     div_grad_error_norm = norm2( &
!       host_field%data(1:dims(1), 1:dims(2), 1:dims(3))) &
!       / product(dims)

!     call backend%allocator%release_block(f_device)
!     call backend%allocator%release_block(f_reference)
!     call host_allocator%release_block(host_field)

!     div_grad_passed = (div_grad_error_norm <= ERROR_TOLERANCE)

!     if (mesh%par%is_root()) then
!       write(stderr,'(6X,A,ES12.4,2X,A)') &
!         'Poisson L2: ', poisson_error_norm, &
!         merge('PASS', 'FAIL', poisson_passed)
!       write(stderr,'(6X,A,ES12.4,2X,A)') &
!         'DivGrad L2: ', div_grad_error_norm, &
!         merge('PASS', 'FAIL', div_grad_passed)
!     end if

!     test_passed     = poisson_passed .and. div_grad_passed
!     poisson_err_out = poisson_error_norm
!     divgrad_err_out = div_grad_error_norm

!   end subroutine run_cos_dir_test

!   ! ================================================================
!   ! Print one summary row
!   ! ================================================================
!   subroutine print_row(n, dir, p_err, dg_err, passed, xfail)
!     integer,          intent(in) :: n
!     character(len=1), intent(in) :: dir
!     real(dp),         intent(in) :: p_err, dg_err
!     logical,          intent(in) :: passed, xfail

!     character(len=4)  :: result_str, expected_str
!     character(len=2)  :: verdict_str
!     character(len=18) :: label

!     result_str   = merge('PASS', 'FAIL', passed)
!     expected_str = merge('FAIL', 'PASS', xfail)
!     verdict_str  = merge('OK', '!!', passed .neqv. xfail)
!     write(label, '(A,I0,A,A,A)') 'cos(', n, '*pi*', dir, ')  '

!     write(stderr,'(2X,I4,2X,A18,ES14.4,ES14.4,2X,A4,4X,A4,4X,A)') &
!       n, label, p_err, dg_err, result_str, expected_str, trim(verdict_str)
!   end subroutine print_row

!   ! ================================================================
!   ! Build mesh+backend for 65x65x64, run cos and sin 3D tests
!   ! ================================================================
!   subroutine run_3d_test(nproc, nrank, n_cases, xfail, &
!                           results, poisson_errs, divgrad_errs)
!     integer, intent(in)  :: nproc, nrank, n_cases
!     logical, intent(in)  :: xfail(n_cases)
!     logical, intent(out) :: results(n_cases)
!     real(dp), intent(out) :: poisson_errs(n_cases), divgrad_errs(n_cases)

!     class(base_backend_t), pointer :: backend
!     class(allocator_t),    pointer :: allocator
!     type(allocator_t),     pointer :: host_allocator

!     type(mesh_t),          target  :: mesh
!     type(dirps_t),         pointer :: xdirps, ydirps, zdirps
!     type(vector_calculus_t)        :: vector_calculus

! #ifdef CUDA
!     type(cuda_backend_t),   target :: cuda_backend
!     type(cuda_allocator_t), target :: cuda_allocator
! #else
!     type(omp_backend_t),    target :: omp_backend
! #endif
!     type(allocator_t),      target :: omp_allocator

!     integer :: nproc_dir(3), i
!     logical :: use_2decomp

!     ! 65 x 65 x 64: 64 intervals in every direction
!     integer, parameter :: DIMS_3D(3) = [513, 513, 512]

!     nproc_dir = [1, 1, nproc]

! #ifdef CUDA
!     use_2decomp = .false.
! #else
!     use_2decomp = .true.
! #endif

!     mesh = mesh_t(DIMS_3D, nproc_dir, [1.0_dp, 1.0_dp, 1.0_dp], &
!                   ['neumann', 'neumann'], &   ! BC_x
!                   ['neumann', 'neumann'], &   ! BC_y
!                   ['periodic ', 'periodic '], &   ! BC_z
!                   use_2decomp=use_2decomp)

! #ifdef CUDA
!     cuda_allocator = cuda_allocator_t(mesh%get_dims(VERT), SZ)
!     allocator      => cuda_allocator
!     omp_allocator  = allocator_t(mesh%get_dims(VERT), SZ)
!     host_allocator => omp_allocator
!     cuda_backend   = cuda_backend_t(mesh, allocator)
!     backend        => cuda_backend
! #else
!     omp_allocator  = allocator_t(mesh%get_dims(VERT), SZ)
!     allocator      => omp_allocator
!     host_allocator => omp_allocator
!     omp_backend    = omp_backend_t(mesh, allocator)
!     backend        => omp_backend
! #endif

!     allocate(xdirps, ydirps, zdirps)
!     xdirps%dir = DIR_X
!     ydirps%dir = DIR_Y
!     zdirps%dir = DIR_Z
!     call allocate_tdsops(xdirps, backend, mesh, 'compact6', 'compact6', &
!                          'classic', 'compact6')
!     call allocate_tdsops(ydirps, backend, mesh, 'compact6', 'compact6', &
!                          'classic', 'compact6')
!     call allocate_tdsops(zdirps, backend, mesh, 'compact6', 'compact6', &
!                          'classic', 'compact6')

!     vector_calculus = vector_calculus_t(backend)
!     call backend%init_poisson_fft(mesh, xdirps, ydirps, zdirps)

!     do i = 1, n_cases
!       call run_cos3d_test(backend, host_allocator, mesh, &
!                           xdirps, ydirps, zdirps, vector_calculus, &
!                           i, results(i), poisson_errs(i), divgrad_errs(i))
!     end do

!   end subroutine run_3d_test

!   ! ================================================================
!   ! Run one 3D Poisson test (Checks 1 & 2)
!   !
!   ! test_id = TEST_COS_XYZ:  f = cos(32*pi*x)*cos(32*pi*y)*cos(32*pi*z)
!   ! test_id = TEST_SIN_XYZ:  f = cos(32*pi*x)*cos(32*pi*y)*sin(32*pi*z)
!   !
!   ! Analytical solution: p = -f / (3*(32*pi)^2)
!   ! ================================================================
!   subroutine run_cos3d_test(backend, host_allocator, mesh, &
!                              xdirps, ydirps, zdirps, vector_calculus, &
!                              test_id, test_passed, &
!                              poisson_err_out, divgrad_err_out)
!     class(base_backend_t), pointer, intent(in)  :: backend
!     type(allocator_t),     pointer, intent(in)  :: host_allocator
!     type(mesh_t),                   intent(in)  :: mesh
!     type(dirps_t),         pointer, intent(in)  :: xdirps, ydirps, zdirps
!     type(vector_calculus_t),        intent(in)  :: vector_calculus
!     integer,                        intent(in)  :: test_id
!     logical,                        intent(out) :: test_passed
!     real(dp),                       intent(out) :: poisson_err_out
!     real(dp),                       intent(out) :: divgrad_err_out

!     class(field_t), pointer :: f_device, f_reference, f_result
!     class(field_t), pointer :: host_field, host_analytical, temp
!     class(field_t), pointer :: dpdx, dpdy, dpdz, gradient_input

!     integer  :: dims(3), ii, jj, kk
!     real(dp) :: coords(3), n_pi, n_pi_sq_3, fval
!     real(dp) :: poisson_error_norm, div_grad_error_norm
!     logical  :: poisson_passed, div_grad_passed
!     logical  :: use_sin_z

!     integer, parameter :: N_WAVE = 32

!     n_pi      = real(N_WAVE, dp) * pi
!     n_pi_sq_3 = 3.0_dp * n_pi * n_pi   ! eigenvalue: 3*(n*pi)^2
!     dims      = mesh%get_dims(CELL)
!     use_sin_z = (test_id == TEST_SIN_XYZ)

!     if (mesh%par%is_root()) then
!       if (use_sin_z) then
!         write(stderr,'(4X,A)') &
!           'cos(32*pi*x)*cos(32*pi*y)*sin(32*pi*z) ---'
!       else
!         write(stderr,'(4X,A)') &
!           'cos(32*pi*x)*cos(32*pi*y)*cos(32*pi*z) ---'
!       end if
!     end if

!     ! ------------------------------------------------------------------
!     ! Fill RHS
!     ! ------------------------------------------------------------------
!     f_device    => backend%allocator%get_block(DIR_C, CELL)
!     f_reference => backend%allocator%get_block(DIR_X)
!     host_field  => host_allocator%get_block(DIR_C)

!     do kk = 1, dims(3)
!       do jj = 1, dims(2)
!         do ii = 1, dims(1)
!           coords = mesh%get_coordinates(ii, jj, kk, CELL)
!           if (use_sin_z) then
!             fval = cos(n_pi*coords(1)) * cos(n_pi*coords(2)) &
!                  * sin(n_pi*coords(3))
!           else
!             fval = cos(n_pi*coords(1)) * cos(n_pi*coords(2)) &
!                  * cos(n_pi*coords(3))
!           end if
!           host_field%data(ii, jj, kk) = fval
!         end do
!       end do
!     end do

!     call backend%set_field_data(f_device, host_field%data, DIR_C)
!     call f_device%set_data_loc(CELL)
!     call host_allocator%release_block(host_field)

!     call backend%reorder(f_reference, f_device, RDR_C2X)

!     ! ------------------------------------------------------------------
!     ! Solve Poisson: nabla^2 p = f
!     ! ------------------------------------------------------------------
!     temp => backend%allocator%get_block(DIR_C)
!     call backend%poisson_fft%solve_poisson(f_device, temp)
!     call backend%allocator%release_block(temp)

!     ! ------------------------------------------------------------------
!     ! Check 1: p_numerical vs p_analytical = -f / (3*(n*pi)^2)
!     ! ------------------------------------------------------------------
!     host_field => host_allocator%get_block(DIR_C)
!     call backend%get_field_data(host_field%data, f_device)

!     ! Remove arbitrary additive constant
!     host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) = &
!       host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) &
!       - host_field%data(1, 1, 1)

!     host_analytical => host_allocator%get_block(DIR_C)
!     do kk = 1, dims(3)
!       do jj = 1, dims(2)
!         do ii = 1, dims(1)
!           coords = mesh%get_coordinates(ii, jj, kk, CELL)
!           if (use_sin_z) then
!             fval = cos(n_pi*coords(1)) * cos(n_pi*coords(2)) &
!                  * sin(n_pi*coords(3))
!           else
!             fval = cos(n_pi*coords(1)) * cos(n_pi*coords(2)) &
!                  * cos(n_pi*coords(3))
!           end if
!           host_analytical%data(ii, jj, kk) = -fval / n_pi_sq_3
!         end do
!       end do
!     end do

!     host_analytical%data(1:dims(1), 1:dims(2), 1:dims(3)) = &
!       host_analytical%data(1:dims(1), 1:dims(2), 1:dims(3)) &
!       - host_analytical%data(1, 1, 1)

!     host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) = &
!       host_field%data(1:dims(1), 1:dims(2), 1:dims(3)) &
!       - host_analytical%data(1:dims(1), 1:dims(2), 1:dims(3))

!     poisson_error_norm = norm2( &
!       host_field%data(1:dims(1), 1:dims(2), 1:dims(3))) &
!       / product(dims)

!     call host_allocator%release_block(host_analytical)
!     call host_allocator%release_block(host_field)

!     poisson_passed = (poisson_error_norm <= ERROR_TOLERANCE)

!     ! ------------------------------------------------------------------
!     ! Check 2: div(grad(p)) recovers original RHS f
!     ! ------------------------------------------------------------------
!     gradient_input => backend%allocator%get_block(DIR_Z)
!     call backend%reorder(gradient_input, f_device, RDR_C2Z)
!     call backend%allocator%release_block(f_device)

!     dpdx => backend%allocator%get_block(DIR_X)
!     dpdy => backend%allocator%get_block(DIR_X)
!     dpdz => backend%allocator%get_block(DIR_X)

!     call vector_calculus%gradient_c2v( &
!       dpdx, dpdy, dpdz, gradient_input, &
!       xdirps%stagder_p2v, xdirps%interpl_p2v, &
!       ydirps%stagder_p2v, ydirps%interpl_p2v, &
!       zdirps%stagder_p2v, zdirps%interpl_p2v)

!     call backend%allocator%release_block(gradient_input)

!     f_result => backend%allocator%get_block(DIR_Z)

!     call vector_calculus%divergence_v2c( &
!       f_result, dpdx, dpdy, dpdz, &
!       xdirps%stagder_v2p, xdirps%interpl_v2p, &
!       ydirps%stagder_v2p, ydirps%interpl_v2p, &
!       zdirps%stagder_v2p, zdirps%interpl_v2p)

!     call backend%allocator%release_block(dpdx)
!     call backend%allocator%release_block(dpdy)
!     call backend%allocator%release_block(dpdz)

!     f_device => backend%allocator%get_block(DIR_X)
!     call backend%reorder(f_device, f_result, RDR_Z2X)
!     call backend%allocator%release_block(f_result)

!     call backend%vecadd(-1.0_dp, f_reference, 1.0_dp, f_device)

!     host_field => host_allocator%get_block(DIR_C)
!     call backend%get_field_data(host_field%data, f_device)

!     div_grad_error_norm = norm2( &
!       host_field%data(1:dims(1), 1:dims(2), 1:dims(3))) &
!       / product(dims)

!     call backend%allocator%release_block(f_device)
!     call backend%allocator%release_block(f_reference)
!     call host_allocator%release_block(host_field)

!     div_grad_passed = (div_grad_error_norm <= ERROR_TOLERANCE)

!     if (mesh%par%is_root()) then
!       write(stderr,'(6X,A,ES12.4,2X,A)') &
!         'Poisson L2: ', poisson_error_norm, &
!         merge('PASS', 'FAIL', poisson_passed)
!       write(stderr,'(6X,A,ES12.4,2X,A)') &
!         'DivGrad L2: ', div_grad_error_norm, &
!         merge('PASS', 'FAIL', div_grad_passed)
!     end if

!     test_passed     = poisson_passed .and. div_grad_passed
!     poisson_err_out = poisson_error_norm
!     divgrad_err_out = div_grad_error_norm

!   end subroutine run_cos3d_test

!   ! ================================================================
!   ! Print one 3D summary row
!   ! ================================================================
!   subroutine print_row_3d(label, p_err, dg_err, passed, xfail)
!     character(len=*), intent(in) :: label
!     real(dp),         intent(in) :: p_err, dg_err
!     logical,          intent(in) :: passed, xfail

!     character(len=4) :: result_str, expected_str
!     character(len=2) :: verdict_str

!     result_str   = merge('PASS', 'FAIL', passed)
!     expected_str = merge('FAIL', 'PASS', xfail)
!     verdict_str  = merge('OK', '!!', passed .neqv. xfail)

!     write(stderr,'(2X,A29,ES14.4,ES14.4,2X,A4,4X,A4,4X,A)') &
!       label, p_err, dg_err, result_str, expected_str, trim(verdict_str)
!   end subroutine print_row_3d


! end program test_poisson_nyquist_sweep

! program test_1d_stability
!   !! 1D convection-diffusion stability test for compact6 + Dirichlet BC.
!   !!
!   !! Solve:   du/dt + U * du/dx = nu * d2u/dx2
!   !!          u(0) = 1, u(L) = 1
!   !!          u(x, 0) = 1 + small noise
!   !!
!   !! Integrate with AB3 and track max|u - 1| over time.
!   !!
!   !! Purpose: isolate whether the compact6 Dirichlet boundary stencil is
!   !! stable under advection-diffusion time integration with noisy IC, in
!   !! the absence of pressure coupling and 3D effects.
!   !!
!   !! Expected outcomes:
!   !!   max|u - 1| decays monotonically  =>  closure is stable
!   !!   max|u - 1| grows exponentially    =>  closure + AB3 unstable on its own
!   !!   max|u - 1| bounded but non-zero  =>  closure has steady-state error

!   use iso_fortran_env
!   use cudafor

!   use m_common,          only: dp, pi, BC_DIRICHLET
!   use m_cuda_common,     only: SZ
!   use m_cuda_tdsops,     only: cuda_tdsops_t, cuda_tdsops_init
!   use m_cuda_exec_thom,  only: exec_thom_tds_compact

!   implicit none

!   ! -------- problem parameters (match 3D cylinder run as closely as poss) --
!   integer,  parameter :: N       = 65              ! grid points
!   real(dp), parameter :: L       = 20.0_dp          ! domain length
!   real(dp), parameter :: U       = 1.0_dp           ! advection velocity
!   real(dp), parameter :: Re      = 300.0_dp         ! Reynolds number
!   real(dp), parameter :: nu      = U * L / Re       ! kinematic viscosity
!   real(dp), parameter :: dt      = 0.0075_dp        ! time step
!   integer,  parameter :: n_iters = 2000             ! total iterations
!   integer,  parameter :: n_print = 50               ! print interval

!   real(dp), parameter :: noise_amp = 0.25_dp        ! IC noise amplitude
!   integer,  parameter :: n_block   = 64             ! CUDA blocks (dummy)

!   ! AB3 coefficients
!   real(dp), parameter :: ab3_a =  23.0_dp / 12.0_dp
!   real(dp), parameter :: ab3_b = -16.0_dp / 12.0_dp
!   real(dp), parameter :: ab3_c =   5.0_dp / 12.0_dp

!   ! -------- fields --------
!   real(dp) :: dx
!   real(dp), allocatable :: u_host(:, :, :), du_host(:, :, :)
!   real(dp), allocatable :: rhs_host(:, :, :)
!   real(dp), allocatable :: rhs_nm1(:, :, :), rhs_nm2(:, :, :)
!   real(dp), allocatable :: u_analytic(:)

!   real(dp), device, allocatable :: u_dev(:, :, :)
!   real(dp), device, allocatable :: du_dev(:, :, :)

!   type(cuda_tdsops_t) :: op1, op2
!   type(dim3) :: blocks, threads

!   integer  :: i, j, k, iter, ndevs, devnum, ierr
!   real(dp) :: umax, umin, umax_err, dt_cfl, dt_diff
!   real(dp) :: r

!   ! -------- CUDA setup --------
!   ierr = cudaGetDeviceCount(ndevs)
!   ierr = cudaSetDevice(0)
!   ierr = cudaGetDevice(devnum)

!   ! -------- grid --------
!   dx = L / real(N - 1, dp)

!   ! -------- stability estimates --------
!   dt_cfl  = dx / U
!   dt_diff = 0.5_dp * dx * dx / nu

!   write (*, '(A)') '============================================='
!   write (*, '(A)') ' 1D convection-diffusion stability test'
!   write (*, '(A)') '============================================='
!   write (*, '(A,I0)')        '  N          = ', N
!   write (*, '(A,ES12.4)')    '  L          = ', L
!   write (*, '(A,ES12.4)')    '  dx         = ', dx
!   write (*, '(A,ES12.4)')    '  U          = ', U
!   write (*, '(A,ES12.4)')    '  nu         = ', nu
!   write (*, '(A,ES12.4)')    '  Re         = ', Re
!   write (*, '(A,ES12.4)')    '  dt         = ', dt
!   write (*, '(A,ES12.4,A)')  '  CFL (U)    = ', dt / dt_cfl,  ' (< 1 stable)'
!   write (*, '(A,ES12.4,A)')  '  diff lim   = ', dt / dt_diff, ' (< 0.5 stable)'
!   write (*, '(A,I0)')        '  n_iters    = ', n_iters
!   write (*, '(A,ES12.4)')    '  noise amp  = ', noise_amp
!   write (*, '(A)')           '  BC         = Dirichlet(1) / Dirichlet(1)'
!   write (*, '(A)') '---------------------------------------------'

!   ! -------- allocate --------
!   allocate (u_host  (SZ, N, n_block))
!   allocate (du_host (SZ, N, n_block))
!   allocate (rhs_host(SZ, N, n_block))
!   allocate (rhs_nm1 (SZ, N, n_block))
!   allocate (rhs_nm2 (SZ, N, n_block))
!   allocate (u_analytic(N))
!   allocate (u_dev (SZ, N, n_block))
!   allocate (du_dev(SZ, N, n_block))

!   blocks  = dim3(n_block, 1, 1)
!   threads = dim3(SZ, 1, 1)

!   ! -------- initial condition: u = 1 + noise (reproducible) --------
!   call random_seed()
!   do j = 1, N
!     call random_number(r)                ! r in [0, 1)
!     u_host(:, j, :) = 1.0_dp + noise_amp * (2.0_dp * r - 1.0_dp)
!   end do
!   ! Enforce Dirichlet IC at boundaries
!   u_host(:, 1, :) = 1.0_dp
!   u_host(:, N, :) = 1.0_dp
!   u_analytic      = 1.0_dp

!   rhs_nm1 = 0.0_dp
!   rhs_nm2 = 0.0_dp

!   ! -------- build compact6 Dirichlet operators --------
!   op1 = cuda_tdsops_init(N, dx, operation='first-deriv',  scheme='compact6', &
!                          bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET)
!   op2 = cuda_tdsops_init(N, dx, operation='second-deriv', scheme='compact6', &
!                          bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET)

!   ! -------- initial monitor --------
!   call report(0, 0.0_dp, u_host, u_analytic, N, n_block)

!   ! -------- time loop --------
!   do iter = 1, n_iters

!     ! ---- copy to device and compute du/dx ----
!     u_dev  = u_host
!     du_dev = 0.0_dp
!     call exec_thom_tds_compact(du_dev, u_dev, op1, blocks, threads)
!     du_host = du_dev
!     ! Start RHS with -U * du/dx
!     rhs_host = -U * du_host

!     ! ---- compute d2u/dx2 ----
!     du_dev = 0.0_dp
!     call exec_thom_tds_compact(du_dev, u_dev, op2, blocks, threads)
!     du_host = du_dev
!     ! Add nu * d2u/dx2 to RHS
!     rhs_host = rhs_host + nu * du_host

!     ! ---- AB3 step (bootstrap: use Euler on iter=1, AB2 on iter=2) ----
!     if (iter == 1) then
!       u_host = u_host + dt * rhs_host
!     else if (iter == 2) then
!       u_host = u_host + dt * (1.5_dp * rhs_host - 0.5_dp * rhs_nm1)
!     else
!       u_host = u_host + dt * ( ab3_a * rhs_host &
!                              + ab3_b * rhs_nm1  &
!                              + ab3_c * rhs_nm2)
!     end if

!     ! ---- enforce Dirichlet BCs at i=1 and i=N ----
!     u_host(:, 1, :) = 1.0_dp
!     u_host(:, N, :) = 1.0_dp

!     ! ---- shift RHS history for AB3 ----
!     rhs_nm2 = rhs_nm1
!     rhs_nm1 = rhs_host

!     ! ---- monitor ----
!     if (mod(iter, n_print) == 0 .or. iter == 1 .or. iter == 2) then
!       call report(iter, iter * dt, u_host, u_analytic, N, n_block)
!     end if

!     ! ---- bailout on blow-up ----
!     umax = maxval(u_host)
!     umin = minval(u_host)
!     if (umax > 1.0e3_dp .or. umin < -1.0e3_dp .or. umax /= umax) then
!       write (*, '(A,I0)') 'BLOW-UP at iter = ', iter
!       exit
!     end if
!   end do

!   write (*, '(A)') '---------------------------------------------'
!   write (*, '(A)') ' Done.'

! contains

!   subroutine report(iter, t, u, u_analytic, n, nblk)
!     integer,  intent(in) :: iter, n, nblk
!     real(dp), intent(in) :: t
!     real(dp), intent(in) :: u(:, :, :)
!     real(dp), intent(in) :: u_analytic(:)
!     real(dp) :: umax, umin, err_max
!     integer  :: j, jmax

!     umax = maxval(u(1, :, 1))
!     umin = minval(u(1, :, 1))
!     err_max = 0.0_dp
!     jmax = 1
!     do j = 1, n
!       if (abs(u(1, j, 1) - u_analytic(j)) > err_max) then
!         err_max = abs(u(1, j, 1) - u_analytic(j))
!         jmax = j
!       end if
!     end do

!     write (*, '(A,I6,A,F10.5,A,ES12.4,A,ES12.4,A,ES12.4,A,I4)') &
!       ' iter=', iter, '  t=', t, &
!       '  umin=', umin, '  umax=', umax, &
!       '  |u-1|_inf=', err_max, '  argmax j=', jmax
!   end subroutine report

! end program test_1d_stability

! program test_2d_stability
!   !! 2D scalar transport stability test.
!   !!
!   !! Solve:   du/dt + U * du/dx = nu * (d2u/dx2 + d2u/dy2)
!   !!          u(0,y) = 1  (x-Dirichlet inflow)
!   !!          u(L,y) = 1  (x-Dirichlet outflow, simplified)
!   !!          du/dy = 0 at y=0 and y=Ly  (y-Neumann, both walls)
!   !!          u(x,y,0) = 1 + small noise
!   !!
!   !! Integrate with AB3 and track max|u - 1| over time.
!   !!
!   !! Purpose: test whether the combination of Dirichlet-x + Neumann-y in
!   !! the derivative stencils is stable under advection-diffusion with
!   !! noisy IC. No pressure, no nonlinear coupling, no 3D. This isolates
!   !! whether the 2D BC combination alone is the source of the 3D blow-up.
!   !!
!   !! Expected outcomes (same interpretation as 1D test):
!   !!   decay   -> 2D BC combination is stable; problem needs pressure or v/w
!   !!   grow    -> 2D BC combination is itself unstable; isolated bug
!   !!   bounded -> steady-state BC error, not the dominant blow-up mechanism
!   !!
!   !! For u tangential to y-walls (du/dy = 0 at wall), we use sym=.true.
!   !! for the y derivatives (same as transeq dispatches for u in y-direction).

!   use iso_fortran_env
!   use cudafor

!   use m_common,          only: dp, pi, BC_DIRICHLET, BC_NEUMANN
!   use m_cuda_common,     only: SZ
!   use m_cuda_tdsops,     only: cuda_tdsops_t, cuda_tdsops_init
!   use m_cuda_exec_thom,  only: exec_thom_tds_compact

!   implicit none

!   ! -------- problem parameters (match 3D cylinder run) --------
!   integer,  parameter :: NX      = 65             ! x-grid points (Dirichlet)
!   integer,  parameter :: NY      = 65             ! y-grid points (Neumann)
!   real(dp), parameter :: LX      = 20.0_dp
!   real(dp), parameter :: LY      = 12.0_dp
!   real(dp), parameter :: U       = 1.0_dp          ! advection velocity (x)
!   real(dp), parameter :: Re      = 300.0_dp
!   real(dp), parameter :: nu      = U * LX / Re
!   real(dp), parameter :: dt      = 0.0075_dp
!   integer,  parameter :: n_iters = 2000
!   integer,  parameter :: n_print = 50

!   real(dp), parameter :: noise_amp = 0.25_dp
!   integer,  parameter :: n_block   = 32    ! needs to be >= NY / SZ roughly

!   ! AB3 coefficients
!   real(dp), parameter :: ab3_a =  23.0_dp / 12.0_dp
!   real(dp), parameter :: ab3_b = -16.0_dp / 12.0_dp
!   real(dp), parameter :: ab3_c =   5.0_dp / 12.0_dp

!   ! -------- fields --------
!   !
!   ! Layout convention to match how transeq_cuda_dist sees the data:
!   !   field_x(SZ, NX, n_block)  with lane index = different y positions
!   !   field_y(SZ, NY, n_block)  with lane index = different x positions
!   !
!   ! To avoid juggling reorders, I keep two copies of u:
!   !   u_x  (SZ, NX, n_block)  -- bundled along x, lane = y within block
!   !   u_y  (SZ, NY, n_block)  -- bundled along y, lane = x within block
!   !
!   ! I manually reorder between them for each directional derivative.
!   ! For simplicity I use n_block = ceil(max(NX,NY)/SZ) and pad unused lanes.

!   real(dp) :: dx, dy
!   real(dp), allocatable :: uf(:, :)             ! (NX, NY) reference layout
!   real(dp), allocatable :: rhs(:, :)
!   real(dp), allocatable :: rhs_nm1(:, :), rhs_nm2(:, :)

!   real(dp), allocatable :: u_x(:, :, :), tmp_x(:, :, :)
!   real(dp), allocatable :: u_y(:, :, :), tmp_y(:, :, :)

!   real(dp), device, allocatable :: u_x_dev(:, :, :), tmp_x_dev(:, :, :)
!   real(dp), device, allocatable :: u_y_dev(:, :, :), tmp_y_dev(:, :, :)

!   type(cuda_tdsops_t) :: op_dx_1st, op_dx_2nd
!   type(cuda_tdsops_t) :: op_dy_2nd_sym         ! Neumann sym for tangential u
!   type(dim3) :: blocks, threads

!   integer  :: i, j, iter, ndevs, devnum, ierr
!   real(dp) :: umax, umin, r

!   ! -------- CUDA setup --------
!   ierr = cudaGetDeviceCount(ndevs)
!   ierr = cudaSetDevice(0)
!   ierr = cudaGetDevice(devnum)

!   ! -------- grid --------
!   dx = LX / real(NX - 1, dp)
!   dy = LY / real(NY - 1, dp)

!   write (*, '(A)') '============================================='
!   write (*, '(A)') ' 2D scalar transport stability test'
!   write (*, '(A)') '============================================='
!   write (*, '(A,I0,A,I0)')  '  Grid       = ', NX, ' x ', NY
!   write (*, '(A,ES12.4,A,ES12.4)') '  L          = ', LX, ' x ', LY
!   write (*, '(A,ES12.4,A,ES12.4)') '  dx, dy     = ', dx, ', ', dy
!   write (*, '(A,ES12.4)')    '  U          = ', U
!   write (*, '(A,ES12.4)')    '  nu         = ', nu
!   write (*, '(A,ES12.4)')    '  Re         = ', Re
!   write (*, '(A,ES12.4)')    '  dt         = ', dt
!   write (*, '(A,ES12.4)')    '  CFL (U)    = ', U * dt / dx
!   write (*, '(A,ES12.4)')    '  diff x     = ', nu * dt / dx**2
!   write (*, '(A,ES12.4)')    '  diff y     = ', nu * dt / dy**2
!   write (*, '(A,I0)')        '  n_iters    = ', n_iters
!   write (*, '(A,ES12.4)')    '  noise amp  = ', noise_amp
!   write (*, '(A)')           '  BC x       = Dirichlet(1) / Dirichlet(1)'
!   write (*, '(A)')           '  BC y       = Neumann(sym) / Neumann(sym)'
!   write (*, '(A)') '---------------------------------------------'

!   ! -------- allocate --------
!   allocate (uf     (NX, NY))
!   allocate (rhs    (NX, NY))
!   allocate (rhs_nm1(NX, NY))
!   allocate (rhs_nm2(NX, NY))

!   allocate (u_x  (SZ, NX, n_block), tmp_x (SZ, NX, n_block))
!   allocate (u_y  (SZ, NY, n_block), tmp_y (SZ, NY, n_block))
!   allocate (u_x_dev  (SZ, NX, n_block), tmp_x_dev (SZ, NX, n_block))
!   allocate (u_y_dev  (SZ, NY, n_block), tmp_y_dev (SZ, NY, n_block))

!   blocks  = dim3(n_block, 1, 1)
!   threads = dim3(SZ, 1, 1)

!   ! -------- initial condition: u = 1 + noise --------
!   call random_seed()
!   do j = 1, NY
!     do i = 1, NX
!       call random_number(r)
!       uf(i, j) = 1.0_dp + noise_amp * (2.0_dp * r - 1.0_dp)
!     end do
!   end do

!   ! Enforce Dirichlet x-BCs in IC
!   uf(1,  :) = 1.0_dp
!   uf(NX, :) = 1.0_dp

!   rhs_nm1 = 0.0_dp
!   rhs_nm2 = 0.0_dp

!   ! -------- build operators --------
!   !
!   ! For u tangential to y-walls (not normal), the sym convention is:
!   !   du/dy uses der1st_sym (sym=.true.)   -- we skip du/dy since U along x
!   !   d2u/dy2 uses der2nd_sym (sym=.true.)
!   !
!   ! For x-Dirichlet, no sym parameter.

!   op_dx_1st = cuda_tdsops_init(NX, dx, operation='first-deriv', &
!                                scheme='compact6', &
!                                bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET)
!   op_dx_2nd = cuda_tdsops_init(NX, dx, operation='second-deriv', &
!                                scheme='compact6', &
!                                bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET)
!   op_dy_2nd_sym = cuda_tdsops_init(NY, dy, operation='second-deriv', &
!                                    scheme='compact6', &
!                                    bc_start=BC_NEUMANN, bc_end=BC_NEUMANN, &
!                                    sym=.true.)

!   ! -------- initial monitor --------
!   call report(0, 0.0_dp, uf, NX, NY)

!   ! -------- time loop --------
!   do iter = 1, n_iters

!     ! ---- pack u into DIR_X bundles ----
!     ! u_x(lane, i, block) with lane+block encoding y index
!     call pack_x(uf, u_x, NX, NY, n_block)

!     ! ---- compute du/dx ----
!     u_x_dev = u_x
!     tmp_x_dev = 0.0_dp
!     call exec_thom_tds_compact(tmp_x_dev, u_x_dev, op_dx_1st, blocks, threads)
!     tmp_x = tmp_x_dev
!     call unpack_x_to_rhs(tmp_x, rhs, NX, NY, n_block, -U)   ! rhs = -U*du/dx

!     ! ---- compute d2u/dx2 ----
!     tmp_x_dev = 0.0_dp
!     call exec_thom_tds_compact(tmp_x_dev, u_x_dev, op_dx_2nd, blocks, threads)
!     tmp_x = tmp_x_dev
!     call unpack_x_add_to_rhs(tmp_x, rhs, NX, NY, n_block, nu)  ! += nu*d2u/dx2

!     ! ---- pack u into DIR_Y bundles ----
!     call pack_y(uf, u_y, NX, NY, n_block)

!     ! ---- compute d2u/dy2 with sym=true ----
!     u_y_dev = u_y
!     tmp_y_dev = 0.0_dp
!     call exec_thom_tds_compact(tmp_y_dev, u_y_dev, op_dy_2nd_sym, blocks, threads)
!     tmp_y = tmp_y_dev
!     call unpack_y_add_to_rhs(tmp_y, rhs, NX, NY, n_block, nu)  ! += nu*d2u/dy2

!     ! ---- AB3 step ----
!     if (iter == 1) then
!       uf = uf + dt * rhs
!     else if (iter == 2) then
!       uf = uf + dt * (1.5_dp * rhs - 0.5_dp * rhs_nm1)
!     else
!       uf = uf + dt * (ab3_a * rhs + ab3_b * rhs_nm1 + ab3_c * rhs_nm2)
!     end if

!     ! ---- enforce Dirichlet x-BCs (Neumann y handled by stencil) ----
!     uf(1,  :) = 1.0_dp
!     uf(NX, :) = 1.0_dp

!     ! ---- shift RHS history ----
!     rhs_nm2 = rhs_nm1
!     rhs_nm1 = rhs

!     ! ---- monitor ----
!     if (mod(iter, n_print) == 0 .or. iter == 1 .or. iter == 2) then
!       call report(iter, iter * dt, uf, NX, NY)
!     end if

!     umax = maxval(uf)
!     umin = minval(uf)
!     if (umax > 1.0e3_dp .or. umin < -1.0e3_dp .or. umax /= umax) then
!       write (*, '(A,I0)') 'BLOW-UP at iter = ', iter
!       exit
!     end if
!   end do

!   write (*, '(A)') '---------------------------------------------'
!   write (*, '(A)') ' Done.'

! contains

!   ! ==========================================================
!   ! Pack 2D array u(NX, NY) into DIR_X bundle u_x(SZ, NX, nb).
!   ! Each (lane, block) pair maps to one y index.  Unused lanes
!   ! are padded with u(i, 1) (arbitrary valid value).
!   ! ==========================================================
!   subroutine pack_x(u_in, u_pack, nx_, ny_, nb)
!     integer,  intent(in)  :: nx_, ny_, nb
!     real(dp), intent(in)  :: u_in(nx_, ny_)
!     real(dp), intent(out) :: u_pack(SZ, nx_, nb)
!     integer :: i, j, lane, blk

!     u_pack = 0.0_dp
!     do j = 1, ny_
!       blk  = (j - 1) / SZ + 1
!       lane = mod(j - 1, SZ) + 1
!       if (blk > nb) cycle
!       do i = 1, nx_
!         u_pack(lane, i, blk) = u_in(i, j)
!       end do
!     end do
!     ! Fill unused lanes with valid-but-unused data (keeps stencils sane)
!     do blk = 1, nb
!       do lane = 1, SZ
!         j = (blk - 1) * SZ + lane
!         if (j > ny_) then
!           do i = 1, nx_
!             u_pack(lane, i, blk) = 1.0_dp
!           end do
!         end if
!       end do
!     end do
!   end subroutine pack_x

!   ! ==========================================================
!   ! Pack 2D array u(NX, NY) into DIR_Y bundle u_y(SZ, NY, nb).
!   ! Each (lane, block) pair maps to one x index.
!   ! ==========================================================
!   subroutine pack_y(u_in, u_pack, nx_, ny_, nb)
!     integer,  intent(in)  :: nx_, ny_, nb
!     real(dp), intent(in)  :: u_in(nx_, ny_)
!     real(dp), intent(out) :: u_pack(SZ, ny_, nb)
!     integer :: i, j, lane, blk

!     u_pack = 0.0_dp
!     do i = 1, nx_
!       blk  = (i - 1) / SZ + 1
!       lane = mod(i - 1, SZ) + 1
!       if (blk > nb) cycle
!       do j = 1, ny_
!         u_pack(lane, j, blk) = u_in(i, j)
!       end do
!     end do
!     ! Pad unused lanes
!     do blk = 1, nb
!       do lane = 1, SZ
!         i = (blk - 1) * SZ + lane
!         if (i > nx_) then
!           do j = 1, ny_
!             u_pack(lane, j, blk) = 1.0_dp
!           end do
!         end if
!       end do
!     end do
!   end subroutine pack_y

!   ! ==========================================================
!   ! Unpack x-bundle back to 2D, writing to rhs with a scale.
!   !   rhs(i, j) = scale * tmp_x(lane(j), i, blk(j))
!   ! ==========================================================
!   subroutine unpack_x_to_rhs(tmp_pack, rhs_out, nx_, ny_, nb, scale)
!     integer,  intent(in)  :: nx_, ny_, nb
!     real(dp), intent(in)  :: tmp_pack(SZ, nx_, nb)
!     real(dp), intent(in)  :: scale
!     real(dp), intent(out) :: rhs_out(nx_, ny_)
!     integer :: i, j, lane, blk

!     do j = 1, ny_
!       blk  = (j - 1) / SZ + 1
!       lane = mod(j - 1, SZ) + 1
!       if (blk > nb) cycle
!       do i = 1, nx_
!         rhs_out(i, j) = scale * tmp_pack(lane, i, blk)
!       end do
!     end do
!   end subroutine unpack_x_to_rhs

!   ! ==========================================================
!   ! Add x-bundle contribution (scaled) to rhs.
!   ! ==========================================================
!   subroutine unpack_x_add_to_rhs(tmp_pack, rhs_out, nx_, ny_, nb, scale)
!     integer,  intent(in)    :: nx_, ny_, nb
!     real(dp), intent(in)    :: tmp_pack(SZ, nx_, nb)
!     real(dp), intent(in)    :: scale
!     real(dp), intent(inout) :: rhs_out(nx_, ny_)
!     integer :: i, j, lane, blk

!     do j = 1, ny_
!       blk  = (j - 1) / SZ + 1
!       lane = mod(j - 1, SZ) + 1
!       if (blk > nb) cycle
!       do i = 1, nx_
!         rhs_out(i, j) = rhs_out(i, j) + scale * tmp_pack(lane, i, blk)
!       end do
!     end do
!   end subroutine unpack_x_add_to_rhs

!   ! ==========================================================
!   ! Add y-bundle contribution (scaled) to rhs.
!   ! ==========================================================
!   subroutine unpack_y_add_to_rhs(tmp_pack, rhs_out, nx_, ny_, nb, scale)
!     integer,  intent(in)    :: nx_, ny_, nb
!     real(dp), intent(in)    :: tmp_pack(SZ, ny_, nb)
!     real(dp), intent(in)    :: scale
!     real(dp), intent(inout) :: rhs_out(nx_, ny_)
!     integer :: i, j, lane, blk

!     do i = 1, nx_
!       blk  = (i - 1) / SZ + 1
!       lane = mod(i - 1, SZ) + 1
!       if (blk > nb) cycle
!       do j = 1, ny_
!         rhs_out(i, j) = rhs_out(i, j) + scale * tmp_pack(lane, j, blk)
!       end do
!     end do
!   end subroutine unpack_y_add_to_rhs

!   ! ==========================================================
!   ! Report monitor
!   ! ==========================================================
!   subroutine report(iter, t, u_in, nx_, ny_)
!     integer,  intent(in) :: iter, nx_, ny_
!     real(dp), intent(in) :: t, u_in(nx_, ny_)
!     real(dp) :: umax, umin, err_max
!     integer  :: i, j, imax, jmax

!     umax = maxval(u_in)
!     umin = minval(u_in)
!     err_max = 0.0_dp
!     imax = 1; jmax = 1
!     do j = 1, ny_
!       do i = 1, nx_
!         if (abs(u_in(i, j) - 1.0_dp) > err_max) then
!           err_max = abs(u_in(i, j) - 1.0_dp)
!           imax = i; jmax = j
!         end if
!       end do
!     end do

!     write (*, '(A,I6,A,F10.5,A,ES12.4,A,ES12.4,A,ES12.4,A,I3,A,I3,A)') &
!       ' iter=', iter, '  t=', t, &
!       '  umin=', umin, '  umax=', umax, &
!       '  |u-1|_inf=', err_max, &
!       '  argmax (i=', imax, ', j=', jmax, ')'
!   end subroutine report

! end program test_2d_stability

! program test_2d_uv_coupled
!   !! 2D coupled u-v advection-diffusion stability test (NO PRESSURE).
!   !!
!   !! Solve the momentum equations WITHOUT the incompressibility constraint:
!   !!    du/dt = -u*du/dx - v*du/dy + nu*(d2u/dx2 + d2u/dy2)
!   !!    dv/dt = -u*dv/dx - v*dv/dy + nu*(d2v/dx2 + d2v/dy2)
!   !!
!   !! BCs (matching the 3D cylinder case, x and y only):
!   !!    u: Dirichlet(1) at x=0 and x=L
!   !!       Neumann at y=0 and y=Ly (tangential, sym=.true. for u in y)
!   !!    v: Dirichlet(0) at x=0 and x=L
!   !!       Neumann at y=0 and y=Ly (normal, sym=.false. for v in y)
!   !!
!   !! IC:
!   !!    u(x,y,0) = 1 + noise
!   !!    v(x,y,0) = 0 + noise
!   !!
!   !! This test stresses nonlinear cross-coupling:
!   !!   d(u)/dy with sym=.true.  (u tangential to y-walls)
!   !!   d(v)/dy with sym=.false. (v normal to y-walls)
!   !!   d(u)/dx and d(v)/dx at x-Dirichlet walls
!   !!   Nonlinear convection: v*du/dy and u*dv/dx etc.
!   !!
!   !! NO pressure projection, NO incompressibility.
!   !!
!   !! Purpose: determine whether the nonlinear cross-coupling with
!   !! different sym dispatches for u and v is stable in 2D without
!   !! pressure coupling. If it blows up -> the bug is in transeq's
!   !! convective cross-term handling. If it stays stable -> the bug
!   !! requires the pressure projection.

!   use iso_fortran_env
!   use cudafor

!   use m_common,          only: dp, pi, BC_DIRICHLET, BC_NEUMANN
!   use m_cuda_common,     only: SZ
!   use m_cuda_tdsops,     only: cuda_tdsops_t, cuda_tdsops_init
!   use m_cuda_exec_thom,  only: exec_thom_tds_compact

!   implicit none

!   ! -------- problem parameters (match 3D cylinder run) --------
!   integer,  parameter :: NX      = 65
!   integer,  parameter :: NY      = 65
!   real(dp), parameter :: LX      = 20.0_dp
!   real(dp), parameter :: LY      = 12.0_dp
!   real(dp), parameter :: U0      = 1.0_dp          ! inflow u-velocity
!   real(dp), parameter :: V0      = 0.0_dp          ! inflow v-velocity
!   real(dp), parameter :: Re      = 300.0_dp
!   real(dp), parameter :: nu      = U0 * LX / Re
!   real(dp), parameter :: dt      = 0.0075_dp
!   integer,  parameter :: n_iters = 2000
!   integer,  parameter :: n_print = 50

!   real(dp), parameter :: noise_amp = 0.25_dp
!   integer,  parameter :: n_block   = 32

!   ! AB3 coefficients
!   real(dp), parameter :: ab3_a =  23.0_dp / 12.0_dp
!   real(dp), parameter :: ab3_b = -16.0_dp / 12.0_dp
!   real(dp), parameter :: ab3_c =   5.0_dp / 12.0_dp

!   real(dp) :: dx, dy
!   real(dp), allocatable :: uf(:, :), vf(:, :)
!   real(dp), allocatable :: rhs_u(:, :), rhs_v(:, :)
!   real(dp), allocatable :: rhs_u_nm1(:, :), rhs_u_nm2(:, :)
!   real(dp), allocatable :: rhs_v_nm1(:, :), rhs_v_nm2(:, :)

!   real(dp), allocatable :: u_x(:, :, :), tmp_x(:, :, :)
!   real(dp), allocatable :: u_y(:, :, :), tmp_y(:, :, :)

!   real(dp), device, allocatable :: u_x_dev(:, :, :), tmp_x_dev(:, :, :)
!   real(dp), device, allocatable :: u_y_dev(:, :, :), tmp_y_dev(:, :, :)

!   ! Operators
!   type(cuda_tdsops_t) :: op_dx_1st, op_dx_2nd
!   type(cuda_tdsops_t) :: op_dy_1st_sym,  op_dy_2nd_sym     ! u: tangential (sym=T)
!   type(cuda_tdsops_t) :: op_dy_1st_nsym, op_dy_2nd_nsym    ! v: normal     (sym=F)
!   type(dim3) :: blocks, threads

!   real(dp), allocatable :: dudx(:, :), dudy(:, :), d2udx2(:, :), d2udy2(:, :)
!   real(dp), allocatable :: dvdx(:, :), dvdy(:, :), d2vdx2(:, :), d2vdy2(:, :)

!   integer  :: i, j, iter, ndevs, devnum, ierr
!   real(dp) :: umax, umin, vmax, vmin, r

!   ierr = cudaGetDeviceCount(ndevs)
!   ierr = cudaSetDevice(0)
!   ierr = cudaGetDevice(devnum)

!   dx = LX / real(NX - 1, dp)
!   dy = LY / real(NY - 1, dp)

!   write (*, '(A)') '============================================='
!   write (*, '(A)') ' 2D coupled u-v advection-diffusion test'
!   write (*, '(A)') '  (NO pressure, NO incompressibility)'
!   write (*, '(A)') '============================================='
!   write (*, '(A,I0,A,I0)')           '  Grid       = ', NX, ' x ', NY
!   write (*, '(A,ES12.4,A,ES12.4)')   '  L          = ', LX, ' x ', LY
!   write (*, '(A,ES12.4,A,ES12.4)')   '  dx, dy     = ', dx, ', ', dy
!   write (*, '(A,ES12.4,A,ES12.4)')   '  U0, V0     = ', U0, ', ', V0
!   write (*, '(A,ES12.4)')            '  nu         = ', nu
!   write (*, '(A,ES12.4)')            '  Re         = ', Re
!   write (*, '(A,ES12.4)')            '  dt         = ', dt
!   write (*, '(A,ES12.4,A,ES12.4)')   '  noise, dt  = ', noise_amp, ', ', dt
!   write (*, '(A)')                   '  u: Dir(1)/Dir(1), Neumann-sym'
!   write (*, '(A)')                   '  v: Dir(0)/Dir(0), Neumann-nsym'
!   write (*, '(A)') '---------------------------------------------'

!   allocate (uf(NX, NY),        vf(NX, NY))
!   allocate (rhs_u(NX, NY),     rhs_v(NX, NY))
!   allocate (rhs_u_nm1(NX, NY), rhs_u_nm2(NX, NY))
!   allocate (rhs_v_nm1(NX, NY), rhs_v_nm2(NX, NY))
!   allocate (dudx(NX, NY), dudy(NX, NY))
!   allocate (d2udx2(NX, NY), d2udy2(NX, NY))
!   allocate (dvdx(NX, NY), dvdy(NX, NY))
!   allocate (d2vdx2(NX, NY), d2vdy2(NX, NY))

!   allocate (u_x(SZ, NX, n_block),     tmp_x(SZ, NX, n_block))
!   allocate (u_y(SZ, NY, n_block),     tmp_y(SZ, NY, n_block))
!   allocate (u_x_dev(SZ, NX, n_block), tmp_x_dev(SZ, NX, n_block))
!   allocate (u_y_dev(SZ, NY, n_block), tmp_y_dev(SZ, NY, n_block))

!   blocks  = dim3(n_block, 1, 1)
!   threads = dim3(SZ, 1, 1)

!   ! -------- initial condition --------
!   call random_seed()
!   do j = 1, NY
!     do i = 1, NX
!       call random_number(r)
!       uf(i, j) = U0 + noise_amp * (2.0_dp * r - 1.0_dp)
!       call random_number(r)
!       vf(i, j) = V0 + noise_amp * (2.0_dp * r - 1.0_dp)
!     end do
!   end do
!   uf(1,  :) = U0; uf(NX, :) = U0
!   vf(1,  :) = V0; vf(NX, :) = V0

!   rhs_u_nm1 = 0.0_dp; rhs_u_nm2 = 0.0_dp
!   rhs_v_nm1 = 0.0_dp; rhs_v_nm2 = 0.0_dp

!   ! -------- operators --------
!   op_dx_1st = cuda_tdsops_init(NX, dx, operation='first-deriv', scheme='compact6', &
!                                bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET)
!   op_dx_2nd = cuda_tdsops_init(NX, dx, operation='second-deriv', scheme='compact6', &
!                                bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET)

!   op_dy_1st_sym = cuda_tdsops_init(NY, dy, operation='first-deriv', scheme='compact6', &
!                                    bc_start=BC_NEUMANN, bc_end=BC_NEUMANN, sym=.true.)
!   op_dy_2nd_sym = cuda_tdsops_init(NY, dy, operation='second-deriv', scheme='compact6', &
!                                    bc_start=BC_NEUMANN, bc_end=BC_NEUMANN, sym=.true.)

!   op_dy_1st_nsym = cuda_tdsops_init(NY, dy, operation='first-deriv', scheme='compact6', &
!                                     bc_start=BC_NEUMANN, bc_end=BC_NEUMANN, sym=.false.)
!   op_dy_2nd_nsym = cuda_tdsops_init(NY, dy, operation='second-deriv', scheme='compact6', &
!                                     bc_start=BC_NEUMANN, bc_end=BC_NEUMANN, sym=.false.)

!   call report(0, 0.0_dp, uf, vf, NX, NY)

!   ! -------- time loop --------
!   do iter = 1, n_iters

!     ! === derivatives of u ===
!     call pack_x(uf, u_x, NX, NY, n_block)
!     u_x_dev = u_x

!     tmp_x_dev = 0.0_dp
!     call exec_thom_tds_compact(tmp_x_dev, u_x_dev, op_dx_1st, blocks, threads)
!     tmp_x = tmp_x_dev
!     call unpack_x(tmp_x, dudx, NX, NY, n_block)

!     tmp_x_dev = 0.0_dp
!     call exec_thom_tds_compact(tmp_x_dev, u_x_dev, op_dx_2nd, blocks, threads)
!     tmp_x = tmp_x_dev
!     call unpack_x(tmp_x, d2udx2, NX, NY, n_block)

!     call pack_y(uf, u_y, NX, NY, n_block)
!     u_y_dev = u_y

!     tmp_y_dev = 0.0_dp
!     call exec_thom_tds_compact(tmp_y_dev, u_y_dev, op_dy_1st_sym, blocks, threads)
!     tmp_y = tmp_y_dev
!     call unpack_y(tmp_y, dudy, NX, NY, n_block)

!     tmp_y_dev = 0.0_dp
!     call exec_thom_tds_compact(tmp_y_dev, u_y_dev, op_dy_2nd_sym, blocks, threads)
!     tmp_y = tmp_y_dev
!     call unpack_y(tmp_y, d2udy2, NX, NY, n_block)

!     ! === derivatives of v ===
!     call pack_x(vf, u_x, NX, NY, n_block)
!     u_x_dev = u_x

!     tmp_x_dev = 0.0_dp
!     call exec_thom_tds_compact(tmp_x_dev, u_x_dev, op_dx_1st, blocks, threads)
!     tmp_x = tmp_x_dev
!     call unpack_x(tmp_x, dvdx, NX, NY, n_block)

!     tmp_x_dev = 0.0_dp
!     call exec_thom_tds_compact(tmp_x_dev, u_x_dev, op_dx_2nd, blocks, threads)
!     tmp_x = tmp_x_dev
!     call unpack_x(tmp_x, d2vdx2, NX, NY, n_block)

!     call pack_y(vf, u_y, NX, NY, n_block)
!     u_y_dev = u_y

!     tmp_y_dev = 0.0_dp
!     call exec_thom_tds_compact(tmp_y_dev, u_y_dev, op_dy_1st_nsym, blocks, threads)
!     tmp_y = tmp_y_dev
!     call unpack_y(tmp_y, dvdy, NX, NY, n_block)

!     tmp_y_dev = 0.0_dp
!     call exec_thom_tds_compact(tmp_y_dev, u_y_dev, op_dy_2nd_nsym, blocks, threads)
!     tmp_y = tmp_y_dev
!     call unpack_y(tmp_y, d2vdy2, NX, NY, n_block)

!     ! === assemble RHS ===
!     do j = 1, NY
!       do i = 1, NX
!         rhs_u(i, j) = -uf(i,j)*dudx(i,j) - vf(i,j)*dudy(i,j) &
!                        + nu*(d2udx2(i,j) + d2udy2(i,j))
!         rhs_v(i, j) = -uf(i,j)*dvdx(i,j) - vf(i,j)*dvdy(i,j) &
!                        + nu*(d2vdx2(i,j) + d2vdy2(i,j))
!       end do
!     end do

!     ! === AB3 step ===
!     if (iter == 1) then
!       uf = uf + dt * rhs_u
!       vf = vf + dt * rhs_v
!     else if (iter == 2) then
!       uf = uf + dt * (1.5_dp*rhs_u - 0.5_dp*rhs_u_nm1)
!       vf = vf + dt * (1.5_dp*rhs_v - 0.5_dp*rhs_v_nm1)
!     else
!       uf = uf + dt * (ab3_a*rhs_u + ab3_b*rhs_u_nm1 + ab3_c*rhs_u_nm2)
!       vf = vf + dt * (ab3_a*rhs_v + ab3_b*rhs_v_nm1 + ab3_c*rhs_v_nm2)
!     end if

!     ! Enforce x-Dirichlet BCs
!     uf(1,  :) = U0; uf(NX, :) = U0
!     vf(1,  :) = V0; vf(NX, :) = V0

!     ! Shift RHS history
!     rhs_u_nm2 = rhs_u_nm1; rhs_u_nm1 = rhs_u
!     rhs_v_nm2 = rhs_v_nm1; rhs_v_nm1 = rhs_v

!     if (mod(iter, n_print) == 0 .or. iter == 1 .or. iter == 2) then
!       call report(iter, iter * dt, uf, vf, NX, NY)
!     end if

!     umax = maxval(uf); umin = minval(uf)
!     vmax = maxval(vf); vmin = minval(vf)
!     if (umax > 1.0e3_dp .or. umin < -1.0e3_dp .or. &
!         vmax > 1.0e3_dp .or. vmin < -1.0e3_dp .or. &
!         umax /= umax .or. vmax /= vmax) then
!       write (*, '(A,I0)') 'BLOW-UP at iter = ', iter
!       exit
!     end if
!   end do

!   write (*, '(A)') '---------------------------------------------'
!   write (*, '(A)') ' Done.'

! contains

!   subroutine pack_x(u_in, u_pack, nx_, ny_, nb)
!     integer,  intent(in)  :: nx_, ny_, nb
!     real(dp), intent(in)  :: u_in(nx_, ny_)
!     real(dp), intent(out) :: u_pack(SZ, nx_, nb)
!     integer :: i, j, lane, blk
!     u_pack = 0.0_dp
!     do j = 1, ny_
!       blk  = (j - 1) / SZ + 1
!       lane = mod(j - 1, SZ) + 1
!       if (blk > nb) cycle
!       do i = 1, nx_
!         u_pack(lane, i, blk) = u_in(i, j)
!       end do
!     end do
!     do blk = 1, nb
!       do lane = 1, SZ
!         j = (blk - 1) * SZ + lane
!         if (j > ny_) u_pack(lane, :, blk) = 1.0_dp
!       end do
!     end do
!   end subroutine pack_x

!   subroutine pack_y(u_in, u_pack, nx_, ny_, nb)
!     integer,  intent(in)  :: nx_, ny_, nb
!     real(dp), intent(in)  :: u_in(nx_, ny_)
!     real(dp), intent(out) :: u_pack(SZ, ny_, nb)
!     integer :: i, j, lane, blk
!     u_pack = 0.0_dp
!     do i = 1, nx_
!       blk  = (i - 1) / SZ + 1
!       lane = mod(i - 1, SZ) + 1
!       if (blk > nb) cycle
!       do j = 1, ny_
!         u_pack(lane, j, blk) = u_in(i, j)
!       end do
!     end do
!     do blk = 1, nb
!       do lane = 1, SZ
!         i = (blk - 1) * SZ + lane
!         if (i > nx_) u_pack(lane, :, blk) = 1.0_dp
!       end do
!     end do
!   end subroutine pack_y

!   subroutine unpack_x(tmp_pack, out, nx_, ny_, nb)
!     integer,  intent(in)  :: nx_, ny_, nb
!     real(dp), intent(in)  :: tmp_pack(SZ, nx_, nb)
!     real(dp), intent(out) :: out(nx_, ny_)
!     integer :: i, j, lane, blk
!     do j = 1, ny_
!       blk  = (j - 1) / SZ + 1
!       lane = mod(j - 1, SZ) + 1
!       if (blk > nb) cycle
!       do i = 1, nx_
!         out(i, j) = tmp_pack(lane, i, blk)
!       end do
!     end do
!   end subroutine unpack_x

!   subroutine unpack_y(tmp_pack, out, nx_, ny_, nb)
!     integer,  intent(in)  :: nx_, ny_, nb
!     real(dp), intent(in)  :: tmp_pack(SZ, ny_, nb)
!     real(dp), intent(out) :: out(nx_, ny_)
!     integer :: i, j, lane, blk
!     do i = 1, nx_
!       blk  = (i - 1) / SZ + 1
!       lane = mod(i - 1, SZ) + 1
!       if (blk > nb) cycle
!       do j = 1, ny_
!         out(i, j) = tmp_pack(lane, j, blk)
!       end do
!     end do
!   end subroutine unpack_y

!   subroutine report(iter, t, u_in, v_in, nx_, ny_)
!     integer,  intent(in) :: iter, nx_, ny_
!     real(dp), intent(in) :: t, u_in(nx_, ny_), v_in(nx_, ny_)
!     real(dp) :: umax, umin, vmax, vmin
!     umax = maxval(u_in); umin = minval(u_in)
!     vmax = maxval(v_in); vmin = minval(v_in)
!     write (*, '(A,I6,A,F9.4,A,ES11.3,A,ES11.3,A,ES11.3,A,ES11.3)') &
!       ' iter=', iter, '  t=', t, &
!       '  umin=', umin, ' umax=', umax, &
!       ' vmin=', vmin, ' vmax=', vmax
!   end subroutine report

! end program test_2d_uv_coupled

program test_2d_uv_skew
  !! 2D coupled u-v SKEW-SYMMETRIC momentum equations (NO PRESSURE).
  !!
  !! Exactly matches transeq's skew-symmetric form:
  !!   du/dt = -0.5*(u*du/dx + d(uu)/dx) - 0.5*(v*du/dy + d(vu)/dy)
  !!           + nu*(d2u/dx2 + d2u/dy2)
  !!   dv/dt = -0.5*(u*dv/dx + d(uv)/dx) - 0.5*(v*dv/dy + d(vv)/dy)
  !!           + nu*(d2v/dx2 + d2v/dy2)
  !!
  !! BCs:
  !!   u: Dirichlet(1) at x=0,L;  Neumann-sym at y=0,Ly  (tangential in y)
  !!   v: Dirichlet(0) at x=0,L;  Neumann-nsym at y=0,Ly (normal in y)
  !!
  !! Operator dispatch (matches transeq):
  !!
  !!   x-direction (Dirichlet at both ends, no sym flag relevant):
  !!     For u (convecting self):
  !!         du/dx  : der1st      (Dirichlet)
  !!         d(uu)/dx: der1st_sym  (Dirichlet -- sym flag doesn't matter)
  !!         d2u/dx2: der2nd      (Dirichlet)
  !!     For v (convected by u):
  !!         dv/dx  : der1st_sym
  !!         d(uv)/dx: der1st
  !!         d2v/dx2: der2nd_sym
  !!
  !!   y-direction (Neumann):
  !!     In y, v is the convecting velocity.  For v (convecting self):
  !!         dv/dy  : der1st      (Neumann, sym=.false.)  [normal, anti-sym]
  !!         d(vv)/dy: der1st_sym (Neumann, sym=.true.)   [self-sq, sym]
  !!         d2v/dy2: der2nd      (Neumann, sym=.false.)
  !!     For u (convected by v):
  !!         du/dy  : der1st_sym  (Neumann, sym=.true.)   [tangent, sym]
  !!         d(vu)/dy: der1st     (Neumann, sym=.false.)  [normal-tang]
  !!         d2u/dy2: der2nd_sym  (Neumann, sym=.true.)
  !!
  !! IC: u = 1 + noise, v = 0 + noise
  !! Time: AB3
  !!
  !! Purpose: if skew-symmetric form in 2D WITHOUT pressure is stable,
  !! then the 3D blow-up needs pressure coupling; if unstable, there is
  !! a genuine bug in the operator-dispatch logic or a stencil coefficient.

  use iso_fortran_env
  use cudafor

  use m_common,          only: dp, pi, BC_DIRICHLET, BC_NEUMANN
  use m_cuda_common,     only: SZ
  use m_cuda_tdsops,     only: cuda_tdsops_t, cuda_tdsops_init
  use m_cuda_exec_thom,  only: exec_thom_tds_compact

  implicit none

  ! -------- parameters --------
  integer,  parameter :: NX      = 65
  integer,  parameter :: NY      = 65
  real(dp), parameter :: LX      = 20.0_dp
  real(dp), parameter :: LY      = 12.0_dp
  real(dp), parameter :: U0      = 1.0_dp
  real(dp), parameter :: V0      = 0.0_dp
  real(dp), parameter :: Re      = 300.0_dp
  real(dp), parameter :: nu      = U0 * LX / Re
  real(dp), parameter :: dt      = 0.0075_dp
  integer,  parameter :: n_iters = 2000
  integer,  parameter :: n_print = 50

  real(dp), parameter :: noise_amp = 0.25_dp
  integer,  parameter :: n_block   = 32

  real(dp), parameter :: ab3_a =  23.0_dp / 12.0_dp
  real(dp), parameter :: ab3_b = -16.0_dp / 12.0_dp
  real(dp), parameter :: ab3_c =   5.0_dp / 12.0_dp

  real(dp) :: dx, dy
  real(dp), allocatable :: uf(:, :), vf(:, :)
  real(dp), allocatable :: rhs_u(:, :), rhs_v(:, :)
  real(dp), allocatable :: rhs_u_nm1(:, :), rhs_u_nm2(:, :)
  real(dp), allocatable :: rhs_v_nm1(:, :), rhs_v_nm2(:, :)
  real(dp), allocatable :: uu(:, :), vv(:, :), uv(:, :)

  real(dp), allocatable :: buf_x(:, :, :), tmp_x(:, :, :)
  real(dp), allocatable :: buf_y(:, :, :), tmp_y(:, :, :)
  real(dp), device, allocatable :: buf_x_dev(:, :, :), tmp_x_dev(:, :, :)
  real(dp), device, allocatable :: buf_y_dev(:, :, :), tmp_y_dev(:, :, :)

  ! Operators
  ! x-direction (Dirichlet)
  type(cuda_tdsops_t) :: op_dx_1st, op_dx_1st_sym
  type(cuda_tdsops_t) :: op_dx_2nd, op_dx_2nd_sym
  ! y-direction (Neumann)
  type(cuda_tdsops_t) :: op_dy_1st,     op_dy_1st_sym     ! nsym, sym
  type(cuda_tdsops_t) :: op_dy_2nd,     op_dy_2nd_sym

  type(dim3) :: blocks, threads

  real(dp), allocatable :: tmpa(:, :), tmpb(:, :)
  integer :: i, j, iter, ndevs, devnum, ierr
  real(dp) :: umax, umin, vmax, vmin, r

  ierr = cudaGetDeviceCount(ndevs)
  ierr = cudaSetDevice(0)
  ierr = cudaGetDevice(devnum)

  dx = LX / real(NX - 1, dp)
  dy = LY / real(NY - 1, dp)

  write (*, '(A)') '============================================='
  write (*, '(A)') ' 2D u-v SKEW-SYMMETRIC momentum test'
  write (*, '(A)') '  (NO pressure, NO incompressibility)'
  write (*, '(A)') '  Matches transeq dispatch exactly'
  write (*, '(A)') '============================================='
  write (*, '(A,I0,A,I0)')           '  Grid       = ', NX, ' x ', NY
  write (*, '(A,ES12.4,A,ES12.4)')   '  L          = ', LX, ' x ', LY
  write (*, '(A,ES12.4,A,ES12.4)')   '  dx, dy     = ', dx, ', ', dy
  write (*, '(A,ES12.4,A,ES12.4)')   '  U0, V0     = ', U0, ', ', V0
  write (*, '(A,ES12.4)')            '  nu, Re     = ', nu
  write (*, '(A,ES12.4)')            '  Re         = ', Re
  write (*, '(A,ES12.4)')            '  dt         = ', dt
  write (*, '(A,ES12.4)')            '  noise      = ', noise_amp
  write (*, '(A)') '---------------------------------------------'

  allocate (uf(NX, NY),        vf(NX, NY))
  allocate (uu(NX, NY), vv(NX, NY), uv(NX, NY))
  allocate (rhs_u(NX, NY),     rhs_v(NX, NY))
  allocate (rhs_u_nm1(NX, NY), rhs_u_nm2(NX, NY))
  allocate (rhs_v_nm1(NX, NY), rhs_v_nm2(NX, NY))
  allocate (tmpa(NX, NY), tmpb(NX, NY))

  allocate (buf_x(SZ, NX, n_block), tmp_x(SZ, NX, n_block))
  allocate (buf_y(SZ, NY, n_block), tmp_y(SZ, NY, n_block))
  allocate (buf_x_dev(SZ, NX, n_block), tmp_x_dev(SZ, NX, n_block))
  allocate (buf_y_dev(SZ, NY, n_block), tmp_y_dev(SZ, NY, n_block))

  blocks  = dim3(n_block, 1, 1)
  threads = dim3(SZ, 1, 1)

  ! -------- IC --------
  call random_seed()
  do j = 1, NY
    do i = 1, NX
      call random_number(r)
      uf(i, j) = U0 + noise_amp * (2.0_dp * r - 1.0_dp)
      call random_number(r)
      vf(i, j) = V0 + noise_amp * (2.0_dp * r - 1.0_dp)
    end do
  end do
  uf(1,  :) = U0; uf(NX, :) = U0
  vf(1,  :) = V0; vf(NX, :) = V0

  rhs_u_nm1 = 0.0_dp; rhs_u_nm2 = 0.0_dp
  rhs_v_nm1 = 0.0_dp; rhs_v_nm2 = 0.0_dp

  ! -------- operators --------
  ! x-Dirichlet: sym flag is irrelevant for Dirichlet branch, but we pass
  ! the appropriate one for clarity
  op_dx_1st = cuda_tdsops_init(NX, dx, operation='first-deriv', scheme='compact6', &
                               bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET)
  op_dx_1st_sym = cuda_tdsops_init(NX, dx, operation='first-deriv', scheme='compact6', &
                                   bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET, &
                                   sym=.true.)
  op_dx_2nd = cuda_tdsops_init(NX, dx, operation='second-deriv', scheme='compact6', &
                               bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET)
  op_dx_2nd_sym = cuda_tdsops_init(NX, dx, operation='second-deriv', scheme='compact6', &
                                   bc_start=BC_DIRICHLET, bc_end=BC_DIRICHLET, &
                                   sym=.true.)

  ! y-Neumann
  op_dy_1st = cuda_tdsops_init(NY, dy, operation='first-deriv', scheme='compact6', &
                               bc_start=BC_NEUMANN, bc_end=BC_NEUMANN, sym=.false.)
  op_dy_1st_sym = cuda_tdsops_init(NY, dy, operation='first-deriv', scheme='compact6', &
                                   bc_start=BC_NEUMANN, bc_end=BC_NEUMANN, sym=.true.)
  op_dy_2nd = cuda_tdsops_init(NY, dy, operation='second-deriv', scheme='compact6', &
                               bc_start=BC_NEUMANN, bc_end=BC_NEUMANN, sym=.false.)
  op_dy_2nd_sym = cuda_tdsops_init(NY, dy, operation='second-deriv', scheme='compact6', &
                                   bc_start=BC_NEUMANN, bc_end=BC_NEUMANN, sym=.true.)

  call report(0, 0.0_dp, uf, vf, NX, NY)

  ! -------- time loop --------
  do iter = 1, n_iters

    ! Precompute products
    do j = 1, NY
      do i = 1, NX
        uu(i, j) = uf(i,j) * uf(i,j)
        vv(i, j) = vf(i,j) * vf(i,j)
        uv(i, j) = uf(i,j) * vf(i,j)
      end do
    end do

    ! Initialize RHS to zero
    rhs_u = 0.0_dp
    rhs_v = 0.0_dp

    ! =============================================================
    ! x-direction contributions
    ! =============================================================

    ! ---- For u-equation: d(u)/dx with der1st (Dirichlet) ----
    call deriv_x(uf, tmpa, op_dx_1st)
    ! ---- For u-equation: d(uu)/dx with der1st_sym ----
    call deriv_x(uu, tmpb, op_dx_1st_sym)
    ! -0.5 * (u * du/dx + d(uu)/dx)
    do j = 1, NY
      do i = 1, NX
        rhs_u(i,j) = rhs_u(i,j) - 0.5_dp * (uf(i,j)*tmpa(i,j) + tmpb(i,j))
      end do
    end do
    ! ---- For u-equation: nu * d2u/dx2 with der2nd ----
    call deriv_x(uf, tmpa, op_dx_2nd)
    do j = 1, NY
      do i = 1, NX
        rhs_u(i,j) = rhs_u(i,j) + nu * tmpa(i,j)
      end do
    end do

    ! ---- For v-equation: d(v)/dx with der1st_sym ----
    call deriv_x(vf, tmpa, op_dx_1st_sym)
    ! ---- For v-equation: d(uv)/dx with der1st ----
    call deriv_x(uv, tmpb, op_dx_1st)
    ! -0.5 * (u * dv/dx + d(uv)/dx)
    do j = 1, NY
      do i = 1, NX
        rhs_v(i,j) = rhs_v(i,j) - 0.5_dp * (uf(i,j)*tmpa(i,j) + tmpb(i,j))
      end do
    end do
    ! ---- For v-equation: nu * d2v/dx2 with der2nd_sym ----
    call deriv_x(vf, tmpa, op_dx_2nd_sym)
    do j = 1, NY
      do i = 1, NX
        rhs_v(i,j) = rhs_v(i,j) + nu * tmpa(i,j)
      end do
    end do

    ! =============================================================
    ! y-direction contributions
    ! =============================================================

    ! ---- For v-equation: d(v)/dy with der1st (sym=F) ----
    call deriv_y(vf, tmpa, op_dy_1st)
    ! ---- For v-equation: d(vv)/dy with der1st_sym (sym=T) ----
    call deriv_y(vv, tmpb, op_dy_1st_sym)
    ! -0.5 * (v * dv/dy + d(vv)/dy)
    do j = 1, NY
      do i = 1, NX
        rhs_v(i,j) = rhs_v(i,j) - 0.5_dp * (vf(i,j)*tmpa(i,j) + tmpb(i,j))
      end do
    end do
    ! ---- For v-equation: nu * d2v/dy2 with der2nd (sym=F) ----
    call deriv_y(vf, tmpa, op_dy_2nd)
    do j = 1, NY
      do i = 1, NX
        rhs_v(i,j) = rhs_v(i,j) + nu * tmpa(i,j)
      end do
    end do

    ! ---- For u-equation: d(u)/dy with der1st_sym (sym=T) ----
    call deriv_y(uf, tmpa, op_dy_1st_sym)
    ! ---- For u-equation: d(vu)/dy with der1st (sym=F) ----
    call deriv_y(uv, tmpb, op_dy_1st)
    ! -0.5 * (v * du/dy + d(vu)/dy)
    do j = 1, NY
      do i = 1, NX
        rhs_u(i,j) = rhs_u(i,j) - 0.5_dp * (vf(i,j)*tmpa(i,j) + tmpb(i,j))
      end do
    end do
    ! ---- For u-equation: nu * d2u/dy2 with der2nd_sym (sym=T) ----
    call deriv_y(uf, tmpa, op_dy_2nd_sym)
    do j = 1, NY
      do i = 1, NX
        rhs_u(i,j) = rhs_u(i,j) + nu * tmpa(i,j)
      end do
    end do

    ! =============================================================
    ! AB3 step
    ! =============================================================
    if (iter == 1) then
      uf = uf + dt * rhs_u
      vf = vf + dt * rhs_v
    else if (iter == 2) then
      uf = uf + dt * (1.5_dp*rhs_u - 0.5_dp*rhs_u_nm1)
      vf = vf + dt * (1.5_dp*rhs_v - 0.5_dp*rhs_v_nm1)
    else
      uf = uf + dt * (ab3_a*rhs_u + ab3_b*rhs_u_nm1 + ab3_c*rhs_u_nm2)
      vf = vf + dt * (ab3_a*rhs_v + ab3_b*rhs_v_nm1 + ab3_c*rhs_v_nm2)
    end if

    uf(1,  :) = U0; uf(NX, :) = U0
    vf(1,  :) = V0; vf(NX, :) = V0

    rhs_u_nm2 = rhs_u_nm1; rhs_u_nm1 = rhs_u
    rhs_v_nm2 = rhs_v_nm1; rhs_v_nm1 = rhs_v

    if (mod(iter, n_print) == 0 .or. iter == 1 .or. iter == 2) then
      call report(iter, iter * dt, uf, vf, NX, NY)
    end if

    umax = maxval(uf); umin = minval(uf)
    vmax = maxval(vf); vmin = minval(vf)
    if (umax > 1.0e3_dp .or. umin < -1.0e3_dp .or. &
        vmax > 1.0e3_dp .or. vmin < -1.0e3_dp .or. &
        umax /= umax .or. vmax /= vmax) then
      write (*, '(A,I0)') 'BLOW-UP at iter = ', iter
      exit
    end if
  end do

  write (*, '(A)') '---------------------------------------------'
  write (*, '(A)') ' Done.'

contains

  subroutine deriv_x(field, out, op)
    real(dp), intent(in)  :: field(NX, NY)
    real(dp), intent(out) :: out(NX, NY)
    type(cuda_tdsops_t), intent(in) :: op
    call pack_x(field, buf_x, NX, NY, n_block)
    buf_x_dev = buf_x
    tmp_x_dev = 0.0_dp
    call exec_thom_tds_compact(tmp_x_dev, buf_x_dev, op, blocks, threads)
    tmp_x = tmp_x_dev
    call unpack_x(tmp_x, out, NX, NY, n_block)
  end subroutine

  subroutine deriv_y(field, out, op)
    real(dp), intent(in)  :: field(NX, NY)
    real(dp), intent(out) :: out(NX, NY)
    type(cuda_tdsops_t), intent(in) :: op
    call pack_y(field, buf_y, NX, NY, n_block)
    buf_y_dev = buf_y
    tmp_y_dev = 0.0_dp
    call exec_thom_tds_compact(tmp_y_dev, buf_y_dev, op, blocks, threads)
    tmp_y = tmp_y_dev
    call unpack_y(tmp_y, out, NX, NY, n_block)
  end subroutine

  subroutine pack_x(u_in, u_pack, nx_, ny_, nb)
    integer,  intent(in)  :: nx_, ny_, nb
    real(dp), intent(in)  :: u_in(nx_, ny_)
    real(dp), intent(out) :: u_pack(SZ, nx_, nb)
    integer :: i, j, lane, blk
    u_pack = 0.0_dp
    do j = 1, ny_
      blk  = (j - 1) / SZ + 1
      lane = mod(j - 1, SZ) + 1
      if (blk > nb) cycle
      do i = 1, nx_
        u_pack(lane, i, blk) = u_in(i, j)
      end do
    end do
    do blk = 1, nb
      do lane = 1, SZ
        j = (blk - 1) * SZ + lane
        if (j > ny_) u_pack(lane, :, blk) = 1.0_dp
      end do
    end do
  end subroutine

  subroutine pack_y(u_in, u_pack, nx_, ny_, nb)
    integer,  intent(in)  :: nx_, ny_, nb
    real(dp), intent(in)  :: u_in(nx_, ny_)
    real(dp), intent(out) :: u_pack(SZ, ny_, nb)
    integer :: i, j, lane, blk
    u_pack = 0.0_dp
    do i = 1, nx_
      blk  = (i - 1) / SZ + 1
      lane = mod(i - 1, SZ) + 1
      if (blk > nb) cycle
      do j = 1, ny_
        u_pack(lane, j, blk) = u_in(i, j)
      end do
    end do
    do blk = 1, nb
      do lane = 1, SZ
        i = (blk - 1) * SZ + lane
        if (i > nx_) u_pack(lane, :, blk) = 1.0_dp
      end do
    end do
  end subroutine

  subroutine unpack_x(tmp_pack, out, nx_, ny_, nb)
    integer,  intent(in)  :: nx_, ny_, nb
    real(dp), intent(in)  :: tmp_pack(SZ, nx_, nb)
    real(dp), intent(out) :: out(nx_, ny_)
    integer :: i, j, lane, blk
    do j = 1, ny_
      blk  = (j - 1) / SZ + 1
      lane = mod(j - 1, SZ) + 1
      if (blk > nb) cycle
      do i = 1, nx_
        out(i, j) = tmp_pack(lane, i, blk)
      end do
    end do
  end subroutine

  subroutine unpack_y(tmp_pack, out, nx_, ny_, nb)
    integer,  intent(in)  :: nx_, ny_, nb
    real(dp), intent(in)  :: tmp_pack(SZ, ny_, nb)
    real(dp), intent(out) :: out(nx_, ny_)
    integer :: i, j, lane, blk
    do i = 1, nx_
      blk  = (i - 1) / SZ + 1
      lane = mod(i - 1, SZ) + 1
      if (blk > nb) cycle
      do j = 1, ny_
        out(i, j) = tmp_pack(lane, j, blk)
      end do
    end do
  end subroutine

  subroutine report(iter, t, u_in, v_in, nx_, ny_)
    integer,  intent(in) :: iter, nx_, ny_
    real(dp), intent(in) :: t, u_in(nx_, ny_), v_in(nx_, ny_)
    real(dp) :: umax, umin, vmax, vmin
    umax = maxval(u_in); umin = minval(u_in)
    vmax = maxval(v_in); vmin = minval(v_in)
    write (*, '(A,I6,A,F9.4,A,ES11.3,A,ES11.3,A,ES11.3,A,ES11.3)') &
      ' iter=', iter, '  t=', t, &
      '  umin=', umin, ' umax=', umax, &
      ' vmin=', vmin, ' vmax=', vmax
  end subroutine report

end program test_2d_uv_skew