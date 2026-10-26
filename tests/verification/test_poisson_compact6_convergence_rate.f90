program test_poisson_compact6_convergence_rate
  !! Compact6 derivative convergence rate test for all BC types.
  !!
  !! All test scenarios are defined in the scenario table near the top of
  !! the main program.  Each scenario specifies:
  !!   - label           : human-readable name for output
  !!   - func_type       : FUNC_COS or FUNC_SIN
  !!   - wavenumber      : k in cos(k*x) or sin(k*x)
  !!   - test_kind       : TEST_DERIV, TEST_STAGDER_V2P, TEST_STAGDER_P2V,
  !!                        TEST_INTERPL_V2P, or TEST_INTERPL_P2V
  !!   - sym             : Neumann symmetry flag for the input field
  !!   - skip_periodic   : skip for periodic BCs
  !!   - skip_dirichlet  : skip for Dirichlet BCs
  !!
  !! Constraints on test functions:
  !!   Periodic BCs:  wavenumber must be a multiple of 2*pi (function
  !!                  must be periodic on [0,1)).
  !!   Neumann BCs:   stagder_v2p (sym=.false.) => odd input (sin)
  !!                  stagder_p2v (sym=.true.)  => even input (cos)
  !!                  interpl_v2p (sym=.true.)  => even input (cos)
  !!                  interpl_p2v (sym=.true.)  => even input (cos)
  !!   Dirichlet BCs: staggered/interpl not supported.
  !!
  !! Staggered operator sym conventions:
  !!   stagder  v2p: sym=.false.    p2v: sym=.true.
  !!   interpl  v2p: sym=.true.     p2v: sym=.true.
  !!
  !! Expected convergence rates when doubling resolution:
  !!   Periodic:   6
  !!   Neumann:    6
  !!   Dirichlet:  5
  !!
  !! Pattern follows test_poisson_bc: using the context to allow for backend-agnostic behaviour.

  use iso_fortran_env

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_backend_runtime, only: backend_runtime_t, backend_is_cuda, tdsops_slot_t
  use m_tdsops_setup, only: allocate_tdsops
  use m_tdsops, only: dirps_t
  use m_test_utils, only: initialise_mpi, finalise_test
  use m_mesh, only: mesh_t
  use m_vector_calculus, only: vector_calculus_t

  use m_common, only: dp, pi, BC_PERIODIC, BC_NEUMANN, BC_DIRICHLET, dp, pi, MPI_X3D2_DP, DIR_C, DIR_X, DIR_Y, DIR_Z, &
                      CELL, RDR_C2Z, RDR_C2X, RDR_Z2X

  implicit none

  ! ------------------------------------------------------------------
  ! Test scenario kind constants
  ! ------------------------------------------------------------------
  integer, parameter :: TEST_DERIV       = 1  ! collocated 1st + 2nd derivative
  integer, parameter :: TEST_STAGDER_V2P = 2  ! staggered derivative v2p only
  integer, parameter :: TEST_STAGDER_P2V = 3  ! staggered derivative p2v only
  integer, parameter :: TEST_INTERPL_V2P = 4  ! staggered interpolation v2p only
  integer, parameter :: TEST_INTERPL_P2V = 5  ! staggered interpolation p2v only

  ! Function type constants
  integer, parameter :: FUNC_COS = 1
  integer, parameter :: FUNC_SIN = 2

  ! ------------------------------------------------------------------
  ! Scenario descriptor
  ! ------------------------------------------------------------------
  type :: test_scenario_t
    character(len=30) :: label          ! e.g. 'cos(3*pi*x) deriv'
    integer           :: func_type      ! FUNC_COS or FUNC_SIN
    real(dp)          :: wavenumber      ! k in cos(k*x) or sin(k*x)
    integer           :: test_kind      ! TEST_DERIV, TEST_STAGDER_*, TEST_INTERPL_*
    logical           :: sym            ! Neumann sym flag for input field
    logical           :: skip_periodic  ! skip this scenario for periodic BCs
    logical           :: skip_dirichlet ! skip this scenario for Dirichlet BCs
  end type test_scenario_t

  ! ------------------------------------------------------------------
  ! The context of the test. This contains the backend runtime and associated allocator for the tests. This also contains the given TDSOps slot.
  ! ------------------------------------------------------------------
  type :: test_context_t
    type(backend_runtime_t) :: runtime
    class(base_backend_t), pointer :: backend => null()
    type(allocator_t), pointer :: host_allocator => null()
    type(mesh_t) :: mesh
    type(dirps_t), pointer :: xdirps => null(), ydirps => null(), zdirps => null()
    type(vector_calculus_t) :: vector_calculus
    type(tdsops_slot_t) :: tdsops_slot
  end type test_context_t

  ! ------------------------------------------------------------------
  ! *** SCENARIO TABLE -- add/remove/reorder test cases here ***
  !
  ! Periodicity rule: k must be m*2*pi for integer m.
  !   3*pi is NOT periodic on [0,1) => skip_periodic=.true.
  !   2*pi, 4*pi ARE periodic       => skip_periodic=.false.
  !
  ! Parity rules for Neumann BCs:
  !   stagder_v2p  (sym=.false.) => use odd  input (sin)
  !   stagder_p2v  (sym=.true.)  => use even input (cos)
  !   interpl_v2p  (sym=.true.)  => use even input (cos)
  !   interpl_p2v  (sym=.true.)  => use even input (cos)
  ! ------------------------------------------------------------------
  integer, parameter :: n_scenarios = 4

    ! ------------------------------------------------------------------
  ! Define all test scenarios (single place to add new tests)
  ! ------------------------------------------------------------------

  type(test_scenario_t), parameter :: scenarios(*) = [&
    test_scenario_t( &
      label='sin(2*pi*x) stagder_v2p', func_type=FUNC_SIN, &
      wavenumber=2._dp*pi, test_kind=TEST_STAGDER_V2P, &
      sym=.false., skip_periodic=.true., skip_dirichlet=.true.), &
    ! test_scenario_t( &
    !   label='sin(3*pi*x) stagder_v2p', func_type=FUNC_SIN, &
    !   wavenumber=3._dp*pi, test_kind=TEST_STAGDER_V2P, &
    !   sym=.false., skip_periodic=.true., skip_dirichlet=.true.), &
    ! test_scenario_t( &
    !   label='cos(10*pi*x) stagder_p2v', func_type=FUNC_COS, &
    !   wavenumber=10._dp*pi, test_kind=TEST_STAGDER_P2V, &
    !   sym=.true., skip_periodic=.true., skip_dirichlet=.true.), &
    test_scenario_t( &
      label='cos(3*pi*x) stagder_p2v', func_type=FUNC_COS, &
      wavenumber=3._dp*pi, test_kind=TEST_STAGDER_P2V, &
      sym=.true., skip_periodic=.true., skip_dirichlet=.true.) &
      ]

  ! BC configuration
  integer, parameter :: n_bcs = 3
  integer :: bc_types(n_bcs)
  real(dp) :: expected_rates(n_bcs)
  character(len=10) :: bc_names(n_bcs)

  ! Grid configuration
  integer, parameter :: n_grids = 6
  integer :: grid_sizes(n_grids)
  integer, parameter :: n_block = 64

  ! Loop indices
  integer :: ibc, iscen, ierr, ndevs, devnum
  real(dp) :: expected_rate

  ! Convergence error arrays
  real(dp) :: errors(n_grids)
  real(dp) :: err_a(n_grids), err_b(n_grids)

  logical :: allpass = .true.

  ! ------------------------------------------------------------------
  ! Define BC types
  ! ------------------------------------------------------------------
  bc_types = [BC_PERIODIC, BC_DIRICHLET, BC_NEUMANN]
  bc_names = ['periodic  ', 'dirichlet ', 'neumann   ']
  expected_rates = [6._dp, 5._dp, 6._dp]

  ! ------------------------------------------------------------------
  ! Main loop: BCs x Scenarios
  ! ------------------------------------------------------------------
  do ibc = 1, n_bcs
    expected_rate = expected_rates(ibc)

    print *, ''
    print *, '=== Testing ', trim(bc_names(ibc)), &
             ' (expected rate:', expected_rate, ') ==='

    ! Grid sizes
    if (bc_types(ibc) == BC_PERIODIC) then
      grid_sizes = [16, 32, 64, 128, 256, 512]
    else
      grid_sizes = [17, 33, 65, 129, 257, 513]
    end if

    do iscen = 1, n_scenarios
      ! Skip scenarios incompatible with this BC type
      if (scenarios(iscen)%skip_periodic .and. &
          bc_types(ibc) == BC_PERIODIC) cycle
      if (scenarios(iscen)%skip_dirichlet .and. &
          bc_types(ibc) == BC_DIRICHLET) cycle

      print *, ''
      print *, '--- ', trim(scenarios(iscen)%label), ' ---'

      ! Run the convergence study for this scenario
      select case (scenarios(iscen)%test_kind)
      case (TEST_DERIV)
        call run_deriv_convergence( &
          scenarios(iscen), bc_types(ibc), &
          grid_sizes, n_grids, n_block, err_a, err_b)
        call print_convergence_pair( &
          scenarios(iscen)%label, grid_sizes, &
          err_a, err_b, expected_rate, allpass, &
          col_a='der1st', col_b='der2nd')

      case (TEST_STAGDER_V2P)
        call run_stag_v2p_convergence( &
          scenarios(iscen), bc_types(ibc), &
          grid_sizes, n_grids, n_block, &
          'stag-deriv', 'compact6', .false., errors)
        call print_convergence_single( &
          scenarios(iscen)%label, grid_sizes, &
          errors, expected_rate, allpass)

      case (TEST_STAGDER_P2V)
        call run_stag_p2v_convergence( &
          scenarios(iscen), bc_types(ibc), &
          grid_sizes, n_grids, n_block, &
          'stag-deriv', 'compact6', .true., errors)
        call print_convergence_single( &
          scenarios(iscen)%label, grid_sizes, &
          errors, expected_rate, allpass)

      case (TEST_INTERPL_V2P)
        call run_stag_v2p_convergence( &
          scenarios(iscen), bc_types(ibc), &
          grid_sizes, n_grids, n_block, &
          'interpolate', 'classic', .true., errors)
        call print_convergence_single( &
          scenarios(iscen)%label, grid_sizes, &
          errors, expected_rate, allpass)

      case (TEST_INTERPL_P2V)
        call run_stag_p2v_convergence( &
          scenarios(iscen), bc_types(ibc), &
          grid_sizes, n_grids, n_block, &
          'interpolate', 'classic', .true., errors)
        call print_convergence_single( &
          scenarios(iscen)%label, grid_sizes, &
          errors, expected_rate, allpass)
      end select
    end do ! iscen

    ! Write convergence data to file
    call write_convergence_file(bc_names(ibc), bc_types(ibc), &
                                scenarios, n_scenarios, &
                                grid_sizes, n_grids, n_block)
  end do ! ibc

  ! Final verdict
  if (allpass) then
    print *, 'ALL TESTS PASSED SUCCESSFULLY.'
  else
    error stop 'SOME TESTS FAILED.'
  end if

contains

  ! ================================================================
  ! Set up the context of the case, and the associated runtime and allocators.
  ! ================================================================

  type(test_context_t) function setup_ctx(scenario, grid_size, bc_type)
    type(test_scenario_t), intent(in) :: scenario
    integer, intent(in) :: grid_size, bc_type
    type(test_context_t), intent(out) :: ctx
    type(backend_runtime_t), target :: runtime
    logical :: use_2decomp

    use_2decomp = .not. backend_is_cuda
    bc_name = bc_names(bc_type)

    ctx%mesh = mesh_t(grid_size, nproc_dir, 1.0_dp, bc_names, use_2decomp=use_2decomp)
    ctx%runtime => runtime
    call ctx%runtime%init(ctx%mesh)
    ctx%backend => ctx%runtime%backend
    ctx%host_allocator = ctx%runtime%host_allocator

  end function setup_ctx
    

  ! ================================================================
  ! Evaluate test function f(x), f'(x), f''(x)
  ! ================================================================
  pure real(dp) function eval_f(func_type, k, x)
    integer, intent(in) :: func_type
    real(dp), intent(in) :: k, x

    select case (func_type)
    case (FUNC_COS); eval_f = cos(k * x)
    case (FUNC_SIN); eval_f = sin(k * x)
    case default;    eval_f = 0._dp
    end select

  end function eval_f

  pure real(dp) function eval_df(func_type, k, x)
    integer, intent(in) :: func_type
    real(dp), intent(in) :: k, x

    select case (func_type)
    case (FUNC_COS); eval_df = -k * sin(k * x)
    case (FUNC_SIN); eval_df = k * cos(k * x)
    case default;     eval_df = 0._dp
    end select
  end function eval_df

  pure real(dp) function eval_d2f(func_type, k, x)
    integer, intent(in) :: func_type
    real(dp), intent(in) :: k, x

    select case (func_type)
    case (FUNC_COS); eval_d2f = -k * k * cos(k * x)
    case (FUNC_SIN); eval_d2f = -k * k * sin(k * x)
    case default;     eval_d2f = 0._dp
    end select
  end function eval_d2f

  ! ================================================================
  ! Compute dx for a given grid size and BC type
  ! ================================================================
  pure real(dp) function compute_dx(n, bc_type)
    integer, intent(in) :: n, bc_type
    if (bc_type == BC_PERIODIC) then
      compute_dx = 1._dp / n
    else
      compute_dx = 1._dp / (n - 1)
    end if
  end function compute_dx

  ! ================================================================
  ! Compute staggered input/output sizes for a given direction
  ! ================================================================
  pure subroutine get_stag_sizes(bc_type, n, direction, n_in, n_out)
    integer, intent(in) :: bc_type, n
    character(*), intent(in) :: direction
    integer, intent(out) :: n_in, n_out

    if (bc_type == BC_PERIODIC) then
      n_in  = n
      n_out = n
    else
      if (direction == 'v2p') then
        n_in  = n       ! v-grid: n vertices
        n_out = n - 1   ! p-grid: n-1 cells
      else  ! p2v
        n_in  = n - 1   ! p-grid: n-1 cells
        n_out = n       ! v-grid: n vertices
      end if
    end if
  end subroutine get_stag_sizes

  ! ================================================================
  ! Collocated derivative convergence study (der1st + der2nd)
  ! ================================================================
  subroutine run_deriv_convergence(scen, bc_type, grid_sizes, n_grids, &
                                   n_block, err_1st, err_2nd)
    type(test_scenario_t), intent(in) :: scen
    integer, intent(in) :: bc_type, grid_sizes(:), n_grids, n_block
    real(dp), intent(out) :: err_1st(:), err_2nd(:)

    integer :: ig, n, j
    real(dp) :: dx, x, k
    real(dp), allocatable :: u(:,:,:), du(:,:,:), ref(:,:,:)
    class(field_t), pointer :: u_field, du_field
    type(test_context_t) :: ctx, ctx2

    k = scen%wavenumber

    do ig = 1, n_grids
      n = grid_sizes(ig)
      dx = compute_dx(n, bc_type)

      ctx = setup_ctx(scen, n, bc_type)

      allocate(u(SZ, n, n_block), du(SZ, n, n_block), ref(SZ, n, n_block))

      ! Fill input field
      do j = 1, n
        x = (j - 1) * dx
        u(:, j, :) = eval_f(scen%func_type, k, x)
      end do

      ! --- First derivative ---
      ref = eval_df(scen%func_type, k, x)

      if (bc_type == BC_NEUMANN) then
        call ctx%backend%allocate_tdsops(ctx%xdirps, ctx%backend, ctx%mesh operation='first-deriv', &
                              scheme='compact6', &
                              bc_start=bc_type, bc_end=bc_type, &
                              sym=scen%sym)
      else
        call ctx%backend%allocate_tdsops(n, dx, operation='first-deriv', &
                              scheme='compact6', &
                              bc_start=bc_type, bc_end=bc_type)
      end if

      call ctx%backend%set_field_data(u_field, u, DIR_X)
      call ctx%backend%tds_solve(du_field, u_field, ctx%tdsops_slot%op)
      call ctx%backend%get_field_data(du, du_field, DIR_X)

      err_1st(ig) = norm2(du - ref) / sqrt(real(SZ * n * n_block, dp))

      ! --- Second derivative ---
      do j = 1, n
        x = (j - 1) * dx
        ref(:, j, :) = eval_d2f(scen%func_type, k, x)
      end do

      ctx2 = setup_ctx(scen, n, bc_type)


      if (bc_type == BC_NEUMANN) then
        call ctx2%backend%allocate_tdsops(n, dx, operation='second-deriv', &
                              scheme='compact6', &
                              bc_start=bc_type, bc_end=bc_type, &
                              sym=scen%sym)
      else
        call ctx2%backend%allocate_tdsops(n, dx, operation='second-deriv', &
                              scheme='compact6', &
                              bc_start=bc_type, bc_end=bc_type)
      end if

      call ctx2%backend%set_field_data(u_field, u, DIR_X)
      call ctx2%backend%tds_solve(du_field, u_field, ctx%tdsops_slot%op)
      call ctx2%backend%get_field_data(du, du_field, DIR_X)

      err_2nd(ig) = norm2(du - ref) / sqrt(real(SZ * n * n_block, dp))

      ! Debug output for smallest Neumann grid
      if (bc_type == BC_NEUMANN .and. ig == 1) then
        print *, '  [DEBUG] N =', n
        print *, '  du(1, 1:5, 1)  =', du(1, 1:5, 1)
        print *, '  ref(1, 1:5, 1) =', ref(1, 1:5, 1)
        print *, '  du(1, n-4:n, 1)  =', du(1, n-4:n, 1)
        print *, '  ref(1, n-4:n, 1) =', ref(1, n-4:n, 1)
        print *, '  err at j=1:', abs(du(1,1,1) - ref(1,1,1))
        print *, '  err at j=n:', abs(du(1,n,1) - ref(1,n,1))
      end if

      deallocate(u, du, ref)
    end do
  end subroutine run_deriv_convergence

  ! ================================================================
  ! Staggered v2p convergence study
  ! ================================================================
subroutine run_stag_v2p_convergence(scen, bc_type, grid_sizes, &
                                      n_grids, n_block, operation, &
                                      scheme, sym_v2p, err_v2p)
    type(test_scenario_t), intent(in) :: scen
    integer, intent(in) :: bc_type, grid_sizes(:), n_grids, n_block
    character(*), intent(in) :: operation, scheme
    logical, intent(in) :: sym_v2p
    real(dp), intent(out) :: err_v2p(:)

    integer :: ig, n, n_in, n_out
    real(dp) :: k
    logical :: is_interpl
    class(field_t), pointer :: u_stag_field, du_stag_field
    type(test_context_t) :: ctx, ctx2

    k = scen%wavenumber
    is_interpl = (operation == 'interpolate')

    do ig = 1, n_grids
      n = grid_sizes(ig)

      ctx = setup_ctx(scen, n, bc_type)

      call get_stag_sizes(bc_type, n, 'v2p', n_in, n_out)

      block
        integer :: j, n_buf
        real(dp) :: dx, x
        real(dp), allocatable :: u_stag(:,:,:), du_stag(:,:,:), ref(:,:,:)

        dx = compute_dx(n, bc_type)
        n_buf = max(n_in, n_out)

        ! Allocate host AND device buffers at the common size n_buf so the
        ! kernel's out-of-bounds staggered reads/writes land in valid memory
        ! and the D<->H copies use matching descriptors.
        allocate(u_stag(SZ, n_buf, n_block), du_stag(SZ, n_buf, n_block), &
                 ref(SZ, n_buf, n_block))

        u_stag = 0._dp; du_stag = 0._dp; ref = 0._dp

        do j = 1, n_in
          x = (j - 1) * dx
          u_stag(:, j, :) = eval_f(scen%func_type, k, x)
        end do
        do j = 1, n_out
          x = (real(j, dp) - 0.5_dp) * dx
          if (is_interpl) then
            ref(:, j, :) = eval_f(scen%func_type, k, x)
          else
            ref(:, j, :) = eval_df(scen%func_type, k, x)
          end if
        end do

        call ctx%backend%alloc_tdsops(n_out, dx, operation=operation, &
                              scheme=scheme, &
                              bc_start=bc_type, bc_end=bc_type, &
                              from_to='v2p', sym=sym_v2p)

        call ctx%backend%set_field_data(u_stag_field, u_stag, DIR_X)
        call ctx%backend%tds_solve(du_stag_field, u_stag_field, ctx%tdsops_slot%op)
        call ctx%backend%get_field_data(du_stag, du_stag_field, DIR_X)

        ! Only the first n_out slices carry valid output.
        err_v2p(ig) = norm2(du_stag(:, 1:n_out, :) - ref(:, 1:n_out, :)) / &
                      sqrt(real(SZ * n_out * n_block, dp))

        if (bc_type == BC_NEUMANN .and. ig == 1) then
          call print_stag_debug(trim(scen%label), n, n_out, &
                                du_stag(:, 1:n_out, :), &
                                ref(:, 1:n_out, :), n_block)
        end if

        deallocate(u_stag, du_stag, ref)
      end block
    end do
  end subroutine run_stag_v2p_convergence

  ! ================================================================
  ! Staggered p2v convergence study
  ! ================================================================
 subroutine run_stag_p2v_convergence(scen, bc_type, grid_sizes, &
                                      n_grids, n_block, operation, &
                                      scheme, sym_p2v, err_p2v)
    type(test_scenario_t), intent(in) :: scen
    integer, intent(in) :: bc_type, grid_sizes(:), n_grids, n_block
    character(*), intent(in) :: operation, scheme
    logical, intent(in) :: sym_p2v
    real(dp), intent(out) :: err_p2v(:)

    integer :: ig, n, n_in, n_out
    real(dp) :: k
    logical :: is_interpl
    type(test_context_t) :: ctx
    class(field_t), pointer :: u_stag_field, du_stag_field

    k = scen%wavenumber
    is_interpl = (operation == 'interpolate')

    do ig = 1, n_grids
      n = grid_sizes(ig)

      ctx = setup_ctx(scen, n, bc_type)

      call get_stag_sizes(bc_type, n, 'p2v', n_in, n_out)

      block
        integer :: j, n_buf
        real(dp) :: dx, x
        real(dp), allocatable :: u_stag(:,:,:), du_stag(:,:,:), ref(:,:,:)

        dx = compute_dx(n, bc_type)
        n_buf = max(n_in, n_out)

        allocate(u_stag(SZ, n_buf, n_block), du_stag(SZ, n_buf, n_block), &
                 ref(SZ, n_buf, n_block))

        u_stag = 0._dp; du_stag = 0._dp; ref = 0._dp

        do j = 1, n_in
          x = (real(j, dp) - 0.5_dp) * dx
          u_stag(:, j, :) = eval_f(scen%func_type, k, x)
        end do
        do j = 1, n_out
          x = (j - 1) * dx
          if (is_interpl) then
            ref(:, j, :) = eval_f(scen%func_type, k, x)
          else
            ref(:, j, :) = eval_df(scen%func_type, k, x)
          end if
        end do

        call ctx%backend%alloc_tdsops(n_out, dx, operation=operation, &
                              scheme=scheme, &
                              bc_start=bc_type, bc_end=bc_type, &
                              from_to='p2v', sym=sym_p2v)

        call ctx%backend%set_field_data(u_stag_field, u_stag, DIR_X)
        call ctx%backend%tds_solve(du_stag_field, u_stag_field, ctx%tdsops_slot%op)
        call ctx%backend%get_field_data(du_stag, du_stag_field, DIR_X)


        err_p2v(ig) = norm2(du_stag(:, 1:n_out, :) - ref(:, 1:n_out, :)) / &
                      sqrt(real(SZ * n_out * n_block, dp))

        if (bc_type == BC_NEUMANN .and. ig == 1) then
          call print_stag_debug(trim(scen%label), n, n_out, &
                                du_stag(:, 1:n_out, :), &
                                ref(:, 1:n_out, :), n_block)
        end if

        deallocate(u_stag, du_stag, ref)
      end block
    end do
  end subroutine run_stag_p2v_convergence

  ! ================================================================
  ! Debug output helper for staggered tests
  ! ================================================================
  subroutine print_stag_debug(label, n, n_out, du_stag, ref, n_block)
    character(*), intent(in) :: label
    integer, intent(in) :: n, n_out, n_block
    real(dp), intent(in) :: du_stag(:,:,:), ref(:,:,:)

    integer :: j

    print *, '  [DEBUG ', trim(label), '] N =', n, ' n_out =', n_out
    do j = 1, n_out
      print '(A,I4,A,ES22.14,A,ES22.14,A,ES12.4)', &
        '  j=', j, '  du=', du_stag(1,j,1), &
        '  ref=', ref(1,j,1), &
        '  err=', abs(du_stag(1,j,1) - ref(1,j,1))
    end do
    print *, '  Per-lane L2 errors (block=1):'
    do j = 1, SZ
      print '(A,I4,A,ES12.4)', '  lane=', j, '  L2=', &
        norm2(du_stag(j,:,1) - ref(j,:,1)) / sqrt(real(n_out, dp))
    end do
  end subroutine print_stag_debug

  ! ! ================================================================
  ! ! Convergence table printer for paired tests (der1st + der2nd)
  ! ! ================================================================
  subroutine print_convergence_pair(label, grid_sizes, err_a, err_b, &
  !                                   expected_rate, allpass, col_a, col_b)
    character(*), intent(in) :: label, col_a, col_b
    integer, intent(in) :: grid_sizes(:)
    real(dp), intent(in) :: err_a(:), err_b(:)
    real(dp), intent(in) :: expected_rate
    logical, intent(inout) :: allpass

    integer :: ig, ng
    real(dp) :: rate_a, rate_b
    character(len=80) :: header

    ng = size(grid_sizes)

    print *, ''
    write(header, '(4X,A6,2A18,2A12)') &
      'N', col_a//' L2', col_b//' L2', 'rate_a', 'rate_b'
    print '(A,A,A)', ' Convergence study ', trim(label), ' (L2):'
    print '(A)', trim(header)

    do ig = 1, ng
      if (ig == 1) then
        print '(4X,I6,2ES18.6)', grid_sizes(ig), err_a(ig), err_b(ig)
      else
        rate_a = log(err_a(ig-1)/err_a(ig)) / log(2._dp)
        rate_b = log(err_b(ig-1)/err_b(ig)) / log(2._dp)
        print '(4X,I6,2ES18.6,2F12.2)', &
          grid_sizes(ig), err_a(ig), err_b(ig), rate_a, rate_b
      end if
    end do

    rate_a = log(err_a(ng-1)/err_a(ng)) / log(2._dp)
    rate_b = log(err_b(ng-1)/err_b(ng)) / log(2._dp)

    print *, 'rate_'//trim(col_a), rate_a
    if (abs(rate_a - expected_rate) > 0.25_dp) then
      allpass = .false.
      print *, 'Check '//trim(label)//' '//trim(col_a)// &
               ' convergence rate... failed'
    else
      print *, 'Check '//trim(label)//' '//trim(col_a)// &
               ' convergence rate... passed'
    end if

    print *, 'rate_'//trim(col_b), rate_b
    if (abs(rate_b - expected_rate) > 0.25_dp) then
      allpass = .false.
      print *, 'Check '//trim(label)//' '//trim(col_b)// &
               ' convergence rate... failed'
    else
      print *, 'Check '//trim(label)//' '//trim(col_b)// &
               ' convergence rate... passed'
    end if
  end subroutine print_convergence_pair

  ! ================================================================
  ! Convergence table printer for single-direction tests
  ! ================================================================
  subroutine print_convergence_single(label, grid_sizes, errors, &
                                    expected_rate, allpass)
  character(*), intent(in) :: label
  integer, intent(in) :: grid_sizes(:)
  real(dp), intent(in) :: errors(:)
  real(dp), intent(in) :: expected_rate
  logical, intent(inout) :: allpass

  integer :: ig, ng
  real(dp) :: rate, best_rate
  real(dp), parameter :: roundoff_floor = 1.0e-12_dp

  ng = size(grid_sizes)

  print *, ''
  print '(A,A,A)', ' Convergence study ', trim(label), ' (L2):'
  print '(4X,A6,A18,A12)', 'N', 'L2 error', 'rate'

  best_rate = -huge(1._dp)
  do ig = 1, ng
    if (ig == 1) then
      print '(4X,I6,ES18.6)', grid_sizes(ig), errors(ig)
    else
      rate = log(errors(ig-1)/errors(ig)) / log(2._dp)
      print '(4X,I6,ES18.6,F12.2)', grid_sizes(ig), errors(ig), rate
      ! Only consider rates computed above the round-off floor.
      if (errors(ig) > roundoff_floor) best_rate = max(best_rate, rate)
    end if
  end do

  print *, 'best rate (above round-off) ', best_rate
  if (abs(best_rate - expected_rate) > 0.25_dp) then
    allpass = .false.
    print *, 'Check '//trim(label)//' convergence rate... failed'
  else
    print *, 'Check '//trim(label)//' convergence rate... passed'
  end if
end subroutine print_convergence_single

  ! ================================================================
  ! Write convergence data file for one BC type
  ! ================================================================
  subroutine write_convergence_file(bc_name, bc_type, scenarios, &
                                    n_scenarios, grid_sizes, n_grids, &
                                    n_block)
    character(*), intent(in) :: bc_name
    integer, intent(in) :: bc_type, n_scenarios, n_grids, n_block
    type(test_scenario_t), intent(in) :: scenarios(:)
    integer, intent(in) :: grid_sizes(:)

    integer :: ig, iscen, col
    real(dp) :: dx, ea, eb
    real(dp) :: row(1 + 2*n_scenarios)

    open(unit=20, file='convergence_'//trim(bc_name)//'.dat', &
         status='replace')

    do ig = 1, n_grids
      dx = compute_dx(grid_sizes(ig), bc_type)
      row(1) = dx
      col = 2

      do iscen = 1, n_scenarios
        if (scenarios(iscen)%skip_periodic .and. &
            bc_type == BC_PERIODIC) cycle
        if (scenarios(iscen)%skip_dirichlet .and. &
            bc_type == BC_DIRICHLET) cycle

        call run_scenario_single_grid( &
          scenarios(iscen), bc_type, grid_sizes(ig), n_block, ea, eb)

        if (scenarios(iscen)%test_kind == TEST_DERIV) then
          row(col)   = ea
          row(col+1) = eb
          col = col + 2
        else
          row(col) = ea
          col = col + 1
        end if
      end do

      write(20, '(999ES20.10)') row(1:col-1)
    end do

    close(20)
  end subroutine write_convergence_file

  ! ================================================================
  ! Run a single scenario at a single grid size
  ! ================================================================
  subroutine run_scenario_single_grid(scen, bc_type, n, n_block, ea, eb)
    type(test_scenario_t), intent(in) :: scen
    integer, intent(in) :: bc_type, n, n_block
    real(dp), intent(out) :: ea, eb

    eb = 0._dp

    select case (scen%test_kind)
    case (TEST_DERIV)
      call run_deriv_single(scen, bc_type, n, n_block, ea, eb)
    case (TEST_STAGDER_V2P)
      call run_stag_v2p_single(scen, bc_type, n, n_block, &
                               'stag-deriv', 'compact6', .false., ea)
    case (TEST_STAGDER_P2V)
      call run_stag_p2v_single(scen, bc_type, n, n_block, &
                               'stag-deriv', 'compact6', .true., ea)
    case (TEST_INTERPL_V2P)
      call run_stag_v2p_single(scen, bc_type, n, n_block, &
                               'interpolate', 'classic', .true., ea)
    case (TEST_INTERPL_P2V)
      call run_stag_p2v_single(scen, bc_type, n, n_block, &
                               'interpolate', 'classic', .true., ea)
    end select
  end subroutine run_scenario_single_grid

  ! ================================================================
  ! Single-grid collocated derivative (for file output)
  ! ================================================================
  subroutine run_deriv_single(scen, bc_type, n, n_block, err_1st, err_2nd)
    type(test_scenario_t), intent(in) :: scen
    integer, intent(in) :: bc_type, n, n_block
    real(dp), intent(out) :: err_1st, err_2nd

    integer :: j
    real(dp) :: dx, x, k
    real(dp), allocatable :: u(:,:,:), du(:,:,:), ref(:,:,:)
    class(field_t), pointer :: u_field, du_field
    type(test_context_t) :: ctx, ctx2

    k = scen%wavenumber
    dx = compute_dx(n, bc_type)

    ctx = setup_ctx(scen, n, bc_type)

    allocate(u(SZ, n, n_block), du(SZ, n, n_block), ref(SZ, n, n_block))

    do j = 1, n
      x = (j - 1) * dx
      u(:, j, :) = eval_f(scen%func_type, k, x)
    end do

    do j = 1, n
      x = (j - 1) * dx
      ref(:, j, :) = eval_df(scen%func_type, k, x)
    end do
    if (bc_type == BC_NEUMANN) then
      call ctx%backend%alloc_tdsops(n, dx, operation='first-deriv', &
                            scheme='compact6', &
                            bc_start=bc_type, bc_end=bc_type, &
                            sym=scen%sym)
    else
      call ctx%backend%alloc_tdsops(n, dx, operation='first-deriv', &
                            scheme='compact6', &
                            bc_start=bc_type, bc_end=bc_type)
    end if

    call ctx%backend%set_field_data(u_field, u, DIR_X)
    call ctx%backend%tds_solve(du_field, u_field, ctx%tdsops_slot%op)
    call ctx%backend%get_field_data(du, du_field, DIR_X)

    err_1st = norm2(du - ref) / sqrt(real(SZ * n * n_block, dp))

    ctx2 = setup_ctx(scen, n, bc_type)

    do j = 1, n
      x = (j - 1) * dx
      ref(:, j, :) = eval_d2f(scen%func_type, k, x)
    end do
    if (bc_type == BC_NEUMANN) then
      call ctx2%backend%alloc_tdsops(n, dx, operation='second-deriv', &
                            scheme='compact6', &
                            bc_start=bc_type, bc_end=bc_type, &
                            sym=scen%sym)
    else
      call ctx2%backend%alloc_tdsops(n, dx, operation='second-deriv', &
                            scheme='compact6', &
                            bc_start=bc_type, bc_end=bc_type)
    end if

    call ctx2%backend%set_field_data(u_field, u, DIR_X)
    call ctx2%backend%tds_solve(du_field, u_field, ctx%tdsops_slot%op)
    call ctx2%backend%get_field_data(du, du_field, DIR_X)

    err_2nd = norm2(du - ref) / sqrt(real(SZ * n * n_block, dp))

    deallocate(u, du, ref)
  end subroutine run_deriv_single

  ! ================================================================
  ! Single-grid staggered v2p (for file output)
  ! ================================================================
  subroutine run_stag_v2p_single(scen, bc_type, n, n_block, &
                                  operation, scheme, sym_v2p, err)
    type(test_scenario_t), intent(in) :: scen
    integer, intent(in) :: bc_type, n, n_block
    character(*), intent(in) :: operation, scheme
    logical, intent(in) :: sym_v2p
    real(dp), intent(out) :: err

    integer :: j, n_in, n_out
    real(dp) :: dx, x, k
    logical :: is_interpl
    real(dp), allocatable :: u_stag(:,:,:), du_stag(:,:,:), ref(:,:,:)
    type(test_context_t), pointer :: ctx

    k = scen%wavenumber
    dx = compute_dx(n, bc_type)
    is_interpl = (operation == 'interpolate')

    ctx = setup_ctx(scen, n, bc_type)

    call get_stag_sizes(bc_type, n, 'v2p', n_in, n_out)

    allocate(u_stag(SZ, n_in, n_block), du_stag(SZ, n_out, n_block), &
             ref(SZ, n_out, n_block))

    do j = 1, n_in
      x = (j - 1) * dx
      u_stag(:, j, :) = eval_f(scen%func_type, k, x)
    end do
    do j = 1, n_out
      x = (real(j, dp) - 0.5_dp) * dx
      if (is_interpl) then
        ref(:, j, :) = eval_f(scen%func_type, k, x)
      else
        ref(:, j, :) = eval_df(scen%func_type, k, x)
      end if
    end do

    call cxt%backend%alloc_tdsops(n_out, dx, operation=operation, scheme=scheme, &
                          bc_start=bc_type, bc_end=bc_type, &
                          from_to='v2p', sym=sym_v2p)

    call ctx%backend%set_field_data(u_stag_fi    err = norm2(du_stag - ref) / sqrt(real(SZ * n_out * n_block, dp))eld, DIR_X)
    call ctx%backend%tds_solve(du_stag_field, u_stag_field, ctx%tdsops_slot%op)
    call ctx%backend%get_field_data(du_stag, du_stag_field, DIR_X)
    
    err = norm2(du_stag - ref) / sqrt(real(SZ * n_out * n_block, dp))

    deallocate(u_stag, du_stag, ref)
  end subroutine run_stag_v2p_single

  ! =============================================    err = norm2(du_stag - ref) / sqrt(real(SZ * n_out * n_block, dp))===================
  ! Single-grid staggered p2v (for file output)
  ! ================================================================
  subroutine run_stag_p2v_single(scen, bc_type, n, n_block, &
                                  operation, scheme, sym_p2v, err)
    type(test_scenario_t), intent(in) :: scen
    integer, intent(in) :: bc_type, n, n_block
    character(*), intent(in) :: operation, scheme
    logical, intent(in) :: sym_p2v
    real(dp), intent(out) :: err

    integer :: j, n_in, n_out
    real(dp) :: dx, x, k
    logical :: is_interpl
    real(dp), allocatable :: u_stag(:,:,:), du_stag(:,:,:), ref(:,:,:)
    class(field_t), pointer :: u_stag_field, du_stag_field
    type(test_context_t), pointer :: ctx

    k = scen%wavenumber
    dx = compute_dx(n, bc_type)
    is_interpl = (operation == 'interpolate')

    ctx = setup_ctx(scen, n, bc_type)

    call get_stag_sizes(bc_type, n, 'p2v', n_in, n_out)

    allocate(u_stag(SZ, n_in, n_block), du_stag(SZ, n_out, n_block), &
             ref(SZ, n_out, n_block))

    do j = 1, n_in
      x = (real(j, dp) - 0.5_dp) * dx
      u_stag(:, j, :) = eval_f(scen%func_type, k, x)
    end do
    do j = 1, n_out
      x = (j - 1) * dx
      if (is_interpl) then
        ref(:, j, :) = eval_f(scen%func_type, k, x)
      else
        ref(:, j, :) = eval_df(scen%func_type, k, x)
      end if
    end do

    call ctx%backend%alloc_tdsops(n, dx, operation=operation, scheme=scheme, bc_start=bc_type, bc_end=bc_type, &
     from_to='p2v', sym=sym_p2v)

    call ctx%backend%set_field_data(u_stag_field, u_stag, DIR_X)
    call ctx%backend%tds_solve(du_stag_field, u_stag_field, ctx%tdsops_slot%op)
    call ctx%backend%get_field_data(du_stag, du_stag_field, DIR_X)
    
    err = norm2(du_stag - ref) / sqrt(real(SZ * n_out * n_block, dp))

    deallocate(u_stag, du_stag, ref)
  end subroutine run_stag_p2v_single

end program test_poisson_compact6_convergence_rate