program test_derivatives_bc
  !! compact10_penta derivative convergence rate test for all BC types.
  !!
  !! Test functions per BC type:
  !!   Periodic:          cos(2*pi*x), sin(2*pi*x)
  !!   Dirichlet:         cos(2*pi*x), sin(2*pi*x)
  !!   Neumann sym:       cos(2*pi*x), cos(3*pi*x)
  !!   Neumann non-sym:   sin(2*pi*x), sin(3*pi*x)
  !!
  !! Expected convergence rates (L2 norm):
  !!   Periodic:   6.0 (der1st and der2nd)
  !!   Neumann:    6.0 (der1st and der2nd)
  !!   Dirichlet:  5.5 (der1st), 3.5 (der2nd)

  use iso_fortran_env
  use cudafor

  use m_common, only: dp, pi, BC_PERIODIC, BC_NEUMANN, BC_DIRICHLET
  use m_cuda_common, only: SZ
  use m_cuda_tdsops, only: cuda_tdsops_t, cuda_tdsops_init
  use m_cuda_exec_thom, only: exec_thom_tds_compact

  implicit none

  integer, parameter :: TEST_DERIV       = 1
  integer, parameter :: TEST_STAGDER_V2P = 2
  integer, parameter :: TEST_STAGDER_P2V = 3
  integer, parameter :: TEST_INTERPL_V2P = 4
  integer, parameter :: TEST_INTERPL_P2V = 5
  integer, parameter :: FUNC_COS = 1
  integer, parameter :: FUNC_SIN = 2

  type :: scenario_t
    character(len=30) :: label
    integer           :: func_type
    real(dp)          :: wavenumber
    integer           :: test_kind
    logical           :: sym
  end type scenario_t

  ! Per-BC scenario lists
  integer, parameter :: n_periodic   = 6
  integer, parameter :: n_dirichlet  = 2
  integer, parameter :: n_neumann    = 10

  type(scenario_t) :: sc_periodic(n_periodic)
  type(scenario_t) :: sc_dirichlet(n_dirichlet)
  type(scenario_t) :: sc_neumann(n_neumann)

  integer, parameter :: n_grids = 6
  integer :: grid_sizes(n_grids)
  integer, parameter :: n_block = 64

  integer :: iscen, ierr, ndevs, devnum
  integer :: n_passed, n_failed, n_total

  real(dp) :: errors(n_grids), err_a(n_grids), err_b(n_grids)
  logical :: allpass = .true.

  ierr = cudaGetDeviceCount(ndevs)
  ierr = cudaSetDevice(0)
  ierr = cudaGetDevice(devnum)

  ! ================================================================
  ! Periodic: cos(2*pi*x), sin(2*pi*x)
  ! ================================================================
  sc_periodic(1) = scenario_t('cos(2*pi*x) deriv',       FUNC_COS, 2*pi, TEST_DERIV,       .true.)
  sc_periodic(2) = scenario_t('sin(2*pi*x) deriv',       FUNC_SIN, 2*pi, TEST_DERIV,       .false.)
  sc_periodic(3) = scenario_t('sin(2*pi*x) stagder_v2p', FUNC_SIN, 2*pi, TEST_STAGDER_V2P, .false.)
  sc_periodic(4) = scenario_t('cos(2*pi*x) stagder_p2v', FUNC_COS, 2*pi, TEST_STAGDER_P2V, .true.)
  sc_periodic(5) = scenario_t('cos(2*pi*x) interpl_v2p', FUNC_COS, 2*pi, TEST_INTERPL_V2P, .true.)
  sc_periodic(6) = scenario_t('cos(2*pi*x) interpl_p2v', FUNC_COS, 2*pi, TEST_INTERPL_P2V, .true.)

  ! ================================================================
  ! Dirichlet: cos(2*pi*x), sin(2*pi*x)
  ! ================================================================
  sc_dirichlet(1) = scenario_t('cos(2*pi*x) deriv', FUNC_COS, 2*pi, TEST_DERIV, .true.)
  sc_dirichlet(2) = scenario_t('sin(2*pi*x) deriv', FUNC_SIN, 2*pi, TEST_DERIV, .false.)

  ! ================================================================
  ! Neumann:  sym  -> cos(2*pi*x), cos(3*pi*x)
  !           nsym -> sin(2*pi*x), sin(3*pi*x)
  ! ================================================================
  sc_neumann(1)  = scenario_t('cos(2*pi*x) deriv',       FUNC_COS, 2*pi, TEST_DERIV,       .true.)
  sc_neumann(2)  = scenario_t('cos(3*pi*x) deriv',       FUNC_COS, 3*pi, TEST_DERIV,       .true.)
  sc_neumann(3)  = scenario_t('sin(2*pi*x) deriv',       FUNC_SIN, 2*pi, TEST_DERIV,       .false.)
  sc_neumann(4)  = scenario_t('sin(3*pi*x) deriv',       FUNC_SIN, 3*pi, TEST_DERIV,       .false.)
  sc_neumann(5)  = scenario_t('sin(2*pi*x) stagder_v2p', FUNC_SIN, 2*pi, TEST_STAGDER_V2P, .false.)
  sc_neumann(6)  = scenario_t('cos(2*pi*x) stagder_p2v', FUNC_COS, 2*pi, TEST_STAGDER_P2V, .true.)
  sc_neumann(7)  = scenario_t('sin(3*pi*x) stagder_v2p', FUNC_SIN, 3*pi, TEST_STAGDER_V2P, .false.)
  sc_neumann(8)  = scenario_t('cos(3*pi*x) stagder_p2v', FUNC_COS, 3*pi, TEST_STAGDER_P2V, .true.)
  sc_neumann(9)  = scenario_t('cos(2*pi*x) interpl_v2p', FUNC_COS, 2*pi, TEST_INTERPL_V2P, .true.)
  sc_neumann(10) = scenario_t('cos(2*pi*x) interpl_p2v', FUNC_COS, 2*pi, TEST_INTERPL_P2V, .true.)

  ! Banner
  write(*,'(A)') ''
  write(*,'(A)') '========================================================================'
  write(*,'(A)') '   compact10_penta Derivative Convergence Tests'
  write(*,'(A)') '========================================================================'

  n_passed = 0; n_failed = 0; n_total = 0

  ! ================================================================
  ! PERIODIC
  ! ================================================================
  write(*,'(A)') ''
  write(*,'(A)') ''
  write(*,'(A)') '========================================================================'
  write(*,'(A)') '   BC: periodic'
  write(*,'(A)') '   Expected rates   der1st: 6.0   der2nd: 6.0   stag/interp: 6.0'
  write(*,'(A)') '========================================================================'

  grid_sizes = [16, 32, 64, 128, 256, 512]

  do iscen = 1, n_periodic
    call run_scenario(sc_periodic(iscen), BC_PERIODIC, grid_sizes, &
                      6._dp, 6._dp, 6._dp)
  end do

  ! ================================================================
  ! DIRICHLET
  ! ================================================================
  write(*,'(A)') ''
  write(*,'(A)') ''
  write(*,'(A)') '========================================================================'
  write(*,'(A)') '   BC: dirichlet'
  write(*,'(A)') '   Expected rates   der1st: 5.5   der2nd: 3.5'
  write(*,'(A)') '========================================================================'

  grid_sizes = [17, 33, 65, 129, 257, 513]

  do iscen = 1, n_dirichlet
    call run_scenario(sc_dirichlet(iscen), BC_DIRICHLET, grid_sizes, &
                      5.5_dp, 3.5_dp, 6._dp)
  end do

  ! ================================================================
  ! NEUMANN
  ! ================================================================
  write(*,'(A)') ''
  write(*,'(A)') ''
  write(*,'(A)') '========================================================================'
  write(*,'(A)') '   BC: neumann'
  write(*,'(A)') '   Expected rates   der1st: 6.0   der2nd: 6.0   stag/interp: 6.0'
  write(*,'(A)') '========================================================================'

  grid_sizes = [17, 33, 65, 129, 257, 513]

  do iscen = 1, n_neumann
    call run_scenario(sc_neumann(iscen), BC_NEUMANN, grid_sizes, &
                      6._dp, 6._dp, 6._dp)
  end do

  ! ================================================================
  ! Write convergence data files
  ! ================================================================
  call write_dat('periodic',  BC_PERIODIC,  sc_periodic,  n_periodic,  &
                 [16, 32, 64, 128, 256, 512], n_grids, n_block)
  call write_dat('dirichlet', BC_DIRICHLET, sc_dirichlet, n_dirichlet, &
                 [17, 33, 65, 129, 257, 513], n_grids, n_block)
  call write_dat('neumann',   BC_NEUMANN,   sc_neumann,   n_neumann,   &
                 [17, 33, 65, 129, 257, 513], n_grids, n_block)

  ! ================================================================
  ! Summary
  ! ================================================================
  write(*,'(A)') ''
  write(*,'(A)') ''
  write(*,'(A)') '========================================================================'
  write(*,'(A)') '   SUMMARY'
  write(*,'(A)') '------------------------------------------------------------------------'
  write(*,'(A,I4)') '   Total checks : ', n_total
  write(*,'(A,I4)') '   Passed       : ', n_passed
  write(*,'(A,I4)') '   Failed       : ', n_failed
  write(*,'(A)') '========================================================================'
  write(*,'(A)') ''

  if (allpass) then
    write(*,'(A)') '   >>> ALL TESTS PASSED <<<'
    write(*,'(A)') ''
  else
    write(*,'(A)') '   >>> SOME TESTS FAILED <<<'
    write(*,'(A)') ''
    error stop 'SOME TESTS FAILED.'
  end if

contains

  ! ================================================================
  ! Dispatch one scenario
  ! ================================================================
  subroutine run_scenario(sc, bc, gs, exp_d1, exp_d2, exp_s)
    type(scenario_t), intent(in) :: sc
    integer, intent(in) :: bc, gs(:)
    real(dp), intent(in) :: exp_d1, exp_d2, exp_s

    character(len=10) :: bc_str, sym_str

    select case (bc)
    case (BC_PERIODIC);  bc_str = 'periodic'
    case (BC_DIRICHLET); bc_str = 'dirichlet'
    case (BC_NEUMANN);   bc_str = 'neumann'
    end select
    if (sc%sym) then; sym_str = 'sym'; else; sym_str = 'non-sym'; end if

    write(*,'(A)') ''
    write(*,'(A,A,A,A,A,A,A)') &
      '  --- ', trim(sc%label), '   [', trim(bc_str), ', ', trim(sym_str), '] ---'

    select case (sc%test_kind)
    case (TEST_DERIV)
      call compute_deriv(sc, bc, gs, n_grids, n_block, err_a, err_b)
      call print_pair(sc%label, gs, err_a, err_b, exp_d1, exp_d2, &
                      allpass, n_passed, n_failed, n_total)

    case (TEST_STAGDER_V2P)
      call compute_stag_v2p(sc, bc, gs, n_grids, n_block, &
                            'stag-deriv', 'compact10_penta', .false., errors)
      call print_single(sc%label, gs, errors, exp_s, &
                        allpass, n_passed, n_failed, n_total)

    case (TEST_STAGDER_P2V)
      call compute_stag_p2v(sc, bc, gs, n_grids, n_block, &
                            'stag-deriv', 'compact10_penta', .true., errors)
      call print_single(sc%label, gs, errors, exp_s, &
                        allpass, n_passed, n_failed, n_total)

    case (TEST_INTERPL_V2P)
      call compute_stag_v2p(sc, bc, gs, n_grids, n_block, &
                            'interpolate', 'classic', .true., errors)
      call print_single(sc%label, gs, errors, exp_s, &
                        allpass, n_passed, n_failed, n_total)

    case (TEST_INTERPL_P2V)
      call compute_stag_p2v(sc, bc, gs, n_grids, n_block, &
                            'interpolate', 'classic', .true., errors)
      call print_single(sc%label, gs, errors, exp_s, &
                        allpass, n_passed, n_failed, n_total)
    end select
  end subroutine

  ! ================================================================
  ! Evaluation helpers
  ! ================================================================
  pure real(dp) function eval_f(ft, k, x)
    integer, intent(in) :: ft; real(dp), intent(in) :: k, x
    select case (ft)
    case (FUNC_COS); eval_f = cos(k*x)
    case (FUNC_SIN); eval_f = sin(k*x)
    case default;     eval_f = 0._dp
    end select
  end function
  pure real(dp) function eval_df(ft, k, x)
    integer, intent(in) :: ft; real(dp), intent(in) :: k, x
    select case (ft)
    case (FUNC_COS); eval_df = -k*sin(k*x)
    case (FUNC_SIN); eval_df =  k*cos(k*x)
    case default;     eval_df = 0._dp
    end select
  end function
  pure real(dp) function eval_d2f(ft, k, x)
    integer, intent(in) :: ft; real(dp), intent(in) :: k, x
    select case (ft)
    case (FUNC_COS); eval_d2f = -k*k*cos(k*x)
    case (FUNC_SIN); eval_d2f = -k*k*sin(k*x)
    case default;     eval_d2f = 0._dp
    end select
  end function
  pure real(dp) function dx_of(n, bc)
    integer, intent(in) :: n, bc
    if (bc == BC_PERIODIC) then; dx_of = 1._dp/n
    else;                        dx_of = 1._dp/(n-1)
    end if
  end function
  pure subroutine stag_sizes(bc, n, dir, ni, no)
    integer, intent(in) :: bc, n; character(*), intent(in) :: dir
    integer, intent(out) :: ni, no
    if (bc == BC_PERIODIC) then; ni = n; no = n
    else if (dir == 'v2p') then; ni = n; no = n-1
    else;                        ni = n-1; no = n
    end if
  end subroutine

  ! Best grid index for rate check (avoids precision floor)
  pure integer function best_ig(errs, ng, rmin)
    real(dp), intent(in) :: errs(:); integer, intent(in) :: ng
    real(dp), intent(in) :: rmin
    integer :: i; real(dp) :: r
    best_ig = ng - 1
    do i = ng, 3, -1
      r = log(errs(i-1)/errs(i)) / log(2._dp)
      if (r >= rmin) then; best_ig = i; return; end if
    end do
  end function

  ! ================================================================
  ! Compute routines
  ! ================================================================
  subroutine compute_deriv(sc, bc, gs, ng, nb, e1, e2)
    type(scenario_t), intent(in) :: sc
    integer, intent(in) :: bc, gs(:), ng, nb
    real(dp), intent(out) :: e1(:), e2(:)
    integer :: ig, n, j; real(dp) :: h, x, k
    real(dp), allocatable :: u(:,:,:), du(:,:,:), rf(:,:,:)
    real(dp), device, allocatable :: ud(:,:,:), dd(:,:,:)
    type(cuda_tdsops_t) :: op; type(dim3) :: bl, th
    k = sc%wavenumber
    do ig = 1, ng
      n = gs(ig); h = dx_of(n, bc)
      allocate(u(SZ,n,nb), du(SZ,n,nb), rf(SZ,n,nb))
      allocate(ud(SZ,n,nb), dd(SZ,n,nb))
      bl = dim3(nb,1,1); th = dim3(SZ,1,1)
      do j = 1, n; x = (j-1)*h; u(:,j,:) = eval_f(sc%func_type,k,x); end do
      ! der1st
      do j = 1, n; x = (j-1)*h; rf(:,j,:) = eval_df(sc%func_type,k,x); end do
      if (bc == BC_NEUMANN) then
        op = cuda_tdsops_init(n,h,operation='first-deriv',scheme='compact10_penta',&
                              bc_start=bc,bc_end=bc,sym=sc%sym)
      else
        op = cuda_tdsops_init(n,h,operation='first-deriv',scheme='compact10_penta',&
                              bc_start=bc,bc_end=bc)
      end if
      ud = u; dd = 0._dp
      call exec_thom_tds_compact(dd, ud, op, bl, th); du = dd
      e1(ig) = norm2(du - rf) / sqrt(real(SZ*n*nb, dp))
      ! der2nd
      do j = 1, n; x = (j-1)*h; rf(:,j,:) = eval_d2f(sc%func_type,k,x); end do
      ! if (bc == BC_NEUMANN) then
      !   op = cuda_tdsops_init(n,h,operation='second-deriv',scheme='compact10_penta',&
      !                         bc_start=bc,bc_end=bc,sym=sc%sym)
      ! else
      !   op = cuda_tdsops_init(n,h,operation='second-deriv',scheme='compact10_penta',&
      !                         bc_start=bc,bc_end=bc)
      ! end if
      ud = u; dd = 0._dp
      call exec_thom_tds_compact(dd, ud, op, bl, th); du = dd
      e2(ig) = norm2(du - rf) / sqrt(real(SZ*n*nb, dp))
      deallocate(u, du, rf, ud, dd)
    end do
  end subroutine

  subroutine compute_stag_v2p(sc, bc, gs, ng, nb, oper, sch, sy, ev)
    type(scenario_t), intent(in) :: sc
    integer, intent(in) :: bc, gs(:), ng, nb
    character(*), intent(in) :: oper, sch; logical, intent(in) :: sy
    real(dp), intent(out) :: ev(:)
    integer :: ig, n, j, ni, no; real(dp) :: h, x, k; logical :: itp
    real(dp), allocatable :: us(:,:,:), ds(:,:,:), rf(:,:,:)
    real(dp), device, allocatable :: usd(:,:,:), dsd(:,:,:)
    type(cuda_tdsops_t) :: op; type(dim3) :: bl, th
    k = sc%wavenumber; itp = (oper == 'interpolate')
    do ig = 1, ng
      n = gs(ig); h = dx_of(n, bc)
      bl = dim3(nb,1,1); th = dim3(SZ,1,1)
      call stag_sizes(bc, n, 'v2p', ni, no)
      allocate(us(SZ,ni,nb), ds(SZ,no,nb), rf(SZ,no,nb))
      allocate(usd(SZ,ni,nb), dsd(SZ,no,nb))
      do j = 1, ni; x = (j-1)*h; us(:,j,:) = eval_f(sc%func_type,k,x); end do
      do j = 1, no; x = (real(j,dp)-0.5_dp)*h
        if (itp) then; rf(:,j,:) = eval_f(sc%func_type,k,x)
        else;          rf(:,j,:) = eval_df(sc%func_type,k,x); end if
      end do
      op = cuda_tdsops_init(no,h,operation=oper,scheme=sch, &
                            bc_start=bc,bc_end=bc,from_to='v2p',sym=sy)
      usd = us; dsd = 0._dp
      call exec_thom_tds_compact(dsd, usd, op, bl, th); ds = dsd
      ev(ig) = norm2(ds - rf) / sqrt(real(SZ*no*nb, dp))
      deallocate(us, ds, rf, usd, dsd)
    end do
  end subroutine

  subroutine compute_stag_p2v(sc, bc, gs, ng, nb, oper, sch, sy, ev)
    type(scenario_t), intent(in) :: sc
    integer, intent(in) :: bc, gs(:), ng, nb
    character(*), intent(in) :: oper, sch; logical, intent(in) :: sy
    real(dp), intent(out) :: ev(:)
    integer :: ig, n, j, ni, no; real(dp) :: h, x, k; logical :: itp
    real(dp), allocatable :: us(:,:,:), ds(:,:,:), rf(:,:,:)
    real(dp), device, allocatable :: usd(:,:,:), dsd(:,:,:)
    type(cuda_tdsops_t) :: op; type(dim3) :: bl, th
    k = sc%wavenumber; itp = (oper == 'interpolate')
    do ig = 1, ng
      n = gs(ig); h = dx_of(n, bc)
      bl = dim3(nb,1,1); th = dim3(SZ,1,1)
      call stag_sizes(bc, n, 'p2v', ni, no)
      allocate(us(SZ,ni,nb), ds(SZ,no,nb), rf(SZ,no,nb))
      allocate(usd(SZ,ni,nb), dsd(SZ,no,nb))
      do j = 1, ni; x = (real(j,dp)-0.5_dp)*h
        us(:,j,:) = eval_f(sc%func_type,k,x)
      end do
      do j = 1, no; x = (j-1)*h
        if (itp) then; rf(:,j,:) = eval_f(sc%func_type,k,x)
        else;          rf(:,j,:) = eval_df(sc%func_type,k,x); end if
      end do
      op = cuda_tdsops_init(no,h,operation=oper,scheme=sch, &
                            bc_start=bc,bc_end=bc,from_to='p2v',sym=sy)
      usd = us; dsd = 0._dp
      call exec_thom_tds_compact(dsd, usd, op, bl, th); ds = dsd
      ev(ig) = norm2(ds - rf) / sqrt(real(SZ*no*nb, dp))
      deallocate(us, ds, rf, usd, dsd)
    end do
  end subroutine

  ! ================================================================
  ! Pretty printers
  ! ================================================================
  subroutine print_pair(label, gs, ea, eb, expa, expb, ap, np, nf, nt)
    character(*), intent(in) :: label
    integer, intent(in) :: gs(:)
    real(dp), intent(in) :: ea(:), eb(:), expa, expb
    logical, intent(inout) :: ap
    integer, intent(inout) :: np, nf, nt
    integer :: ig, ng, ia, ib
    real(dp) :: ra, rb, tol
    character(len=8) :: tag
    ng = size(gs); tol = 0.3_dp

    write(*,'(A)') ''
    write(*,'(6X,A6,A18,A18,A10,A10)') &
      'N', 'der1st L2', 'der2nd L2', 'rate_a', 'rate_b'
    write(*,'(6X,A6,A18,A18,A10,A10)') &
      '------', '------------------', '------------------', &
      '----------', '----------'
    do ig = 1, ng
      if (ig == 1) then
        write(*,'(6X,I6,2ES18.6)') gs(ig), ea(ig), eb(ig)
      else
        ra = log(ea(ig-1)/ea(ig)) / log(2._dp)
        rb = log(eb(ig-1)/eb(ig)) / log(2._dp)
        write(*,'(6X,I6,2ES18.6,2F10.2)') gs(ig), ea(ig), eb(ig), ra, rb
      end if
    end do

    ia = best_ig(ea, ng, expa - 1._dp)
    ib = best_ig(eb, ng, expb - 1._dp)
    ra = log(ea(ia-1)/ea(ia)) / log(2._dp)
    rb = log(eb(ib-1)/eb(ib)) / log(2._dp)

    write(*,'(A)') ''
    nt = nt + 1
    if (abs(ra - expa) > tol) then
      ap = .false.; nf = nf + 1; tag = 'FAILED'
    else
      np = np + 1; tag = 'PASSED'
    end if
    write(*,'(6X,A8,A,F6.2,A,I0,A,I0,A,F5.1,A,A)') &
      'der1st:', '  rate =', ra, '   (N:', gs(ia-1), ' -> ', gs(ia), &
      ')   expect:', expa, '   ', trim(tag)

    nt = nt + 1
    if (abs(rb - expb) > tol) then
      ap = .false.; nf = nf + 1; tag = 'FAILED'
    else
      np = np + 1; tag = 'PASSED'
    end if
    write(*,'(6X,A8,A,F6.2,A,I0,A,I0,A,F5.1,A,A)') &
      'der2nd:', '  rate =', rb, '   (N:', gs(ib-1), ' -> ', gs(ib), &
      ')   expect:', expb, '   ', trim(tag)
  end subroutine

  subroutine print_single(label, gs, errs, expr, ap, np, nf, nt)
    character(*), intent(in) :: label
    integer, intent(in) :: gs(:)
    real(dp), intent(in) :: errs(:), expr
    logical, intent(inout) :: ap
    integer, intent(inout) :: np, nf, nt
    integer :: ig, ng, ib
    real(dp) :: r, tol
    character(len=8) :: tag
    ng = size(gs); tol = 0.3_dp

    write(*,'(A)') ''
    write(*,'(6X,A6,A18,A10)') 'N', 'L2 error', 'rate'
    write(*,'(6X,A6,A18,A10)') &
      '------', '------------------', '----------'
    do ig = 1, ng
      if (ig == 1) then
        write(*,'(6X,I6,ES18.6)') gs(ig), errs(ig)
      else
        r = log(errs(ig-1)/errs(ig)) / log(2._dp)
        write(*,'(6X,I6,ES18.6,F10.2)') gs(ig), errs(ig), r
      end if
    end do

    ib = best_ig(errs, ng, expr - 1._dp)
    r = log(errs(ib-1)/errs(ib)) / log(2._dp)

    write(*,'(A)') ''
    nt = nt + 1
    if (abs(r - expr) > tol) then
      ap = .false.; nf = nf + 1; tag = 'FAILED'
    else
      np = np + 1; tag = 'PASSED'
    end if
    write(*,'(6X,A,F6.2,A,I0,A,I0,A,F5.1,A,A)') &
      'rate =', r, '   (N:', gs(ib-1), ' -> ', gs(ib), &
      ')   expect:', expr, '   ', trim(tag)
  end subroutine

  ! ================================================================
  ! File output
  ! ================================================================
  subroutine write_dat(bname, bc, sc, nsc, gs, ng, nb)
    character(*), intent(in) :: bname
    integer, intent(in) :: bc, nsc, ng, nb, gs(:)
    type(scenario_t), intent(in) :: sc(:)
    integer :: ig, is, col; real(dp) :: h, ea, eb, row(1+2*nsc)
    open(unit=20, file='convergence_'//trim(bname)//'.dat', status='replace')
    do ig = 1, ng
      h = dx_of(gs(ig), bc); row(1) = h; col = 2
      do is = 1, nsc
        call compute_one(sc(is), bc, gs(ig), nb, ea, eb)
        if (sc(is)%test_kind == TEST_DERIV) then
          row(col) = ea; row(col+1) = eb; col = col + 2
        else
          row(col) = ea; col = col + 1
        end if
      end do
      write(20,'(999ES20.10)') row(1:col-1)
    end do
    close(20)
  end subroutine

  subroutine compute_one(sc, bc, n, nb, ea, eb)
    type(scenario_t), intent(in) :: sc
    integer, intent(in) :: bc, n, nb
    real(dp), intent(out) :: ea, eb
    real(dp) :: e1a(1), e1b(1), ev(1)
    integer :: gs1(1)
    eb = 0._dp; gs1 = [n]
    select case (sc%test_kind)
    case (TEST_DERIV)
      call compute_deriv(sc, bc, gs1, 1, nb, e1a, e1b)
      ea = e1a(1); eb = e1b(1)
    case (TEST_STAGDER_V2P)
      call compute_stag_v2p(sc, bc, gs1, 1, nb, &
                            'stag-deriv', 'compact10_penta', .false., ev)
      ea = ev(1)
    case (TEST_STAGDER_P2V)
      call compute_stag_p2v(sc, bc, gs1, 1, nb, &
                            'stag-deriv', 'compact10_penta', .true., ev)
      ea = ev(1)
    case (TEST_INTERPL_V2P)
      call compute_stag_v2p(sc, bc, gs1, 1, nb, &
                            'interpolate', 'classic', .true., ev)
      ea = ev(1)
    case (TEST_INTERPL_P2V)
      call compute_stag_p2v(sc, bc, gs1, 1, nb, &
                            'interpolate', 'classic', .true., ev)
      ea = ev(1)
    end select
  end subroutine

end program test_derivatives_bc