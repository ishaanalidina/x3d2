program repro_penta_periodic_cos
  use cudafor
  use m_common,      only: dp, pi, BC_PERIODIC
  use m_cuda_common, only: SZ
  use m_cuda_tdsops, only: cuda_tdsops_t, cuda_tdsops_init
  use m_cuda_exec_dist, only: exec_dist_penta_periodic
  implicit none

  integer, parameter :: n_halo  = 4
  integer, parameter :: n_block = 1
  integer, parameter :: n_grids = 6
  integer :: grid_sizes(n_grids), ig, n, j
  real(dp) :: dx, err, prev_err, rate

  real(dp), allocatable         :: u(:,:,:), du(:,:,:), ref(:,:,:)
  real(dp), device, allocatable :: u_dev(:,:,:), du_dev(:,:,:)
  real(dp), device, allocatable :: u_s(:,:,:), u_e(:,:,:)
  type(cuda_tdsops_t) :: op
  type(dim3) :: blocks, threads

  grid_sizes = [16, 32, 64, 128, 256, 512]

  print '(a)', 'BC_PERIODIC: f = cos(2*pi*x), f'' = -2*pi*sin(2*pi*x)'
  print '(4X,A6,A18,A12)', 'N', 'L2 error', 'rate'

  prev_err = 0._dp
  do ig = 1, n_grids
    n  = grid_sizes(ig)
    dx = 1._dp / n

    allocate(u(SZ,n,n_block), du(SZ,n,n_block), ref(SZ,n,n_block))
    allocate(u_dev(SZ,n,n_block), du_dev(SZ,n,n_block))
    allocate(u_s(SZ,n_halo,n_block), u_e(SZ,n_halo,n_block))

    do j = 1, n
      u(:,j,:)   =          cos(2._dp*pi*(j-1)*dx)
      ref(:,j,:) = -2._dp*pi*sin(2._dp*pi*(j-1)*dx)
    end do

    op = cuda_tdsops_init(n, dx, operation='first-deriv', &
                          scheme='compact10_penta', &
                          bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)

    u_dev  = u
    du_dev = 0._dp

    ! Periodic wrap: fill all 4 halo cells each side.
    u_s(:,4,:) = u_dev(:,n,  :)
    u_s(:,3,:) = u_dev(:,n-1,:)
    u_s(:,2,:) = u_dev(:,n-2,:)
    u_s(:,1,:) = u_dev(:,n-3,:)

    u_e(:,1,:) = u_dev(:,1,:)
    u_e(:,2,:) = u_dev(:,2,:)
    u_e(:,3,:) = u_dev(:,3,:)
    u_e(:,4,:) = u_dev(:,4,:)

    blocks  = dim3(n_block, 1, 1)
    threads = dim3(SZ, 1, 1)
    call exec_dist_penta_periodic(du_dev, u_dev, u_s, u_e, op, blocks, threads)
    du = du_dev

    err = norm2(du - ref) / sqrt(real(SZ*n*n_block, dp))

    if (ig == 1) then
      print '(4X,I6,ES18.6)', n, err
    else
      rate = log(prev_err/err) / log(2._dp)
      print '(4X,I6,ES18.6,F12.2)', n, err, rate
    end if
    prev_err = err

    deallocate(u, du, ref, u_dev, du_dev, u_s, u_e)
  end do

end program repro_penta_periodic_cos