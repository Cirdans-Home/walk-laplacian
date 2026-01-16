!
!
!                             AMG4PSBLAS version 1.0
!    Algebraic Multigrid Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.7)
!
!    (C) Copyright 2021
!
!        Salvatore Filippone
!        Pasqua D'Ambra
!        Fabio Durastante
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the AMG4PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AMG4PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!
!
program amg_dkatz
  use psb_base_mod
  use psb_ext_mod
  use psb_linsolve_mod
  use psb_util_mod
#ifdef PSB_HAVE_CUDA
  use psb_cuda_mod
#endif
  use amg_prec_mod
  use data_input
  implicit none


  ! input parameters

  character(len=40) :: kmethd, mtrx_file, rhs_file, guess_file, sol_file, part
  character(len=2)  :: filefmt

  ! Krylov solver data
  type solverdata
    character(len=40)  :: kmethd      ! Krylov solver
    integer(psb_ipk_)  :: istopc      ! stopping criterion
    integer(psb_ipk_)  :: itmax       ! maximum number of iterations
    integer(psb_ipk_)  :: itrace      ! tracing
    integer(psb_ipk_)  :: irst        ! restart
    real(psb_dpk_)     :: eps         ! stopping tolerance
  end type solverdata
  type(solverdata)       :: s_choice

  ! preconditioner data
  type precdata

    ! preconditioner type
    character(len=40)  :: descr       ! verbose description of the prec
    character(len=10)  :: ptype       ! preconditioner type

    integer(psb_ipk_)  :: outer_sweeps ! number of outer sweeps: sweeps for 1-level,
                                       ! AMG cycles for ML
    ! general AMG data
    character(len=16)  :: mlcycle      ! AMG cycle type
    integer(psb_ipk_)  :: maxlevs     ! maximum number of levels in AMG preconditioner

    ! AMG aggregation
    character(len=16)  :: aggr_prol    ! aggregation type: SMOOTHED, NONSMOOTHED
    character(len=16)  :: par_aggr_alg    ! parallel aggregation algorithm: DEC, SYMDEC
    character(len=32)  :: aggr_type   ! Type of aggregation SOC1, SOC2, MATCHBOXP
    integer(psb_ipk_)  :: aggr_size   ! Requested size of the aggregates for MATCHBOXP
    character(len=16)  :: aggr_ord    ! ordering for aggregation: NATURAL, DEGREE
    character(len=16)  :: aggr_filter ! filtering: FILTER, NO_FILTER
    real(psb_dpk_)     :: mncrratio  ! minimum aggregation ratio
    real(psb_dpk_), allocatable :: athresv(:) ! smoothed aggregation threshold vector
    integer(psb_ipk_)  :: thrvsz      ! size of threshold vector
    real(psb_dpk_)     :: athres      ! smoothed aggregation threshold
    integer(psb_ipk_)  :: csize       ! minimum size of coarsest matrix

    ! AMG smoother or pre-smoother; also 1-lev preconditioner
    character(len=16)  :: smther      ! (pre-)smoother type: BJAC, AS
    integer(psb_ipk_)  :: jsweeps     ! (pre-)smoother / 1-lev prec. sweeps
    integer(psb_ipk_)  :: novr        ! number of overlap layers
    character(len=16)  :: restr       ! restriction over application of AS
    character(len=16)  :: prol        ! prolongation over application of AS
    character(len=16)  :: solve       ! local subsolver type: ILU, MILU, ILUT,
                                      ! UMF, MUMPS, SLU, FWGS, BWGS, JAC
    integer(psb_ipk_)  :: fill        ! fill-in for incomplete LU factorization
    real(psb_dpk_)     :: thr         ! threshold for ILUT factorization

    ! AMG post-smoother; ignored by 1-lev preconditioner
    character(len=16)  :: smther2     ! post-smoother type: BJAC, AS
    integer(psb_ipk_)  :: jsweeps2    ! post-smoother sweeps
    integer(psb_ipk_)  :: novr2       ! number of overlap layers
    character(len=16)  :: restr2      ! restriction  over application of AS
    character(len=16)  :: prol2       ! prolongation over application of AS
    character(len=16)  :: solve2      ! local subsolver type: ILU, MILU, ILUT,
                                      ! UMF, MUMPS, SLU, FWGS, BWGS, JAC
    integer(psb_ipk_)  :: fill2       ! fill-in for incomplete LU factorization
    real(psb_dpk_)     :: thr2        ! threshold for ILUT factorization

    ! coarsest-level solver
    character(len=16)  :: cmat        ! coarsest matrix layout: REPL, DIST
    character(len=16)  :: csolve      ! coarsest-lev solver: BJAC, SLUDIST (distr.
                                      ! mat.); UMF, MUMPS, SLU, ILU, ILUT, MILU
                                      ! (repl. mat.)
    character(len=16)  :: csbsolve    ! coarsest-lev local subsolver: ILU, ILUT,
                                      ! MILU, UMF, MUMPS, SLU
    integer(psb_ipk_)  :: cfill       ! fill-in for incomplete LU factorization
    real(psb_dpk_)     :: cthres      ! threshold for ILUT factorization
    integer(psb_ipk_)  :: cjswp       ! sweeps for GS or JAC coarsest-lev subsolver

  end type precdata
  type(precdata)       :: p_choice

  ! sparse matrices
  type(psb_dspmat_type)  :: a
  type(psb_ldspmat_type) :: aux_a
  ! To work with different matrix formats we use mold objects
  type(psb_d_ell_sparse_mat), target    :: aell
  type(psb_d_csr_sparse_mat), target   :: acsr
  type(psb_d_coo_sparse_mat), target   :: acoo
  type(psb_d_hll_sparse_mat), target    :: ahll
  type(psb_d_hdia_sparse_mat), target  :: ahdia
  type(psb_d_dns_sparse_mat), target   :: adns
  class(psb_d_base_sparse_mat), pointer :: amold => acsr
  class(psb_d_base_sparse_mat), pointer :: amold2 => acsr

  type(psb_d_base_vect_type), target   :: dvect
  class(psb_d_base_vect_type), pointer :: vmold => dvect
  type(psb_i_base_vect_type), target   :: ivect
  class(psb_i_base_vect_type), pointer :: imold => ivect

  ! preconditioner data
  type(amg_dprec_type)  :: prec
  ! dense matrices
  real(psb_dpk_), allocatable, target ::  aux_b(:,:), d(:), aux_g(:,:), aux_x(:,:)
  real(psb_dpk_), allocatable , save  :: x_col_glob(:), r_col_glob(:)
  real(psb_dpk_), pointer  :: b_col_glob(:), ref_col_glob(:), guess_col_glob(:)
  type(psb_d_vect_type)    :: b_col, x_col, r_col, ref_col

  ! communications data structure
  type(psb_desc_type):: desc_a

  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: iam, np
  integer(psb_lpk_) :: lnp
  ! solver paramters
  integer(psb_ipk_) :: iter, ircode, nlv
  integer(psb_epk_) :: amatsize, precsize, descsize
  real(psb_dpk_)    :: err

  character(len=5)  :: afmt
  character(len=20) :: name, renum, ch_err
  integer(psb_ipk_), parameter :: iunit=12
  integer(psb_ipk_) :: iparm(20)

  ! Defining variables
  real(psb_dpk_) :: alpha,mu

  ! other variables
  integer(psb_ipk_)  :: i, info, j, k, m_problem, nl
  integer(psb_lpk_)  :: lbw, ubw, prf
  real(psb_dpk_)     :: t0, t1, t2, t3, tprec, thier, tslv
  real(psb_dpk_)     :: resmx, resmxp, xdiffn2, xdiffni, xni, xn2
  integer(psb_ipk_)  :: nrhs, nv
  integer(psb_ipk_), allocatable :: ivg(:), ipv(:)
  integer(psb_lpk_), allocatable :: perm(:), glob_indexes(:)
  logical   :: have_guess=.false., have_ref=.false.

#ifdef PSB_HAVE_CUDA
  type(psb_d_cuda_elg_sparse_mat), target   :: aelg
  type(psb_d_cuda_csrg_sparse_mat), target  :: acsrg
  type(psb_d_cuda_hlg_sparse_mat), target   :: ahlg
  type(psb_d_cuda_hdiag_sparse_mat), target :: ahdiag
  type(psb_d_vect_cuda), target         :: dvgpu
  type(psb_i_vect_cuda), target         :: ivgpu
#endif

  call psb_init(ctxt)
  call psb_info(ctxt,iam,np)
#ifdef PSB_HAVE_CUDA
  call psb_cuda_init(ctxt,iam)
#endif

  if (iam < 0) then
    ! This should not happen, but just in case
    call psb_exit(ctxt)
    stop
  endif


  name='amg_dkatz'
  if(psb_get_errstatus() /= 0) goto 9999
  info=psb_success_
  call psb_set_errverbosity(itwo)
  !
  ! Hello world
  !
  if (iam == psb_root_) then
    write(psb_out_unit,*) ' '
    write(psb_out_unit,*) 'Using PSBLAS version: ',psb_version_string_
    write(psb_out_unit,*) 'Using AMG4PSBLAS version: ',amg_version_string_
#ifdef PSB_HAVE_CUDA
    write(psb_out_unit,*) 'This is the CUDA Capable Version'
#endif
    write(psb_out_unit,*) 'This is the ',trim(name),' test program'
    write(psb_out_unit,*) ' '
  end if
  !
  ! get parameters
  !
  call get_parms(ctxt,mtrx_file,rhs_file,guess_file,sol_file,filefmt, &
       & part,afmt,s_choice,p_choice,alpha,mu)

#ifdef PSB_HAVE_CUDA
  select case(psb_toupper(afmt))
  case('ELG')
    amold => aelg
    vmold  => dvgpu
    imold  => ivgpu
  case('HLG')
    call psi_set_hksz(32)
    amold => ahlg
    vmold  => dvgpu
    imold  => ivgpu
  case('HDIAG')
    amold => ahdiag
    vmold  => dvgpu
    imold  => ivgpu
  case('CSRG')
    amold => acsrg
    vmold  => dvgpu
    imold  => ivgpu
  case('HDIA')
    amold => ahdia
  case('CSR')
    amold => acsr
  case('DNS')
    amold => adns
  case('ELL')
    amold => aell
  case('HLL')
    call psi_set_hksz(32)
    amold => ahll
  case default
    write(*,*) 'Unknown format defaulting to HLG'
    amold => ahlg
    vmold  => dvgpu
    imold  => ivgpu
  end select
#else
  select case(psb_toupper(afmt))
  case('ELL')
    amold => aell
  case('HLL')
    call psi_set_hksz(32)
    amold => ahll
  case('HDIA')
    amold => ahdia
  case('CSR')
    amold => acsr
  case('DNS')
    amold => adns
  case default
    write(*,*) 'Unknown format defaulting to CSR'
    amold => acsr
  end select
#endif

  call psb_barrier(ctxt)
  t1 = psb_wtime()
  ! read the input matrix to be processed and (possibly) the rhs,
  ! the initial guess and the reference solution
  nrhs = 1

  if (iam == psb_root_) then
    select case(psb_toupper(filefmt))

    case('MM')
      ! For Matrix Market we have an input file for the matrix
      ! and (optional) separate files for the rhs, the initial guess
      ! and the reference solution
      call mm_mat_read(aux_a,info,iunit=iunit,filename=mtrx_file)
      if ((info == psb_success_).and.(rhs_file /= 'NONE')) &
           & call mm_array_read(aux_b,info,iunit=iunit,filename=rhs_file)
      if ((info == psb_success_).and.(guess_file /= 'NONE')) then
        call mm_array_read(aux_g,info,iunit=iunit,filename=guess_file)
        have_guess = .true.
      end if
      if ((info == psb_success_).and.(sol_file /= 'NONE')) then
        call mm_array_read(aux_x,info,iunit=iunit,filename=sol_file)
        have_ref = .true.
      end if

    case ('HB')
      ! For Harwell-Boeing we have a single file which may or may not
      ! contain rhs, initial guess and reference solution.
      call hb_read(aux_a,info,iunit=iunit,b=aux_b,&
           & g=aux_g,x=aux_x,filename=mtrx_file)
      have_guess = allocated(aux_g)
      have_ref   = allocated(aux_x)

    case default
      info = -1
      write(psb_err_unit,*) 'Wrong choice for fileformat ', filefmt
    end select

    if (info /= psb_success_) then
      write(psb_err_unit,*) 'Error while reading input matrix '
      call psb_abort(ctxt)
    end if

    m_problem = aux_a%get_nrows()
    call psb_bcast(ctxt,m_problem)
    call psb_bcast(ctxt,have_guess)
    call psb_bcast(ctxt,have_ref)

    ! At this point aux_b may still be unallocated
    if (psb_size(aux_b,dim=ione) == m_problem) then
      ! if any rhs were present, broadcast the first one
      write(psb_err_unit,'("Ok, got an rhs ")')
      b_col_glob =>aux_b(:,1)
    else
      write(psb_out_unit,'("Generating an rhs...")')
      call psb_realloc(m_problem,1,aux_b,ircode)
      if (ircode /= 0) then
        call psb_errpush(psb_err_alloc_dealloc_,name)
        goto 9999
      endif
      b_col_glob => aux_b(:,1)
      do i=1, m_problem
        b_col_glob(i) = 1.d0
      enddo
    endif

    if ((have_guess).and.(psb_size(aux_g,dim=ione) == m_problem)) then
      ! if any initial guess were present, broadcast the first one
      write(psb_err_unit,'("Ok, got an initial guess ")')
      guess_col_glob =>aux_g(:,1)
    else
      write(psb_out_unit,'("Generating an initial guess...")')
      call psb_realloc(m_problem,1,aux_g,ircode)
      if (ircode /= 0) then
        call psb_errpush(psb_err_alloc_dealloc_,name)
        goto 9999
      endif
      guess_col_glob => aux_g(:,1)
      do i=1, m_problem
        guess_col_glob(i) = 0.d0
      enddo
    endif

    if ((have_ref).and.(psb_size(aux_x,dim=ione) == m_problem)) then
      ! if any reference were present, broadcast the first one
      write(psb_err_unit,'("Ok, got a reference solution ")')
      ref_col_glob =>aux_x(:,1)
    else
      write(psb_out_unit,'("No reference solution...")')
    endif

    ! clean zeros in the input matrix
    call aux_a%clean_zeros(info)

  else
    call psb_bcast(ctxt,m_problem)
    call psb_bcast(ctxt,have_guess)
    call psb_bcast(ctxt,have_ref)
  end if


  !
  ! Renumbering (NONE for the moment)
  !
  if (iam==psb_root_) then
    renum='NONE'

    call psb_cmp_bwpf(aux_a,lbw,ubw,prf,info)
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,*) 'Bandwidth and profile: ',lbw,ubw,prf
    write(psb_out_unit,*) 'Renumbering algorithm: ',psb_toupper(renum)

    if (trim(psb_toupper(renum))/='NONE') then
      call psb_mat_renum(renum,aux_a,info,perm=perm)
      if (info /= 0) then
        write(psb_err_unit,*) 'Error from RENUM',info
        goto 9999
      end if
      call psb_gelp('N',perm(1:m_problem),b_col_glob(1:m_problem),info)
      call psb_cmp_bwpf(aux_a,lbw,ubw,prf,info)
      write(psb_out_unit,*) 'Bandwidth and profile (renumberd):',lbw,ubw,prf
    end if

    write(psb_out_unit,'(" ")')

  end if

  !
  ! switch over different partition types
  !
  select case (psb_toupper(part))
  case('BLOCK')
    call psb_barrier(ctxt)
    if (iam == psb_root_) write(psb_out_unit,'("Partition type: block")')
    call psb_dmatdist_zerodiag(aux_a, a,  ctxt, desc_a,info,fmt=afmt,parts=part_block,mold=amold)
  case('GRAPH')
    if (iam == psb_root_) then
      write(psb_out_unit,'("Partition type: graph")')
      write(psb_out_unit,'(" ")')
      !      write(psb_err_unit,'("Build type: graph")')
      call aux_a%cscnv(info,type='csr')
      lnp = np
      call build_mtpart(aux_a,lnp)

    endif
    call distr_mtpart(psb_root_,ctxt)
    call getv_mtpart(ivg)
    call psb_dmatdist_zerodiag(aux_a, a, ctxt,desc_a,info,fmt=afmt,vg=ivg,mold=amold)
  case default
    if (iam == psb_root_) write(psb_out_unit,'("Partition type: block")')
    call psb_dmatdist_zerodiag(aux_a, a,  ctxt, desc_a,info,fmt=afmt,parts=part_block,mold=amold)
  end select

  !
  ! Build the matrix I − αA − α^2(I − D)
  !
  d = a%rowsum(info) ! local row sum
  do i=1,desc_a%get_local_rows()
    d(i) = 1.d0 -alpha*alpha*(1.d0 - d(i))
  end do
  call psb_scal(-alpha,a,info)
  glob_indexes = desc_a%get_global_indices(owned=.true.)
  nl = desc_a%get_local_rows()
  call desc_a%indxmap%set_state(psb_desc_bld_)
  call psb_spins(nl,glob_indexes,glob_indexes,d,a,desc_a,info)

  call psb_barrier(ctxt)
  t0 = psb_wtime()
  call psb_cdasb(desc_a,info,mold=imold)
  t1 = psb_wtime()
  if(info /= psb_success_)then
    info=psb_err_from_subroutine_
    ch_err='psb_cdasb'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  call psb_barrier(ctxt)
  t2 = psb_wtime()
  call psb_spasb(a,desc_a,info,mold=amold)
  t3 = psb_wtime()
  if(info /= psb_success_)then
    info=psb_err_from_subroutine_
    ch_err='psb_spasb'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (iam == psb_root_) then
    write(psb_out_unit,*) 'descriptor assembly: ',t1-t0
    write(psb_out_unit,*) 'sparse matrix assembly: ',t3-t2
  end if

  !
  ! Scatter rhs, initial guess and reference solution
  !
  call psb_geall(b_col,desc_a,info)
  call psb_geall(x_col,desc_a,info)
  if (have_ref) call psb_geall(ref_col,desc_a,info)

  if (iam == psb_root_) write(psb_out_unit,'("Scatter rhs")')
  call psb_scatter(b_col_glob,b_col,desc_a,info)
  if (iam == psb_root_) write(psb_out_unit,'("Scatter initial guess")')
  call psb_scatter(guess_col_glob,x_col,desc_a,info)
  if (have_ref) then
    if (iam == psb_root_) write(psb_out_unit,'("Scatter reference solution")')
    call psb_scatter(ref_col_glob,ref_col,desc_a,info)
  end if
  call psb_geasb(x_col,desc_a,info,mold=vmold)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_geasb x_col'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  call psb_geasb(b_col,desc_a,info,mold=vmold)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_geasb b_col'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  t2 = psb_wtime() - t1
  call psb_amx(ctxt, t2)

  if (iam == psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("Time to read and partition matrix, rhs(, guess, ref sol) : ",es12.5)')t2
    write(psb_out_unit,'(" ")')
  end if

  !
  ! initialize the preconditioner
  !
  call prec%init(ctxt,p_choice%ptype,info)
  select case(trim(psb_toupper(p_choice%ptype)))
  case ('NONE','NOPREC')
    ! Do nothing, keep defaults

  case ('JACOBI','GS','FWGS','FBGS')
    ! 1-level sweeps from "outer_sweeps"
    call prec%set('smoother_sweeps', p_choice%jsweeps, info)

  case ('BJAC')
    call prec%set('smoother_sweeps', p_choice%jsweeps, info)
    call prec%set('sub_solve',       p_choice%solve,   info)
    if (psb_toupper(p_choice%solve)=='MUMPS') &
         & call prec%set('mumps_loc_glob','local_solver',info)
    call prec%set('sub_fillin',      p_choice%fill,    info)
    call prec%set('sub_iluthrs',     p_choice%thr,     info)

  case('AS')
    call prec%set('smoother_sweeps', p_choice%jsweeps, info)
    call prec%set('sub_ovr',         p_choice%novr,    info)
    call prec%set('sub_restr',       p_choice%restr,   info)
    call prec%set('sub_prol',        p_choice%prol,    info)
    call prec%set('sub_solve',       p_choice%solve,   info)
    if (psb_toupper(p_choice%solve)=='MUMPS') &
         & call prec%set('mumps_loc_glob','local_solver',info)
    call prec%set('sub_fillin',      p_choice%fill,    info)
    call prec%set('sub_iluthrs',     p_choice%thr,     info)

  case ('ML')
    ! multilevel preconditioner

    call prec%set('ml_cycle',        p_choice%mlcycle,    info)
    call prec%set('outer_sweeps',    p_choice%outer_sweeps,info)
    if (p_choice%csize>0)&
         & call prec%set('min_coarse_size', p_choice%csize,      info)
    if (p_choice%mncrratio>1)&
         & call prec%set('min_cr_ratio',   p_choice%mncrratio, info)
    if (p_choice%maxlevs>0)&
         & call prec%set('max_levs',    p_choice%maxlevs,    info)
    if (p_choice%athres >= dzero) &
         & call prec%set('aggr_thresh',     p_choice%athres,  info)
    if (p_choice%thrvsz>0) then
      do k=1,min(p_choice%thrvsz,size(prec%precv)-1)
        call prec%set('aggr_thresh',     p_choice%athresv(k),  info,ilev=(k+1))
      end do
    end if

    call prec%set('aggr_prol',       p_choice%aggr_prol,   info)
    call prec%set('par_aggr_alg',    p_choice%par_aggr_alg,   info)
    call prec%set('aggr_type',       p_choice%aggr_type, info)
    call prec%set('aggr_size',       p_choice%aggr_size, info)

    call prec%set('aggr_ord',        p_choice%aggr_ord,   info)
    call prec%set('aggr_filter',     p_choice%aggr_filter,info)


    call prec%set('smoother_type',   p_choice%smther,     info)
    call prec%set('smoother_sweeps', p_choice%jsweeps,    info)

    select case (psb_toupper(p_choice%smther))
    case ('GS','BWGS','FBGS','JACOBI','L1-JACOBI')
      ! do nothing
    case default
      call prec%set('sub_ovr',         p_choice%novr,       info)
      call prec%set('sub_restr',       p_choice%restr,      info)
      call prec%set('sub_prol',        p_choice%prol,       info)
      call prec%set('sub_solve',       p_choice%solve,      info)
      if (psb_toupper(p_choice%solve)=='MUMPS') &
           & call prec%set('mumps_loc_glob','local_solver',info)
      call prec%set('sub_fillin',      p_choice%fill,       info)
      call prec%set('sub_iluthrs',     p_choice%thr,        info)
    end select

    if (psb_toupper(p_choice%smther2) /= 'NONE') then
      call prec%set('smoother_type',   p_choice%smther2,   info,pos='post')
      call prec%set('smoother_sweeps', p_choice%jsweeps2,  info,pos='post')
      select case (psb_toupper(p_choice%smther2))
      case ('GS','BWGS','FBGS','JACOBI','L1-JACOBI')
        ! do nothing
      case default
        call prec%set('sub_ovr',         p_choice%novr2,     info,pos='post')
        call prec%set('sub_restr',       p_choice%restr2,    info,pos='post')
        call prec%set('sub_prol',        p_choice%prol2,     info,pos='post')
        call prec%set('sub_solve',       p_choice%solve2,    info,pos='post')
        if (psb_toupper(p_choice%solve2)=='MUMPS') &
             & call prec%set('mumps_loc_glob','local_solver',info)
        call prec%set('sub_fillin',      p_choice%fill2,     info,pos='post')
        call prec%set('sub_iluthrs',     p_choice%thr2,      info,pos='post')
      end select
    end if

    call prec%set('coarse_solve',    p_choice%csolve,    info)
    if (psb_toupper(p_choice%csolve) == 'BJAC') &
         &  call prec%set('coarse_subsolve', p_choice%csbsolve,  info)
    call prec%set('coarse_mat',      p_choice%cmat,      info)
    call prec%set('coarse_fillin',   p_choice%cfill,     info)
    call prec%set('coarse_iluthrs',  p_choice%cthres,    info)
    call prec%set('coarse_sweeps',   p_choice%cjswp,     info)

  end select

  ! build the preconditioner
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call prec%hierarchy_build(a,desc_a,info)
  thier = psb_wtime()-t1
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='amg_hierarchy_bld')
    goto 9999
  end if
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call prec%smoothers_build(a,desc_a,info,amold=amold,vmold=vmold,imold=imold)
  tprec = psb_wtime()-t1
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='amg_smoothers_bld')
    goto 9999
  end if

  call psb_amx(ctxt, thier)
  call psb_amx(ctxt, tprec)

  if(iam == psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("Preconditioner: ",a)') trim(p_choice%descr)
    write(psb_out_unit,'("Preconditioner time: ",es12.5)')thier+tprec
    write(psb_out_unit,'(" ")')
  end if

  !
  ! iterative method parameters
  !
  call prec%allocate_wrk(info,vmold=vmold)
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call psb_krylov(s_choice%kmethd,a,prec,b_col,x_col,s_choice%eps,&
       & desc_a,info,itmax=s_choice%itmax,iter=iter,err=err,itrace=s_choice%itrace,&
       & istop=s_choice%istopc,irst=s_choice%irst)
  call psb_barrier(ctxt)
  tslv = psb_wtime() - t1

  call psb_amx(ctxt,tslv)

  ! compute residual norms
  call psb_geall(r_col,desc_a,info)
  call r_col%zero()
  call psb_geasb(r_col,desc_a,info,mold=vmold)
  call psb_geaxpby(done,b_col,dzero,r_col,desc_a,info)
  call psb_spmm(-done,a,x_col,done,r_col,desc_a,info)
  resmx  = psb_genrm2(r_col,desc_a,info)
  resmxp = psb_geamax(r_col,desc_a,info)

  ! compute error in solution
  if (have_ref) then
    call psb_geaxpby(-done,x_col,done,ref_col,desc_a,info)
    xdiffn2  = psb_genrm2(ref_col,desc_a,info)
    xdiffni  = psb_geamax(ref_col,desc_a,info)
    xn2      = psb_genrm2(ref_col,desc_a,info)
    xni      = psb_geamax(ref_col,desc_a,info)
  end if

  amatsize = a%sizeof()
  descsize = desc_a%sizeof()
  precsize = prec%sizeof()
  call psb_sum(ctxt,amatsize)
  call psb_sum(ctxt,descsize)
  call psb_sum(ctxt,precsize)
  call prec%descr(info,iout=psb_out_unit)
  if (iam == psb_root_) then
    write(psb_out_unit,'("Matrix: ",a)')mtrx_file
    write(psb_out_unit,'("Computed solution on ",i8," processors")')np
    write(psb_out_unit,'("Linear system size                 : ",i12)') desc_a%get_global_rows()
    write(psb_out_unit,'("Value of alpha                     : ",es12.5)') alpha
    write(psb_out_unit,'("Value of mu                        : ",es12.5)') mu
    write(psb_out_unit,'("Krylov method                      : ",a)') trim(s_choice%kmethd)
    write(psb_out_unit,'("Preconditioner                     : ",a)') trim(p_choice%descr)
    write(psb_out_unit,'("Iterations to convergence          : ",i12)')iter
    write(psb_out_unit,'("Relative error estimate on exit    : ",es12.5)') err
    write(psb_out_unit,'("Number of levels in hierarchy      : ",i12)') prec%get_nlevs()
    write(psb_out_unit,'("Time to build hierarchy            : ",es12.5)')thier
    write(psb_out_unit,'("Time to build smoothers            : ",es12.5)')tprec
    write(psb_out_unit,'("Total time for preconditioner      : ",es12.5)')tprec+thier
    write(psb_out_unit,'("Time to solve system               : ",es12.5)')tslv
    write(psb_out_unit,'("Time per iteration                 : ",es12.5)')tslv/iter
    write(psb_out_unit,'("Total time                         : ",es12.5)')tslv+tprec+thier
    write(psb_out_unit,'("Residual 2-norm                    : ",es12.5)')resmx
    write(psb_out_unit,'("Residual inf-norm                  : ",es12.5)')resmxp
    write(psb_out_unit,'("Total memory occupation for A      : ",i12)')amatsize
    write(psb_out_unit,'("Total memory occupation for DESC_A : ",i12)')descsize
    write(psb_out_unit,'("Total memory occupation for PREC   : ",i12)')precsize
    write(psb_out_unit,'("Storage format for A               : ",a  )')a%get_fmt()
    write(psb_out_unit,'("Storage format for DESC_A          : ",a  )')desc_a%get_fmt()
    if (have_ref) then
      write(psb_out_unit,'(" ")')
      write(psb_out_unit,'(2x,a10,9x,a8,4x,a20,5x,a8)') &
        & '||X-XREF||','||XREF||','||X-XREF||/||XREF||','(2-norm)'
      write(psb_out_unit,'(1x,3(e12.6,6x))') xdiffn2,xn2,xdiffn2/xn2
      write(psb_out_unit,'(" ")')
      write(psb_out_unit,'(2x,a10,9x,a8,4x,a20,4x,a10)') &
        & '||X-XREF||','||XREF||','||X-XREF||/||XREF||','(inf-norm)'
      write(psb_out_unit,'(1x,3(e12.6,6x))') xdiffni,xni,xdiffni/xni
    end if

  end if

!   call psb_gather(x_col_glob,x_col,desc_a,info,root=psb_root_)
!   if (info == psb_success_) &
!        & call psb_gather(r_col_glob,r_col,desc_a,info,root=psb_root_)
!   if (info /= psb_success_) goto 9999
!   if (iam == psb_root_) then
!     write(psb_err_unit,'(" ")')
!     write(psb_err_unit,'("Saving x on file")')
!     write(20,*) 'Matrix: ',mtrx_file
!     write(20,*) 'Krylov method:',trim(s_choice%kmethd)
!     write(20,*) 'Preconditioner:',trim(p_choice%descr)
!     write(20,*) 'Computed solution on ',np,' processors.'
!     write(20,*) 'Iterations to convergence: ',iter
!     write(20,*) 'Error estimate (infinity norm) on exit:', &
!          & ' ||r||/||b|| (inf-norm) = ',err
!     write(20,'(" Residual 2-norm 2         : ",es12.5)')resmx
!     write(20,'(" Residual inf-norm         : ",es12.5)')resmxp
!     write(20,'(a8,4(2x,a20))') 'I','X(I)','R(I)','B(I)'
!     do i=1,m_problem
!       write(20,998) i,x_col_glob(i),r_col_glob(i),b_col_glob(i)
!     enddo
!   end if
! 998 format(i8,4(2x,g20.14))
! 993 format(i6,4(1x,e12.6))


  call psb_gefree(b_col,desc_a,info)
  call psb_gefree(x_col,desc_a,info)
  call psb_gefree(r_col,desc_a,info)
  if (have_ref) call psb_gefree(ref_col,desc_a,info)
  call psb_spfree(a, desc_a,info)
  call prec%free(info)
  call psb_cdfree(desc_a,info)
if(info /= psb_success_) then
  info=psb_err_from_subroutine_
  ch_err='free routine'
  call psb_errpush(info,name,a_err=ch_err)
  goto 9999
end if

#ifdef PSB_HAVE_CUDA
  call psb_cuda_exit()
#endif
  call psb_exit(ctxt)
  stop

9999 continue
  call psb_error(ctxt)

contains
  !
  ! get iteration parameters from standard input
  !
  subroutine get_parms(ctxt,mtrx,rhs,guess,sol,filefmt,part,afmt,solve,prec,alpha,mu)

    implicit none

    type(psb_ctxt_type) :: ctxt
    character(len=*)    :: mtrx, rhs, guess, sol, filefmt, afmt, part
    type(solverdata)    :: solve
    type(precdata)      :: prec
    integer(psb_ipk_)   :: iam, nm, np, inp_unit
    real(psb_dpk_)      :: alpha, mu
    character(len=1024)   :: filename

    call psb_info(ctxt,iam,np)

    if (iam == psb_root_) then
      ! read input data
      if (command_argument_count()>0) then
        call get_command_argument(1,filename)
        inp_unit = 30
        open(inp_unit,file=filename,action='read',iostat=info)
        if (info /= 0) then
          write(psb_err_unit,*) 'Could not open file ',filename,' for input'
          call psb_abort(ctxt)
          stop
        else
          write(psb_err_unit,*) 'Opened file ',trim(filename),' for input'
        end if
      else
        inp_unit=psb_inp_unit
      end if
      !
      ! input files
      call read_data(mtrx,inp_unit)            ! matrix file
      call read_data(rhs,inp_unit)             ! rhs file
      call read_data(guess,inp_unit)           ! starting guess file
      call read_data(sol,inp_unit)             ! solution file (for comparison)
      call read_data(filefmt,inp_unit)         ! format of files
      call read_data(afmt,inp_unit)            ! matrix storage format
      call read_data(part,inp_unit)            ! partition type
      call read_data(alpha,inp_unit)           ! alpha value
      call read_data(mu,inp_unit)              ! mu
      ! Krylov solver data
      call read_data(solve%kmethd,inp_unit)    ! Krylov solver
      call read_data(solve%istopc,inp_unit)    ! stopping criterion
      call read_data(solve%itmax,inp_unit)     ! max num iterations
      call read_data(solve%itrace,inp_unit)    ! tracing
      call read_data(solve%irst,inp_unit)      ! restart
      call read_data(solve%eps,inp_unit)       ! tolerance
      ! preconditioner type
      call read_data(prec%descr,inp_unit)      ! verbose description of the prec
      call read_data(prec%ptype,inp_unit)      ! preconditioner type
      ! First smoother / 1-lev preconditioner
      call read_data(prec%smther,inp_unit)     ! smoother type
      call read_data(prec%jsweeps,inp_unit)    ! (pre-)smoother / 1-lev prec sweeps
      call read_data(prec%novr,inp_unit)       ! number of overlap layers
      call read_data(prec%restr,inp_unit)      ! restriction  over application of AS
      call read_data(prec%prol,inp_unit)       ! prolongation over application of AS
      call read_data(prec%solve,inp_unit)      ! local subsolver
      call read_data(prec%fill,inp_unit)       ! fill-in for incomplete LU
      call read_data(prec%thr,inp_unit)        ! threshold for ILUT
      ! Second smoother/ AMG post-smoother (if NONE ignored in main)
      call read_data(prec%smther2,inp_unit)     ! smoother type
      call read_data(prec%jsweeps2,inp_unit)    ! (post-)smoother sweeps
      call read_data(prec%novr2,inp_unit)       ! number of overlap layers
      call read_data(prec%restr2,inp_unit)      ! restriction  over application of AS
      call read_data(prec%prol2,inp_unit)       ! prolongation over application of AS
      call read_data(prec%solve2,inp_unit)      ! local subsolver
      call read_data(prec%fill2,inp_unit)       ! fill-in for incomplete LU
      call read_data(prec%thr2,inp_unit)        ! threshold for ILUT
      ! general AMG data
      call read_data(prec%mlcycle,inp_unit)     ! AMG cycle type
      call read_data(prec%outer_sweeps,inp_unit) ! number of 1lev/outer sweeps
      call read_data(prec%maxlevs,inp_unit)    ! max number of levels in AMG prec
      call read_data(prec%csize,inp_unit)       ! min size coarsest mat
      ! aggregation
      call read_data(prec%aggr_prol,inp_unit)    ! aggregation type
      call read_data(prec%par_aggr_alg,inp_unit)    ! parallel aggregation alg
      call read_data(prec%aggr_type,inp_unit)   ! type of aggregation
      call read_data(prec%aggr_size,inp_unit) ! Requested size of the aggregates for MATCHBOXP
      call read_data(prec%aggr_ord,inp_unit)    ! ordering for aggregation
      call read_data(prec%mncrratio,inp_unit)  ! minimum aggregation ratio
      call read_data(prec%aggr_filter,inp_unit) ! filtering
      call read_data(prec%athres,inp_unit)      ! smoothed aggr thresh
      call read_data(prec%thrvsz,inp_unit)      ! size of aggr thresh vector
      if (prec%thrvsz > 0) then
        call psb_realloc(prec%thrvsz,prec%athresv,info)
        call read_data(prec%athresv,inp_unit)   ! aggr thresh vector
      else
        read(inp_unit,*)                        ! dummy read to skip a record
      end if
      ! coasest-level solver
      call read_data(prec%csolve,inp_unit)      ! coarsest-lev solver
      call read_data(prec%csbsolve,inp_unit)    ! coarsest-lev subsolver
      call read_data(prec%cmat,inp_unit)        ! coarsest mat layout
      call read_data(prec%cfill,inp_unit)       ! fill-in for incompl LU
      call read_data(prec%cthres,inp_unit)      ! Threshold for ILUT
      call read_data(prec%cjswp,inp_unit)       ! sweeps for GS/JAC subsolver
      if (inp_unit /= psb_inp_unit) then
        close(inp_unit)
      end if
    end if

    call psb_bcast(ctxt,mtrx)
    call psb_bcast(ctxt,rhs)
    call psb_bcast(ctxt,guess)
    call psb_bcast(ctxt,sol)
    call psb_bcast(ctxt,filefmt)
    call psb_bcast(ctxt,afmt)
    call psb_bcast(ctxt,part)
    call psb_bcast(ctxt,alpha)
    call psb_bcast(ctxt,mu)

    call psb_bcast(ctxt,solve%kmethd)
    call psb_bcast(ctxt,solve%istopc)
    call psb_bcast(ctxt,solve%itmax)
    call psb_bcast(ctxt,solve%itrace)
    call psb_bcast(ctxt,solve%irst)
    call psb_bcast(ctxt,solve%eps)

    call psb_bcast(ctxt,prec%descr)
    call psb_bcast(ctxt,prec%ptype)

    ! broadcast first (pre-)smoother / 1-lev prec data
    call psb_bcast(ctxt,prec%smther)
    call psb_bcast(ctxt,prec%jsweeps)
    call psb_bcast(ctxt,prec%novr)
    call psb_bcast(ctxt,prec%restr)
    call psb_bcast(ctxt,prec%prol)
    call psb_bcast(ctxt,prec%solve)
    call psb_bcast(ctxt,prec%fill)
    call psb_bcast(ctxt,prec%thr)
    ! broadcast second (post-)smoother
    call psb_bcast(ctxt,prec%smther2)
    call psb_bcast(ctxt,prec%jsweeps2)
    call psb_bcast(ctxt,prec%novr2)
    call psb_bcast(ctxt,prec%restr2)
    call psb_bcast(ctxt,prec%prol2)
    call psb_bcast(ctxt,prec%solve2)
    call psb_bcast(ctxt,prec%fill2)
    call psb_bcast(ctxt,prec%thr2)

    ! broadcast AMG parameters
    call psb_bcast(ctxt,prec%mlcycle)
    call psb_bcast(ctxt,prec%outer_sweeps)
    call psb_bcast(ctxt,prec%maxlevs)

    call psb_bcast(ctxt,prec%aggr_prol)
    call psb_bcast(ctxt,prec%par_aggr_alg)
    call psb_bcast(ctxt,prec%aggr_type)
    call psb_bcast(ctxt,prec%aggr_size)
    call psb_bcast(ctxt,prec%aggr_ord)
    call psb_bcast(ctxt,prec%aggr_filter)
    call psb_bcast(ctxt,prec%mncrratio)
    call psb_bcast(ctxt,prec%thrvsz)
    if (prec%thrvsz > 0) then
      if (iam /= psb_root_) call psb_realloc(prec%thrvsz,prec%athresv,info)
      call psb_bcast(ctxt,prec%athresv)
    end if
    call psb_bcast(ctxt,prec%athres)

    call psb_bcast(ctxt,prec%csize)
    call psb_bcast(ctxt,prec%cmat)
    call psb_bcast(ctxt,prec%csolve)
    call psb_bcast(ctxt,prec%csbsolve)
    call psb_bcast(ctxt,prec%cfill)
    call psb_bcast(ctxt,prec%cthres)
    call psb_bcast(ctxt,prec%cjswp)


  end subroutine get_parms

  subroutine psb_dmatdist_zerodiag(a_glob, a, ctxt, desc_a,&
     & info, parts, vg, vsz, inroot,fmt,mold)
     !
     ! an utility subroutine to distribute a matrix among processors
     ! according to a user defined data distribution, using
     ! sparse matrix subroutines.
     !
     !  type(psb_dspmat)                       :: a_glob
     !     on entry: this contains the global sparse matrix as follows:
     !
     !  type(psb_dspmat_type)                            :: a
     !     on exit : this will contain the local sparse matrix.
     !
     !       interface parts
     !         !   .....user passed subroutine.....
     !         subroutine parts(global_indx,n,np,pv,nv)
     !           implicit none
     !           integer(psb_ipk_), intent(in)  :: global_indx, n, np
     !           integer(psb_ipk_), intent(out) :: nv
     !           integer(psb_ipk_), intent(out) :: pv(*)
     !
     !       end subroutine parts
     !       end interface
     !     on entry:  subroutine providing user defined data distribution.
     !        for each global_indx the subroutine should return
     !        the list  pv of all processes owning the row with
     !        that index; the list will contain nv entries.
     !        usually nv=1; if nv >1 then we have an overlap in the data
     !        distribution.
     !
     !  integer(psb_ipk_) :: ctxt
     !     on entry: the PSBLAS parallel environment context.
     !
     !  type (desc_type)                  :: desc_a
     !     on exit : the updated array descriptor
     !
     !  integer(psb_ipk_), optional    :: inroot
     !     on entry: specifies processor holding a_glob. default: 0
     !     on exit : unchanged.
     !
     use psb_base_mod
     use psb_mat_mod
     implicit none

     ! parameters
     type(psb_ldspmat_type)      :: a_glob
     type(psb_ctxt_type) :: ctxt
     type(psb_dspmat_type)      :: a
     type(psb_desc_type)        :: desc_a
     integer(psb_ipk_), intent(out)       :: info
     integer(psb_ipk_), optional       :: inroot
     character(len=*), optional :: fmt
     class(psb_d_base_sparse_mat), optional :: mold
     procedure(psb_parts), optional  :: parts
     integer(psb_ipk_), optional     :: vg(:)
     integer(psb_ipk_), optional     :: vsz(:)

     ! local variables
     logical           :: use_parts, use_vg, use_vsz
     integer(psb_ipk_) :: np, iam, np_sharing, root, iproc
     integer(psb_ipk_) :: err_act, il, inz
     integer(psb_lpk_) :: k_count, liwork,  nnzero, nrhs,&
          & i, ll, nz, isize, nnr, err
     integer(psb_lpk_) :: i_count, j_count, nrow, ncol, ig, lastigp
     integer(psb_ipk_), allocatable  :: iwork(:), iwrk2(:)
     integer(psb_lpk_), allocatable  :: irow(:),icol(:)
     real(psb_dpk_), allocatable    :: val(:)
     integer(psb_ipk_), parameter    :: nb=30
     real(psb_dpk_)              :: t0, t1, t2, t3, t4, t5
     character(len=20)           :: name, ch_err

     info = psb_success_
     err  = 0
     name = 'psb_d_mat_dist'
     call psb_erractionsave(err_act)

     ! executable statements
     if (present(inroot)) then
       root = inroot
     else
       root = psb_root_
     end if
     call psb_info(ctxt, iam, np)
     if (iam == root) then
       nrow = a_glob%get_nrows()
       ncol = a_glob%get_ncols()
       if (nrow /= ncol) then
         write(psb_err_unit,*) 'a rectangular matrix ? ',nrow,ncol
         info=-1
         call psb_errpush(info,name)
         goto 9999
       endif
       nnzero = a_glob%get_nzeros()
       nrhs   = 1
     endif

     use_parts = present(parts)
     use_vg    = present(vg)
     use_vsz   = present(vsz)
     if (count((/ use_parts, use_vg, use_vsz /)) /= 1) then
       info=psb_err_no_optional_arg_
       call psb_errpush(info,name,a_err=" vg, vsz, parts")
       goto 9999
     endif

     ! broadcast informations to other processors
     call psb_bcast(ctxt,nrow, root)
     call psb_bcast(ctxt,ncol, root)
     call psb_bcast(ctxt,nnzero, root)
     call psb_bcast(ctxt,nrhs, root)
     liwork = max(np, nrow + ncol)
     allocate(iwork(liwork), iwrk2(np),stat = info)
     if (info /= psb_success_) then
       info=psb_err_alloc_request_
       call psb_errpush(info,name,l_err=(/liwork/),a_err='integer')
       goto 9999
     endif
     if (iam == root) then
       write (*, fmt = *) 'start matdist',root, size(iwork),&
            &nrow, ncol, nnzero,nrhs
     endif
     if (use_parts) then
       call psb_cdall(ctxt,desc_a,info,mg=nrow,parts=parts)
     else if (use_vg) then
       call psb_cdall(ctxt,desc_a,info,vg=vg)
     else if (use_vsz) then
       call psb_cdall(ctxt,desc_a,info,nl=vsz(iam+1))
     else
       info = -1
     end if
     if(info /= psb_success_) then
       info=psb_err_from_subroutine_
       ch_err='psb_cdall'
       call psb_errpush(info,name,a_err=ch_err)
       goto 9999
     end if
     inz = ((nnzero+np-1)/np)
     call psb_spall(a,desc_a,info,nnz=inz,dupl=psb_dupl_add_)
     if(info /= psb_success_) then
       info=psb_err_from_subroutine_
       ch_err='psb_spall'
       call psb_errpush(info,name,a_err=ch_err)
       goto 9999
     end if

     isize = 3*nb*max(((nnzero+nrow)/nrow),nb)
     allocate(val(isize),irow(isize),icol(isize),stat=info)
     if(info /= psb_success_) then
       info=psb_err_from_subroutine_
       ch_err='Allocate'
       call psb_errpush(info,name,a_err=ch_err)
       goto 9999
     end if

     i_count   = 1
     if (use_vsz) then
       iproc = 0
       lastigp = vsz(iproc+1)
     end if
     do while (i_count <= nrow)

       if (use_parts) then
         call parts(i_count,nrow,np,iwork, np_sharing)
         !
         ! np_sharing allows for overlap in the data distribution.
         ! If an index is overlapped, then we have to send its row
         ! to multiple processes. NOTE: we are assuming the output
         ! from PARTS is coherent, otherwise a deadlock is the most
         ! likely outcome.
         !
         j_count = i_count
         if (np_sharing == 1) then
           iproc   = iwork(1)
           do
             j_count = j_count + 1
             if (j_count-i_count >= nb) exit
             if (j_count > nrow) exit
             call parts(j_count,nrow,np,iwrk2, np_sharing)
             if (np_sharing /= 1 ) exit
             if (iwrk2(1) /= iproc ) exit
           end do
         end if
       else if (use_vg) then
         np_sharing = 1
         j_count = i_count
         iproc   = vg(i_count)
         iwork(1:np_sharing) = iproc
         do
           j_count = j_count + 1
           if (j_count-i_count >= nb) exit
           if (j_count > nrow) exit
           if (vg(j_count) /= iproc ) exit
         end do
       else if (use_vsz) then
         np_sharing = 1
         j_count = i_count
         iwork(1:np_sharing) = iproc
         do
           j_count = j_count + 1
           if (j_count-i_count >= nb) exit
           if (j_count > nrow) exit
           if (j_count > lastigp) exit
         end do
       end if

       ! now we should insert rows i_count..j_count-1
       nnr = j_count - i_count

       if (iam == root) then

         ll = 0
         do i= i_count, j_count-1
           call a_glob%csget(i,i,nz,&
                & irow,icol,val,info,nzin=ll,append=.true.)
           if (info /= psb_success_) then
             if (nz >min(size(irow(ll+1:)),size(icol(ll+1:)),size(val(ll+1:)))) then
               write(psb_err_unit,*) 'Allocation failure? This should not happen!'
             end if
             call psb_errpush(info,name,a_err=ch_err)
             goto 9999
           end if
           ll = ll + nz
         end do

         do k_count = 1, np_sharing
           iproc = iwork(k_count)

           if (iproc == iam) then
             il = ll
             call psb_spins(il,irow,icol,val,a,desc_a,info)
             if(info /= psb_success_) then
               info=psb_err_from_subroutine_
               ch_err='psb_spins'
               call psb_errpush(info,name,a_err=ch_err)
               goto 9999
             end if
           else
             call psb_snd(ctxt,nnr,iproc)
             call psb_snd(ctxt,ll,iproc)
             call psb_snd(ctxt,irow(1:ll),iproc)
             call psb_snd(ctxt,icol(1:ll),iproc)
             call psb_snd(ctxt,val(1:ll),iproc)
             call psb_rcv(ctxt,ll,iproc)
           endif
         end do
       else if (iam /= root) then

         do k_count = 1, np_sharing
           iproc = iwork(k_count)
           if (iproc == iam) then
             call psb_rcv(ctxt,nnr,root)
             call psb_rcv(ctxt,ll,root)
             if (ll > size(irow)) then
               write(psb_err_unit,*) iam,'need to reallocate ',ll
               deallocate(val,irow,icol)
               allocate(val(ll),irow(ll),icol(ll),stat=info)
               if(info /= psb_success_) then
                 info=psb_err_from_subroutine_
                 ch_err='Allocate'
                 call psb_errpush(info,name,a_err=ch_err)
                 goto 9999
               end if

             endif
             call psb_rcv(ctxt,irow(1:ll),root)
             call psb_rcv(ctxt,icol(1:ll),root)
             call psb_rcv(ctxt,val(1:ll),root)
             call psb_snd(ctxt,ll,root)
             il = ll
             call psb_spins(il,irow,icol,val,a,desc_a,info)
             if(info /= psb_success_) then
               info=psb_err_from_subroutine_
               ch_err='psspins'
               call psb_errpush(info,name,a_err=ch_err)
               goto 9999
             end if
           endif
         end do
       endif
       i_count = j_count
       if ((use_vsz).and.(j_count <= nrow)) then
         if (j_count > lastigp) then
           iproc = iproc + 1
           lastigp = lastigp + vsz(iproc+1)
         end if
       end if
     end do

     deallocate(val,irow,icol,iwork,iwrk2,stat=info)
     if(info /= psb_success_)then
       info=psb_err_from_subroutine_
       ch_err='deallocate'
       call psb_errpush(info,name,a_err=ch_err)
       goto 9999
     end if

     if (iam == root) write (*, fmt = *) 'end matdist'

     call psb_erractionrestore(err_act)
     return

   9999 call psb_error_handler(ctxt,err_act)

     return

end subroutine psb_dmatdist_zerodiag

end program amg_dkatz
