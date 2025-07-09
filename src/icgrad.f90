!       interface to convolution by 2dcomp lib
!       Q. Wang 2017 Mar 9

!
!subroutine icifft(data, nside, param )
!        use decomp_2d
!        use decomp_2d_fft
!
!        use MPI
!
!        implicit none
!
!        integer :: i, j, k, l,m,n, nhalf
!        integer, dimension(3) :: nside
!        integer, dimension(3) :: fft_start, fft_end, fft_size
!        real(cputype), dimension(xsize(1)*xsize(2)*xsize(3)) :: data
!        real(cputype),dimension(3) :: param
!        real(cputype), parameter :: M_PI = 3.1415926
!        real*8 k2, nor, kv, kx, ky, kz
!        real*8 smooth, box, ismth2, gf, pref, mag, delta, phase, smth
!        real*8 delta_re, delta_im
!        real*8 fx, fy, fz, ff
!        real(mytype), allocatable, dimension(:,:,:) :: in2
!        complex(mytype), allocatable, dimension(:,:,:) :: out
!        complex(mytype) tmps
!
!        call decomp_2d_fft_init
!
!        call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)
!
!        allocate (out(fft_start(1):fft_end(1), &
!                 fft_start(2):fft_end(2), fft_start(3):fft_end(3)))
!
!        nhalf = nside(1)/2
!
!        smooth = param(1)
!        box = param(2)     
!        nor = 2*M_PI/(box)
!       ! nor = 2*M_PI/(nside(1))
!
!
!        do k=fft_start(3), fft_end(3), 1
!           do j=fft_start(2), fft_end(2), 1
!              do i=fft_start(1), fft_end(1), 1
!                 out(i,j,k) = 0
!              end do
!           end do
!        end do
!
!        do k=fft_start(3),fft_end(3),1
!           n = k - 1 
!           if (n>nhalf) then
!              n = n-nside(3)
!           endif
!
!           do j=fft_start(2),fft_end(2),1
!              m = j - 1
!              if (m>nhalf) then
!                 m = m-nside(2)
!              endif  
!
!              do i=fft_start(1),fft_end(1),1
!                 l = i - 1
!                 if (l>nhalf) then
!                    l = l-nside(1)
!                 endif  
!
!                 kx = REAL(l)*nor
!                 ky = REAL(m)*nor
!                 kz = REAL(n)*nor
!         
!                 k2 = kx*kx + ky*ky + kz*kz
!
!                call ps_ic_lmn(l,m,n, delta_re, delta_im)
!        ! gf = 1 / ( k2  *nor *nor)
!        gf = 1 / ( k2 )
!
!        if (l==0.and.m==0.and.n==0) then
!        gf = 0
!        endif
!
!!        if (k2==0) then
!!        gf = 0
!!        endif
!
!        if (l>nhalf.and.m>nhalf.and.n>nhalf) then
!        gf = 0
!        endif
!
!        fx=1
!        fy=1
!        fz=1
!        if (i>1) then
!        fx=l*M_PI/nside(1)
!        fx=sin(fx)/fx
!        endif
!
!        if (j>1) then
!        fy=m*M_PI/nside(2)
!        fy=sin(fy)/fy
!        endif
!        if (k>1) then
!        fz=n*M_PI/nside(3)
!        fz=sin(fz)/fz
!        endif
!        ff = 1/(fx*fy*fz)
!        smth = ff*ff
!
!        gf = gf*smth
!
!        out(i,j,k) = cmplx(-delta_re*gf, delta_im*gf)
!
!
!        end do
!        end do
!        end do
!
!        ! compute c2r transform
!
!        allocate (in2(xstart(1):xend(1), &
!                 xstart(2):xend(2),xstart(3):xend(3)))
!
!        call decomp_2d_fft_3d(out,in2) 
!
!        do k=0, xsize(3)-1, 1
!           do j=0, xsize(2)-1, 1
!              do i=0, xsize(1)-1, 1
!                 data((i*xsize(2)+j)*xsize(3)+k+1) = &
!                     in2(i+xstart(1), j+xstart(2), k+xstart(3))
!              end do
!           end do
!        end do
!
!        deallocate(in2,out)
!
!        call decomp_2d_fft_finalize
!
!        return
!end subroutine icifft
!!
subroutine icgrad(data, nside, param )
        use decomp_2d
        use decomp_2d_fft

        use MPI
        use data_prec

        implicit none

        integer :: i, j, k,dd, l,m,n,ll,mm,nn,ii,jj,kk, nhalf
        integer :: idx
        integer, dimension(3) :: nside
        integer, dimension(3) :: fft_start, fft_end, fft_size
        real(cputype), dimension(xsize(1)*xsize(2)*xsize(3)) :: data
        real(cputype),dimension(3) :: param
        real(cputype), parameter :: M_PI = 3.1415926
        real*8 k2, nor, kv, kx, ky, kz
        real*8 smooth, box, ismth2, gf, pref, mag, delta, phase
        real*8 delta_re, delta_im
        real*8 fx, fy, fz, ff, smth
        real(mytype), allocatable, dimension(:,:,:) :: in2
        complex(mytype), allocatable, dimension(:,:,:) :: out
        complex(mytype) tmps

        dd = param(3) 
        call decomp2d_init 
        call decomp_2d_fft_init
#ifdef MEMLOG
        idx = 16
        call mock_mymalloc(xsize(1) * xsize(2) * xsize(3) * SIZEOF(mytype) * 2 * 2, idx)
        idx = 17
        call mock_mymalloc(xsize(1) * xsize(2) * xsize(3) * SIZEOF(mytype) * 2, idx)
        call mock_myfree(0, idx)
#endif

        call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)

        allocate (out(fft_start(1):fft_end(1), &
                                fft_start(2):fft_end(2), fft_start(3):fft_end(3)))
#ifdef MEMLOG
        idx = 14
        call mock_mymalloc(fft_size(1) * fft_size(2) * fft_size(3) * SIZEOF(mytype) * 2, idx)
#endif

        ! Green function
        nhalf = nside(1)/2
        smooth = param(1)
        box = param(2)     
        nor = 2*M_PI/(box)


        do k=fft_start(3),fft_end(3),1
        do j=fft_start(2),fft_end(2),1
        do i=fft_start(1),fft_end(1),1
        out(i,j,k) = 0
        end do
        end do
        end do

        do k=fft_start(3),fft_end(3),1
        n = k - 1 
        if (n>nhalf) then
        n = n-nside(3)
        endif

        do j=fft_start(2),fft_end(2),1
        m = j - 1
        if (m>nhalf) then
        m = m-nside(2)
        endif  


        do i=fft_start(1),fft_end(1),1
        l = i - 1
        if (l>nhalf) then
        l = l-nside(1)
        endif  

        kx = REAL(l*l)
        ky = REAL(m*m)
        kz = REAL(n*n)
        k2 = kx + ky + kz

!        call ps_ic_k(k2*nor*nor, delta_re, delta_im)

        call ps_ic_lmn(l, m, n, delta_re, delta_im)

        fx=1
        fy=1
        fz=1

#ifdef IC_CORRECT_CIC 
        if (l.NE.0) then
        fx=REAL(l)*M_PI/nside(1)
        fx=sin(fx)/fx
        endif

        if (m.NE.0) then
        fy=REAL(m)*M_PI/nside(2)
        fy=sin(fy)/fy
        endif

        if (n.NE.0) then
        fz=REAL(n)*M_PI/nside(3)
        fz=sin(fz)/fz
        endif
#endif

        ff = 1/(fx*fy*fz)
        smth = ff*ff


        if (dd.EQ.1) then
        kv=REAL(l*nor)
        endif

        if (dd.EQ.2) then
        kv=REAL(m*nor)
        endif

        if (dd.EQ.3) then
        kv=REAL(n*nor)
        endif


        gf = kv / ( k2  *nor *nor)

        gf = gf * smth

        if (l==0.and.m==0.and.n==0) then
        gf = 0
        endif
        out(i,j,k) = cmplx(-delta_re*gf, delta_im*gf)

        end do
        end do
        end do

        ! compute c2r transform
        allocate(in2(xstart(1):xend(1), &
                                xstart(2):xend(2),xstart(3):xend(3)))
#ifdef MEMLOG
        idx = 15
        call mock_mymalloc(xsize(1) * xsize(2) * xsize(3) * SIZEOF(mytype), idx)
#endif

        call decomp_2d_fft_3d(out,in2) 

        do k=0,xsize(3)-1,1
        do j=0,xsize(2)-1,1
        do i=0,xsize(1)-1,1
        data((i*xsize(2)+j)*xsize(3)+k+1) = &
                in2(i+xstart(1), j+xstart(2), k+xstart(3))
        end do
        end do
        end do

        deallocate(in2,out)
#ifdef MEMLOG
        idx = 15
        call mock_myfree(in2, idx)
        idx = 14
        call mock_myfree(out, idx)
#endif
        call decomp_2d_fft_finalize
#ifdef MEMLOG
        idx = 16
        call mock_myfree(0, idx)
#endif
        call decomp2d_free 
        return
end subroutine icgrad

