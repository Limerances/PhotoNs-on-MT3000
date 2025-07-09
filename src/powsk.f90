

subroutine powsk(data, nside, param, ps, psn)
        use decomp_2d
        use decomp_2d_fft

        use MPI
        use data_prec

        implicit none

        integer :: c,i, j, k, l,m,n, nhalf, idxk
        integer :: idx
        integer, dimension(3) :: nside
        integer, dimension(3) :: fft_start, fft_end, fft_size
        integer, dimension(xsize(1)) :: psn
        real(cputype), dimension(xsize(1)) :: ps
        real(cputype), dimension(xsize(1)*xsize(2)*xsize(3)) :: data
        real(cputype),dimension(2) :: param
        
        real(cputype), parameter :: M_PI = 3.1415926
        real(mytype) k2, k2i, dltr, dlti
        real(mytype) fx, fy, fz, ff, smth
        real(mytype) smooth, box, ismth2, gf, pref
        real(mytype), allocatable, dimension(:,:,:) :: in
        complex(mytype), allocatable, dimension(:,:,:) :: out
        
        call decomp2d_init 
        call decomp_2d_fft_init
#ifdef MEMLOG
        idx = 257
        call mock_mymalloc(xsize(1) * xsize(2) * xsize(3) * SIZEOF(mytype) * 2 * 2, idx)
        idx = 256
        call mock_mymalloc(xsize(1) * xsize(2) * xsize(3) * SIZEOF(mytype) * 2, idx)
        call mock_myfree(0, idx)
#endif

        allocate(in(xstart(1):xend(1), &
                                xstart(2):xend(2),xstart(3):xend(3)))
#ifdef MEMLOG
        idx = 258
        call mock_mymalloc(xsize(1) * xsize(2) * xsize(3) * SIZEOF(mytype), idx)
#endif


        do k=0,xsize(3)-1,1
        do j=0,xsize(2)-1,1
        do i=0,xsize(1)-1,1
                in( i + xstart(1), j + xstart(2), k + xstart(3) )=&
                       data(((i)*xsize(2)+j)*xsize(3)+k+1)
        end do
        end do
        end do

        call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)

        allocate (out(fft_start(1):fft_end(1), &
                                fft_start(2):fft_end(2), fft_start(3):fft_end(3)))
#ifdef MEMLOG
        idx = 259
        call mock_mymalloc(fft_size(1) * fft_size(2) * fft_size(3) * SIZEOF(mytype) * 2, idx)
#endif


        ! compute r2c transform
        call decomp_2d_fft_3d(in,out)

        ! Green function
        nhalf = nside(1)/2
        smooth = param(1)
        box = param(2)


        ismth2 = 2*M_PI*smooth/box
        ismth2 = ismth2 * ismth2
        pref =  box * box / ( M_PI * nside(1) * nside(2) * nside(3) )


        do k=fft_start(3),fft_end(3),1
        n = k - 1 
        if (n>nhalf) then
                n = n - nside(3)
        endif
            do j=fft_start(2),fft_end(2),1
            m = j - 1
            if (m>nhalf) then
                m = m - nside(2)
            endif
                do, i=fft_start(1),fft_end(1),1
                l = i - 1
                if (l>nhalf) then
                        l = l - nside(1)
                endif

        fx=1
        fy=1
        fz=1

#ifdef PK_CORRECT_CIC 
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



                k2 = REAL(n*n + m*m + l*l)

                idxk = floor( sqrt(k2) ) + 1
                
                if (idxk < xsize(1) ) then
                    dltr = REALPART(out(i,j,k))
                    dlti = IMAGPART(out(i,j,k))
!                write(*,*) dltr, dlti, idxk
!                    ps(idxk) = ps(idxk) + out(i,j,k)*out(i,j,k)

#ifdef PK_CORRECT_CIC 
                       dltr = dltr*smth                 
                       dlti = dlti*smth                 
#endif

                    ps(idxk) = ps(idxk) + dltr*dltr + dlti*dlti
                    psn(idxk)= psn(idxk) + 1
                endif

                end do
            end do
        end do

        deallocate(in,out)
#ifdef MEMLOG
        idx = 258
        call mock_myfree(in, idx)
        idx = 259
        call mock_myfree(out, idx)
#endif
        call decomp_2d_fft_finalize
#ifdef MEMLOG
        idx = 257
        call mock_myfree(0, idx)
#endif
        call decomp2d_free 
        return

end subroutine powsk

