! re2c.f90
program re2c
    implicit none
    integer :: rc
    real    :: re

    do
        read (*, *, iostat=rc) re
        if (rc /= 0) exit
        print '(f0.2)', re * 5 / 4
    end do
end program re2c
