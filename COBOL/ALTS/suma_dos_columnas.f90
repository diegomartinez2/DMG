program suma_dos_columnas
    implicit none
    real :: num_col1, num_col2
    real :: total_col1 = 0.0
    real :: total_col2 = 0.0
    integer :: io_status
    integer, parameter :: unit = 10

    open(unit=unit, file='DATOS.DAT', status='old', action='read', iostat=io_status)
    if (io_status /= 0) then
        print *, 'Error al abrir DATOS.DAT'
        stop
    end if

    do
        read(unit, *, iostat=io_status) num_col1, num_col2
        if (io_status /= 0) exit
        total_col1 = total_col1 + num_col1
        total_col2 = total_col2 + num_col2
    end do

    close(unit)

    print '(A, F15.2)', 'TOTAL COLUMNA 1: ', total_col1
    print '(A, F15.2)', 'TOTAL COLUMNA 2: ', total_col2

end program suma_dos_columnas
