PROGRAM Parallel_Hello_World
USE OMP_LIB

INTEGER :: Suma


    Suma = 0;
    !$OMP PARALLEL DO DEFAULT (PRIVATE), REDUCTION(+:Suma)

    DO i=1,1000000000
        Suma = Suma + i
    END DO
    !$OMP END PARALLEL DO


PRINT *, "Total Sum: ", Suma

END PROGRAM Parallel_Hello_World
