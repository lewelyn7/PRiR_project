EXECS=mpi_hello_world
MPICC?=mpicc

all: ${EXECS}

mpi_parallel_qselect: mpi_parallel_qselect.c
    ${MPICC} -o mpi_parallel_qselect mpi_parallel_qselect.c

clean:
    rm ${EXECS}