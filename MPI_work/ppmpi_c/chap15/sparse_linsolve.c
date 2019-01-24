/* sparse_linsolve.c
 *
 * Uses iterative methods in PETSc to solve a random 
 * sparse linear system Ax = b.
 *
 * Input:
 *    n:  order of the linear system
 *    diagonal:  value of diagonal entries in matrix
 *    prob:  probability that an off-diagonal entry
 *        is nonzero
 *    initial_dist:  the initial distribution of the
 *        coefficient matrix:  0 indicates that 
 *        the entire matrix initially resides on
 *        process 0, 1 indicates that it will be
 *        distributed among the processes.
 *
 * Output:
 *    Information from PETSc:  solver, tolerances,
 *        type of preconditioning, and storage 
 *        information for the coefficient matrix 
 *    error:  2-norm of error in solution
 *    iterations:  number of iterations before 
 *        solver terminated
 *
 * To compile:
 *    See PETSc makefile.
 *
 * To run:
 *    mpirun -np p linsolve [petsc options]
 *
 * Petsc Options:
 *    Can specify solver and preconditioner.  Options have 
 *    the form
 *        -option_name option_value
 *    Some valid option names and values follow.
 *        1.  ksp_method:  the Krylov Subspace method used 
 *            by the solver.  Possible values:  richardson, 
 *            chebychev, cg, gmres, tcqmr, bcgs, cgs, 
 *            tfqmr, cr, lsqr, preonly.
 *        2.  ksp_rtol:  decrease in residual norm relative
 *            to size of initial residual for convergence.
 *            Values are doubles.
 *        3.  ksp_atol:  absolute size of the residual norm
 *            for convergence.  Values are doubles.
 *        4.  ksp_max_it:  maximum number of iterations
 *            before terminating solve.  Values are ints. 
 *        5.  ksp_monitor: display residual norms. 
 *            (No value.)
 *        6.  pc_method:  the preconditioning method.
 *            Some values are none, jacobi, bjacobi, sor,
 *            icc, ilu, lu.  (lu simply uses LU 
 *            factorization to solve the system.  This is
 *            only available when run with one process.)
 *        7.  pc_sor_omega:  the relaxation factor
 *            for sor.  Value is a double.
 *        8.  pc_sor_its:  the number of inner iterations
 *            to use in the sor.  Value is an int.
 *        9.  pc_bjacobi_blocks:  the number of blocks to
 *            to use in block jacobi.  Value is an int.
 *
 * Algorithm:
 *    1.  Initialize MPI and PETSc.
 *    2.  Build derived datatype for input data.
 *    3a. Process 0 read and broadcast input data.
 *    3b. Processes != 0 receive input data.
 *    4.  if (initial_dist == 0)
 *            process 0 initializes the matrix
 *        else
 *            each process initializes its rows
 *    5.  All processes begin nonblocking distribution of 
 *        coefficient matrix using MatAssemblyBegin.
 *    6.  Use VecCreate and VecDuplicate to create 
 *        vectors b, x, and exact.
 *    7.  Set entries in exact to 1 using VecSet.
 *    8.  Complete distribution of matrix with
 *        MatAssemblyEnd.
 *    9.  Use MatMult to set rhs b = A*exact.
 *    10. Create solver context using SLESCreate.
 *    11. Associate matrix and preconditioning matrix
 *        with solver using SLESSetOperators.
 *    12. Get command line options specifying solver
 *        options such as type of solver and type of
 *        preconditioning.  Use SLESSetFromOptions.
 *    13  Solve the system with SLESSolve.
 *    14. Print information on solver and matrix
 *        using SLESView.
 *    15. Compute norm of error, ||x - exact||_2, using
 *        VecAXPY and VecNorm.
 *    16. Print error norm and number of iterations
 *        using MPIU_printf.
 *    17. Free storage used by PETSc.
 *    18. Shut down PETSc and MPI.
 *
 * Notes:
 *    1.  Our PETSc matrices and vectors are distributed 
 *        by block panels.
 *    2.  It is not necessary to initialize A on process 0:
 *        any entries can be inserted on any process; 
 *        MatAssemblyBegin/End will correctly distribute
 *        them among the processes.
 *    3.  Solver context consists of such things as a
 *        communicator for solver communication and data
 *        structures needed for solver operations such
 *        as matrix-vector multiply.
 *    4.  PETSc provides an extensive error-handling and
 *        traceback facility which we have not illustrated
 *        in this example to make source more readable.
 *
 * See Chap 15, pp. 350 & ff, in PPMPI.
 */
#include <stdio.h>
#include <stdlib.h> /* Needed for drand48 and srand48 */
#include "sles.h"   /* Includes headers needed by PETSc */
#include "mpi.h"

/* Required in all PETSc programs */
static char help[] = "Solve a random sparse linear system";

void Get_input(int my_rank, int* n_ptr, 
               double* diagonal_ptr, double* prob_ptr,
               int* initial_dist_ptr);

void Build_input_datatype(MPI_Datatype* input_datatype_ptr, 
    int* n_ptr, double* diagonal_ptr, double* prob_ptr,
    int* initial_dist_ptr);

void Initialize_matrix(int my_rank, Mat A, int n,
                       double diagonal, double prob,
                       int initial_dist);

void Allocate(int my_rank, char* name, void* list, 
              int size, int datatype);

/*======================================================*/
int main(int argc,char **argv) {
    Vec       x;             /* computed solution       */
    Vec       b;             /* right-hand side         */
    Vec       exact;         /* exact solution          */
    int       p;
    int       my_rank;
    int       n;             /* order of system         */
    double    diagonal;      /* diagonal matrix entries */
    double    prob;          /* probability of an off-  */
                             /* diagonal non-zero       */
    int       initial_dist;  /* =0 matrix initialized   */
                             /* on process 0. =1 dist-  */
                             /* tributed initialization */
    Mat       A;             /* coefficient matrix      */
    double    one = 1.0;
    double    minus_one = -1.0;
    SLES      sles;          /* context for solver      */
    int       iterations;    /* number of iterations    */
                             /* used by solver          */
    double    error;         /* 2-norm of error in sol- */
                             /* ution                   */

    MPI_Init(&argc, &argv);
    PetscInitialize(&argc,&argv,0,0,help);

    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&p);  

    Get_input(my_rank, &n, &diagonal, &prob, &initial_dist);

    MatCreate(MPI_COMM_WORLD, n, n, &A);

    Initialize_matrix(my_rank, A, n, diagonal, prob,
                      initial_dist);

    MatAssemblyBegin(A, FINAL_ASSEMBLY);
    
    /* Create and set vectors */
    VecCreate(MPI_COMM_WORLD,n,&x);
    VecDuplicate(x, &exact);
    VecDuplicate(x, &b);
    VecSet(&one, exact);

    MatAssemblyEnd(A, FINAL_ASSEMBLY);

    MatMult(A, exact, b);

    /* Set up solver context */
    SLESCreate(MPI_COMM_WORLD, &sles);

    /* Identify coefficient matrix and preconditioning */
    /* matrix                                          */
    SLESSetOperators(sles, A, A, 0);

    /* Use solver and preconditioner specified */
    /* by command line options                 */
    SLESSetFromOptions(sles);

    /* Now solve the system */
    SLESSolve(sles, b, x, &iterations);

#ifdef DEBUG
    VecView(x, STDOUT_VIEWER_WORLD);
    MatView(A, STDOUT_VIEWER_WORLD);
#endif

    SLESView(sles, STDOUT_VIEWER_WORLD);

    /* Check solution */
    VecAXPY(&minus_one, exact, x);
    VecNorm(x, NORM_2, &error);
    if (error >= 1.0e-12)
        MPIU_printf(MPI_COMM_WORLD, 
            "Norm of error %g, Iterations %d\n",
            error, iterations);
    else
        MPIU_printf(MPI_COMM_WORLD, 
            "Norm of error < 1.0e-12, Iterations %d\n",
            iterations);

    /* Pause before quitting */
    {
        int anything;
        if (my_rank == 0) {
            printf("Enter a number to continue\n");
            scanf("%d",&anything);
        }
        MPI_Bcast(&anything, 1, MPI_INT, 0, MPI_COMM_WORLD);
        printf("Process %d > anything = %d\n", my_rank, anything);
    }

    /* Free work space */
    VecDestroy(x); 
    VecDestroy(exact);
    VecDestroy(b);
    MatDestroy(A);
    SLESDestroy(sles);

    PetscFinalize();
    MPI_Finalize();
}  /* main */


/*========================================================
 *
 * Process 0 read and broadcast input data
 */
void Get_input(int my_rank, int* n_ptr, 
               double* diagonal_ptr, double* prob_ptr,
               int* initial_dist_ptr) {

    MPI_Datatype input_datatype;

    Build_input_datatype(&input_datatype, n_ptr, 
            diagonal_ptr, prob_ptr, initial_dist_ptr);

    if (my_rank == 0) 
        scanf("%d %lf %lf %d", n_ptr, diagonal_ptr, 
              prob_ptr, initial_dist_ptr);

    MPI_Bcast(n_ptr, 1, input_datatype, 0, MPI_COMM_WORLD);

}  /* Get_input */


/*========================================================
 *
 * Build derived datatype for distributing input data
 */
void Build_input_datatype(MPI_Datatype* input_datatype_ptr, 
    int* n_ptr, double* diagonal_ptr, double* prob_ptr,
    int* initial_dist_ptr){

    int           array_of_block_lengths[4];
    MPI_Aint      array_of_displacements[4];
    MPI_Datatype  array_of_types[4]; 
    MPI_Aint      base_address; 
    MPI_Aint      temp_address;
    int           i;

    for (i = 0; i < 4; i++) {
        array_of_block_lengths[i] = 1;
    }
    array_of_types[0] = MPI_INT;
    array_of_types[1] = MPI_DOUBLE;
    array_of_types[2] = MPI_DOUBLE;
    array_of_types[3] = MPI_INT;
 
    /* Compute displacements from n */
    array_of_displacements[0] = 0; 
    MPI_Address(n_ptr, &base_address);
    MPI_Address(diagonal_ptr, &temp_address);
    array_of_displacements[1] = 
        temp_address - base_address;
    MPI_Address(prob_ptr, &temp_address);
    array_of_displacements[2] = 
        temp_address - base_address;
    MPI_Address(initial_dist_ptr, &temp_address);
    array_of_displacements[3] = 
        temp_address - base_address;

    MPI_Type_struct(4, array_of_block_lengths, 
        array_of_displacements, array_of_types, 
        input_datatype_ptr);
    MPI_Type_commit(input_datatype_ptr);

}  /* Build_input_datatype */


/*========================================================
 *
 * Assign values to matrix entries:  Diagonal entries
 *    get value in diagonal.  Off-diagonals get zero
 *    or a random value in the range (-1,1).
 *    initial_dist = 0: entire matrix initialized by 
 *        process 0.
 *    initial_dist = 1: each processes initializes its
 *        rows.
 */
void Initialize_matrix(int my_rank, Mat A, int n, 
                  double diagonal, double prob,
                  int initial_dist) {
    int*      columns;  /* temporary storage for col */
    double*   temp_row; /* indices and row entries   */
    int       nonzero_count;
    int       my_min_row;
    int       my_max_row;
    int       i, j;

    if (initial_dist == 0) { 
        if (my_rank == 0) {
            my_min_row = 0;
            my_max_row = n;
        } else {
            return;
        }
    } else {
        MatGetOwnershipRange(A, &my_min_row, &my_max_row);
    }
    printf("Process %d: my_min_row = %d, my_max_row = %d\n",
           my_rank, my_min_row, my_max_row);
    fflush(stdout);
        
    /* Allocate temporary storage */
    Allocate(my_rank, "columns", &columns, n, 0);
    Allocate(my_rank, "temp_row", &temp_row, n, 1);

    /* Seed random number generator */
    srand48((long) (my_rank*my_rank));
    for (i = my_min_row; i < my_max_row; i++) {
        nonzero_count = 0;
        for (j = 0; j < n; j++) {
            if (i == j) {
                temp_row[nonzero_count] = diagonal;
                columns[nonzero_count] = j;
                nonzero_count++;
            } else if (drand48() <= prob) {
                temp_row[nonzero_count] = 
                    2.0*drand48()-1.0;
                columns[nonzero_count] = j;
                nonzero_count++;
            }
        }
        /* Insert entries in a single row (row i) into */
        /* the matrix                                  */
        MatSetValues(A, 1, &i, nonzero_count, columns, 
                     temp_row, INSERT_VALUES);
    }
}  /* Initialize_matrix */


/*========================================================
 * 
 * Allocate a list of ints or doubles.  On error exit.  
 * datatype = 0 => int,  datatype = 1 => double.
 */
void Allocate(int my_rank, char* name, void* list, 
              int size, int datatype) {
    int error = 0;

    if (datatype == 0) {
        *((int**)list) = (int*) malloc(size*sizeof(int));
        if (*((int**)list) == (int*) NULL) error = 1;
    } else {
        *((double**)list) = 
            (double*) malloc(size*sizeof(double));
        if (*((double**)list) == (double*) NULL) error = 1;
    }

    if (error) {
        fprintf(stderr, 
            "Process %d > Malloc failed for %s!\n", 
            my_rank, name); 
        fprintf(stderr, "Process %d > size = %d\n", 
            my_rank, size);  
        fprintf(stderr, "Process %d > Quitting.\n", 
            my_rank);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
}  /* Allocate */
