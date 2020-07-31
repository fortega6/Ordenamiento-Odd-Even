//mpicc odd_even.c -o oe
//mpirun -np 4 ./oe <entrada.in

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int cmpfunc (const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}

int odd(int);
void MAX(int *my, int *recv, int *temp, int local_n);
void MIN(int *my, int *recv, int *temp, int local_n);
void asignacion(int *L, int *temp, int local_n);

int main(int argc, char* argv[])
{
	double t1, t2;  //Tiempos
	int *G;		// Vector globar
	int proc_id;    // id del proceso
	int nproc;      // numero de procesos
	int n;          // tamaño del vector
    int dimension;  // cantidad real de elementos a ordenar
    int i;
								// recv: vector que recive el proceso						
	int *temp, *L, *recv; // L: es el vector local de cada proceso
    int local_n;				// Tamaño del vector local
    MPI_Status status;

    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &nproc );
    MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );
    	
	if(proc_id==0) 
	{
		/*scanf("%d", &n); G = malloc(n * sizeof(int)); G[0] = 6; printf("%d\n", G[0]);*/
		int factor_correccion;
		int n_parche; //Se define un nuevo N para hacer una division exacta, para los casos inpares
		scanf("%d", &n);
		printf("%d\n", n);
		dimension = n;
		factor_correccion = (n % nproc);
		
		if(factor_correccion > 0)
		{
			n_parche = n + (nproc - factor_correccion);
			G = malloc(n_parche * sizeof(int));
			printf("PARCHE %d %d\n",n_parche,(n % nproc));
			for(i = n; i < n_parche;i++)
			{
				G[i] = n_parche;
			}
			
			for(i = 0; i < n; i++)
			{
				scanf("%d", &G[i]);
			}
			n = n_parche;
		}
		else
		{
			G = malloc(n * sizeof(int));
			for(i = 0; i < n; i++)
			{
				scanf("%d", &G[i]);
			}
		}
	}	// El proceso 0 hace la lectura inicial
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);	// Se les pasa el vector glogal a los demas procesos
	MPI_Barrier(MPI_COMM_WORLD);
	if(proc_id != 0) G = malloc(n * sizeof(int));
	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Bcast(G, n, MPI_INT, 0, MPI_COMM_WORLD);	// Se les pasa el vector glogal a los demas procesos
	local_n = n/nproc;    // Se divide el tamaño del vector entre el numero de procesos
	t1 = MPI_Wtime();
	
	L = malloc(local_n * sizeof(int));
	temp = malloc(local_n * sizeof(int));
	recv = malloc(local_n * sizeof(int));
	
	MPI_Scatter (G, local_n, MPI_INT, L, local_n, MPI_INT, 0, MPI_COMM_WORLD);
	
	qsort(L, local_n, sizeof(int), cmpfunc); // Ordenar vector local de cada proceso
	
	for(i = 0; i < nproc/2+1; i++)
    {
        if(odd(proc_id))
        {
            if(proc_id < (nproc-1))
            {
                MPI_Recv( recv, local_n, MPI_INTEGER, proc_id+1, 1, MPI_COMM_WORLD, &status);
                MAX(L, recv, temp, local_n);
                MPI_Send( temp, local_n, MPI_INTEGER, proc_id+1, 1, MPI_COMM_WORLD);
                MIN(L, recv, temp, local_n);
                asignacion(L, temp, local_n); // A <- Min(temp, A)
            }
            MPI_Send( L, local_n, MPI_INTEGER, proc_id-1, 2, MPI_COMM_WORLD);
            MPI_Recv( L, local_n, MPI_INTEGER, proc_id-1, 2, MPI_COMM_WORLD, &status);
        }
        else
        {
            if(proc_id > 0)
            {
                MPI_Send( L, local_n, MPI_INTEGER, proc_id-1, 1, MPI_COMM_WORLD);
                MPI_Recv( L, local_n, MPI_INTEGER, proc_id-1, 1, MPI_COMM_WORLD, &status);
            }

            if(proc_id < (nproc-1))
            {
                MPI_Recv( recv, local_n, MPI_INTEGER, proc_id+1, 2, MPI_COMM_WORLD, &status);
                MAX(L, recv, temp, local_n);
                MPI_Send( temp, local_n, MPI_INTEGER, proc_id+1, 2, MPI_COMM_WORLD);
                MIN(L, recv, temp, local_n);
                asignacion(L, temp, local_n); // A <- Min(temp, A)
            }
        }
    }
    
    MPI_Gather (L, local_n, MPI_INT, G, local_n, MPI_INT, 0, MPI_COMM_WORLD);
    
	t2 = MPI_Wtime();
	
	/*if(proc_id==0)
    {
		printf("\nVector ordenado\n");
		for(i = 0; i < dimension; i++)
		{
			printf("%d ", G[i]);
		}
		printf("\n");
	}*/
	if(proc_id==0)
		printf("Tiempo %lf\n", t2-t1);
		
	MPI_Finalize();
}

int odd(int proc_id)
{
    return(proc_id%2 != 0);
}

void MAX(int *my, int *recv, int *temp, int local_n)
{
    int m_i, r_i, t_i;
	m_i = r_i = t_i = local_n-1;
	
    for(t_i = local_n-1; t_i >= 0; t_i--)
    {
        if(my[m_i] >= recv[r_i])
        {
            temp[t_i] = my[m_i];
            m_i--;
        }
        else
        {
            temp[t_i] = recv[r_i];
            r_i--;
        }
    }
}

void MIN(int *my, int *recv, int *temp, int local_n)
{
    int m_i, r_i, t_i;
	m_i = r_i = t_i = 0;

    for(t_i = 0; t_i < local_n; t_i++)
    {
        if(my[m_i] <= recv[r_i])
        {
            temp[t_i] = my[m_i];
            m_i++;
        }
        else
        {
            temp[t_i] = recv[r_i];
            r_i++;
        }
    }
}

void asignacion(int *L, int *temp, int local_n)
{
    int i;
    for(i = 0; i < local_n; i++)
    {
        L[i] = temp[i];
    }
}
