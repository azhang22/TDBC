//+++++++++++++++++++++++++++filename: main.cpp +++++++++++++++++++++++++++++++//
#include <stdlib.h>
#include <string.h>
#include "airpore.h"
#include "math.h"
#include "mpi.h"
#include <time.h>

int main(int argc, char *argv[])
{
	int mpi_rank,mpi_size;
	clock_t begin,end;
	double time_spent;
	begin = clock();
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	char infile[100]=" ";
	strcpy_s(infile,"input.txt");
	airpore *AirPore;
printf("%d,%d\n",mpi_rank,1);
	AirPore=new airpore(infile);
printf("%d,%d\n",mpi_rank,2);
//	goto loop;
	AirPore->get_output(mpi_rank,mpi_size);
printf("%d,%d\n",mpi_rank,3);
//	loop:
	if (mpi_rank==0)AirPore->get_FFT_y1(2);// for the long-distance sound propagatio
	delete AirPore;
	MPI_Finalize();
		
	end = clock();
	time_spent = ((double) (end - begin)) / CLOCKS_PER_SEC;
	if (mpi_rank == 0) printf("Total time used (wall time) : %f sec \n",time_spent);

	return 0;
}
