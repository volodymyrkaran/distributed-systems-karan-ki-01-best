#include <iostream>
#include <cstdlib>
#include <time.h> 
#include <mpi.h>
using namespace std;

int main(int argc, char** argv) {
    double* x = new double[PointsNum];
    double* y = new double[PointsNum];

    int ProcRank, ProcNum;
    int PointsNum = 1000000, Count = 0;
    MPI_Status Status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);

    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    int* ReceivedPointsNum = new int[ProcNum];
    double* ReceivedPi = new double[ProcNum];
    double* ReceivedTime = new double[ProcNum];
    double* ReceivedStatus = new double[ProcNum];

    srand(time(NULL));
    if (ProcRank == 0){
        for (int i = 0; i < PointsNum; i++){
            x[i] = (double)rand() / RAND_MAX;
            y[i] = (double)rand() / RAND_MAX; 
        }
    }

    MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(y, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int k;
    if (ProcNum == 1) {
        k = PointsNum / ProcNum;
    }
    else {
        k = PointsNum / (ProcNum - 1);
    }

    int points_in_set = k * ProcRank;
    if (ProcRank == ProcNum - 1) points_in_set = N;

    if (ProcRank != 0) {
        double start = MPI_Wtime();

        for (int i = 0; i < points_in_set; i++) {
            double z = sqrt(x[i] * x[i] + y[i] * y[i]);
            if (z <= 1){
                Count++;
            }
        }
        double pi = 4 * ((double)Count / (double)points_in_set);

        double end = MPI_Wtime();
        double seconds = (end - start);

        MPI_Send(&pi, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        MPI_Send(&seconds, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
        MPI_Send(&points_in_set, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
    }
    else if (ProcRank == 0) {
        for (int i = 1; i < ProcNum; i++) {
            MPI_Recv(&ReceivedPi[i], 1, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&ReceivedTime[i], 1, MPI_DOUBLE, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&ReceivedPointsNum[i], 1, MPI_INT, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &Status);
            ReceivedStatus[i] = Status.MPI_SOURCE;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (ProcRank == 0) {
        for (int i = 1; i < ProcNum; i++) {
            cout << "status:" << ReceivedStatus[i] << " pi=" << ReceivedPi[i] << " time:" << ReceivedTime[i]
                << " pnumber:" << ReceivedPointsNum[i] << endl;
        }
    }
    MPI_Finalize();
    delete[] x;
    delete[] y;
    delete[] ReceivedPointsNum;
    delete[] ReceivedPi;
    delete[] ReceivedTime;
    delete[] ReceivedStatus;

    return 0;
}