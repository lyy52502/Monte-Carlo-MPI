#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#define R_ssa 15
#define X_LEN 7
#define T_FINAL 100.0


int P[R_ssa][X_LEN] = {
    { 1,  0,  0,  0,  0,  0,  0},
    {-1,  0,  0,  0,  0,  0,  0},
    {-1,  0,  1,  0,  0,  0,  0},
    { 0,  1,  0,  0,  0,  0,  0},
    { 0, -1,  0,  0,  0,  0,  0},
    { 0, -1,  0,  1,  0,  0,  0},
    { 0,  0, -1,  0,  0,  0,  0},
    { 0,  0, -1,  0,  1,  0,  0},
    { 0,  0,  0, -1,  0,  0,  0},
    { 0,  0,  0, -1,  0,  1,  0},
    { 0,  0,  0,  0, -1,  0,  0},
    { 0,  0,  0,  0, -1,  0,  1},
    { 0,  0,  0,  0,  0, -1,  0},
    { 1,  0,  0,  0,  0,  0, -1},
    { 0,  0,  0,  0,  0,  0, -1}
};

void prop(int *x, double *w) {
    const double LAMBDA_H = 20;
    const double LAMBDA_M = 0.5;
    const double B = 0.075;
    const double BETA_H = 0.3;
    const double BETA_M = 0.5;
    const double MU_H = 0.015;
    const double MU_M = 0.02;
    const double DELTA_H = 0.05;
    const double DELTA_M = 0.15;
    const double ALFA_H = 0.6;
    const double ALFA_M = 0.6;
    const double R = 0.05;
    const double OMEGA = 0.02;
    const double NU_H = 0.5;
    const double NU_M = 0.15;

    w[0]  = LAMBDA_H;
    w[1]  = MU_H * x[0];
    w[2]  = (B * BETA_H * x[0] * x[5]) / (1 + NU_H * x[5]);
    w[3]  = LAMBDA_M;
    w[4]  = MU_M * x[1];
    w[5]  = (B * BETA_M * x[1] * x[4]) / (1 + NU_M * x[4]);
    w[6]  = MU_H * x[2];
    w[7]  = ALFA_H * x[2];
    w[8]  = MU_M * x[3];
    w[9]  = ALFA_M * x[3];
    w[10] = (MU_H + DELTA_H) * x[4];
    w[11] = R * x[4];
    w[12] = (MU_M + DELTA_M) * x[5];
    w[13] = OMEGA * x[6];
    w[14] = MU_H * x[6];
}

void SSA(double T, int *x0, int *x_final) {
    double t = 0.0;
    int x[X_LEN];
    double w[R_ssa];
    double a0, tau;
    int i, r;
    double u1, u2, sum;

    for (i = 0; i < X_LEN; i++) x[i] = x0[i];

    while (t < T) {
        prop(x, w);
        a0 = 0.0;
        for (i = 0; i < R_ssa; i++) a0 += w[i];
        if (a0 == 0.0) break;

        u1 = (double) rand() / RAND_MAX;
        u2 = (double) rand() / RAND_MAX;
        tau = -log(u1) / a0;
        t += tau;
        if (t >= T) break;

        sum = 0.0;
        for (r = 0; r < R_ssa; r++) {
            sum += w[r];
            if (sum >= u2 * a0) break;
        }

        for (i = 0; i < X_LEN; i++)
            x[i] += P[r][i];
    }

    for (i = 0; i < X_LEN; i++)
        x_final[i] = x[i];
}

void build_histogram(int *values, int N, int bins) {
    int min_val = values[0];
    int max_val = values[0];
    for (int i = 1; i < N; i++) {
        if (values[i] < min_val) min_val = values[i];
        if (values[i] > max_val) max_val = values[i];
    }

    int *hist = calloc(bins, sizeof(int));
    double bin_width = (double)(max_val - min_val) / bins;

    for (int i = 0; i < N; i++) {
        int bin_index = (int)((values[i] - min_val) / bin_width);
        if (bin_index == bins) bin_index = bins - 1;
        hist[bin_index]++;
    }

    printf("Histogram bins and counts:\n");
    for (int i = 0; i < bins; i++) {
        double bin_start = min_val + i * bin_width;
        double bin_end = bin_start + bin_width;
        printf("Bin %2d: [%6.2f - %6.2f) Count: %d\n", i, bin_start, bin_end, hist[i]);
    }

    free(hist);
}
void run_simulation(double *result) {
    int x0[X_LEN] = {900, 900, 30, 330, 50, 270, 20};  
    int x_final[X_LEN];
    SSA(T_FINAL, x0, x_final);
    for (int i = 0; i < X_LEN; i++) {
        result[i] = (double)x_final[i];
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    srand(time(NULL) + rank);
    if (argc < 2) {
        if (rank == 0) {
            printf("Usage: %s <total_runs>\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    int total_runs = atoi(argv[1]);

    
    double start_time = MPI_Wtime();
    
    /*Here, I used the dynamic loading balancing, so each woker gets a new tasks
    only when it finishes the previous one, I learn this way from the chagpt and the course of HPP*/
    double *all_results = NULL;
    if (rank == 0) {
        all_results = malloc(total_runs * X_LEN * sizeof(double));
        if (all_results == NULL) {
            fprintf(stderr, "Master: Failed to allocate memory for results.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
     if (size == 1) {
        // If only one process, run all simulations sequentially on rank 0
        if (rank == 0) {
            printf("Running in single-process mode (rank 0 performs all %d simulations).\n", total_runs);
            for (int i = 0; i < total_runs; ++i) {
                double result[X_LEN];
                run_simulation(result);
                memcpy(&all_results[i * X_LEN], result, X_LEN * sizeof(double));
            }
        }
    } else {
        // Master-worker with dynamic load balancing
        if (rank == 0) {
            int num_sent = 0;
            int num_results = 0;

            // 1. Master sends initial tasks to all workers (excluding itself for now)
            for (int i = 1; i < size && num_sent < total_runs; i++) {
                MPI_Send(&num_sent, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                num_sent++;
            }

            // 2. Master continues to receive results and send new tasks.
            while (num_results < total_runs) {
                double result[X_LEN];
                MPI_Status status;

                // 3. Receive results from any worker
                MPI_Recv(result, X_LEN, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                int sender = status.MPI_SOURCE;
                memcpy(&all_results[num_results * X_LEN], result, X_LEN * sizeof(double));
                num_results++;

                // 4. Send next task if available
                if (num_sent < total_runs) {
                    MPI_Send(&num_sent, 1, MPI_INT, sender, 0, MPI_COMM_WORLD);
                    num_sent++;
                } else {
                    
                    int stop_signal = -1;
                    MPI_Send(&stop_signal, 1, MPI_INT, sender, 0, MPI_COMM_WORLD);
                }
            }

            // After all results are collected, send stop signals to any remaining idle workers
            for (int i = 1; i < size; ++i) {
                int stop_signal = -1;
                MPI_Send(&stop_signal, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }

        } else { 
            int task_index;
            while (1) {
                MPI_Recv(&task_index, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (task_index == -1) break;

                double result[X_LEN];
                run_simulation(result);
                MPI_Send(result, X_LEN, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
        }
    }
    if (rank == 0) {
        FILE *fp = fopen("histogram_data.csv", "w");
        if (fp != NULL) {
            for (int i = 0; i < total_runs; i++) {
                double *x = &all_results[i * X_LEN];
                fprintf(fp, "%.2f,%.2f,%.2f\n", x[0], x[1], x[2]);
            }
            fclose(fp);
        } else {
            printf("Error opening file for writing.\n");
        }

        int *S_values = malloc(total_runs * sizeof(int));
        for (int i = 0; i < total_runs; i++) {
            S_values[i] = (int) all_results[i * X_LEN];
        }
        build_histogram(S_values, total_runs, 20);
        free(S_values);
        free(all_results); 

        double end_time = MPI_Wtime();
        printf("Total simulation time: %.3f seconds\n", end_time - start_time);
    }

    MPI_Finalize();
    return 0;
}