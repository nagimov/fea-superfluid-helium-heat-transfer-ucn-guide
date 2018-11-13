#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "params.h"
#include "thermal1.h"
#include <pthread.h>
#define N_THREADS 1

double dt = 0.0000005;  // timestep of integration
const int t_steps = 200000000;  // number of integration steps
const int update_thermals = 1000;  // timesteps before thermal properties recalculated
const int update_output_file = 10000;  // timesteps before output file is updated on disk
const int update_log = 1;  // timesteps before output to log
const double he3_temp = 0.7;  // temperature of helium-3 to be set at each iteration
const double t_eps_max = 0.01;  // maximal allowed divergence of F*dT
double one_third = 1.0/3.0;
struct mesh_op {
    int x_current;
    int y_current;
    int x_proc;
    int y_proc;
    double dL;
    double dA;
    int op_type;
};
struct mesh_op mesh_ops[OPS_MESH];

double mesh_T[L_MESH][H_MESH];
double (*mesh_T_p)[L_MESH][H_MESH];  // pointer to array of doubles, if without brackets - array of pointers!
double mesh_DT_0[L_MESH][H_MESH];
double mesh_DT[L_MESH][H_MESH];
double mesh_DT_prev_iter[L_MESH][H_MESH];
double mesh_DT_prev_time[L_MESH][H_MESH];
double mesh_DT_temp[L_MESH][H_MESH];

double mesh_v[L_MESH][H_MESH];
double mesh_m[L_MESH][H_MESH];
double mesh_Q[L_MESH][H_MESH];
double mesh_C[L_MESH][H_MESH];

double mesh_F[L_MESH][H_MESH];
double mesh_F_prev[L_MESH][H_MESH];
double mesh_F_prev_time[L_MESH][H_MESH];

double mesh_FS[L_MESH][H_MESH];
double mesh_FS_prev[L_MESH][H_MESH];
double mesh_FS_prev_time[L_MESH][H_MESH];

double mesh_DF[L_MESH][H_MESH];
double mesh_DFS[L_MESH][H_MESH];

double mesh_F_temp[N_THREADS][L_MESH][H_MESH];

double he4_k_kap[L_MESH][H_MESH];
double he4_f_inv[L_MESH][H_MESH];
double he4_C[L_MESH][H_MESH];
double he3_k_kap[L_MESH][H_MESH];
double cu_C[L_MESH][H_MESH];

double q_deposit, q_he4_cu, q_cu_he3;

double T_eps;

// imported from thermal.h
// const double d_cu = ...
// const double cp_cu = ...

// structure for bounds for array crunching
struct bounds {
    int lo;
    int hi;
    int thread_num;
};
struct bounds B[N_THREADS];

int SYNC_THREADS[N_THREADS];
int END_THREADS = 0;
int FS;

void LoadCSV(char filename[256], double data[L_MESH][H_MESH]) {
    int rowIndex = 0;
    char line[64000];
    char* token = NULL;
    FILE* fp = fopen(filename, "r");
    while (fgets( line, sizeof(line), fp) != NULL && rowIndex < H_MESH) {
        int colIndex = 0;
        for (token = strtok( line, ","); token != NULL && colIndex < L_MESH; token = strtok(NULL, ",")) {
            data[colIndex++][rowIndex] = atof(token);
        }
        rowIndex++;
    }
    fclose(fp);
    printf("file %s loaded, first lines content:\n", filename);
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 7; ++j)
            printf("%10.6lf", data[i][j]);
        printf(" ... ");
        putchar('\n');
    }
}

void LoadOpsCSV(char filename[256], struct mesh_op mesh_ops[OPS_MESH]) {
    int rowIndex = 0;
    char line[64000];
    char* token = NULL;
    FILE* fp = fopen(filename, "r");
    while (fgets( line, sizeof(line), fp) != NULL && rowIndex < OPS_MESH) {
        token = strtok( line, ",");
        mesh_ops[rowIndex].x_current = atoi(token);  // x_current
        token = strtok( NULL, ",");
        mesh_ops[rowIndex].y_current = atoi(token);  // y_current
        token = strtok( NULL, ",");
        mesh_ops[rowIndex].x_proc = atoi(token);  // x_proc
        token = strtok( NULL, ",");
        mesh_ops[rowIndex].y_proc = atoi(token);  // y_proc
        token = strtok( NULL, ",");
        mesh_ops[rowIndex].dL = atof(token);  // dL
        token = strtok( NULL, ",");
        mesh_ops[rowIndex].dA = atof(token);  // dA
        token = strtok( NULL, ",");
        mesh_ops[rowIndex].op_type = atoi(token);  // op_type
        rowIndex++;
    }
    fclose(fp);
    printf("file %s loaded, first lines content:\n", filename);
    for (int i = 0; i < 2; ++i) {
            printf("%3d", mesh_ops[i].x_current);
            printf("%3d", mesh_ops[i].y_current);
            printf("%3d", mesh_ops[i].x_proc);
            printf("%3d", mesh_ops[i].y_proc);
            printf("%12.8lf", mesh_ops[i].dL);
            printf("%12.8lf", mesh_ops[i].dA);
            printf("%3d", mesh_ops[i].op_type);
        putchar('\n');
    }
}

void SaveCSV(char filename[256], double data[L_MESH][H_MESH]) {
    FILE *fp;
    int i, j;
    fp = fopen(filename, "w+");
    for (j = 0; j < H_MESH; j++) {
        fprintf(fp, "%16.12lf", data[0][j]);
        for(i = 1; i < L_MESH; i++) {
            fprintf(fp, ",%16.12lf", data[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void process_cu_cu(int thread_num, int op_num) {  // case 22
	struct mesh_op mesh_op_current = mesh_ops[op_num];
    int x_current = mesh_op_current.x_current;
    int y_current = mesh_op_current.y_current;
    int x_proc = mesh_op_current.x_proc;
    int y_proc = mesh_op_current.y_proc;
    double dA = mesh_op_current.dA;
    double dL = mesh_op_current.dL;
    double T1 = (*mesh_T_p)[x_current][y_current];
    double T2 = (*mesh_T_p)[x_proc][y_proc];
    double dq = k_cu * dA / dL * (T1 - T2);  // flux from 1 to 2 (1 is warmer)
    mesh_F_temp[thread_num][x_current][y_current] -= dq;
    mesh_F_temp[thread_num][x_proc][y_proc] += dq;
}

void set_temp_3he() {  // case 33
    // constant temperature of helium-3 is assumed
    // therefore, no thermal conductivity within helium-3
    struct mesh_op mesh_op_current;
    int x_current, y_current, x_proc, y_proc;
    int h;
    for (h = 0; h < OPS_MESH; ++h) {
        mesh_op_current = mesh_ops[h];
        if (mesh_op_current.op_type == 33) {
            x_current = mesh_op_current.x_current;
            y_current = mesh_op_current.y_current;
            x_proc = mesh_op_current.x_proc;
            y_proc = mesh_op_current.y_proc;
            mesh_T[x_current][y_current] = he3_temp;
            mesh_T[x_proc][y_proc] = he3_temp;
        }
    }
}

void process_4he_4he(int thread_num, int op_num) {  // case 44
	struct mesh_op mesh_op_current = mesh_ops[op_num];
    int x_current = mesh_op_current.x_current;
    int y_current = mesh_op_current.y_current;
    int x_proc = mesh_op_current.x_proc;
    int y_proc = mesh_op_current.y_proc;
    double dA = mesh_op_current.dA;
    double dL = mesh_op_current.dL;
    double T1 = (*mesh_T_p)[x_current][y_current];
    double T2 = (*mesh_T_p)[x_proc][y_proc];
    double f_inv_avg = 0.5 * (he4_f_inv[x_current][y_current] + he4_f_inv[x_proc][y_proc]);
    double dq_cubed = dA*dA*dA * f_inv_avg * (T1 - T2) / dL;  // flux from 1 to 2 (1 is warmer)
    //dq_cubed = 0;
    mesh_F_temp[thread_num][x_current][y_current] -= dq_cubed;
    mesh_F_temp[thread_num][x_proc][y_proc] += dq_cubed;
}

void process_cu_4he(int thread_num, int op_num) {  // case 24
	struct mesh_op mesh_op_current = mesh_ops[op_num];
    int x_current = mesh_op_current.x_current;
    int y_current = mesh_op_current.y_current;
    int x_proc = mesh_op_current.x_proc;
    int y_proc = mesh_op_current.y_proc;
    double dA = mesh_op_current.dA;
    double T_cu = (*mesh_T_p)[x_current][y_current];
    double T_he = (*mesh_T_p)[x_proc][y_proc];
    double k_kap = he4_k_kap[x_proc][y_proc];
    double dq = k_kap * dA * (T_cu - T_he);  // flux from 1 to 2 (1 is warmer copper)
    //dq = 0;
    mesh_F_temp[thread_num][x_current][y_current] -= dq;
    mesh_F_temp[thread_num][x_proc][y_proc] += dq;
    q_he4_cu -= dq;
}

void process_4he_cu(int thread_num, int op_num) {  // case 42
	struct mesh_op mesh_op_current = mesh_ops[op_num];
    int x_current = mesh_op_current.x_current;
    int y_current = mesh_op_current.y_current;
    int x_proc = mesh_op_current.x_proc;
    int y_proc = mesh_op_current.y_proc;
    double dA = mesh_op_current.dA;
    double T_he = (*mesh_T_p)[x_current][y_current];
    double T_cu = (*mesh_T_p)[x_proc][y_proc];
    double k_kap = he4_k_kap[x_current][y_current];
    double dq = k_kap * dA * (T_he - T_cu);  // flux from 1 to 2 (1 is warmer helium)
    //dq = 0;
    mesh_F_temp[thread_num][x_current][y_current] -= dq;
    mesh_F_temp[thread_num][x_proc][y_proc] += dq;
    q_he4_cu += dq;
}

void process_cu_3he(int thread_num, int op_num) {  // case 23
	struct mesh_op mesh_op_current = mesh_ops[op_num];
    int x_current = mesh_op_current.x_current;
    int y_current = mesh_op_current.y_current;
    int x_proc = mesh_op_current.x_proc;
    int y_proc = mesh_op_current.y_proc;
    double dA = mesh_op_current.dA;
    double T_cu = (*mesh_T_p)[x_current][y_current];
    double T_he = (*mesh_T_p)[x_proc][y_proc];
    double k_kap = he3_k_kap[x_proc][y_proc];
    double dq = k_kap * dA * (T_cu - T_he);  // flux from 1 to 2 (1 is warmer copper)
    mesh_F_temp[thread_num][x_current][y_current] -= dq;
    mesh_F_temp[thread_num][x_proc][y_proc] += dq;
    q_cu_he3 += dq;
}

void process_3he_cu(int thread_num, int op_num) {  // case 32
	struct mesh_op mesh_op_current = mesh_ops[op_num];
    int x_current = mesh_op_current.x_current;
    int y_current = mesh_op_current.y_current;
    int x_proc = mesh_op_current.x_proc;
    int y_proc = mesh_op_current.y_proc;
    double dA = mesh_op_current.dA;
    double T_he = (*mesh_T_p)[x_current][y_current];
    double k_kap = he3_k_kap[x_current][y_current];
    double T_cu = (*mesh_T_p)[x_proc][y_proc];
    double dq = k_kap * dA * (T_he - T_cu);  // flux from 1 to 2 (1 is warmer helium)
    mesh_F_temp[thread_num][x_current][y_current] -= dq;
    mesh_F_temp[thread_num][x_proc][y_proc] += dq;
    q_cu_he3 -= dq;
}

void zero_threaded_mesh(double mesh[N_THREADS][L_MESH][H_MESH]) {
	int th, i, j;
	for (th = 0; th < N_THREADS; ++th) {
		for (i = 0; i < L_MESH; ++i) {
			for (j = 0; j < H_MESH; ++j) {
				mesh[th][i][j] = 0.0;
			}
		}
	}
}

void zero_unthreaded_mesh(double mesh[L_MESH][H_MESH]) {
	int i, j;
	for (i = 0; i < L_MESH; ++i) {
		for (j = 0; j < H_MESH; ++j) {
			mesh[i][j] = 0.0;
		}
	}
}

void join_threaded_mesh(double threaded[N_THREADS][L_MESH][H_MESH], double joined[L_MESH][H_MESH]) {
    zero_unthreaded_mesh(joined);
	double joiner;
	int th, i, j;
	for (i = 0; i < L_MESH; ++i) {
		for (j = 0; j < H_MESH; ++j) {
			joiner = 0;
			for (th = 0; th < N_THREADS; ++th) {
				joiner += threaded[th][i][j];
			}
			joined[i][j] = joiner;
		}
	}
}

void plus_equal_mesh(double mesh[L_MESH][H_MESH], double plus[L_MESH][H_MESH]) {
    int i, j;
    for (i = 0; i < L_MESH; ++i) {
        for (j = 0; j < H_MESH; ++j) {
            mesh[i][j] += plus[i][j];
        }
    }
}

double mesh_eps(double mesh1[L_MESH][H_MESH], double mesh2[L_MESH][H_MESH]) {
	double eps = 0.0;
	int i, j;
	for (i = 0; i < L_MESH; ++i) {
		for (j = 0; j < H_MESH; ++j) {
		    if ((mesh_m[i][j] != 1) && (mesh_m[i][j] != 3)) {
                eps += fabs(mesh1[i][j] - mesh2[i][j]);
            }
		}
	}
    return eps;
}

void assign_mesh(double mesh1[L_MESH][H_MESH], double mesh2[L_MESH][H_MESH]) {
	int i, j;
	for (i = 0; i < L_MESH; ++i) {
		for (j = 0; j < H_MESH; ++j) {
            mesh1[i][j] = mesh2[i][j];
		}
	}
}




static void *crunching(void *arg) {
    // multiple threads run mesh_ops operations according to assigned low/high bounds
    struct bounds *data = (struct bounds *)arg;
    int lo = (*data).lo;
    int hi = (*data).hi;
    int thread_num = (*data).thread_num;
    printf("worker %d started for bounds [%d %d] \n", thread_num, lo, hi); //todo remove

    int i;

    while (END_THREADS != 1) {  // END_THREADS tells threads to terminate
        if (SYNC_THREADS[thread_num] == 1) {  // SYNC_THREADS allows threads to start crunching
            //printf("worker %d working... \n", thread_num );
            for (i = lo; i <= hi; ++i) {
                switch(mesh_ops[i].op_type) {
                    case 22:
                        if (FS == 0) process_cu_cu(thread_num, i);
                        break;
                    /* process_3he_3he resets helium-3 temperatures directly bypassing mesh_T_diff matrix,
                    so it has to called outside of threads */
                    //    process_3he_3he(thread_num, mesh_ops[i], dt, mesh_T);
                    //case 33:
                    //    break;
                    case 44:
                        if (FS == 1) process_4he_4he(thread_num, i);
                        break;
                    case 24:
                        if (FS == 0) process_cu_4he(thread_num, i);
                        break;
                    case 42:
                        if (FS == 0) process_4he_cu(thread_num, i);
                        break;
                    case 23:
                        if (FS == 0) process_cu_3he(thread_num, i);
                        break;
                    case 32:
                        if (FS == 0) process_3he_cu(thread_num, i);
                        break;
                }
            }
            SYNC_THREADS[thread_num] = 0;  // thread disables itself until SYNC_THREADS is back to 1
            //printf("worker %d stopped... \n", thread_num );
        }
    }
    return 0;
}

void calc_DT() {  // calculate mesh_DT_temp
// DT = Q + F + DF + (FS + DFS) ^ 1/3)
    int i, j;
    int istherenan = 0;
    double sign = 1.0;
    double unpow, powed;
    zero_unthreaded_mesh(mesh_DT_temp);
    for (i = 0; i < L_MESH; ++i) {
        for (j = 0; j < H_MESH; ++j) {
            if ((mesh_m[i][j] != 1) && (mesh_m[i][j] != 3)) {
                unpow = 0.5 * mesh_FS[i][j] + 0.5 * mesh_FS_prev_time[i][j] - mesh_DFS[i][j];
                if (unpow < 0) {
                    powed = -pow(-unpow, one_third);
                } else {
                    powed = pow(unpow, one_third);
                }
                mesh_DT_temp[i][j] = (mesh_Q[i][j] + (0.5 * mesh_F[i][j] + 0.5 * mesh_F_prev_time[i][j] - mesh_DF[i][j]) + powed ) * dt / mesh_C[i][j];
                if (mesh_DT_temp[i][j] > 1) {
                    printf("%d %d too big 0! \n", i, j);
                    SaveCSV("mesh_T_post.csv", mesh_T);
                }
                if isnan(mesh_DT_temp[i][j]) istherenan = 1;
            }
        }
    }
    if (istherenan) printf("there was a nan in calc_DT!!!");
}





int main() {
    int t_steps_digits = floor (log10 (abs (t_steps))) + 1; // for fancy logging

    // read input files
    LoadOpsCSV("mesh_ops.csv", mesh_ops);
    LoadCSV("mesh_v.csv", mesh_v);
    LoadCSV("mesh_m.csv", mesh_m);
    LoadCSV("mesh_Q.csv", mesh_Q);
    LoadCSV("mesh_T.csv", mesh_T);

    int h, i, j, nan, th, relaxation_steps;
    int m_current;
    int t = 0;
    double T_current, T_max;
    double t_diff;
    double he4_cp, he4_d, v;
    double q_join;

    // calculate total heat deposition
    q_deposit = 0;
    for (i = 0; i < L_MESH; ++i) {
    	for (j = 0; j < H_MESH; ++j) {
        	q_deposit += mesh_Q[i][j];
    	}
    }
    printf("q_deposit = %f; \n", q_deposit);

    // zero matrices
    zero_unthreaded_mesh(mesh_C);

    // set initial constant thermals for copper, since they are constant
    for (i = 0; i < L_MESH; ++i) {
        for (j = 0; j < H_MESH; ++j) {
            v = mesh_v[i][j];
            m_current = mesh_m[i][j];
            if (m_current == 2) {
        		mesh_C[i][j] = v * cp_cu * d_cu;
            }
        }
    }

    // split workload between workers
    int worker_length = OPS_MESH / N_THREADS;
    for (i = 0; i < N_THREADS; ++i) {
        B[i].thread_num = i;
        B[i].lo = i * worker_length;
        if (i == N_THREADS - 1) {
            B[i].hi = OPS_MESH;
        } else {
            B[i].hi = i * worker_length + worker_length - 1;
        }
    }
    // pointer to parameters to be passed to worker
    struct bounds **data = malloc(N_THREADS * sizeof(struct bounds*));
    for (i = 0; i < N_THREADS; i++) {
        data[i] = malloc(sizeof(struct bounds));
        data[i]->lo = B[i].lo;
        data[i]->hi = B[i].hi;
        data[i]->thread_num = B[i].thread_num;
    }
    // create thread objects
    pthread_t threads[N_THREADS];

    nan = 0;  // flag defining whether there are nans in mesh_T

    // disallow thread to process heat transfer
    for (th = 0; th < N_THREADS; ++th) {
        SYNC_THREADS[th] = 0;
    }
    // launch workers
    for(th = 0; th < N_THREADS; th++) {
        pthread_create(&threads[th], NULL, crunching, data[th]);
    }

    zero_unthreaded_mesh(mesh_F);
    zero_unthreaded_mesh(mesh_FS);
    zero_unthreaded_mesh(mesh_DT);

    assign_mesh(mesh_DT_prev_time, mesh_DT);  // zero'ed

    // main iteration loop
    for (t = 0; t < t_steps; ++t) {
        q_he4_cu = 0;
        q_cu_he3 = 0;
        // if enough timesteps elapsed, update thermal properties
        if (t % update_thermals == 0) {
            printf("updating thermal properties...\n");
            for (i = 0; i < L_MESH; ++i) {
                for (j = 0; j < H_MESH; ++j) {
                    T_current = mesh_T[i][j];
                    v = mesh_v[i][j];
                    m_current = mesh_m[i][j];
                    switch(m_current) {
                    	case 3:
                    		he3_k_kap[i][j] = k_kap_3he(T_current);
                    		break;
                    	case 4:
                    		he4_cp = cp_4he(T_current);
                    		he4_d = d_4he(T_current);
                    		mesh_C[i][j] = v * he4_cp * he4_d; // maybe /dt
                    		he4_k_kap[i][j] = k_kap_4he(T_current);
                    		he4_f_inv[i][j] = f_inv(T_current);
                    		break;
                    }
                }
            }
        }
        // if enough timesteps elapsed, save output to disk
        if ((t + 1) % update_output_file == 0) {  // t + 1 here so t = 0 never triggers
            printf("saving output file...\n");
            SaveCSV("mesh_T_post.csv", mesh_T);
        }
        // reset he-3 temperature
        set_temp_3he();

        // heat transfer

// F = K*T (all except helium-4)  --  mesh_F
        FS = 0;
        // set pointers to temp and heat flux
        mesh_T_p = &mesh_T;

        zero_threaded_mesh(mesh_F_temp);
        // allow thread to process heat transfer and wait for threads to complete
        for (th = 0; th < N_THREADS; ++th) { SYNC_THREADS[th] = 1; }
                // ...here threads started by pthread_create calculate heat transfer...
        for (th = 0; th < N_THREADS; th++) { while (SYNC_THREADS[th] != 0) {} }

    	// join heat fluxes
    	assign_mesh(mesh_F_prev_time, mesh_F);
    	join_threaded_mesh(mesh_F_temp, mesh_F);

// FS = KS*T (helium-4 only)  --  mesh_FS
        FS = 1;
        // set pointers to temp and heat flux
        mesh_T_p = &mesh_T;

        zero_threaded_mesh(mesh_F_temp);
        // allow thread to process heat transfer and wait for threads to complete
        for (th = 0; th < N_THREADS; ++th) { SYNC_THREADS[th] = 1; }
                // ...here threads started by pthread_create calculate heat transfer...
        for (th = 0; th < N_THREADS; th++) { while (SYNC_THREADS[th] != 0) {} }

    	// join heat fluxes
    	assign_mesh(mesh_FS_prev_time, mesh_FS);
    	join_threaded_mesh(mesh_F_temp, mesh_FS);

// DT = Q + F + DF + (FS + DFS) ^ 1/3)
// initially DF = 0 and DFS = 0

        zero_unthreaded_mesh(mesh_DFS);
        zero_unthreaded_mesh(mesh_DF);

    	// calculate mesh_DT_0
        calc_DT();
        assign_mesh(mesh_DT_0, mesh_DT_temp);

        T_eps = 200.0;
        assign_mesh(mesh_DT_prev_iter, mesh_DT_0);
        relaxation_steps = 0;
        while (T_eps > t_eps_max) {
            relaxation_steps += 1;

    // DF = K*DT (all except helium-4)  --  mesh_DF
            FS = 0;
            // set pointers
            mesh_T_p = &mesh_DT_prev_iter;

            // allow thread to process heat transfer and wait for threads to complete
            for (th = 0; th < N_THREADS; ++th) { SYNC_THREADS[th] = 1; }
                    // ...here threads started by pthread_create calculate heat transfer...
            for (th = 0; th < N_THREADS; ++th) { while (SYNC_THREADS[th] != 0) {} }

            // join heat fluxes
            join_threaded_mesh(mesh_F_temp, mesh_DF);

    // DFS = KS*DT (helium-4 only)  --  mesh_DFS
            FS = 1;
            // set pointers
            mesh_T_p = &mesh_DT_prev_iter;
            // zero mesh_F
            zero_threaded_mesh(mesh_F_temp);

            // allow thread to process heat transfer and wait for threads to complete
            for (th = 0; th < N_THREADS; ++th) { SYNC_THREADS[th] = 1; }
                    // ...here threads started by pthread_create calculate heat transfer...
            for (th = 0; th < N_THREADS; ++th) { while (SYNC_THREADS[th] != 0) {} }

            // join heat fluxes
            join_threaded_mesh(mesh_F_temp, mesh_DFS);


    // DT = Q + F + DF + (FS + DFS) ^ 1/3)
    // here DF and DFS are not zero

            // recalculate mesh_t_step
            calc_DT();
            assign_mesh(mesh_DT, mesh_DT_temp);

            T_eps = mesh_eps(mesh_DT, mesh_DT_prev_iter);  // todo
            assign_mesh(mesh_DT_prev_iter, mesh_DT);
            //printf("relaxation... T_eps = %f \n", T_eps);
        }

        //double mesh_T_max = 0;
        //int mesh_T_max_i, mesh_T_max_j;
        for (i = 0; i < L_MESH; ++i) {
            for (j = 0; j < H_MESH; ++j) {
                mesh_T[i][j] += 1.0 * mesh_DT[i][j];// + 0.5 * mesh_DT_prev_time[i][j];
                //if (mesh_DT[i][j] > mesh_T_max) {
                //    mesh_T_max = mesh_DT[i][j];
                //    mesh_T_max_i = i;
                //    mesh_T_max_j = j;
                //}
            }
        }

        assign_mesh(mesh_DT_prev_time, mesh_DT);

        // find max temp and nans, if present
        T_max = mesh_T[0][0];
        for (i = 0; i < L_MESH; ++i) {
            for (j = 0; j < H_MESH; ++j) {
                T_current = mesh_T[i][j];
                if (isnan(T_current)) {
                    T_max = 1.0;
                } else {
                    if (T_current > T_max) {
                        T_max = T_current;
                    }
                }
            }
        }
        // logging
        if (t % update_log == 0) {
            fflush(stdout);  // need to flush output, so it can be read from python
            // otherwise windows buffers the output and subprocess call hangs
            printf("step %*d out of %*d; ", t_steps_digits, t, t_steps_digits, t_steps);
            printf("T_max = %8.6lf; ", T_max);
            printf("q_he4_cu = %8.6lf; ", q_he4_cu);
            printf("q_cu_he3 = %8.6lf; ", q_cu_he3);
            printf("T_eps = %8.6lf; ", T_eps);
            printf("dt = %10.8lf; ", dt);
            printf("rel.steps = %d; ", relaxation_steps);
            putchar('\n');
        }
    }
    // end main loop
    // join workers
    END_THREADS = 1;
    for(th = 0; th < N_THREADS; th++) {
        pthread_join(threads[th], NULL);
    }
    // write output file
    SaveCSV("mesh_T_post.csv", mesh_T);
    return 0;
}
