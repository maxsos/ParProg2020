//
// Created by - on 03.10.2020.
//
#include <vector>
#include <fstream>
#include <omp.h>

using namespace std;

int create_threads(int threads_count) {
    int res = 0;
    #pragma omp parallel num_threads(threads_count)
    {
        res++;
    }
    return res;
}


int main() {
    omp_set_nested(1);
    printf("[NESTED] %d\n", omp_get_nested());

    const int NUM_THREADS_ = 20;
    const int NUM_CREATED_THREADS_ = 20;

    double start, end, begin, total;
    vector<vector<double>> data;
    vector<double> temp;

    for(int num_thread = 1; num_thread < NUM_THREADS_; ++num_thread) {
        temp.clear();
        begin = omp_get_wtime();
        for (int num_created_thread = 1; num_created_thread < NUM_CREATED_THREADS_; ++num_created_thread) {
            start = omp_get_wtime();
            int res = 0;
            #pragma omp parallel num_threads(num_thread)
            {
                for (int i = 0; i < 1000; ++i) {
                    res += create_threads(num_created_thread);
                }
            }
            end = omp_get_wtime();
            temp.push_back(end - start);
        }
        data.push_back(temp);
        total += (end - begin)/60;
        printf("finished for %d for %lf min\n", num_thread, (end - begin)/60);
        printf("total %lf minutes\n", total);
    }

    printf("Completed calculation\n");

    ofstream of("D:\\CLionProjects\\parrprog\\ParProg2020\\OpenMP1\\data3.txt");

    for (int i = 0; i < NUM_THREADS_ - 1; ++i) {
        for (int j = 0; j < NUM_CREATED_THREADS_ - 1; ++j) {
            of << data[i][j];
            if (j != NUM_CREATED_THREADS_ - 2)
                of << " ";
        }
        of << endl;
    }

    of.close();
    printf("FINISHED\n");
    return 0;
}