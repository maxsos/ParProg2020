#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

double acceleration(double t)
{
  return sin(t);
}

void calc(double* trace, uint32_t traceSize, double t0, double dt, double y0, double y1, int rank, int size)
{
 
  MPI_Bcast(&traceSize, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
  MPI_Bcast(&t0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  size_t begin_int = traceSize / size * rank + ((traceSize % size > rank) ? rank : (traceSize % size));     
  size_t end_int = traceSize / size * (rank + 1) + ((traceSize % size > rank) ? rank + 1 : traceSize % size);
  size_t len = end_int - begin_int;

  double* data = new double[len];
  double start_v = 0, start_y = y0, end_v, end_y, b_star;

  t0 = t0 + dt * begin_int;
  
  // Sighting shot
  data[0] = start_y;
  data[1] = start_y + dt*start_v;
  for (uint32_t i = 2; i < len; i++)
  {
    data[i] = dt*dt*acceleration(t0 + (i - 1)*dt) + 2*data[i - 1] - data[i - 2];
  }

  /// Syncronize the information
  end_y = data[len - 1];
  end_v = (data[len - 1] - data[len - 2])/dt;

  if (size > 1) {
    double v = 0;
    if (rank > 0) {
      double prev_y, prev_v;
      MPI_Recv(&prev_y, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, NULL);
      MPI_Recv(&prev_v, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, NULL);

      start_y = prev_y;
      start_v = prev_v;

      end_y += prev_y + prev_v * dt * len;
      end_v += prev_v;
    }

    if (rank != size - 1) {
      MPI_Send(&end_y, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
      MPI_Send(&end_v, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
    }

    if (rank == size - 1) {
      MPI_Send(&end_y, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    if (rank == 0) {
      MPI_Recv(&b_star, 1, MPI_DOUBLE, size - 1, 0, MPI_COMM_WORLD, NULL);
      v = (y1 - b_star) / (dt * traceSize);
    }

    MPI_Bcast(&v, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    start_v += v;
    start_y += v * dt * begin_int;
  } else
  {
    start_v = (y1 - data[len - 1]) / (traceSize * dt);
  }



  data[0] = start_y;
  data[1] = start_y + dt * start_v;
  for (uint32_t i = 2; i < len; i++)
  {
    data[i] = dt*dt*acceleration(t0 + (i - 1)*dt) + 2*data[i - 1] - data[i - 2];
  }

  

  if (rank == 0) {
    int* rcounts = new int[size];
    int* displs = new int[size];
    for (int i = 0; i < size; i++) {
      begin_int = traceSize / size * i + ((traceSize % size > i) ? i : traceSize % size);
      end_int = traceSize / size * (i + 1) + ((traceSize % size > i) ? i + 1 : traceSize % size);
      displs[i] = begin_int;
      rcounts[i] = (end_int - begin_int);
    }

    MPI_Gatherv(
      data, len, MPI_DOUBLE,
      trace, rcounts, displs, MPI_DOUBLE,
      0, MPI_COMM_WORLD
    );
    delete[] rcounts;
    delete[] displs;
  } else
  {
    MPI_Gatherv(
      data, len, MPI_DOUBLE,
      nullptr, nullptr, nullptr, MPI_DOUBLE,
      0, MPI_COMM_WORLD
      );
  }
  
  delete[] data;
}


int main(int argc, char** argv)
{
  int rank = 0, size = 0, status = 0;
  uint32_t traceSize = 0;
  double t0 = 0, t1 = 0, dt = 0, y0 = 0, y1 = 0;
  double* trace = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0)
  {
    // Check arguments
    if (argc != 3)
    {
      std::cout << "[Error] Usage <inputfile> <output file>\n";
      status = 1;
      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      return 1;
    }

    // Prepare input file
    std::ifstream input(argv[1]);
    if (!input.is_open())
    {
      std::cout << "[Error] Can't open " << argv[1] << " for write\n";
      status = 1;
      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      return 1;
    }

    // Read arguments from input
    input >> t0 >> t1 >> dt >> y0 >> y1;
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    traceSize = (t1 - t0)/dt;
    trace = new double[traceSize];

    input.close();
  } else {
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (status != 0)
    {
      return 1;
    }
  }

  calc(trace, traceSize, t0, dt, y0, y1, rank, size);

  if (rank == 0)
  {
    // Prepare output file
    std::ofstream output(argv[2]);
    if (!output.is_open())
    {
      std::cout << "[Error] Can't open " << argv[2] << " for read\n";
      delete trace;
      return 1;
    }

    for (uint32_t i = 0; i < traceSize; i++)
    {
      output << "\n" << trace[i];
    }
    output << std::endl;
    output.close();
    delete trace;
  }

  MPI_Finalize();
  return 0;
}
