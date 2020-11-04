#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

void calc(double* arr, uint32_t ySize, uint32_t xSize, int rank, int size)
{
  if (size == 1)
  {
    for (uint32_t y = 4; y < ySize; y++)
    {
      for (uint32_t x = 0; x < xSize; x++)
      {
        arr[y*xSize + x] = sin(arr[(y - 4)*xSize + x]);
      }
    }
    return;
  }

  MPI_Bcast(&ySize, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&xSize, 1, MPI_INT, 0, MPI_COMM_WORLD);


  if (rank == 0) {
    size_t begin_int[size + 1];
    for (size_t i = 0; i < size + 1; i++)
    {
      begin_int[i] = xSize / size * i + ((xSize % size > i) ? i : xSize % size);
    }
    
    for (size_t i = 1; i < size; ++i) {
      for (size_t j = begin_int[i]; j < begin_int[i+1]; ++j)
      {
        double array[] = {arr[j],
                    arr[xSize + j],
                    arr[2*xSize + j],
                    arr[3*xSize + j]};
        MPI_Send(&array, 4, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
      }
    }
    
    for (size_t i = begin_int[0]; i < begin_int[1]; ++i)
      for (size_t j = 4; j < ySize; ++j) 
      {
        arr[j*xSize + i] = sin(arr[(j - 4)*xSize + i]);
      }
    double buf[ySize];

    for (size_t i = 1; i < size; ++i) {
      for (size_t j = begin_int[i]; j < begin_int[i+1]; ++j) 
      {
        MPI_Recv(&buf, ySize, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for(size_t t = 0; t < ySize; ++t)
          arr[t*xSize + j] = buf[t];
      }
    }

  } else {
    size_t begin_int = xSize / size * rank + ((xSize % size > rank) ? rank : (xSize % size));     
    size_t end_int = xSize / size * (rank + 1) + ((xSize % size > rank) ? rank + 1 : xSize % size);

    size_t len = end_int - begin_int;

    double* data = (double *)malloc(ySize * len * sizeof(double));

    for(size_t i = 0; i < len; ++i) 
    {
      int staus;
      double buff[] = {0, 0, 0, 0};
      MPI_Recv(&buff, 4, MPI_DOUBLE, 
               0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      data[i * ySize] = buff[0];
      data[i * ySize + 1] = buff[1];
      data[i * ySize + 2] = buff[2];
      data[i * ySize + 3] = buff[3];
    }
  
    for(size_t x = 0; x < len; ++x)
      for (size_t y = 4; y < ySize; ++y)
      {
        data[x*ySize + y] = sin(data[x*ySize + y - 4]);
      }
  
    for(size_t i = 0; i < len; ++i) 
    {
      MPI_Send((data + i*ySize), ySize, MPI_DOUBLE, 
               0, 0, MPI_COMM_WORLD);          
    }

    free(data);
  }
}

int main(int argc, char** argv)
{
  int rank = 0, size = 0, buf = 0;
  uint32_t ySize = 0, xSize = 0;
  double* arr = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0)
  {
    // Check arguments
    if (argc != 3)
    {
      std::cout << "[Error] Usage <inputfile> <output file>\n";
      buf = 1;
      MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
      return 1;
    }

    // Prepare input file
    std::ifstream input(argv[1]);
    if (!input.is_open())
    {
      std::cout << "[Error] Can't open " << argv[1] << " for write\n";
      buf = 1;
      MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
      return 1;
    }

    // Read arguments from input
    input >> ySize >> xSize;
    MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);

    arr = new double[ySize * xSize];

    for (uint32_t y = 0; y < ySize; y++)
    {
     for (uint32_t x = 0; x < xSize; x++)
      {
        input >> arr[y*xSize + x];
      }
    }
    input.close();
  } else {
    MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (buf != 0)
    {
      return 1;
    }
  }

  calc(arr, ySize, xSize, rank, size);

  if (rank == 0)
  {
    // Prepare output file
    std::ofstream output(argv[2]);
    if (!output.is_open())
    {
      std::cout << "[Error] Can't open " << argv[2] << " for read\n";
      delete arr;
      return 1;
    }
    for (uint32_t y = 0; y < ySize; y++)
    {
      for (uint32_t x = 0; x < xSize; x++)
      {
        output << " " << arr[y*xSize + x];
      }
      output << std::endl;
    }
    output.close();
    delete arr;
  }

  MPI_Finalize();
  return 0;
}
