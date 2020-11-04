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
    for (uint32_t y = 0; y < ySize - 1; y++)
    {
      for (uint32_t x = 3; x < xSize; x++)
      {
        arr[y*xSize + x] = sin(0.00001*arr[(y + 1)*xSize + x - 3]);
      }
    }
    return;
  }
  MPI_Bcast(&ySize, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&xSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int* sendcounts, *rcounts, *displs;
  double* recvbuf;
  size_t end_int;
  size_t begin_int;
  sendcounts = (int *)malloc(size * sizeof(int));
  displs = (int *)malloc(size * sizeof(int));
  rcounts = (int *)malloc(size * sizeof(int));

  for (size_t i = 0; i < size; i++)
  {
    begin_int = ySize / size * i + ((ySize % size > i) ? i : ySize % size);
    end_int = ySize / size * (i + 1) + ((ySize % size > i) ? i + 1 : ySize % size);
    displs[i] = begin_int *  xSize;
    sendcounts[i] = (end_int - begin_int + ((i+1 < size) ? 1 : 0)) * xSize;
    rcounts[i] = (end_int - begin_int)* xSize;
  }



  if (rank == 0) {    
    // for (size_t i = 0; i < size; i++) {
    //   std::cout << "[" << i << "]" 
    //   << "ready to send " << sendcounts[i] << std::endl;
    // }
    // for (size_t i = 0; i < size; i++) {
    //   std::cout << "[" << i << "]" 
    //   << "rcounts " << rcounts[i] << std::endl;
    // }
    begin_int = ySize / size * 1 + ((ySize % size > 1) ? 1 : ySize % size);

    recvbuf = (double *)malloc(sendcounts[0] * sizeof(double));

    MPI_Scatterv(arr, sendcounts, displs, MPI_DOUBLE, 
                 recvbuf, sendcounts[0], MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

    
    for (size_t y = 0; y < begin_int; ++y)
      for (size_t x = 3; x < xSize; ++x) 
      {
        recvbuf[y*xSize + x] = sin(0.00001*recvbuf[(y + 1)*xSize + x - 3]);
      }

    MPI_Gatherv(recvbuf, rcounts[0], MPI_DOUBLE,
                arr, rcounts, displs, MPI_DOUBLE, 
                0, MPI_COMM_WORLD);
     
    free(recvbuf);

  } else {
    size_t begin_int = ySize / size * rank + ((ySize % size > rank) ? rank : (ySize % size));     
    end_int = ySize / size * (rank + 1) + ((ySize % size > rank) ? rank + 1 : ySize % size);

    size_t len = end_int - begin_int + ((rank + 1 < size) ? 1 : 0);

    // std::cout << "[" << rank << "]" << "ready to get " << len * xSize << std::endl;
    

    double* recvbuf = (double *)malloc(xSize * len * sizeof(double));

    MPI_Scatterv(nullptr, nullptr, nullptr, MPI_DOUBLE, 
                 recvbuf, len * xSize, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

    for (size_t y = 0; y < len - 1; ++y)
      for (size_t x = 3; x < xSize; ++x) 
      {
        recvbuf[y*xSize + x] = sin(0.00001*recvbuf[(y + 1)*xSize + x - 3]);
      }

    if (rank + 1 < size)
      --len;

    MPI_Gatherv(recvbuf, len * xSize, MPI_DOUBLE,
                nullptr, nullptr, nullptr, MPI_DOUBLE, 
                0, MPI_COMM_WORLD);          

    free(recvbuf);
  }
  

    free(sendcounts);
    free(displs);    
    free(rcounts);
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
