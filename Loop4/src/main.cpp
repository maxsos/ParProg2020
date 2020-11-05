#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

void calc(double* arr, uint32_t zSize, uint32_t ySize, uint32_t xSize, int rank, int size)
{
  if (size >= 1) {
    for (uint32_t z = 1; z < zSize; z++) {
      for (uint32_t y = 0; y < ySize - 1; y++) {
        for (uint32_t x = 0; x < xSize - 1; x++) {
          arr[z*ySize*xSize + y*xSize + x] = sin(arr[(z - 1)*ySize*xSize + (y + 1)*xSize + x + 1]);
        }
      }
    }
    return;
  }

  MPI_Bcast(&zSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ySize, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&xSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0) {   
    int *sendcounts = (int *)malloc(size * sizeof(int)),
        *rcounts = (int *)malloc(size * sizeof(int)),
        *displs = (int *)malloc(size * sizeof(int));
    double* recvbuf;
    size_t end_int, begin_int;
    
    
    for (size_t i = 0; i < size; i++)
    {
      begin_int = zSize / size * i + ((zSize % size > i) ? i : zSize % size);
      end_int = zSize / size * (i + 1) + ((zSize % size > i) ? i + 1 : zSize % size);
      displs[i] = begin_int *  xSize * ySize;
      sendcounts[i] = (end_int - begin_int + ((i+1 < size) ? 1 : 0)) * xSize * ySize;
      rcounts[i] = (end_int - begin_int)* xSize * ySize;
    } 
    begin_int = zSize / size * 1 + ((zSize % size > 1) ? 1 : zSize % size);

    recvbuf = (double *)malloc(sendcounts[0] * sizeof(double));

    MPI_Scatterv(arr, sendcounts, displs, MPI_DOUBLE, 
                 recvbuf, sendcounts[0], MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

    
    // for (size_t y = 0; y < begin_int; ++y)
    //   for (size_t x = 3; x < xSize; ++x) 
    //   {
    for (uint32_t z = 1; z < begin_int; z++) {
      for (uint32_t y = 0; y < ySize - 1; y++) {
        for (uint32_t x = 0; x < xSize - 1; x++) {
         recvbuf[z*ySize*xSize + y*xSize + x] = 
                sin(
                  recvbuf[(z - 1)*ySize*xSize + (y + 1)*xSize + x + 1]
                );
        }
      }
    }

    MPI_Gatherv(recvbuf, rcounts[0], MPI_DOUBLE,
                arr, rcounts, displs, MPI_DOUBLE, 
                0, MPI_COMM_WORLD);


    free(recvbuf);
    free(sendcounts);
    free(displs);    
    free(rcounts);

  } else {
    size_t begin_int = zSize / size * rank + ((zSize % size > rank) ? rank : (zSize % size));     
    size_t end_int = zSize / size * (rank + 1) + ((zSize % size > rank) ? rank + 1 : zSize % size);

    size_t len = end_int - begin_int + ((rank + 1 < size) ? 1 : 0);

    double* recvbuf = (double *)malloc(ySize * xSize * len * sizeof(double));

    MPI_Scatterv(nullptr, nullptr, nullptr, MPI_DOUBLE, 
                 recvbuf, len * xSize * ySize, MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

    for (uint32_t z = 1; z < len; z++) {
      for (uint32_t y = 0; y < ySize - 1; y++) {
        for (uint32_t x = 0; x < xSize - 1; x++) {
         recvbuf[z*ySize*xSize + y*xSize + x] = 
                sin(
                  recvbuf[(z - 1)*ySize*xSize + (y + 1)*xSize + x + 1]
                );
        }
      }
    }

    if (rank + 1 < size)
      --len;

    MPI_Gatherv(recvbuf, len * xSize * ySize, MPI_DOUBLE,
                nullptr, nullptr, nullptr, MPI_DOUBLE, 
                0, MPI_COMM_WORLD);          

    free(recvbuf);
  }


}

int main(int argc, char** argv)
{
  int rank = 0, size = 0, buf = 0;
  uint32_t zSize = 0, ySize = 0, xSize = 0;
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
    input >> zSize >> ySize >> xSize;
    MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);

    arr = new double[zSize * ySize * xSize];
    for (uint32_t z = 0; z < zSize; z++) {
      for (uint32_t y = 0; y < ySize; y++) {
        for (uint32_t x = 0; x < xSize; x++) {
          input >> arr[z*ySize*xSize + y*xSize + x];
        }
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

  calc(arr, zSize, ySize, xSize, rank, size);

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

    for (uint32_t z = 0; z < zSize; z++) {
      for (uint32_t y = 0; y < ySize; y++) {
        for (uint32_t x = 0; x < xSize; x++) {
          output << " " << arr[z*ySize*xSize + y*xSize + x];
        }
        output << std::endl;
      }
      output << std::endl;
    }
    output.close();
    delete arr;
  }

  MPI_Finalize();
  return 0;
}
