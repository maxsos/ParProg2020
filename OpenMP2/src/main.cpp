#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>

double calc(uint32_t x_last, uint32_t num_threads_)
{
    double sum = 0;
    #pragma omp parrallel for num_threads(num_threads_) reduction(+:sum)
    for (double i = x_last; i != 0; --i)
        sum += 1/i;
    return sum;
}

int main(int argc, char** argv)
{
  // Check arguments
  if (argc != 3)
  {
    std::cout << "[Error] Usage <inputfile> <output file>\n";
    return 1;
  }

  // Prepare input file
  std::ifstream input(argv[1]);
  if (!input.is_open())
  {
    std::cout << "[Error] Can't open " << argv[1] << " for write\n";
    return 1;
  }

  // Prepare output file
  std::ofstream output(argv[2]);
  if (!output.is_open())
  {
    std::cout << "[Error] Can't open " << argv[2] << " for read\n";
    input.close();
    return 1;
  }

  // Read arguments from input
  uint32_t x_last = 0, num_threads = 0;
  input >> x_last >> num_threads;
  // Calculation
  double res = calc(x_last, num_threads);

  // Write result
  output << std::setprecision(15) << res << std::endl;
  // Prepare to exit
  output.close();
  input.close();
  return 0;
}
