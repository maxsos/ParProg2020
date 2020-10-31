#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include <cmath>

double func(double x)
{
  return sin(x);
}

double calc(double x0, double x1, double dx, uint32_t num_threads_)
{
  double sum = dx * func(x0)/2;
  int count = (x1 - x0) / dx;
  int i = 0;
  omp_set_num_threads(num_threads_);

  #pragma omp parallel for num_threads(num_threads_) reduction(+:sum)
  for (i = 1; i < count; ++i) {
    if (i != count - 1)
      sum += func(x0 + dx * i) * dx;
    else 
      sum += func(x0 + dx * i) * dx/2;
  }

  if (count * dx + x0 < x1) {
    sum += (func(x0 + (count-1) * dx) + func(x0 + count * dx)) * dx/2;
  }

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
  double x0 = 0.0, x1 =0.0, dx = 0.0;
  uint32_t num_threads = 0;
  input >> x0 >> x1 >> dx >> num_threads;

  // Calculation
  double res = calc(x0, x1, dx, num_threads);

  // Write result
  output << std::setprecision(13) << res << std::endl;
  // Prepare to exit
  output.close();
  input.close();
  return 0;
}
