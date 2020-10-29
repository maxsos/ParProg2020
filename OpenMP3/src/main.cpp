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
  // Without parallel compution
  // double sum = 0.;
  // for (double i = x0; i < x1; i += dx) {
  //     sum += func(i);
  // }
  // return sum * dx;

  // Way of using many storage
  // double sum = 0.;
  // int count = (x1 - x0) / dx;
  // double* data = new double[count ? count : 1];

  // #pragma omp parallel for num_threads(num_threads_) 
  // for (int i = 0; i <= count; ++i) {
  //   data[i] = sin(x0 + dx * i) * dx;
  // }

  // for (int i = 0; i <= count; i++)
  //   sum += data[i];

  // delete[] data;
  // // std::cout.precision(14);
  // // std::cout << "Sum = " << sum << std::endl;
  // return sum;  
  
  // Way of using less storage than previous oneS
  /*
  double sum = 0.;
  int count = (x1 - x0) / dx;
  omp_set_num_threads(num_threads_);
  double data[num_threads_] = {0};

  #pragma omp parallel for num_threads(num_threads_) 
  for (int i = 0; i <= count; ++i) {
    data[omp_get_thread_num()] += func(x0 + dx * i) * dx;
  }

  for (uint32_t i = 0; i < num_threads_; i++)
    sum += data[i];

  return sum;  
  // */
//*
  double sum = 0.;
  int count = (x1 - x0) / dx;
  int i = 0;
  omp_set_num_threads(num_threads_);

  #pragma omp parallel for num_threads(num_threads_) reduction(+:sum)
  for (i = 0; i <= count; ++i) {
    sum += func(x0 + dx * i) * dx;
  }

  if(i * dx < x1) {
    sum += func()
  }

  return sum;  
// */
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
