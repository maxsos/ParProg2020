#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>

double calc(uint32_t x_last, uint32_t num_threads_)
{
  double sum = 0.;
  double fact = 1.;
  omp_set_num_threads(num_threads_);
  if (x_last == 1 || x_last ==2)
    --x_last;

  double index[num_threads_ + 1], subsum[num_threads_ + 1] =  {0};
  #pragma omp parallel for num_threads(num_threads_)
    for (uint32_t i = 0; i <= num_threads_; ++i)
      index[i] = 1;
  

  // Init coefs 
  #pragma omp parallel for num_threads(num_threads_) ordered
  for (uint32_t i = 0; i <= x_last; i++) {
    index[omp_get_thread_num() + 1] /= i ? i : 1;
    subsum[omp_get_thread_num()] += index[omp_get_thread_num() + 1];
  }

  for (uint32_t i = 0; i <= num_threads_; ++i) {
    sum += fact * subsum[i];
    fact *= index[i + 1];
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
  uint32_t x_last = 0, num_threads = 0;
  input >> x_last >> num_threads;

  // Calculation
  double res = calc(x_last, num_threads);

  // Write result
  output << std::setprecision(16) << res << std::endl;
  // Prepare to exit
  output.close();
  input.close();
  return 0;
}
