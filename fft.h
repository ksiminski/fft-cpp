#ifndef KSI_FFT_H
#define KSI_FFT_H

#include <complex>
#include <vector>

namespace ksi
{
   // Fast Fourier Transform
   class fft
   {
      private:
         bool is_power_of_two(const std::size_t n);
         std::size_t next_power_of_two(const std::size_t n);
         std::size_t resize_if_needed (std::vector<std::complex<double>> & a);
         void bit_reverse(std::vector<std::complex<double>> & a);
         void butterfly (std::vector<std::complex<double>> & a, const bool invert);
         void fast_fourier_transform (std::vector<std::complex<double>> & a, const bool invert);

      public:
         void iterative_fft(std::vector<std::complex<double>> & a);
         void iterative_inverted_fft(std::vector<std::complex<double>> & a);
   };
}

#endif // KSI_FFT_H
