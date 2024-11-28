
#include <complex>
#include <vector>

#include "fft.h"

bool ksi::fft::is_power_of_two(const std::size_t n)
{
   return (n & (n - 1)) == 0;
}

std::size_t ksi::fft::next_power_of_two(const std::size_t n)
{
   std::size_t power = 1;
   while (power < n)
      power <<= 1;
   return power;
}

std::size_t ksi::fft::resize_if_needed (std::vector<std::complex<double>> & a)
{
   std::size_t n = a.size();
   if (not is_power_of_two(n))
   {
      std::size_t power = next_power_of_two(n);
      a.resize(power);
   }
   return a.size();
}

void ksi::fft::bit_reverse(std::vector<std::complex<double>> & a)
{
   std::size_t n = a.size();
   for (std::size_t i = 1, j = 0; i < n; i++) 
   {
      int bit = n >> 1;
      for (; j & bit; bit >>= 1) 
         j ^= bit;
      j ^= bit;
      if (i < j) 
         std::swap(a[i], a[j]);
   }
}

void ksi::fft::butterfly (std::vector<std::complex<double>> & a, const bool invert)
{
   std::size_t n = a.size();
   for (std::size_t len = 2; len <= n; len *= 2) 
   {
      double ang = 2 * M_PI / len * (invert ? -1 : 1);
      std::complex<double> wlen(cos(ang), sin(ang));
      for (std::size_t i = 0; i < n; i += len) 
      {
         std::complex<double> w(1);
         for (std::size_t j = 0; j < len / 2; j++) 
         {
            std::complex<double> u = a[i + j], v = a[i + j + len / 2] * w;
            a[i + j] = u + v;
            a[i + j + len / 2] = u - v;
            w *= wlen;
         }
      }
   }
   if (invert) 
      for (auto &x : a) 
         x /= n;

}

void ksi::fft::fast_fourier_transform (std::vector<std::complex<double>> & a, const bool invert) 
{
   resize_if_needed(a);
   bit_reverse(a);
   butterfly(a, invert);
}

void ksi::fft::iterative_fft(std::vector<std::complex<double>> & a) 
{
   fast_fourier_transform(a, false);
}

void ksi::fft::iterative_inverted_fft(std::vector<std::complex<double>> & a) 
{
   fast_fourier_transform(a, true);
}


