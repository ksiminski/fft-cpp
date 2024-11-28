
#include <complex>
#include <iostream>
#include <sstream>
#include <vector>

#include "fft.h"

template <typename T>
std::ostream & operator << (std::ostream & sos, const std::vector<T> & data)
{
   sos << "[ ";
   for (const auto & d : data)
      sos << d << " ";
   return sos << "]";
}

std::string polynomial_to_string(const std::vector<std::complex<double>> & a)
{
   std::stringstream sos;
   for (std::size_t i = 0; i < a.size(); i++)
   {
      if (a[i].real() != 0)
      {
         if (i > 0)
            sos << " + ";
         sos << a[i].real();
         if (i > 0)
            sos << "x^" <<i;
      }
   }
   return sos.str();
}

int main()
{

   std::cout << "Fast Fourier Transform" << std::endl;

   std::cout << "---------" << std::endl;
   std::cout << "Example 1: FFT and IFFT" << std::endl;
   {
      // simple example of FFT
      std::vector<std::complex<double>> data = {1, 1, 1, 1, 0, 0, 0, 0};
      std::cout << data << std::endl;

      ksi::fft f;
      f.iterative_fft(data);  // fast Fourier transform
      std::cout << data << std::endl;

      f.iterative_inverted_fft(data); // inverted fast Fourier transform 
      std::cout << data << std::endl;

   }

   std::cout << "---------" << std::endl;
   std::cout << "Example 2: FFT and IFFT" << std::endl;
   {
      // zero padding if the size of the input data is not a power of 2
      std::vector<std::complex<double>> data = {1, 1, 1, 1, 0, 0, 0, 0, 0, 0};
      std::cout << data << std::endl;

      ksi::fft f;
      f.iterative_fft(data);  // fast Fourier transform
      std::cout << data << std::endl;

      f.iterative_inverted_fft(data); // inverted fast Fourier transform 
      std::cout << data << std::endl;

   }
   {
      std::cout << "---------" << std::endl;
      std::cout << "Example 3: polynomial multiplication" << std::endl;

      std::vector<std::complex<double>> A = {1, 2, 0, 3};   // A(x) = 1 + 2x + 3x^3
      std::vector<std::complex<double>> B = {4, 0, 5, 0, 2};    // B(x) = 4 + 5x^2 + 2x^4
      std::cout << "A : " << polynomial_to_string(A) << std::endl;
      std::cout << "B : " << polynomial_to_string(B) << std::endl;

      std::size_t size = A.size() + B.size() - 1;

      std::size_t n = 1;
      while (n < size)
         n <<= 1;

      A.resize(n);
      B.resize(n);

      std::vector<std::complex<double>> C(n);

      ksi::fft f;
      f.iterative_fft(A);  // fast Fourier transform
      f.iterative_fft(B);  // fast Fourier transform

      for (std::size_t i = 0; i < n; i++)
         C[i] = A[i] * B[i];

      f.iterative_inverted_fft(C); // inverted fast Fourier transform

      std::cout << "C = A * B : " << polynomial_to_string(C) << std::endl;

   }
   return 0;
}

