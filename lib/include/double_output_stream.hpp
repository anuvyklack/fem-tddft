#ifndef DOUBLE_OUTPUT_STREAM_HEADER
#define DOUBLE_OUTPUT_STREAM_HEADER

#include <iostream>

struct double_ostream
{
  double_ostream (std::ostream& out1, std::ostream& out2)
    : out1(out1), out2(out2)
  {}

  std::ostream& out1;
  std::ostream& out2;
};



/// General '<<' operator suitable for most data types.
template <typename T>
double_ostream& operator<<(double_ostream& os, T const& t)
{
   os.out1 << t;
   os.out2 << t;
   return os;
}



/// Allow for 'std::endl' to be used with a 'OutputStream'.
inline double_ostream& operator<<(double_ostream& os,
                                  std::ostream&(*f)(std::ostream&))
{
   os.out1 << f;
   os.out2 << f;
   return os;
}


#endif // DOUBLE_OUTPUT_STREAM_HEADER
// vim: ts=2 sts=2 sw=2

