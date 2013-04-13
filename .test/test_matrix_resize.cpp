
#include "linalg.h"
#include <sstream>
#include <limits>

int main()
{
   DoubleMatrix m1(2,2);
   m1(1,1) = 1.0;
   m1(1,2) = 2.0;
   m1(2,1) = 3.0;
   m1(2,2) = 4.0;

   // increase number of columns
   {
      DoubleMatrix m2(m1);
      m2.setCols(4);

      for (int r = 1; r <= 2; ++r) {
         for (int c = 1; c <= 2; ++c) {
            if (m1(r,c) != m2(r,c)) {
               cout << "Error: m1(" << r << "," << c << ") != m2("
                    << r << "," << c << "): " << m1(r,c) << " != " << m2(r,c)
                    << endl;
               cout << "m1 = " << m1 << "m2 = " << m2 << endl;
               return 1;
            }
         }
      }
      for (int r = 1; r <= 2; ++r) {
         for (int c = 3; c <= 4; ++c) {
            if (std::fabs(m2(r,c)) > std::numeric_limits<double>::epsilon()) {
               cout << "Error: m2(" << r << "," << c << ") != 0.0: "
                    << m2(r,c) << " != " << 0.0 << endl;
               cout << "m2 = " << m2 << endl;
               return 1;
            }
         }
      }
   }

   // decrease number of columns
   {
      DoubleMatrix m3(m1);
      m3.setCols(1);

      for (int r = 1; r <= 2; ++r) {
         for (int c = 1; c <= 1; ++c) {
            if (m1(r,c) != m3(r,c)) {
               cout << "Error: m1(" << r << "," << c << ") != m3("
                    << r << "," << c << "): " << m1(r,c) << " != " << m3(r,c)
                    << endl;
               cout << "m1 = " << m1 << "m3 = " << m3 << endl;
               return 1;
            }
         }
      }
   }

   // increase number of rows
   {
      DoubleMatrix m4(m1);
      m4.setRows(4);

      for (int r = 1; r <= 2; ++r) {
         for (int c = 1; c <= 2; ++c) {
            if (m1(r,c) != m4(r,c)) {
               cout << "Error: m1(" << r << "," << c << ") != m4("
                    << r << "," << c << "): " << m1(r,c) << " != " << m4(r,c)
                    << endl;
               cout << "m1 = " << m1 << "m4 = " << m4 << endl;
               return 1;
            }
         }
      }
      for (int r = 3; r <= 4; ++r) {
         for (int c = 1; c <= 2; ++c) {
            if (std::fabs(m4(r,c)) > std::numeric_limits<double>::epsilon()) {
               cout << "Error: m4(" << r << "," << c << ") != 0.0: "
                    << m4(r,c) << " != " << 0.0 << endl;
               cout << "m4 = " << m4 << endl;
               return 1;
            }
         }
      }
   }

   // decrease number of rows
   {
      DoubleMatrix m5(m1);
      m5.setRows(1);

      for (int r = 1; r <= 1; ++r) {
         for (int c = 1; c <= 2; ++c) {
            if (m1(r,c) != m5(r,c)) {
               cout << "Error: m1(" << r << "," << c << ") != m5("
                    << r << "," << c << "): " << m1(r,c) << " != " << m5(r,c)
                    << endl;
               cout << "m1 = " << m1 << "m5 = " << m5 << endl;
               return 1;
            }
         }
      }
   }

   return 0;
}
