
/** \file sigma.h
   - Project:     Markov
   - Author:      Ben Allanach 

   - Description: calls HERWIG to determine cross sections etc

*/

extern "C" int inithwg_();

extern "C" int readinp_();

extern "C" int inithwg2_();

extern "C" int calcsig_(int * ipp, double * sigmaGlu, double * sigmaUpL, 
			double * sigmaDnL, double * sigmaChL, 
			double * sigmaStL);

extern "C" int benin_();

