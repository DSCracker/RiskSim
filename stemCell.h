#ifndef STEMCELL_H
#define STEMCELL_H

#include <string>
#include <map>

//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_cdf.h>
//#include <gsl/gsl_rng.h>

//#include "strdiff.h"

class stemCell: public Cell
{
protected:
  double m_pDivSym = 0.0; // the prob of sim div if div
  double m_pDivAsym = 0.0; // the prob of asym div if div
  
  double m_pDivAsymSingle = 0.0;
  double m_pDivAsymDouble = 0.0;
  double m_pDivSymDiff = 0.0;
  
  double m_cellNumberDeath = 0.0;
  double m_cellNumberStay = 0.0;
  double m_cellNumberSymDiv = 0.0;
  double m_cellNumberAsymDiv = 0.0;
  
  double m_cellNumberAsymDivSingle = 0.0; // sub prob under m_cellNumberAsymDiv: generating only one daughter stem
  double m_cellNumberAsymDivDouble = 0.0; // sub..., generating one daughter stem and one daughter progenitor
  
  double m_cellNumberSymDiff = 0.0; // generating two daughter progenitor
  int turnoverTime_S = 0;
  
  double *p2 = new double [nsetState];
  double *p3 = new double [nsetState];
  unsigned long long *n2 = new unsigned long long [nsetState];
  unsigned long long *n3 = new unsigned long long [nsetState];
  double *tempNextM = new double [nsetState];
  double *tempNextM_P = new double [nsetState];
  
  
  unsigned int n2_gsl[32] = {0};
  unsigned int n3_gsl[32] = {0};
  
  //gsl_rng * r;
  //const gsl_rng_type * T;
    
  //std::string sdiff = "0";
  
public:
  //void strdiff(int i);
   
  double m_cellNumberNext = 0.0;
  
  double *m_cellNumberNextM = NULL;
  
  double *m_cellNumberNextM_P = NULL; 
  
  stemCell();
    
  void setStemCellParams(double cellNumber , double pDeath , double pDiv , double pM1 , double pM2 , double pM3 , double pM4, double pM5,
  double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff);
  
  void initCellNumberNext();
  
  void updateProbDiv(double pDeath , double pDiv , double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff);
  
  void updateTurnover(int turnoverTime_S_u);
  
  void getCellNumberByBehavior(unsigned int s);
  
 // void getCellNumberByBehavior_mul(unsigned int s);
  
  void getCellNumberByBehavior_poi(unsigned int s);
  
 // void getCellNumberByBehavior_mul_gsl(unsigned int s);
  
  void updateCellNumberNext(); // take m_State, the current cell state as input, output the updated *m_cellNumberNextM and *m_cellNumberNextM_P; 
  
  bool isTimeToDivide();
  
  void memDelete();
  
  
  //void updateCellNumberNext();
};

#endif
