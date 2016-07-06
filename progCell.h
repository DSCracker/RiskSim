#ifndef PROGCELL_H
#define PROGCELL_H
//#include "strdiff.h"
#include <string>
#include <map>


class progCell: public Cell
{
protected:
  int m_layer = 0; // the layer of the current progCell object

  double m_pDivSym= 0.0; // the prob of sim div if div
  double m_pDivAsym= 0.0; // the prob of asym div if div
  
  double m_pDivAsymSingle= 0.0;
  double m_pDivAsymDouble= 0.0;
  double m_pDivSymDiff= 0.0;
  double m_pDivAsymDiff = 0.0;
  
  double m_cellNumberDeath= 0.0;
  double m_cellNumberStay= 0.0;
  double m_cellNumberSymDiv= 0.0;
  double m_cellNumberAsymDiv= 0.0;
  
  double m_cellNumberAsymDivSingle= 0.0; // sub prob under m_cellNumberAsymDiv: generating only one daughter stem
  double m_cellNumberAsymDivDouble= 0.0; // sub..., generating one daughter stem and one daughter progenitor
  
  double m_cellNumberSymDiff= 0.0; // generating two daughter progenitor
  double m_cellNumberAsymDiff = 0.0;

  int turnoverTime_P= 0;
  int turnoverTime_T= 0;
  
  double *p2 = new double [nsetState];
  double *p3 = new double [nsetState];
  unsigned long long *n2 = new unsigned long long [nsetState];
  unsigned long long *n3 = new unsigned long long [nsetState];
  double *tempNextM = new double [nsetState];
  double *tempNextM_NextLayer = new double [nsetState];
  double *tempNextM_T = new double [nsetState];
  
  
  
public:
  // std::string strdiff(std::string sbig, std::string ssub);
  // double m_cellNumberNext= 0.0;
  
  double *m_cellNumberNextM = NULL;
  
  double *m_cellNumberNextM_NextLayer = NULL;  // -+: changed _T to _NextLayer

  double *m_cellNumberNextM_T = NULL;  // -+: add _T back in
  
  progCell();

  void setProgCellParams(double cellNumber , double pDeath , double pDiv , double pM1 , double pM2 , double pM3 , double pM4, double pM5,
  double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff, double pDivAsymDiff);

  void setLayer(int layer);
  
  void initCellNumberNext();
  
  void updateTurnover(int turnoverTime_P_u, int turnoverTime_T_u);
  
  void updateProbDiv(double pDeath , double pDiv , double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff, double pDivAsymDiff);
  
  double getChildrenTermCellNumber();
  
  
  void getCellNumberByBehavior(unsigned int s);
  
  void getCellNumberByBehavior_mul(unsigned int s);
  
  void getCellNumberByBehavior_poi(unsigned int s);
  
  
  void updateCellNumberNext(); // take m_State, the current cell state as input, output the updated *m_cellNumberNextM and *m_cellNumberNextM_T; 
  
  void clearCell(); // set terminal cell number to zero
  
  bool isTimeToDivide();
  
  void memDelete();
  
  //void updateCellNumberNext();
};

#endif
