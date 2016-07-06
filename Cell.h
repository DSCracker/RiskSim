#ifndef CELL_H
#define CELL_H


#include "Params.h"
#include <string>
#include <map>
#include <set>

class Cell: public Params
{
protected:
  // define the set and map of celltypes

  
   
//  char m_cellType; // "S","P","T": stem, progenitor and terminal cells
// cell dynamics variables
  double m_cellNumber; // cell Number
  double m_cellNumberFromNonM;
  double m_pDeath; // death prob at 1 day
  double m_pDiv; // div prob at 1 day
  double m_pStay; // stay still prob at 1 day
//  double m_pMnone; // prob of no mutation during div

// +++
  // public:
  // --
  //double m_pM0, m_pM1, m_pM2, m_pM3, m_pM4, m_pM5, m_pM12, m_pM13, m_pM14, m_pM15, m_pM23, m_pM24, m_pM25,
  //       m_pM34, m_pM35, m_pM45, m_pM123, m_pM124, m_pM125, m_pM134, m_pM135, m_pM145, m_pM234, m_pM235, m_pM245,m_pM345, m_pM1234, m_pM1235, m_pM1245, m_pM1345, m_pM2345, m_pM12345; // mutation probability: pM0, pM1, pM2, mP3
  
  // ++
  double *m_pM = NULL ;
  
  std::string m_State = "000000"; 
     
  int timeCounter = 0;
  
public:
//double *m_pM = NULL ;
  Cell();
  
  void setCellParams(double cellNumber, double pDeath, double pDiv , double pM1 , double pM2 ,double pM3, double mP4, double pM5 );
  
  void updateState(std::string State);
  
  void updateProbDiv(double pDeath , double pDiv );
  
  void updateProbMutation(double pM1 , double pM2 , double pM3, double pM4, double pM5 );
  
  void updateCellNumber(double cellNumber);
  
  void updateCellNumberFromNonM(double cellNumber);
  
  void addCellNumber(double cellNumber);
  
  void addCellNumberFromNonM();
    
  double getCellNumber();
  
  void setTimeCounter(int timeCounter_u);
  
  void addTimeCounter();
  
  // bool isSubstr(std::string s, std::string ssub);
  
 // std::string strdiff(std::string s, std::string ssub);
  
  // bool isEqstr(std::string s1, std::string s2);
  
};

#endif
