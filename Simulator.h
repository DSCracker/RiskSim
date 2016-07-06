#ifndef SIMULATOR_H
#define SIMULATOR_H
#include <iostream>
#include <string>
#include "Params.h"
#include "Cell.h"
#include "stemCell.h"
#include "progCell.h"

class Simulator: public Params
{
protected:


std::string dynamics_effect_option = "M1";
// list all cell objects
    
    Cell cell;
    stemCell *scm = NULL;
    progCell **pcm = NULL;
    std::string MChar = "0";

public:
    Simulator(); // default constructor: assign default value to member cell objects
    
    void SetCellState();
    void SetprogCellLayer(); // ++
    
    void SetEffectOption(std::string dynamics_effect_option);
    
    void SetStem(double cellNumber , double pDeath , double pDiv , double pM1 , double pM2 , double pM3 ,
   double pM4, double pM5,
  double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff,
  double sltAdv, double mfacDeathRate, double mfacMutRate);
    
    void SetProg(double cellNumber , double pDeath , double pDiv , double pM1 , double pM2 , double pM3 ,
   double pM4, double pM5,
  double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff, double pDivAsymDiff,
  double sltAdv, double mfacDeathRate, double mfacMutRate, int il); // ++ include layer
  
  void SetStemCellNumber(double cellNumber);
  void SetProgCellNumber(double cellNumber, int il); // ++ include layer
  
  void SetStemCellM1Number(double cellNumber); 
  void SetStemCellM3Number(double cellNumber);
  
  void SetStemCellMNumber(double cellNumber);
  
  void SetProgCellM1Number(double cellNumber, int il); // ++ include layer
  void SetProgCellM3Number(double cellNumber, int il); // ++ include layer
  
  void SetProgCellMNumber(double cellNumber, int il); // ++ include layer  
  
    
    void SetStemDiv(double pDeath , double pDiv , double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff,
    double sltAdv, double mfacDeathRate);
    
    void SetStemM1Div(double pDeath , double pDiv , double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff,
    double sltAdv, double mfacDeathRate);
    
    void SetStemM3Div(double pDeath , double pDiv , double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff,
    double sltAdv, double mfacDeathRate);
    
    void SetStemMDiv(double pDeath , double pDiv , double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff,
    double sltAdv, double mfacDeathRate);
    
    
    
    void SetProgDiv(double pDeath , double pDiv , double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff, double pDivAsymDiff,
    double sltAdv, double mfacDeathRate, int il); // ++ include layer
    
    void SetProgM1Div(double pDeath , double pDiv , double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff, double pDivAsymDiff,
    double sltAdv, double mfacDeathRate, int il); // ++ include layer
    
    void SetProgM3Div(double pDeath , double pDiv , double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff, double pDivAsymDiff,
    double sltAdv, double mfacDeathRate, int il); // ++ include layer
    
    void SetProgMDiv(double pDeath , double pDiv , double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff, double pDivAsymDiff,
    double sltAdv, double mfacDeathRate, int il); // ++ include layer
    
    
    void SetStemMut(double pM1 , double pM2 , double pM3, double pM4, double pM5, double mfacMutRate);
    void SetProgMut(double pM1 , double pM2 , double pM3, double pM4, double pM5, double mfacMutRate, int il); // include layer
    
        
    // dynamics functions start here
    void getMChar();
    void StemDivToNext(unsigned int s, std::string thisType);
    void StemNumUpdateFromPrev(std::string thisType, std::string prevType);

    void ProgDivToNext(unsigned int s, int thisLayer, std::string thisType);
    void ProgNumUpdateFromPrev(int thisLayer, int prevLayer, std::string thisType, std::string prevType);

    void updateStemCellNumberFromStem(unsigned int s);
    void updateProgCellNumberFromProg(unsigned int s, int il); // ++ include layer
    void updateProgCellNumberFromStem(); // layer should be 0 

    
    
    void updateStemCellM1NumberFromStemM1(unsigned int s);
    void updateStemCellM3NumberFromStemM3(unsigned int s);
    
    void updateStemCellMNumberFromStemM(unsigned int s);
    
    //+++
    void updateStemCellNumberFromStemNonM1(unsigned int s);//check
    void updateStemCellNumberFromStemNonM3(unsigned int s);//check
    void updateStemCellNumberFromStemNonM(unsigned int s);// check
    
    
    void updateProgCellM1NumberFromProgM1(unsigned int s, int il); // ++ include layer
    void updateProgCellM3NumberFromProgM3(unsigned int s, int il); // ++ include layer
    
    void updateProgCellMNumberFromProgM(unsigned int s, int il); // ++ include layer
    
    
    //+++
    void updateProgCellNumberFromProgNonM1(unsigned int s, int il); // check // ++ include layer
    void updateProgCellNumberFromProgNonM3(unsigned int s, int il); // check // ++ include layer
    void updateProgCellNumberFromProgNonM(unsigned int s, int il); // check // ++ include layer
    
    
    void updateProgCellM1NumberFromStemM1();  // layer should be 0
    void updateProgCellM3NumberFromStemM3();  // layer should be 0
    
    void updateProgCellMNumberFromStemM();  // layer should be 0
    
    // +++
    void updateProgCellNumberFromStemNonM1();// check // layer should be 0
    void updateProgCellNumberFromStemNonM3();// check // layer should be 0
    void updateProgCellNumberFromStemNonM();// check // layer should be 0
    
    
    void addCellNumberFromStemNonM1ToStemM1();// check
    void addCellNumberFromStemNonM3ToStemM3();// check
    void addCellNumberFromStemNonMToStemM();// check
    
    void addCellNumberFromProgNonM1ToProgM1(int il);// check // ++ include layer
    void addCellNumberFromProgNonM3ToProgM3(int il);// check // ++ include layer
    void addCellNumberFromProgNonMToProgM(int il);// check // ++ include layer
    
    // result reporting functions start here
    double getTotalStemCellNumber(std::string thisType);
    double getTotalProgCellNumber(std::string thisType); 
    double getTotalProgCellNumber_byLayer(int il, std::string thisType); // ++ include layer 
    double getTotalTermCellNumber(); // layer should be last
    
    double getTotalStemCellM1Number();
    double getTotalStemCellM3Number();
    
    double getTotalStemCellMNumber();
    
    double getTotalProgCellM1Number(); // 
    double getTotalProgCellM3Number(); // 
    double getTotalProgCellMNumber(); // 
    
    double getCancerStemCellNumber();
    double getCancerProgCellNumber();
    double getCancerTermCellNumber();
    
    void clearProgCell(std::string thisType);
    void cellNumberDebugger();
    
    void memFree();
    
};



#endif
