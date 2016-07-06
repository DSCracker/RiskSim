#include "Params.h"
#include "Cell.h"
#include <string>
#include <map>
#include <set>
#include <assert.h> 

// 2/2/16
Cell::Cell(){
    setMap();
    updateState("0");
    updateCellNumber(1);
    updateProbDiv(0, 1);
    updateProbMutation(0, 0, 0, 0, 0);
}


void Cell::setCellParams(double cellNumber , double pDeath , double pDiv , double pM1 , double pM2 ,double pM3, double pM4, double pM5 ){
     updateCellNumber(cellNumber);
     updateProbDiv(pDeath, pDiv);
     updateProbMutation(pM1, pM2, pM3, pM4, pM5);
  }


void Cell::updateState(std::string State){
     m_State = State;
}

void Cell::updateProbDiv(double pDeath , double pDiv ){
     m_pDeath = pDeath;
     m_pDiv = pDiv;
     //m_pStay = 1.0 - pDeath - pDiv;
     m_pStay = 0; 
  }

void Cell::updateProbMutation(double pM1 , double pM2 , double pM3, double pM4, double pM5){
     m_pM = new double [nsetState];
     for(int i = 0; i < nsetState; i++){
         m_pM[i] = 0;
     } 
     
     if(pM1 > 1.0){ pM1 = 1.0;}
     if(pM2 > 1.0){ pM2 = 1.0;}
     if(pM3 > 1.0){ pM3 = 1.0;}
     if(pM4 > 1.0){ pM4 = 1.0;}
     if(pM5 > 1.0){ pM5 = 1.0;}

     m_pM[setState2ind["0"]] = (1.0 - pM1) * (1.0 - pM2) * (1.0 - pM3) * (1 - pM4) * (1 - pM5);
    
     m_pM[setState2ind["1"]] = pM1 * (1 - pM2) * (1 - pM3) * (1 - pM4) * (1 - pM5);
     m_pM[setState2ind["2"]] = (1 - pM1) * pM2 * (1 - pM3) * (1 - pM4) * (1 - pM5);
     m_pM[setState2ind["3"]] = (1 - pM1) * (1 - pM2) * pM3 * (1 - pM4) * (1 - pM5);
     m_pM[setState2ind["4"]] = (1 - pM1) * (1 - pM2) * (1 - pM3) * pM4 * (1 - pM5);
     m_pM[setState2ind["5"]] = (1 - pM1) * (1 - pM2) * (1 - pM3) * (1 - pM4) * pM5;
    
     m_pM[setState2ind["12"]] = pM1 * pM2 * (1 - pM3) * (1 - pM4) * (1 - pM5);
     m_pM[setState2ind["13"]] = pM1 * (1 - pM2) * pM3 * (1 - pM4) * (1 - pM5);
     m_pM[setState2ind["14"]] = pM1 * (1 - pM2) * (1 - pM3) * pM4 * (1 - pM5);
     m_pM[setState2ind["15"]] = pM1 * (1 - pM2) * (1 - pM3) * (1 - pM4) * pM5;
     m_pM[setState2ind["23"]] = (1 - pM1) * pM2 * pM3 * (1 - pM4) * (1 - pM5);
     m_pM[setState2ind["24"]] = (1 - pM1) * pM2 * (1 - pM3) * pM4 * (1 - pM5);
     m_pM[setState2ind["25"]] = (1 - pM1) * pM2 * (1 - pM3) * (1 - pM4) * pM5;
     m_pM[setState2ind["34"]] = (1 - pM1) * (1 - pM2) * pM3 * pM4 * (1 - pM5);
     m_pM[setState2ind["35"]] = (1 - pM1) * (1 - pM2) * pM3 * (1 - pM4) * pM5;
     m_pM[setState2ind["45"]] = (1 - pM1) * (1 - pM2) * (1 - pM3) * pM4 * pM5;
    
     m_pM[setState2ind["123"]] = pM1 * pM2 * pM3 * (1 - pM4) * (1 - pM5);
     m_pM[setState2ind["124"]] = pM1 * pM2 * (1 - pM3) * pM4 * (1 - pM5);
     m_pM[setState2ind["125"]] = pM1 * pM2 * (1 - pM3) * (1 - pM4) * pM5;
     m_pM[setState2ind["134"]] = pM1 * (1 - pM2) * pM3 * pM4 * (1 - pM5);
     m_pM[setState2ind["135"]] = pM1 * (1 - pM2) * pM3 * (1 - pM4) * pM5;
     m_pM[setState2ind["145"]] = pM1 * (1 - pM2) * (1 - pM3) * pM4 * pM5;
     m_pM[setState2ind["234"]] = (1 - pM1) * pM2 * pM3 * pM4 * (1 - pM5);
     m_pM[setState2ind["235"]] = (1 - pM1) * pM2 * pM3 * (1 - pM4) * pM5;
     m_pM[setState2ind["245"]] = (1 - pM1) * pM2 * (1 - pM3) * pM4 * pM5;
     m_pM[setState2ind["345"]] = (1 - pM1) * (1 - pM2) * pM3 * pM4 * pM5;
    
     m_pM[setState2ind["1234"]] = pM1 * pM2 * pM3 * pM4 * (1 - pM5);
     m_pM[setState2ind["1235"]] = pM1 * pM2 * pM3 * (1 - pM4) * pM5;
     m_pM[setState2ind["1245"]] = pM1 * pM2 * (1 - pM3) * pM4 * pM5;
     m_pM[setState2ind["1345"]] = pM1 * (1 - pM2) * pM3 * pM4 * pM5;
     m_pM[setState2ind["2345"]] = (1 - pM1) * pM2 * pM3 * pM4 * pM5;
    
     m_pM[setState2ind["12345"]] = pM1 * pM2 * pM3 * pM4 * pM5;

  }

void Cell::updateCellNumber(double cellNumber){
     m_cellNumber = cellNumber;
  }

void Cell::updateCellNumberFromNonM(double cellNumber){
     m_cellNumberFromNonM = cellNumber;
}


void Cell::addCellNumber(double cellNumber){
     m_cellNumber += cellNumber;
}


void Cell::addCellNumberFromNonM(){
     addCellNumber(m_cellNumberFromNonM);
}

double Cell::getCellNumber(){
     return m_cellNumber;
  }


void Cell::setTimeCounter(int timeCounter_u){
     timeCounter = timeCounter_u;
}
  
void Cell::addTimeCounter(){
     timeCounter++;
}

/*
bool Cell::isSubstr(std::string s, std::string ssub){
    bool issubstr = true;
    for(int i = 0; i < ssub.length(); i++){
        if(s.find(ssub[i]) == std::string::npos){
            issubstr = false;
        }
    }
    return issubstr;
}
*/
/*
std::string Cell::strdiff(std::string s, std::string ssub){
// return the the chars in s but not in ssub (normally ssub is a subset of set s;)
// the order of chars in string for this function does not matter
    //std::cout<<ssub<<" "<<ssub.length()<<" "<<s<<" "<<s.length()<<std::endl;
    assert(ssub.length() < s.length());
    char *ch = new char [s.length() - ssub.length()];
    int counter = 0;
    for(int i = 0; i < s.length(); i++){
        if(ssub.find(s[i]) == std::string::npos){
        // not found
            ch[counter] = s[i];
            counter ++;
        }
    }
    std::string sdiff = ch;
    delete [] ch;
    return sdiff;
}
*/

/*
bool Cell::isEqstr(std::string s1, std::string s2){
 return ((s1.length() == s2.length()) && (isSubstr(s1, s2)));
}
*/



