#include "Params.h"
#include "Cell.h"
#include "progCell.h"

#include "strdiff.h"
#include "isSubstr.h"
#include "isEqstr.h"
#include "strunion.h"


#include <iostream>
#include <string>
#include <random>
//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_cdf.h>
//#include <gsl/gsl_rng.h>
#include <time.h>       /* clock_t, clock, CLOCKS_TER_SEC */
#include <math.h>       /* sqrt */
#include <assert.h>

progCell::progCell(){
    setMap();
    updateState("0");     
    setProgCellParams(1, 0, 1, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0);
    setLayer(0);
}
/*
std::string progCell::strdiff(std::string sbig, std::string ssub){
// return the the chars in s but not in ssub (normally ssub is a subset of set s;)
// the order of chars in string for this function does not matter
    assert(ssub.length() < sbig.length());
    char *ch = new char [sbig.length() - ssub.length()];
    int counter = 0;
    for(int i = 0; i < sbig.length(); i++){
        if(ssub.find(sbig[i]) == std::string::npos){
        // not found
            ch[counter] = sbig[i];
            counter ++;
        }
    }
    std::string strdiff = ch;
    delete [] ch;
    return strdiff;
}

*/
void progCell::setProgCellParams(double cellNumber, double pDeath, double pDiv, double pM1, double pM2, double pM3, double pM4, double pM5,
  double pDivSym, double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff, double pDivAsymDiff){
    initCellNumberNext(); 
    setCellParams(cellNumber, pDeath, pDiv, pM1, pM2, pM3, pM4, pM5);
    m_pDivSym = pDivSym;
    m_pDivAsym = pDivAsym;
    m_pDivAsymSingle = pDivAsymSingle;
    m_pDivAsymDouble = pDivAsymDouble;
    m_pDivSymDiff = pDivSymDiff;
    m_pDivAsymDiff = pDivAsymDiff;
    //assert((m_pDiv == 1) && (m_pDivAsym == 1));
  }

void progCell::setLayer(int layer){
    m_layer = layer;    
}

void progCell::initCellNumberNext(){
    m_cellNumberNextM = new double [nsetState];
    m_cellNumberNextM_NextLayer = new double [nsetState];
    m_cellNumberNextM_T = new double [nsetState];
    for(int i = 0; i < nsetState; i++){
         m_cellNumberNextM[i] = 0;
         m_cellNumberNextM_NextLayer[i] = 0;
         m_cellNumberNextM_T[i] = 0;
     }
}

void progCell::updateProbDiv(double pDeath, double pDiv, double pDivSym, double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff, double pDivAsymDiff){
     m_pDeath = pDeath;
     m_pDiv = pDiv;
     //m_pStay = 1.0 - pDeath - pDiv;
    
     m_pDivSym = pDivSym;
     m_pDivAsym = pDivAsym;
     m_pDivAsymSingle = pDivAsymSingle;
     m_pDivAsymDouble = pDivAsymDouble;
     m_pDivSymDiff = pDivSymDiff;
     m_pDivAsymDiff = pDivAsymDiff;
     m_pStay = 1.0 - m_pDeath - m_pDivSym - m_pDivAsymSingle - m_pDivAsymDouble - m_pDivSymDiff - m_pDivAsymDiff;
     //assert((m_pDiv == 1) && (m_pDivAsym == 1) && (m_pStay == 0));
  }

void progCell::updateTurnover(int turnoverTime_P_u, int turnoverTime_T_u){
     turnoverTime_P = turnoverTime_P_u;
     turnoverTime_T = turnoverTime_T_u;
}

double progCell::getChildrenTermCellNumber(){
    double total = 0;
    for(int i = 0; i < nsetState; i++){
        total += m_cellNumberNextM_T[i];
    }
    return  total;
}

void progCell::getCellNumberByBehavior(unsigned int s){
    getCellNumberByBehavior_poi(s);
    
}


void progCell::getCellNumberByBehavior_poi(unsigned int s){
    assert((m_pDiv == 1) && (m_pDivAsym == 1) && (m_pStay == 0));
    double p1[7] = {m_pDeath, m_pStay, m_pDivSym, m_pDivAsymSingle,
                    m_pDivAsymDouble, m_pDivSymDiff, m_pDivAsymDiff};
    unsigned long long n1[7] = {0};

    for(int i = 0; i < nsetState; i++){
        p2[i] = m_pM[i];
        p3[i] = m_pM[i];
        n2[i] = 0;
        n3[i] = 0;
    }
    /*  comment out to allow more flexible divition pattern from last layer prog cell to terminal cell
    if (m_layer == number_max_layer_P) {
        for(int i = 0; i < nsetState; i++){
         m_cellNumberNextM[i] = 0;
         m_cellNumberNextM_NextLayer[i] = 0;
         m_cellNumberNextM_T[i] = 0;
     }
        m_cellNumberNextM_T[0] = m_cellNumber * (1 - p1[0]);
        return;
    }
    */
    

  // r is the random number generator, s is a user specified seed
  // --- std::default_random_engine generator;
  // --- generator = std::default_random_engine(s);
  // --- std::poisson_distribution<unsigned long long> distribution;
  // cell Numbers that has different cell behaviors at one day
    
    /*
    for(int i = 0; i < 4; i++){
        distribution = std::poisson_distribution<unsigned long long>(m_cellNumber * p1[i]);
        n1[i] = distribution(generator);
        }
    */
    for(int i = 0; i < 7; i++){
        n1[i] = m_cellNumber * p1[i];
        }
    
    m_cellNumberDeath = n1[0];
    m_cellNumberStay = n1[1];
    m_cellNumberSymDiv = n1[2];
    m_cellNumberAsymDivSingle = n1[3];
    m_cellNumberAsymDivDouble = n1[4];
    m_cellNumberSymDiff = n1[5];
    m_cellNumberAsymDiff = n1[6];
  // generate the cellNumberNext (only for progenitor daughters)
  //  m_cellNumberNext = m_cellNumberStay + 2 * m_cellNumberSymDiv + 1 * m_cellNumberAsymDivSingle + 1 * m_cellNumberAsymDivDouble;
  // among the cellNumberNext, what are the Numbers belonging to each cell mutation state?
  
  // mutations on progenitor daughters
    unsigned long long N2 = 2 * m_cellNumberSymDiv + 1 * m_cellNumberAsymDivSingle + 1 * m_cellNumberAsymDivDouble; // only with division is mutation acquisition possible
    // N2 denote the same layer next gen cells
    //double *p2 = new double [nsetState];
    //double *p3 = new double [nsetState];
    //unsigned long long *n2 = new unsigned long long [nsetState];
    //unsigned long long *n3 = new unsigned long long [nsetState];
    
    
    unsigned long long small_total = 0;
   // ---  generator = std::default_random_engine(s);
    for(int i = 1; i < nsetState; i++){
    // +++ 
        //generator = std::default_random_engine(s + i - 1); // --- s+i-1
        distribution = std::poisson_distribution<unsigned long long>(N2 * p2[i]);
        n2[i] = distribution(generator);
        small_total += n2[i];
    }
    assert(N2 >= small_total);
    assert(setState2ind["0"] == 0);
    n2[0] = N2 - small_total;
    m_cellNumberNextM[0] = n2[0] + m_cellNumberStay;
    for(int i = 1; i < nsetState; i++){
        m_cellNumberNextM[i] = n2[i];
    }
    
  
// now start for terminal daughters

unsigned long long N3= 1 * m_cellNumberAsymDivDouble + 2 * m_cellNumberSymDiff + 1 * m_cellNumberAsymDiff; // only with division is mutation acquisition possible
// N3 denote the next layer next gen cells
    
    small_total = 0;
    // --- generator = std::default_random_engine(s); 
    for(int i = 1; i < nsetState; i++){
    // +++ 
        //generator = std::default_random_engine(s + i - 1); // --- s+i-1
        distribution = std::poisson_distribution<unsigned long long>(N3 * p3[i]);
        n3[i] = distribution(generator);
        small_total += n3[i];
    }
    assert(N3 >= small_total);
    n3[0] = N3 - small_total;
    
    m_cellNumberNextM_NextLayer[0] = n3[0];
    
    for(int i = 1; i < nsetState; i++){
        m_cellNumberNextM_NextLayer[i] = n3[i];
    }
    
    //delete [] p2;
    //delete [] p3;
    //delete [] n2;
    //delete [] n3;
    
}
/*
void progCell::getCellNumberByBehavior_mul(unsigned int s){
  std::default_random_engine generator;
  std::discrete_distribution<double> distribution ;
  
  // cell Numbers that has different cell behaviors at one day
    assert((m_pDiv == 1) && (m_pDivAsym == 1) && (m_pStay == 0));
    double p1[6] = {m_pDeath, m_pStay, m_pDivSym, m_pDivAsymSingle,
                    m_pDivAsymDouble, m_pDivSymDiff};
    unsigned long long n1[6] = {0};

    for(int i = 0; i < 6; i++){
        n1[i] = m_cellNumber * p1[i];
        }
    
    m_cellNumberDeath = n1[0];
    m_cellNumberStay = n1[1];
    m_cellNumberSymDiv = n1[2];
    m_cellNumberAsymDivSingle = n1[3];
    m_cellNumberAsymDivDouble = n1[4];
    m_cellNumberSymDiff = n1[5];

  // generate the cellNumberNext (only for progenitor daughters)
    m_cellNumberNext = m_cellNumberStay + 2 * m_cellNumberSymDiv + 1 * m_cellNumberAsymDivSingle + 1 * m_cellNumberAsymDivDouble;
  // among the cellNumberNext, what are the Numbers belonging to each cell mutation state?
  
  // mutations on progenitor daughters
    unsigned long long N2 = 2 * m_cellNumberSymDiv + 1 * m_cellNumberAsymDivSingle + 1 * m_cellNumberAsymDivDouble; // only with division is mutation acquisition possible
    
    //double *p2 = new double [nsetState];
    //double *p3 = new double [nsetState];
    //unsigned long long *n2 = new unsigned long long [nsetState];
    //unsigned long long *n3 = new unsigned long long [nsetState];
    for(int i = 0; i < nsetState; i++){
        p2[i] = m_pM[i];
        p3[i] = m_pM[i];
        n2[i] = 0;
        n3[i] = 0;
    }
    
    std::vector<double> v;
    for(int i = 0; i < nsetState; i++){
        v.push_back(p2[i]);
    }
    
    distribution = std::discrete_distribution<double>(v.begin(), v.end());
    
    for(int i = 1; i < nsetState; i++){
    // +++
        for (int j=0; j<N2; ++j) {
            generator = std::default_random_engine(s+j-1); // --- s+i-1
            int number = distribution(generator);
            ++n2[number];
        }
    }
    m_cellNumberNextM[0] = n2[0] + m_cellNumberStay;
    for(int i = 1; i < nsetState; i++){
        m_cellNumberNextM[i] = n2[i];
    }
    
  
// now start for terminal daughters

unsigned long long N3= 1 * m_cellNumberAsymDivDouble + 2 * m_cellNumberSymDiff; // only with division is mutation acquisition possible
    
    for(int i = 1; i < nsetState; i++){
    // +++
        for (int j=0; j<N3; ++j) {
            generator = std::default_random_engine(s+j-1); // --- s+i-1
            int number = distribution(generator);
            ++n3[number];
        }
    }
    
    for(int i = 0; i < nsetState; i++){
        m_cellNumberNextM_NextLayer[i] = n3[i];
    }
    
    //delete [] p2;
    //delete [] p3;
    //delete [] n2;
    //delete [] n3;
    
}

*/


void progCell::updateCellNumberNext(){
    //double *tempNextM = new double [nsetState];
    //double *tempNextM_NextLayer = new double [nsetState];
    
    for(int i = 0; i < nsetState; i++){
        tempNextM[i] = 0;
        tempNextM_NextLayer[i] = 0;
        tempNextM_T[i] = 0;
    }
    
    std::string mut_state = "";
    std::string next_state = "";
    int i_m_state = 0;
    int i_mut_state = 0;
    int i_next_state = 0; // index for next_state in setState;
    
    i_m_state = setState2ind[m_State];
    for(int i = 0; i < nsetState; i++) // Loop over all mutation state
    {
        mut_state = setState[i];
        //next_state = strunion("1", "2");
        // --- next_state = strunion(m_State, mut_state); // get next state
        
        // +++
        i_mut_state = setState2ind[mut_state];
        next_state = NextStates[i_m_state][i_mut_state];
        
        i_next_state = setState2ind[next_state];
        tempNextM[i_next_state] += m_cellNumberNextM[i];
        tempNextM_NextLayer[i_next_state] += m_cellNumberNextM_NextLayer[i];
        tempNextM_T[i_next_state] += m_cellNumberNextM_T[i];
    }
    
    for(int i = 0; i < nsetState; i++){
        m_cellNumberNextM[i] = tempNextM[i];
        m_cellNumberNextM_NextLayer[i] = tempNextM_NextLayer[i];
        m_cellNumberNextM_T[i] = tempNextM_T[i];
    }
    

    if (m_layer == number_max_layer_P) {
        for(int i = 0; i < nsetState; i++){
         m_cellNumberNextM[i] = 0;
         m_cellNumberNextM_T[i] = m_cellNumberNextM_NextLayer[i];
         m_cellNumberNextM_NextLayer[i] = 0;
        }
    }

    //delete [] tempNextM;
    //delete [] tempNextM_NextLayer;
}







/*
void progCell::updateCellNumberNext(){
    double *tempNextM = new double [nsetState];
    double *tempNextM_NextLayer = new double [nsetState];
    std::string sdiff = "";
    std::string sdiffj = "";
    for(int i = 0; i < nsetState; i++){
        tempNextM[i] = 0;
        tempNextM_NextLayer[i] = 0;
    }
    int istate = setState2ind[m_State];       
    for(int i = 0; i < nsetState; i++){
        if( (setState[i] == "0") || (isSubstr(m_State, setState[i])) ){
            if (i != istate){
                tempNextM[i] = 0;      
                tempNextM_NextLayer[i] = 0;
            }            
            tempNextM[istate] += m_cellNumberNextM[i];      
            tempNextM_NextLayer[istate] += m_cellNumberNextM_NextLayer[i];
        }else if( (!isEqstr(setState[i], m_State) ) && (isSubstr(setState[i], m_State)) ){
        //std::cout<<setState[i]<<" "<<setState[i].length()<<" "<<m_State<<" "<<m_State.length()<<std::endl;
            sdiff = strdiff(setState[i], m_State); // ---
            //std::cout<<setState[i]<<" "<<m_State<<" "<<sdiff<<std::endl;
            //sdiff = "1"; // +++            
            for(int j = 0; j < nsetState; j++){
                if (isEqstr(setState[j], sdiff)){
                   tempNextM[i] += m_cellNumberNextM[j]; 
                   tempNextM_NextLayer[i] += m_cellNumberNextM_NextLayer[j]; 
                
                }else if (isSubstr(setState[j], sdiff)) {
                   sdiffj = strdiff(setState[j], sdiff); // ---
                   //std::cout<<setState[j]<<" "<<sdiff<<std::endl;
                   //sdiffj = "1"; // +++
                   if(isSubstr(m_State, sdiffj)){
                       tempNextM[i] += m_cellNumberNextM[j]; 
                       tempNextM_NextLayer[i] += m_cellNumberNextM_NextLayer[j]; 
                   }
                }
            }         
        }else{
            tempNextM[i] = 0;
            tempNextM_NextLayer[i] = 0;
        }
    }
    
    for(int i = 0; i < nsetState; i++){
        m_cellNumberNextM[i] = tempNextM[i];
        m_cellNumberNextM_NextLayer[i] = tempNextM_NextLayer[i];
    }
    
    delete [] tempNextM;
    delete [] tempNextM_NextLayer;
    
}
*/

bool progCell::isTimeToDivide(){
    return (timeCounter == turnoverTime_P);
}


void progCell::clearCell(){
    for(int i = 0; i < nsetState; i++){
         m_cellNumberNextM[i] = 0;
         m_cellNumberNextM_T[i] = m_cellNumberNextM_NextLayer[i];
         m_cellNumberNextM_NextLayer[i] = 0;
    }
    updateCellNumber(0.0);
}


void progCell::memDelete(){
    delete [] p2;
    delete [] p3;
    delete [] n2;
    delete [] n3;
    delete [] tempNextM;
    delete [] tempNextM_NextLayer;
    delete [] tempNextM_T;
    delete [] m_cellNumberNextM;
    delete [] m_cellNumberNextM_NextLayer;
    delete [] m_cellNumberNextM_T;
}






