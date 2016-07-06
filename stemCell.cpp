#include "Params.h"
#include "Cell.h"
#include "stemCell.h"

#include "strdiff.h"
#include "isSubstr.h"
#include "isEqstr.h"
#include "strunion.h"

//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_cdf.h>
//#include <gsl/gsl_rng.h>


#include <iostream>
#include <string>
#include <random>

#include <map>
//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_cdf.h>
//#include <gsl/gsl_rng.h>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <math.h>       /* sqrt */
#include <assert.h>
//#include <stdlib.h> 


stemCell::stemCell(){
    setMap();
    updateState("0");
    setStemCellParams(1, 0, 1, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0);
}


void stemCell::setStemCellParams(double cellNumber, double pDeath, double pDiv, double pM1, double pM2, double pM3, double pM4, double pM5,
  double pDivSym, double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff){
    initCellNumberNext();
    setCellParams(cellNumber, pDeath, pDiv, pM1, pM2, pM3, pM4, pM5);
    m_pDivSym = pDivSym;
    m_pDivAsym = pDivAsym;
    m_pDivAsymSingle = pDivAsymSingle;
    m_pDivAsymDouble = pDivAsymDouble;
    m_pDivSymDiff = pDivSymDiff;
    //std::cout<<m_pDiv << " " << m_pDivAsym << std::endl;
    //assert((m_pDiv == 1) && (m_pDivAsym == 1));
  }


void stemCell::initCellNumberNext(){
    m_cellNumberNextM = new double [nsetState];
    m_cellNumberNextM_P = new double [nsetState];
    for(int i = 0; i < nsetState; i++){
         m_cellNumberNextM[i] = 0;
         m_cellNumberNextM_P[i] = 0;
     }
}



void stemCell::updateProbDiv(double pDeath, double pDiv, double pDivSym, double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff){
     m_pDeath = pDeath;
     m_pDiv = pDiv;
     //m_pStay = 1.0 - pDeath - pDiv;
     m_pDivSym = pDivSym;
     m_pDivAsym = pDivAsym;
     m_pDivAsymSingle = pDivAsymSingle;
     m_pDivAsymDouble = pDivAsymDouble;
     m_pDivSymDiff = pDivSymDiff;
     m_pStay = 1.0 - m_pDeath - m_pDivSym - m_pDivAsymSingle - m_pDivAsymDouble - m_pDivSymDiff;
     //assert((m_pDiv == 1) && (m_pDivAsym == 1) && (m_pStay == 0));
  }
  
void stemCell::updateTurnover(int turnoverTime_S_u){
    turnoverTime_S = turnoverTime_S_u;
}



void stemCell::getCellNumberByBehavior(unsigned int s){
    //getCellNumberByBehavior_mul(s);
    getCellNumberByBehavior_poi(s);
    //getCellNumberByBehavior_mul_gsl(s);
}



/*
void stemCell::getCellNumberByBehavior_mul(unsigned int s){
  // multinomial distribution


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
    //std::cout<<"m_State "<<m_State<<" m_cellNumber "<<m_cellNumber<<std::endl;
    
    m_cellNumberDeath = n1[0];
    m_cellNumberStay = n1[1];
    m_cellNumberSymDiv = n1[2];
    m_cellNumberAsymDivSingle = n1[3];
    m_cellNumberAsymDivDouble = n1[4];
    m_cellNumberSymDiff = n1[5];
    //std::cout<<"NUMBER:: "<<m_cellNumberAsymDivDouble<<std::endl;
  // generate the cellNumberNext (only for stem daughters)
    m_cellNumberNext = m_cellNumberStay + 2 * m_cellNumberSymDiv + 1 * m_cellNumberAsymDivSingle + 1 * m_cellNumberAsymDivDouble;
  // among the cellNumberNext, what are the Numbers belonging to each cell mutation state?
  // mutations on stem daughters
    unsigned long long N2 = 2 * m_cellNumberSymDiv + 1 * m_cellNumberAsymDivSingle + 1 * m_cellNumberAsymDivDouble; // only with division is mutation acquisition possible

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

    for(int i = 0; i < nsetState; i++){
        //std::cout<<i<<std::endl;
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
    
   
unsigned long long N3= 1 * m_cellNumberAsymDivDouble + 2 * m_cellNumberSymDiff; // only with division is mutation acquisition possible
//std::cout<<"NUMBER:: "<<N3<<std::endl;
//std::cout<<"NUMBER:: "<<N3<<std::endl;

    // ++
    for(int i = 0; i < nsetState; i++){
        for (int j=0; j<N3; ++j) {
            generator = std::default_random_engine(s+j-1); // --- s+i-1
            int number = distribution(generator);
            ++n3[number];
        }
    }
    for(int i = 0; i < nsetState; i++){
        m_cellNumberNextM_P[i] = n3[i];
    }
    
}

*/

void stemCell::getCellNumberByBehavior_poi(unsigned int s){
  // r is the random number generator, s is a user specified seed
  // --- std::default_random_engine generator;
  // --- generator = std::default_random_engine(s);
  // --- std::poisson_distribution<unsigned long long> distribution;

  // cell Numbers that has different cell behaviors at one day
    assert((m_pDiv == 1) && (m_pDivAsym == 1) && (m_pStay == 0));
    double p1[6] = {m_pDeath, m_pStay, m_pDivSym, m_pDivAsymSingle,
                    m_pDivAsymDouble, m_pDivSymDiff};
    unsigned long long n1[6] = {0};

    for(int i = 0; i < 6; i++){
        n1[i] = m_cellNumber * p1[i];
        }
    //std::cout<<"m_State "<<m_State<<" m_cellNumber "<<m_cellNumber<<std::endl;
    
    m_cellNumberDeath = n1[0];
    m_cellNumberStay = n1[1];
    m_cellNumberSymDiv = n1[2];
    m_cellNumberAsymDivSingle = n1[3];
    m_cellNumberAsymDivDouble = n1[4];
    m_cellNumberSymDiff = n1[5];
    //std::cout<<"NUMBER:: "<<m_cellNumberAsymDivDouble<<std::endl;
  // generate the cellNumberNext (only for stem daughters)
    m_cellNumberNext = m_cellNumberStay + 2 * m_cellNumberSymDiv + 1 * m_cellNumberAsymDivSingle + 1 * m_cellNumberAsymDivDouble;
  // among the cellNumberNext, what are the Numbers belonging to each cell mutation state?
    //std::cout<<"m_State "<<m_State<<" m_cellNumberNext "<<m_cellNumberNext<<std::endl;
  // mutations on stem daughters
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

    unsigned long long small_total = 0;
    
    //std::cout<<"setState2ind_5 "<<setState2ind["5"]<<std::endl;
    // +++ 
   // std::cout<<" N2 "<<N2<<std::endl;
    // --- generator = std::default_random_engine(s); // +++
    for(int i = 1; i < nsetState; i++){
    // +++ 
       // std::cout<<" seed "<<s+i-1;
        
        //generator = std::default_random_engine(s + i - 1); // --- s+i-1
        distribution = std::poisson_distribution<unsigned long long>(N2 * p2[i]);
        
        n2[i] = distribution(generator);
        //if(N2 * p2[i] > 1){
     //   std::cout<<"poisson mean "<<i<<" "<<N2 * p2[i] * 1e9<<std::endl;
        //}
        small_total += n2[i];
    }
    //std::cout<<"m_State "<<m_State<<" N2 "<<N2<<" small total "<<small_total<<std::endl;
    // --- 
    assert(N2 >= small_total);
    assert(setState2ind["0"] == 0);
    n2[0] = N2 - small_total;
    m_cellNumberNextM[0] = n2[0] + m_cellNumberStay;
    for(int i = 1; i < nsetState; i++){
        m_cellNumberNextM[i] = n2[i];
    }
    
    // ---- std::cout<<"cell state: "<<m_State<<" cell number "<<m_cellNumber<<" Next M12 "<<m_cellNumberNextM[setState2ind["12"]]<<std::endl;
//    std::cout<<" Next M0 "<<m_cellNumberNextM[0]<<std::endl;
//    std::cout<<" Next M1 "<<m_cellNumberNextM[1]<<std::endl;
//std::cout<<" m_cellNumberNextM of M3 "<<m_cellNumberNextM[setState2ind["3"]]<<std::endl;    
// now start for progenitor daughters

unsigned long long N3= 1 * m_cellNumberAsymDivDouble + 2 * m_cellNumberSymDiff; // only with division is mutation acquisition possible
//std::cout<<"NUMBER:: "<<N3<<std::endl;
//std::cout<<"NUMBER:: "<<N3<<std::endl;

    small_total = 0;
    // ++ 
    // --- generator = std::default_random_engine(s); 
    for(int i = 1; i < nsetState; i++){
    // +++ 
    //generator = std::default_random_engine(s + i - 1); // --- s+i-1
        distribution = std::poisson_distribution<unsigned long long>(N3 * p3[i]);
        n3[i] = distribution(generator);
        small_total += n3[i];
    }
    // --- 
    assert(N3 >= small_total);
    n3[0] = N3 - small_total;
    
    m_cellNumberNextM_P[0] = n3[0];
    
    for(int i = 1; i < nsetState; i++){
        m_cellNumberNextM_P[i] = n3[i];
    }
    
    //delete [] p2;
    //delete [] p3;
    //delete [] n2;
    //delete [] n3;
}



/*

void stemCell::getCellNumberByBehavior_mul_gsl(unsigned int s){
  gsl_rng * r;
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  size_t K = nsetState;
  gsl_rng_set(r, s);


  // cell Numbers that has different cell behaviors at one day
    assert((m_pDiv == 1) && (m_pDivAsym == 1) && (m_pStay == 0));
    double p1[6] = {m_pDeath, m_pStay, m_pDivSym, m_pDivAsymSingle,
                    m_pDivAsymDouble, m_pDivSymDiff};
    unsigned long long n1[6] = {0};

    for(int i = 0; i < 6; i++){
        n1[i] = m_cellNumber * p1[i];
        }
    //std::cout<<"m_State "<<m_State<<" m_cellNumber "<<m_cellNumber<<std::endl;
    
    m_cellNumberDeath = n1[0];
    m_cellNumberStay = n1[1];
    m_cellNumberSymDiv = n1[2];
    m_cellNumberAsymDivSingle = n1[3];
    m_cellNumberAsymDivDouble = n1[4];
    m_cellNumberSymDiff = n1[5];
    //std::cout<<"NUMBER:: "<<m_cellNumberAsymDivDouble<<std::endl;
  // generate the cellNumberNext (only for stem daughters)
    m_cellNumberNext = m_cellNumberStay + 2 * m_cellNumberSymDiv + 1 * m_cellNumberAsymDivSingle + 1 * m_cellNumberAsymDivDouble;
  // among the cellNumberNext, what are the Numbers belonging to each cell mutation state?
    //std::cout<<"m_State "<<m_State<<" m_cellNumberNext "<<m_cellNumberNext<<std::endl;
  // mutations on stem daughters
    unsigned long long N2 = 2 * m_cellNumberSymDiv + 1 * m_cellNumberAsymDivSingle + 1 * m_cellNumberAsymDivDouble; // only with division is mutation acquisition possible

    for(int i = 0; i < nsetState; i++){
        p2[i] = m_pM[i];
        p3[i] = m_pM[i];
        n2[i] = 0;
        n3[i] = 0;
        n2_gsl[i] = 0;
        n3_gsl[i] = 0;
        
    }
   // unsigned int n2_gsl[32] = {0};
   // unsigned int n3_gsl[32] = {0};
    
    gsl_ran_multinomial(r, K, N2, p2, n2_gsl);
    
    m_cellNumberNextM[0] = n2_gsl[0] + m_cellNumberStay;
    for(int i = 1; i < nsetState; i++){
        m_cellNumberNextM[i] = n2_gsl[i];
    }
    


unsigned long long N3= 1 * m_cellNumberAsymDivDouble + 2 * m_cellNumberSymDiff; // only with division is mutation acquisition possible
    gsl_ran_multinomial(r, K, N3, p3, n3_gsl);
    for(int i = 0; i < nsetState; i++){
        m_cellNumberNextM_P[i] = n3_gsl[i];
    }
gsl_rng_free (r);
}

*/





void stemCell::updateCellNumberNext(){

    
    for(int i = 0; i < nsetState; i++){
        tempNextM[i] = 0;
        tempNextM_P[i] = 0;
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
        tempNextM_P[i_next_state] += m_cellNumberNextM_P[i];
    }
    
    for(int i = 0; i < nsetState; i++){
        m_cellNumberNextM[i] = tempNextM[i];
        m_cellNumberNextM_P[i] = tempNextM_P[i];
    }
    


}





bool stemCell::isTimeToDivide(){
    return (timeCounter == turnoverTime_S) ;
}

void stemCell::memDelete(){
    delete [] p2;
    delete [] p3;
    delete [] n2;
    delete [] n3;
    delete [] tempNextM;
    delete [] tempNextM_P;
    delete [] m_cellNumberNextM;
    delete [] m_cellNumberNextM_P;
}












