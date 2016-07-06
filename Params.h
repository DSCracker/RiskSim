#ifndef PARAMS_H
#define PARAMS_H

#include <iostream>
#include <string>
#include <map>
#include <set>
#include <random>
//#include <map>

class Params
{
protected:
    std::string setState[32] = {"0", "1", "2", "3", "4", "5",
                             "12", "13", "14", "15", "23", "24", "25", "34", "35", "45",
                             "123", "124", "125", "134", "135", "145", "234", "235", "245", "345",
                             "1234", "1235", "1245", "1345", "2345", "12345"};
   
    std::string setMut[5] = {"1","2","3","4","5"};
    std::map<std::string, int> setState2ind;
    std::map<std::string, int> setMut2ind;
    
    int nsetState = 0;
    int nsetMut = 0;
    
    std::string NextStates[32][32];
    
    std::default_random_engine generator;
    std::poisson_distribution<unsigned long long> distribution;
    
                            
    // --
    int version;
    // --- --
    std::string tissue_name;
    int nrolls;
    double stemCellNumberHomeo = 0.0;
    double totalCellNumberHomeo = 0.0;

    int number_div_sym = 0;
    int number_div_asym = 0;
    int number_div_P = 0;
    int number_max_layer_P = 0; // +++ 

    int turnoverTime_S = 0;
    int turnoverTime_P = 0;
    int turnoverTime_T = 0;

    double pDeath = 0.0;
    double pDeath_Asym = 0.0;
    const double pDiv = 1.0;
    double pDivSym = 0.0;
    double pDivSym_Asym = 0.0;
    const double pDivAsym = 1.0;
    const double pDivAsym_Asym = 1.0;
    double pDivAsymSingle = 0.0;
    double pDivAsymDouble = 0.0;
    double pDivSymDiff = 0.0;
    double pDivAsymSingle_Asym = 0.0;
    double pDivAsymDouble_Asym = 0.0;
    double pDivSymDiff_Asym = 0.0;

    double pM1 = 0.0;
    double pM2 = 0.0;
    double pM3 = 0.0;
    double pM4 = 0.0;
    double pM5 = 0.0;

    double pDeath_P = 0.0;
    const double pDiv_P = 1.0;
    double pDivSym_P = 0.0; // prob of layer to two layer
    const double pDivAsym_P = 1.0; 
    double pDivAsymSingle_P = 0.0; // prob of layer to one layer
    double pDivAsymDouble_P = 0.0; // prob of layer to one layer and one layer + 1
    double pDivSymDiff_P = 0.0; // prob of layer to two layer + 1
    double pDivAsymDiff_P = 0.0; // prob of layer to one layer + 1

    double pM1_P = 0.0;
    double pM2_P = 0.0;
    double pM3_P = 0.0;
    double pM4_P = 0.0;
    double pM5_P = 0.0;

    double sltAdv = 0.0;
    double mfacDeathRate = 0.0;
    double mfacMutRate = 0.0;
    double mfacTurnoverTime = 1.0;

    double sltAdv_P = 0.0;
    double mfacDeathRate_P = 0.0;
    double mfacMutRate_P = 0.0;
    double mfacTurnoverTime_P = 1.0;
    
    std::string dynamics_effect_option = "M1";
    
    // --- --
    int number_div_total;
    const int total_days = 29200;
    // ---
public:
    Params();
    void setNextStates();
    void setMap();
    void setVersion();
    void setParams();

};



#endif
