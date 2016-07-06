
#include "Params.h"
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include "strunion.h"

Params::Params(){
    setParams();
}
void Params::setMap(){
    nsetState = sizeof(setState) / sizeof(setState[0]);
    nsetMut = sizeof(setMut) / sizeof(setMut[0]); 
    //std::cout<<"nsetState  "<<nsetState<<std::endl;
    //std::cout<<"nsetMut  "<<nsetMut<<std::endl;
    for(int i = 0; i < nsetState; i++){
       setState2ind[setState[i]] = i;
    }
    for(int i = 0; i < nsetMut; i++){
       setMut2ind[setMut[i]] = i;
    }
    //std::cout<<"nsetState  "<<nsetState<<std::endl;
    //std::cout<<"nsetMut  "<<nsetMut<<std::endl;
}



void Params::setNextStates(){
    for(int i = 0; i < 32; i++){
        for(int j = 0; j < 32; j++){
            NextStates[i][j] = strunion(setState[i], setState[j]);
        }
    }
}

void Params::setVersion(){
    std::ifstream fin("../input/version.txt");
    std::string name;
    int var;
    while (fin >> name >> var)
    {
        if(name == "version"){version = var;}
    }
}

void Params::setParams(){
    setMap();
    setNextStates();
    setVersion();
    std::ifstream fin("../input/params"+std::to_string(version)+".txt");
    std::string name;
    double var;
    while (fin >> name >> var)
    {
        if(name == "nrolls"){nrolls = var;}
        else if(name == "stemCellNumberHomeo") {stemCellNumberHomeo = var;}
        else if(name == "totalCellNumberHomeo"){totalCellNumberHomeo = var;}
        else if(name == "number_div_sym"){number_div_sym = var;}
        else if(name == "number_div_asym"){number_div_asym = var;}
        else if(name == "number_div_P"){number_div_P = var;}
        else if(name == "turnoverTime_S"){turnoverTime_S = var;}
        else if(name == "turnoverTime_P"){turnoverTime_P = var;}
        else if(name == "turnoverTime_T"){turnoverTime_T = var;}
        else if(name == "pDeath"){pDeath = var;}
        else if(name == "pDeath_Asym"){pDeath_Asym = var;}
        else if(name == "pDivSym"){pDivSym = var;}
        else if(name == "pDivSym_Asym"){pDivSym_Asym = var;}
        else if(name == "pDivAsymSingle"){pDivAsymSingle = var;}
        else if(name == "pDivAsymDouble"){pDivAsymDouble = var;}
        else if(name == "pDivSymDiff"){pDivSymDiff = var;}
        else if(name == "pDivAsymSingle_Asym"){pDivAsymSingle_Asym = var;}
        else if(name == "pDivAsymDouble_Asym"){pDivAsymDouble_Asym = var;}
        else if(name == "pDivSymDiff_Asym"){pDivSymDiff_Asym = var;}
        else if(name == "pM1"){pM1 = var;}
        else if(name == "pM2"){pM2 = var;}
        else if(name == "pM3"){pM3 = var;}
        else if(name == "pM4"){pM4 = var;}
        else if(name == "pM5"){pM5 = var;}
        else if(name == "pDeath_P"){pDeath_P = var;}
        else if(name == "pDivSym_P"){pDivSym_P = var;}
        else if(name == "pDivAsymSingle_P"){pDivAsymSingle_P = var;}
        else if(name == "pDivAsymDouble_P"){pDivAsymDouble_P = var;}
        else if(name == "pDivSymDiff_P"){pDivSymDiff_P = var;}
        else if(name == "pDivAsymDiff_P"){pDivAsymDiff_P = var;}
        else if(name == "pM1_P"){pM1_P = var;}
        else if(name == "pM2_P"){pM2_P = var;}
        else if(name == "pM3_P"){pM3_P = var;}
        else if(name == "pM4_P"){pM4_P = var;}
        else if(name == "pM5_P"){pM5_P = var;}
        else if(name == "sltAdv"){sltAdv = var;}
        else if(name == "mfacDeathRate"){mfacDeathRate = var;}
        else if(name == "mfacMutRate"){mfacMutRate = var;}
        else if(name == "mfacTurnoverTime"){mfacTurnoverTime = var;}
        else if(name == "sltAdv_P"){sltAdv_P = var;}
        else if(name == "mfacDeathRate_P"){mfacDeathRate_P = var;}
        else if(name == "mfacMutRate_P"){mfacMutRate_P = var;}
        else if(name == "mfacTurnoverTime_P"){mfacTurnoverTime_P = var;}
        else if(name == "dynamics_effect_option"){
            if(var == 1){
                dynamics_effect_option = "M1";
            }else if(var == 13){
                dynamics_effect_option = "M1M3";
            }
        }
        else{tissue_name = name;}
    }
        number_div_total = number_div_sym + number_div_asym;
        number_max_layer_P = number_div_P ;

        turnoverTime_S = (365.0 * 80.0)/ (double)number_div_total;
        turnoverTime_P = (365.0 * 80.0)/ ((double)number_div_total * (double)number_div_P);
        turnoverTime_T = -1;
}


