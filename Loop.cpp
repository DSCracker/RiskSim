#include "Params.h"
#include "Loop.h"
#include <stdexcept>

#include <iostream>
#include <string>
#include <random>
#include <fstream>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <algorithm>    // std::min

Loop::Loop(){
    // set parameter and initialize cells, seed and count time
    setParams();
    sim.SetEffectOption(dynamics_effect_option); // whether the div effect on M1 or M3
    initCountCellNumber(); // init per natural day cell numbers
    initcountCancerGeneration(); // init array to record generation wise cancer incidence
    initcountMCellGeneration();
    initCell(); // init stem and prog cells
    s0 = 0; // random seed init
    s = 0; // random seed
    t = clock(); // time counter
    //std::cout<<" MChar "<<MChar<<std::endl;
}

void Loop::initCountCellNumber(){ // init per natural day cell numbers
    countStemNonMNumber = new double[365 * 80];
    countStemMNumber = new double[365 * 80];
    countStemNumber = new double[365 * 80];
    countProgNonMNumber = new double[365 * 80];
    countProgMNumber = new double[365 * 80];
    countProgNumber = new double[365 * 80];
    countTermNumber = new double[365 * 80];
    //countTempNumber = new double[365 * 80];

    for(int i = 0; i < 365 * 80; ++i){
        countStemNonMNumber[i] = 0;
        countStemMNumber[i] = 0;
        countStemNumber[i] = 0;
        countProgNonMNumber[i] = 0;
        countProgMNumber[i] = 0;
        countProgNumber[i] = 0;
        countTermNumber[i] = 0;
        //countTempNumber[i] = 0;        
    }
} 

void Loop::computeCountCellNumber(double value, int startTp, int endTp, std::string cellType){
    endTp = std::min(endTp, (365 * 80 - 1));
    if(cellType == "stemNonM"){
        for(int i = startTp; i < endTp; i++){
            countStemNonMNumber[i] += value;
        }
    }else if(cellType == "stemM"){
        for(int i = startTp; i < endTp; i++){
            countStemMNumber[i] += value;
        }
    }else if(cellType == "progNonM"){
        for(int i = startTp; i < endTp; i++){
            countProgNonMNumber[i] += value;
        }
    }else if(cellType == "progM"){
        for(int i = startTp; i < endTp; i++){
            countProgMNumber[i] += value;
        }
    }else if(cellType == "term"){
        for(int i = startTp; i < endTp; i++){
            countTermNumber[i] += value;
        }
    }else{
        throw std::invalid_argument( "invalid cellType" );
    }
}

void Loop::initcountCancerGeneration(){
    countCancerGeneration = new int[number_div_total];
    for(int i = 0; i < number_div_total; ++i){
        countCancerGeneration[i] = 0;
    }
}

void Loop::initcountStemNumberGeneration(){
    countStemNumberGeneration = new int[number_div_total];
    for(int i = 0; i < number_div_total; ++i){
        countStemNumberGeneration[i] = 0;
    }
}

void Loop::initcountMCellGeneration(){
    countMCellGeneration = new int[number_div_total];
    for(int i = 0; i < number_div_total; ++i){
        countMCellGeneration[i] = 0;
    }
}


// ++++++ 
void Loop::initCell(){
    sim.SetStem(1.0, pDeath, pDiv, pM1, pM2, pM3, pM4, pM5, pDivSym, pDivAsym, 
    pDivAsymSingle, pDivAsymDouble, pDivSymDiff, sltAdv, mfacDeathRate, mfacMutRate);
    if (number_max_layer_P >= 0){
        for(int il = 0; il <= number_max_layer_P; il++){
            sim.SetProg(0.0, pDeath_P, pDiv_P, pM1_P, pM2_P, pM3_P, pM4_P, pM5_P, pDivSym_P, pDivAsym_P, pDivAsymSingle_P, pDivAsymDouble_P, pDivSymDiff_P, pDivAsymDiff_P,
                sltAdv_P, mfacDeathRate_P, mfacMutRate_P, il);
        }
    }
}

void Loop::updateOneStemGeneration(){
    if(iterStemGen >= number_div_sym){ 
        // asym
        sim.SetStemDiv(pDeath_Asym, pDiv, pDivSym_Asym, pDivAsym_Asym, pDivAsymSingle_Asym, pDivAsymDouble_Asym,
         pDivSymDiff_Asym, sltAdv, mfacDeathRate);
    }else{
        // sym
        sim.SetStemDiv(pDeath, pDiv, pDivSym, pDivAsym, pDivAsymSingle, pDivAsymDouble, pDivSymDiff, sltAdv, mfacDeathRate);
    }
    s = s0 + iterStemGen * 200; // ---+++ 31
    
    // compute and add in stem Non M cell number iterStemGen to iterStemGen + 1

    double value = sim.getTotalStemCellNumber("stemNonM");
    std::cout<<" globel stem iter: "<<iterStemGen<<" stemNonM cell num "<<value<<std::endl;
    int startTp = iterStemGen * turnoverTime_S;
    int endTp = startTp + turnoverTime_S;
    computeCountCellNumber(value, startTp, endTp, "stemNonM");

    // M stem cells propogate for several more iterations from global generation i to i + 1 
    int num_stem_div_left_sym = 0, num_stem_div_left_asym = 0;
    if(iterStemGen <= number_div_sym - 1){
        num_stem_div_left_sym = 1;
        num_stem_div_left_asym = 0;
    }else{
        num_stem_div_left_sym = 0;
        num_stem_div_left_asym = 1 * ceil(1.0/mfacTurnoverTime);
    }
    //std::cout<<" globel stem iter: "<<iterStemGen<<" prog init cell num "<<sim.getTotalProgCellNumber_byLayer(0, "progM")<<std::endl;
    loopStemCellMLifeTime(num_stem_div_left_sym, num_stem_div_left_asym); // loop M stem cell lineage

    // Non-M stem cells propogate for one iteration from global generation i to i + 1
    sim.StemDivToNext(s, "stemNonM");
    // aggregate
    sim.StemNumUpdateFromPrev("stemNonM", "stemNonM");
    sim.StemNumUpdateFromPrev("stemM", "stemNonM");

    // If prog cell included, add to prog cell class prog children of Non-M stem cell    
    if(number_max_layer_P >= 0){
        sim.ProgNumUpdateFromPrev(0, 0, "all", "stemNonM"); 
    }
    // ::: check if cancer cell occur at gen i+1 ::: //
    if(isCancer()){
        isFirstCancerCell = true;
    }

    if(number_max_layer_P >= 0){
        //loopProgCellLifeTime_simple(); // loop prog Cell lineage 
        loopProgCellLifeTime();
    }  
    //std::cout<<" globel stem iter: "<<iterStemGen<<" prog init cell num after clear"<<sim.getTotalProgCellNumber_byLayer(0, "progM")<<std::endl;
    if(isFirstCancerCell){
        genFirstCancerCell = std::min(iterStemGen,genFirstCancerCell);
    }
   // ---  sim.cellNumberDebugger();   
    // ---
}


void Loop::loopStemCellMLifeTime(int num_stem_div_left_sym, int num_stem_div_left_asym){

    for(int ngenM = 0; ngenM < num_stem_div_left_sym + num_stem_div_left_asym; ngenM++){ // sym M stem cell iteration
        iterStemMGen = ngenM;
        // set sym and asym params
        if(ngenM >= num_stem_div_left_sym){
        // asym
        sim.SetStemMDiv(pDeath_Asym, pDiv, pDivSym_Asym, pDivAsym_Asym, pDivAsymSingle_Asym, 
            pDivAsymDouble_Asym, pDivSymDiff_Asym, sltAdv, mfacDeathRate);
        }else{
        // sym
        sim.SetStemMDiv(pDeath, pDiv, pDivSym, pDivAsym, pDivAsymSingle, pDivAsymDouble, 
            pDivSymDiff, sltAdv, mfacDeathRate);
        }
        // s = s0 + (iterStemGen + ngenM) * 31; // ---+++ 31 --- +1
        s = s0 + (iterStemGen) * 31 + ngenM * 10000;

        // compute and add in stem  M cell number iterStemGen + ngenM to iterStemGen + ngenM + 1
        double value = sim.getTotalStemCellNumber("stemM");
        int startTp = iterStemGen * turnoverTime_S + ngenM * turnoverTime_S * 1.0 / ceil(1.0/mfacTurnoverTime);
        int endTp = startTp + turnoverTime_S * 1.0 / ceil(1.0/mfacTurnoverTime); 
        computeCountCellNumber(value, startTp, endTp, "stemM");

        // propogate
        sim.StemDivToNext(s, "stemM");
        // aggregate 
        sim.StemNumUpdateFromPrev("stemM", "stemM");
        //std::cout<<"global stem iter: "<<iterStemGen << " stem M cell num: "<<sim.getTotalStemCellMNumber()<<std::endl;
        //sim.updateStemCellMNumberFromStemM(s);
        if(number_max_layer_P >= 0){
            sim.ProgNumUpdateFromPrev(0, 0, "progM", "stemM"); 
            //sim.updateProgCellMNumberFromStemM();// get initial M prog children from stem
        }
        //std::cout<<" globel stem iter: "<<iterStemGen<<" prog init cell num "<<sim.getTotalProgCellNumber_byLayer(0, "progM")<<std::endl;
        if(isCancer()){
            isFirstCancerCell = true;
        }

        if(number_max_layer_P >= 0){
            if(ngenM < num_stem_div_left_sym + num_stem_div_left_asym - 1){
                loopProgCellMLifeTime(); // loop M prog Cell lineage
            }
        }
    }
}

void Loop::loopProgCellLifeTime(){
// update division params
    pDeath_P = 0.0;
    pDivSym_P = 0.0;
    pDivAsymSingle_P = 0.0;
    pDivAsymDouble_P = 0.0; //prob that one il generates one il and one (il + 1)
    pDivSymDiff_P = 1.0; // prob that
    pDivAsymDiff_P = 0.0; 
    for(int il = 0; il <= number_max_layer_P; il++){
        sim.SetProgDiv(pDeath_P, pDiv_P, pDivSym_P, pDivAsym_P, pDivAsymSingle_P, pDivAsymDouble_P, pDivSymDiff_P, pDivAsymDiff_P,
            sltAdv_P, mfacDeathRate_P, il);
    }
// loop
    int num_prog_div_left = 1 * ceil(1.0 / mfacTurnoverTime_P);
    int count = 0;
    while (sim.getTotalProgCellNumber("all") > 0){

        // compute and add in prog Non M Cell
        double value = sim.getTotalProgCellNumber("progNonM");
        int startTp = iterStemGen * turnoverTime_S + count * turnoverTime_P;
        int endTp = startTp + turnoverTime_P;
        computeCountCellNumber(value, startTp, endTp, "progNonM");
        
        // only operate on current generation prog cells
        int temp = number_max_layer_P;
        if(temp > count){
            temp = count;
        }
              
        // propogate at M cell faster sub-generation
        for(int j = 0; j < num_prog_div_left; j++){ // j is the M prog cell sub-generation 
            if(j < num_prog_div_left - 1){
                pDeath_P = 0.0;
                pDivSym_P = 1.0;
                pDivAsymSingle_P = 0.0;
                pDivAsymDouble_P = 0.0;
                pDivSymDiff_P = 0.0;
                pDivAsymDiff_P = 0.0; 
            }else{
                pDeath_P = 0.0;
                pDivSym_P = 0.0;
                pDivAsymSingle_P = 0.0;
                pDivAsymDouble_P = 0.0;
                pDivSymDiff_P = 1.0;
                pDivAsymDiff_P = 0.0; 
            }                                      
            // compute and add in prog M Cell
            double value = sim.getTotalProgCellNumber("progM");
            int startTp = iterStemGen * turnoverTime_S + count * turnoverTime_P + 
                            j * turnoverTime_P * 1.0 / ceil(1.0 / mfacTurnoverTime_P);
            int endTp = startTp + turnoverTime_P * 1.0 / ceil(1.0 / mfacTurnoverTime_P);
            computeCountCellNumber(value, startTp, endTp, "progM");

            // set M prog cell division rate config and propogate
            for(int il = 0; il <= temp; il++){
                if(il == number_max_layer_P){
                    pDeath_P = 0.0;
                    pDivSym_P = 0.0;
                    pDivAsymSingle_P = 0.0;
                    pDivAsymDouble_P = 0.0;
                    pDivSymDiff_P = 1.0;  
                    pDivAsymDiff_P = 0.0;                
                }
                sim.SetProgMDiv(pDeath_P, pDiv_P, pDivSym_P, pDivAsym_P, pDivAsymSingle_P, pDivAsymDouble_P, 
                pDivSymDiff_P, pDivAsymDiff_P, sltAdv_P, mfacDeathRate_P, il);                    
                sim.ProgDivToNext(s, il, "progM");
                                           
            }

            for(int il = 0; il <= temp; il++){
                sim.ProgNumUpdateFromPrev(il, il, "progM", "progM"); 
                if(il < number_max_layer_P){
                    sim.ProgNumUpdateFromPrev(il+1, il, "progM", "progM");
                }
            } 

            if(isCancer()){
                isFirstCancerCell = true;
            }    

        } 
        //std::cout<<" globel stem iter: "<<iterStemGen<<" globel prog count: "<<count<<" prog cell num "<<sim.getTotalProgCellNumber("all")<<std::endl;
        // propogate prog NonM cells at global generation for 1 generation
        for(int il = 0; il <= temp; il++){  // for this generation all layer cells
            if (il < number_max_layer_P){ // il is the layer 
                pDeath_P = 0.0;
                pDivSym_P = 0.0;
                pDivAsymSingle_P = 0.0;
                pDivAsymDouble_P = 0.0;
                pDivSymDiff_P = 1.0;
                pDivAsymDiff_P = 0.0;            
            }else{
                //std::cout<<"all layer prog cell number: "<<sim.getTotalProgCellNumber("all")<<std::endl;
                pDeath_P = 0.0;
                pDivSym_P = 0.0;
                pDivAsymSingle_P = 0.0;
                pDivAsymDouble_P = 0.0;
                pDivSymDiff_P = 1.0;  
                pDivAsymDiff_P = 0.0;
            }
            // set global division rate config
            sim.SetProgDiv(pDeath_P, pDiv_P, pDivSym_P, pDivAsym_P, pDivAsymSingle_P, pDivAsymDouble_P, pDivSymDiff_P, pDivAsymDiff_P,
                            sltAdv_P, mfacDeathRate_P, il);
            sim.ProgDivToNext(s, il, "progNonM");                                    
        }

        for(int il = 0; il <= temp; il++){
            sim.ProgNumUpdateFromPrev(il, il, "progNonM", "progNonM"); 
            sim.ProgNumUpdateFromPrev(il, il, "progM", "progNonM"); 
            if(il < number_max_layer_P){
                sim.ProgNumUpdateFromPrev(il+1, il, "progNonM", "progNonM");
                sim.ProgNumUpdateFromPrev(il+1, il, "progM", "progNonM");  
            }
        }

        if(isCancer()){
            isFirstCancerCell = true;
        }    
        //std::cout<<" globel stem iter: "<<iterStemGen<<" globel prog count: "<<count<<" prog cell num "<<sim.getTotalProgCellNumber("all")<<std::endl;
        count += 1;
        //std::cout<<" globel stem iter: "<<iterStemGen<<" globel prog count: "<<count<<" prog cell num "<<sim.getTotalProgCellNumber("all")<<std::endl;
    }
    sim.clearProgCell("all");        
}    


void Loop::loopProgCellLifeTime_simple(){
// assert(num_prog_div_left > 0); // -----
    // update division params
    pDeath_P = 0.0;
    pDivSym_P = 0.0;
    pDivAsymSingle_P = 0.0;
    pDivAsymDouble_P = 0.0;
    pDivSymDiff_P = 1.0;  
    pDivAsymDiff_P = 0.0;
    for (int il = 0; il <= number_max_layer_P; il++){
        sim.SetProgDiv(pDeath_P, pDiv_P, pDivSym_P, pDivAsym_P, pDivAsymSingle_P, pDivAsymDouble_P, pDivSymDiff_P, pDivAsymDiff_P,
            sltAdv_P, mfacDeathRate_P, il);
    }
    int count = 0;    
    while (sim.getTotalProgCellNumber("all") > 0){
        //std::cout<<" globel stem iter: "<<iterStemGen<<" globel prog count: "<<count<<" prog cell num "<<sim.getTotalProgCellNumber("all")<<std::endl;
        // --- while condition -- count <= number_max_layer_P
        // --- while condiion -- sim.getTotalProgCellNumber("all") > 0
        int temp = number_max_layer_P;
        if(temp > count){
            temp = count;
        }
        if(count > 0){
            for(int il = 0; il <= temp; il++){
                sim.ProgNumUpdateFromPrev(il, il, "all", "all");  
                if(il > 0){
                    sim.ProgNumUpdateFromPrev(il, il-1, "all", "all");  
                }              
            }
        }

        if(isCancer()){
                isFirstCancerCell = true;
        }
        
        //std::cout<<"number_max_layer_P: "<<number_max_layer_P<<std::endl;
        //std::cout<<"count: "<<count<<std::endl;
        std::cout<<"0 layer prog cell number: "<<sim.getTotalProgCellNumber_byLayer(0,"all")<<std::endl;
        std::cout<<"1 layer prog cell number: "<<sim.getTotalProgCellNumber_byLayer(1, "all")<<std::endl;
        //std::cout<<"2 layer prog cell number: "<<sim.getTotalProgCellNumber_byLayer(2, "all")<<std::endl;
        //std::cout<<"3 layer prog cell number: "<<sim.getTotalProgCellNumber_byLayer(3, "all")<<std::endl;
        //std::cout<<"4 layer prog cell number: "<<sim.getTotalProgCellNumber_byLayer(4, "all")<<std::endl;
        //std::cout<<"all layer prog cell number: "<<sim.getTotalProgCellNumber("all")<<std::endl;
        //std::cout<<"terminal cell number: "<<sim.getTotalTermCellNumber()<<std::endl;
        

        
        for(int il = 0; il <= temp; il++){            
            if (il < number_max_layer_P){
                pDeath_P = 0.0;
                pDivSym_P = 0.0;
                pDivAsymSingle_P = 0.0;
                pDivAsymDouble_P = 0.0;
                pDivSymDiff_P = 1.0;
                pDivAsymDiff_P = 0.0;            
            }else{
                //std::cout<<"all layer prog cell number: "<<sim.getTotalProgCellNumber("all")<<std::endl;
                pDeath_P = 0.0;
                pDivSym_P = 0.0;
                pDivAsymSingle_P = 0.0;
                pDivAsymDouble_P = 0.0;
                pDivSymDiff_P = 1.0;  
                pDivAsymDiff_P = 0.0;
            }
            sim.SetProgDiv(pDeath_P, pDiv_P, pDivSym_P, pDivAsym_P, pDivAsymSingle_P, pDivAsymDouble_P, pDivSymDiff_P, pDivAsymDiff_P,
                            sltAdv_P, mfacDeathRate_P, il);
            sim.ProgDivToNext(s, il, "all");
            //sim.updateProgCellNumberFromProg(s, il);
            //std::cout<<"all layer prog cell number: "<<sim.getTotalProgCellNumber("all")<<std::endl;                    
        }
        
        //std::cout<<"terminal cell number: "<<sim.getTotalTermCellNumber()<<std::endl;
        count += 1;
        //std::cout<<" globel stem iter: "<<iterStemGen<<" globel prog count: "<<count<<" prog cell num "<<sim.getTotalProgCellNumber("all")<<std::endl;
    }
    sim.clearProgCell("all");
    //std::cout<<"prog cell layer count: "<<count<<std::endl;
    //for(int il = 0; il <= number_max_layer_P; il++){
    //    sim.SetProgCellNumber(0.0, il);
    //}
}


void Loop::loopProgCellMLifeTime(){
// update division params
    pDeath_P = 0.0;
    pDivSym_P = 0.0;
    pDivAsymSingle_P = 0.0;
    pDivAsymDouble_P = 0.0; //prob that one il generates one il and one (il + 1)
    pDivSymDiff_P = 1.0; // prob that
    pDivAsymDiff_P = 0.0; 
    for(int il = 0; il <= number_max_layer_P; il++){
        sim.SetProgMDiv(pDeath_P, pDiv_P, pDivSym_P, pDivAsym_P, pDivAsymSingle_P, pDivAsymDouble_P, pDivSymDiff_P, pDivAsymDiff_P,
            sltAdv_P, mfacDeathRate_P, il);
    }
// loop
    int num_prog_div_left = 1 * ceil(1.0 / mfacTurnoverTime_P);
    int count = 0;
    while (sim.getTotalProgCellNumber("progM") > 0){

        //std::cout<<" globel stem iter: "<<iterStemGen<<" globel prog count: "<<count<<" prog cell num "<<sim.getTotalProgCellNumber("progM")<<std::endl;
        // only operate on current generation prog cells
        int temp = number_max_layer_P;
        if(temp > count){
            temp = count;
        }              

        // propogate at M cell faster sub-generation
        for(int j = 0; j < num_prog_div_left; j++){ // j is the M prog cell sub-generation 
            if(j < num_prog_div_left - 1){
                pDeath_P = 0.0;
                pDivSym_P = 1.0;
                pDivAsymSingle_P = 0.0;
                pDivAsymDouble_P = 0.0;
                pDivSymDiff_P = 0.0;
                pDivAsymDiff_P = 0.0; 
            }else{
                pDeath_P = 0.0;
                pDivSym_P = 0.0;
                pDivAsymSingle_P = 0.0;
                pDivAsymDouble_P = 0.0;
                pDivSymDiff_P = 1.0;
                pDivAsymDiff_P = 0.0; 
            }                                     

            // compute and add in prog M Cell
            double value = sim.getTotalProgCellNumber("progM");
            int startTp = iterStemGen * turnoverTime_S + iterStemMGen * turnoverTime_S * 1.0 / ceil(1.0 / mfacTurnoverTime) + 
                            count * turnoverTime_P + 
                            j * turnoverTime_P * 1.0 / ceil(1.0 / mfacTurnoverTime_P);
            int endTp = startTp + turnoverTime_P * 1.0 / ceil(1.0 / mfacTurnoverTime_P);
            computeCountCellNumber(value, startTp, endTp, "progM");

            // set M prog cell division rate config and propogate
            for(int il = 0; il <= temp; il++){
                if(il == number_max_layer_P){
                    pDeath_P = 0.0;
                    pDivSym_P = 0.0;
                    pDivAsymSingle_P = 0.0;
                    pDivAsymDouble_P = 0.0;
                    pDivSymDiff_P = 1.0;  
                    pDivAsymDiff_P = 0.0;                
                }
                sim.SetProgMDiv(pDeath_P, pDiv_P, pDivSym_P, pDivAsym_P, pDivAsymSingle_P, pDivAsymDouble_P, 
                pDivSymDiff_P, pDivAsymDiff_P, sltAdv_P, mfacDeathRate_P, il);                    
                sim.ProgDivToNext(s, il, "progM");
                                           
            }

            for(int il = 0; il <= temp; il++){
                sim.ProgNumUpdateFromPrev(il, il, "progM", "progM"); 
                if(il < number_max_layer_P){
                    sim.ProgNumUpdateFromPrev(il+1, il, "progM", "progM");
                }
            } 

            if(isCancer()){
                isFirstCancerCell = true;
            }    
        } 

        //std::cout<<" globel stem iter: "<<iterStemGen<<" globel prog count: "<<count<<" prog cell num "<<sim.getTotalProgCellNumber("progM")<<std::endl;
        count += 1;
        //std::cout<<" globel stem iter: "<<iterStemGen<<" globel prog count: "<<count<<" prog cell num "<<sim.getTotalProgCellNumber("progM")<<std::endl;
    }
    sim.clearProgCell("progM");     
}

void Loop::loopOneTissueLifetime(){
    while(iterStemGen < number_div_total){
        //if(iterStemGen>1){break;} // debug

        //std::cout<<"iterStemGen: "<< iterStemGen << " isFirstCancerCell "<<isFirstCancerCell<<std::endl;// debug
        //std::cout<<" seed value for stem-gen "<<iterStemGen<<" is "<<s<<std::endl;        
        //std::cout<<"global stem iter: "<<iterStemGen << " stem M cell num: "<<sim.getTotalStemCellMNumber()<<std::endl;
        //std::cout<<"global stem iter: "<<iterStemGen << " prog cell num: "<<sim.getTotalProgCellNumber("all")<<std::endl;
        updateOneStemGeneration();
        iterStemGen++;
    // ---     if(iterStemGen > 10){break;}
        //if(isFirstCancerCell){
        //    std::cout<<"cancer incidence at stem generation: "<< iterStemGen << std::endl;// debug
        //}
        // s++; // update rng seed
    }
}


void Loop::loopRolls(){
    std::cout<<"params version: "<<version<<std::endl;
    for(int iroll = 0; iroll < nrolls; iroll++){
        //if( iroll != 8 ){continue;} // ++++++
        //if(iroll > 0){break;} // debug
        if( iroll != 98){continue;}
        std::cout<<"roll "<<iroll<<" out of "<<nrolls<<std::endl;
        s0 = (iroll%10000000) * number_div_total * 200 ; // ---+++ 31


        resetiterStemGen();
        resetFirstCancer();
        resetFirstMCell();
        initCell();
        loopOneTissueLifetime();


        if(isFirstCancerCell){
            computecountCancer(genFirstCancerCell); // genFirstCancerCell = in this roll, the generation
            std::cout<<" cancer incidence appeared on roll: "<<iroll<<" at gen: "<<genFirstCancerCell<<std::endl;
            //in which the first cancer cell appeared.
        }
        //std::cout<<"genFirstCancerCell "<<genFirstCancerCell<<std::endl; // debug

        if(isFirstMCell){
            computecountMCell(genFirstMCell);
        }


        // debug 2/14/16
        if(genFirstCancerCell < 2){
            dbug_iroll = iroll;
        }

    }
    t = clock() - t;
    sim.memFree();
}


bool Loop::isCancer(){
    bool iscancer = 0;
    if(number_max_layer_P >= 0){
        iscancer = ((sim.getCancerStemCellNumber() > 0)||(sim.getCancerProgCellNumber()>0)||(sim.getCancerTermCellNumber()>0));
    }else{
        iscancer = ((sim.getCancerStemCellNumber() > 0));
    }
    return iscancer;
}

bool Loop::isMCell(){
    return (sim.getTotalStemCellMNumber() > 0);
}

void Loop::computecountCancer(int ngen){
    countCancer += 1;
    countCancerGeneration[ngen] += 1; // need to replace this array with time index instead of generation index
    if(ngen < number_div_sym){
        countCancerSym += 1;
    }else{
        countCancerAsym += 1;
    }
}

void Loop::computecountMCell(int ngen){
    countMCellGeneration[ngen] += 1;
}



void Loop::writeResults(){
    // aggregate cell number to stem/prog/term
    for(int i = 0; i < 365 * 80; i++){
        countStemNumber[i] = countStemNonMNumber[i] + countStemMNumber[i];
        countProgNumber[i] = countProgNonMNumber[i] + countProgMNumber[i];
    }

    // print basic cancer risk results
    std::cout<<"cancer risk: "<<double(countCancer)/double(nrolls)<<std::endl;
    std::cout<<"cancer total incidences : "<<countCancer<<std::endl;
    std::cout<<"cancer total incidences in sym div : "<<countCancerSym<<std::endl;
    std::cout<<"cancer total incidences in asym div : "<<countCancerAsym<<std::endl;
    std::cout<<"Runing time in seconds: "<< ((float)t)/CLOCKS_PER_SEC <<std::endl;
    std::cout<<"bug iroll number "<<dbug_iroll<<std::endl;

    std::ofstream myfile;
    // output per day cell number
    myfile.open("../output/CellNumber"+std::to_string(version)+".csv");
    myfile << "StemNonM" <<","<< "StemM"<< "," << "Stem"<< "," << "ProgNonM"<< "," << "ProgM"<< "," << "Prog" << "\n";
    for(int i = 0; i < 365 * 80; i++){
        myfile << countStemNonMNumber[i] <<"," << countStemMNumber[i] << "," << countStemNumber[i] << "," << 
                countProgNonMNumber[i]<< "," << countProgMNumber[i]<< "," << countProgNumber[i];
        myfile << "\n";
    } 

    myfile.open("../output/CancerRiskSummary"+std::to_string(version)+".csv");
    myfile << "version "<<version<<"\n";
    myfile << "tissue_name "<<tissue_name<<"\n";
    myfile << "nrolls " << nrolls << "\n";
    myfile << "countCancer "<< countCancer << "\n";
    myfile << "countCancerSym "<< countCancerSym << "\n";
    myfile << "countCancerAsym " << countCancerAsym << "\n";
    myfile.close();

    myfile.open("../output/CancerNumberGeneration"+std::to_string(version)+".csv");
    for(int l = 0; l < number_div_total; l++){
        myfile << countCancerGeneration[l];
        myfile << "\n";
    }
    myfile.close();


    myfile.open("../output/MCellNumberGeneration"+std::to_string(version)+".csv");
    for(int l = 0; l < number_div_total; l++){
        myfile << countMCellGeneration[l];
        myfile << "\n";
    }
    myfile.close();
}
