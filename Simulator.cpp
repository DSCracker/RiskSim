#include <iostream>
#include "Params.h"
#include "Cell.h"
#include "stemCell.h"
#include "progCell.h"
#include "Simulator.h"
#include "strdiff.h"
#include "isSubstr.h"
#include "isEqstr.h"
#include <assert.h>

#include <math.h>       /* sqrt */


Simulator::Simulator(){
    setMap();
    //std::cout<<"nsetState  "<<nsetState<<std::endl;
    //std::cout<<"nsetMut  "<<nsetMut<<std::endl;
    //std::cout<<"setState  "<<setState[31]<<std::endl;
    
    scm = new stemCell [nsetState];
    pcm = new progCell * [number_max_layer_P + 1];
    for (int i = 0; i <= number_max_layer_P; i++){
      pcm[i] = new progCell [nsetState];
    }


    SetCellState();
    SetprogCellLayer();
    SetStem(0.0 ,0.0, 1.0 ,0.0 ,0.0 , 0.0 ,0.0, 0.0,0.0 , 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for(int il = 0; il <= number_max_layer_P; il ++){
      SetProg(0.0 ,0.0, 1.0 ,0.0 ,0.0 , 0.0 ,0.0, 0.0,0.0 , 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, il);
    }
}


void Simulator::SetCellState(){
    for(int i = 0; i < nsetState; i++){
        scm[i].updateState(setState[i]);
    }
    for(int i = 0; i <= number_max_layer_P; i++){
      for(int j = 0; j < nsetState; j++){
        pcm[i][j].updateState(setState[j]);
      }
    }
}

void Simulator::SetprogCellLayer(){
    for(int i = 0; i <= number_max_layer_P; i++){
      for(int j = 0; j < nsetState; j++){
        pcm[i][j].setLayer(i);
      }
    }
}

void Simulator::SetEffectOption(std::string dynamics_effect_option_u){
     dynamics_effect_option = dynamics_effect_option_u;
     getMChar();

}

void Simulator::SetStemCellMNumber(double cellNumber){
    if(dynamics_effect_option == "M1"){
        SetStemCellM1Number(cellNumber);
    }else if(dynamics_effect_option == "M1M3"){
        SetStemCellM3Number(cellNumber);
    }
}

void Simulator::SetProgCellMNumber(double cellNumber, int il){
    if(dynamics_effect_option == "M1"){
        SetProgCellM1Number(cellNumber, il);
    }else if(dynamics_effect_option == "M1M3"){
        SetProgCellM3Number(cellNumber, il);
    }
}

void Simulator::updateStemCellMNumberFromStemM(unsigned int s){
    if(dynamics_effect_option == "M1"){
        updateStemCellM1NumberFromStemM1(s);
    }else if(dynamics_effect_option == "M1M3"){
        updateStemCellM3NumberFromStemM3(s);
    }

}


void Simulator::updateStemCellNumberFromStemNonM(unsigned int s){
    if(dynamics_effect_option == "M1"){
        updateStemCellNumberFromStemNonM1(s);
    }else if(dynamics_effect_option == "M1M3"){
        updateStemCellNumberFromStemNonM3(s);
    }

}


void Simulator::updateProgCellMNumberFromProgM(unsigned int s, int il){
    if(dynamics_effect_option == "M1"){
        updateProgCellM1NumberFromProgM1(s, il);
    }else if(dynamics_effect_option == "M1M3"){
        updateProgCellM3NumberFromProgM3(s, il);
    }
}



void Simulator::updateProgCellNumberFromProgNonM(unsigned int s, int il){
    if(dynamics_effect_option == "M1"){
        updateProgCellNumberFromProgNonM1(s, il);
    }else if(dynamics_effect_option == "M1M3"){
        updateProgCellNumberFromProgNonM3(s, il);
    }
}


void Simulator::updateProgCellMNumberFromStemM(){
    if(dynamics_effect_option == "M1"){
        updateProgCellM1NumberFromStemM1();
    }else if(dynamics_effect_option == "M1M3"){
        updateProgCellM3NumberFromStemM3();
    }
}


void Simulator::updateProgCellNumberFromStemNonM(){
    if(dynamics_effect_option == "M1"){
        updateProgCellNumberFromStemNonM1();
    }else if(dynamics_effect_option == "M1M3"){
        updateProgCellNumberFromStemNonM3();
    }
}





double Simulator::getTotalStemCellMNumber(){
    double total = 0.0;
    if(dynamics_effect_option == "M1"){
        total = getTotalStemCellM1Number();
    }else if(dynamics_effect_option == "M1M3"){
        total = getTotalStemCellM3Number();
    }
    return total;
}

double Simulator::getTotalProgCellMNumber(){
    double total = 0.0;
    if(dynamics_effect_option == "M1"){
        total = getTotalProgCellM1Number();
    }else if(dynamics_effect_option == "M1M3"){
        total = getTotalProgCellM3Number();
    }
    return total;
}


void Simulator::SetStem(double cellNumber , double pDeath , double pDiv , double pM1 , double pM2 , double pM3 ,
   double pM4, double pM5,
  double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff,
  double sltAdv, double mfacDeathRate, double mfacMutRate){
    SetStemCellNumber(cellNumber);
    SetStemDiv(pDeath , pDiv , pDivSym , pDivAsym, pDivAsymSingle, pDivAsymDouble, pDivSymDiff, sltAdv, mfacDeathRate);
    SetStemMut(pM1 , pM2 , pM3, pM4, pM5, mfacMutRate);
}

void Simulator::SetProg(double cellNumber , double pDeath , double pDiv , double pM1 , double pM2 , double pM3 ,
   double pM4, double pM5,
  double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff, double pDivAsymDiff,
  double sltAdv, double mfacDeathRate, double mfacMutRate, int il){
    SetProgCellNumber(cellNumber, il);
    SetProgDiv(pDeath , pDiv , pDivSym , pDivAsym, pDivAsymSingle, pDivAsymDouble, pDivSymDiff, pDivAsymDiff, sltAdv, mfacDeathRate, il);
    SetProgMut(pM1 , pM2 , pM3, pM4, pM5, mfacMutRate, il);
}

void Simulator::SetStemCellNumber(double cellNumber){
     for(int i = 0; i < nsetState; i++){
        if(setState[i] == "0"){
            scm[i].updateCellNumber(cellNumber);
        }else{
            scm[i].updateCellNumber(0);
        }
     }
}

void Simulator::SetProgCellNumber(double cellNumber, int il){
     for(int i = 0; i < nsetState; i++){
        pcm[il][i].updateCellNumber(cellNumber);
     }     
}

//+++
void Simulator::addCellNumberFromStemNonM1ToStemM1(){
    for(int i = 0; i < nsetState; i++){
       if(isSubstr(setState[i], "1")){
           scm[i].addCellNumberFromNonM();
       }    
    }     
}
//+++
void Simulator::addCellNumberFromStemNonM3ToStemM3(){
     for(int i = 0; i < nsetState; i++){
       if(isSubstr(setState[i], "3")){
           scm[i].addCellNumberFromNonM();
       }    
     }    
}
//+++
void Simulator::addCellNumberFromStemNonMToStemM(){
    if(dynamics_effect_option == "M1"){
        addCellNumberFromStemNonM1ToStemM1();
    }else if(dynamics_effect_option == "M1M3"){
        addCellNumberFromStemNonM3ToStemM3();
    }
}
//+++
void Simulator::addCellNumberFromProgNonM1ToProgM1(int il){
     for(int i = 0; i < nsetState; i++){
       if(isSubstr(setState[i], "1")){
          pcm[il][i].addCellNumberFromNonM();
       }    
     } 
}
//+++
void Simulator::addCellNumberFromProgNonM3ToProgM3(int il){
     for(int i = 0; i < nsetState; i++){
       if(isSubstr(setState[i], "3")){
          pcm[il][i].addCellNumberFromNonM();
       }    
     } 
}
//+++
void Simulator::addCellNumberFromProgNonMToProgM(int il){
    if(dynamics_effect_option == "M1"){
        addCellNumberFromProgNonM1ToProgM1(il);
    }else if(dynamics_effect_option == "M1M3"){
        addCellNumberFromProgNonM3ToProgM3(il);
    }
}

void Simulator::SetStemCellM1Number(double cellNumber){
     for(int i = 0; i < nsetState; i++){
        if(isSubstr(setState[i], "1")){
           scm[i].updateCellNumber(cellNumber);
        }       
     } 
}

void Simulator::SetStemCellM3Number(double cellNumber){
     for(int i = 0; i < nsetState; i++){
        if(isSubstr(setState[i], "3")){
           scm[i].updateCellNumber(cellNumber);
        }       
     }
}

void Simulator::SetProgCellM1Number(double cellNumber, int il){
     for(int i = 0; i < nsetState; i++){
        if(isSubstr(setState[i], "1")){
            pcm[il][i].updateCellNumber(cellNumber);        
        }       
     } 
}

void Simulator::SetProgCellM3Number(double cellNumber, int il){
     for(int i = 0; i < nsetState; i++){
        if(isSubstr(setState[i], "3")){
            pcm[il][i].updateCellNumber(cellNumber);          
        }       
     }
     
}

void Simulator::SetStemDiv(double pDeath , double pDiv , double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff,double sltAdv, double mfacDeathRate){
    // the div prob effect is on M1 cells  
    //double pDeathM1 = pDeath * pow((1 - sltAdv), mfacDeathRate);
    double pDeathM1 = pDeath * (1 - sltAdv) *  mfacDeathRate;
    //double pDivM1 = 1 - pDeathM1;
    //double pDivM1 = 1;
    double pDivSymM1 = 1.0 - pDeathM1 - pDivAsymSingle - pDivAsymDouble - pDivSymDiff;
    for(int i = 0; i < nsetState; i++){
       if(isSubstr(setState[i], "1")){
          scm[i].updateProbDiv(pDeathM1, pDiv, pDivSymM1, pDivAsym, pDivAsymSingle, pDivAsymDouble, pDivSymDiff);
        }else{
          scm[i].updateProbDiv(pDeath, pDiv, pDivSym, pDivAsym, pDivAsymSingle, pDivAsymDouble, pDivSymDiff);
       }        
    }  
        
}

void Simulator::SetStemM1Div(double pDeath , double pDiv , double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff,double sltAdv, double mfacDeathRate){
    
    //double pDeathM1 = pDeath * pow((1 - sltAdv), mfacDeathRate);
    double pDeathM1 = pDeath * (1 - sltAdv) *  mfacDeathRate;
    //double pDivM1 = 1 - pDeathM1;
    //double pDivM1 = 1;
    double pDivSymM1 = 1.0 - pDeathM1 - pDivAsymSingle - pDivAsymDouble - pDivSymDiff;
    for(int i = 0; i < nsetState; i++){
       if(isSubstr(setState[i], "1")){
          scm[i].updateProbDiv(pDeathM1, pDiv, pDivSymM1, pDivAsym, pDivAsymSingle, pDivAsymDouble, pDivSymDiff);
        }       
    }
       
}

void Simulator::SetStemM3Div(double pDeath , double pDiv , double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff,double sltAdv, double mfacDeathRate){
    
    //double pDeathM1 = pDeath * pow((1 - sltAdv), mfacDeathRate);
    double pDeathM1 = pDeath * (1 - sltAdv) *  mfacDeathRate;
    //double pDivM1 = 1 - pDeathM1;
    //double pDivM1 = 1;
    double pDivSymM1 = 1.0 - pDeathM1 - pDivAsymSingle - pDivAsymDouble - pDivSymDiff;
    for(int i = 0; i < nsetState; i++){
       if(isSubstr(setState[i], "3")){
          if(isSubstr(setState[i], "1")){       
              scm[i].updateProbDiv(pDeathM1, pDiv, pDivSymM1, pDivAsym, pDivAsymSingle, pDivAsymDouble, pDivSymDiff);
          }else{
              scm[i].updateProbDiv(pDeath, pDiv, pDivSym, pDivAsym, pDivAsymSingle, pDivAsymDouble, pDivSymDiff);
          }
        }       
    }      
}


void Simulator::SetStemMDiv(double pDeath , double pDiv , double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff,
    double sltAdv, double mfacDeathRate){
       if(dynamics_effect_option == "M1"){
        SetStemM1Div(pDeath , pDiv , pDivSym , pDivAsym, pDivAsymSingle, pDivAsymDouble, pDivSymDiff, sltAdv,  mfacDeathRate);
        }else if(dynamics_effect_option == "M1M3"){
        SetStemM3Div(pDeath , pDiv , pDivSym , pDivAsym, pDivAsymSingle, pDivAsymDouble, pDivSymDiff, sltAdv,  mfacDeathRate);
        }
    }

void Simulator::SetProgMDiv(double pDeath , double pDiv , double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff, double pDivAsymDiff,
    double sltAdv, double mfacDeathRate, int il){
       if(dynamics_effect_option == "M1"){
        SetProgM1Div(pDeath , pDiv , pDivSym , pDivAsym, pDivAsymSingle, pDivAsymDouble, pDivSymDiff, pDivAsymDiff, sltAdv,  mfacDeathRate, il);
        }else if(dynamics_effect_option == "M1M3"){
        SetProgM3Div(pDeath , pDiv , pDivSym , pDivAsym, pDivAsymSingle, pDivAsymDouble, pDivSymDiff, pDivAsymDiff, sltAdv,  mfacDeathRate, il);
        }
    }


void Simulator::SetProgDiv(double pDeath , double pDiv , double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff, double pDivAsymDiff,
  double sltAdv, double mfacDeathRate, int il){

       //double pDeathM1 = pDeath * pow((1 - sltAdv), mfacDeathRate);
       double pDeathM1 = pDeath * (1 - sltAdv) * mfacDeathRate;
       //double pDivM1 = 1 - pDeathM1;
       double pDivSymM1 = 1.0 - pDeathM1 - pDivAsymSingle - pDivAsymDouble - pDivSymDiff;
       
       for(int i = 0; i < nsetState; i++){
         if(isSubstr(setState[i], "1")){
           pcm[il][i].updateProbDiv(pDeathM1, pDiv, pDivSymM1, pDivAsym, pDivAsymSingle, pDivAsymDouble, pDivSymDiff, pDivAsymDiff);
       }else{
           pcm[il][i].updateProbDiv(pDeath, pDiv, pDivSym, pDivAsym, pDivAsymSingle, pDivAsymDouble, pDivSymDiff, pDivAsymDiff);
       }        
       }        
}
//:::::::pause::::continued::::///
void Simulator::SetProgM1Div(double pDeath , double pDiv , double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff, double pDivAsymDiff,
  double sltAdv, double mfacDeathRate, int il){

       //double pDeathM1 = pDeath * pow((1 - sltAdv), mfacDeathRate);
       double pDeathM1 = pDeath * (1 - sltAdv) * mfacDeathRate;
       //double pDivM1 = 1 - pDeathM1;
       double pDivSymM1 = 1.0 - pDeathM1 - pDivAsymSingle - pDivAsymDouble - pDivSymDiff;
    
       for(int i = 0; i < nsetState; i++){
         if(isSubstr(setState[i], "1")){
           pcm[il][i].updateProbDiv(pDeathM1, pDiv, pDivSymM1, pDivAsym, pDivAsymSingle, pDivAsymDouble, pDivSymDiff, pDivAsymDiff);
       }
       }
}

void Simulator::SetProgM3Div(double pDeath , double pDiv , double pDivSym , double pDivAsym, double pDivAsymSingle, double pDivAsymDouble, double pDivSymDiff, double pDivAsymDiff,
  double sltAdv, double mfacDeathRate, int il){

       //double pDeathM1 = pDeath * pow((1 - sltAdv), mfacDeathRate);
       double pDeathM1 = pDeath * (1 - sltAdv) * mfacDeathRate;
       //double pDivM1 = 1 - pDeathM1;
       double pDivSymM1 = 1.0 - pDeathM1 - pDivAsymSingle - pDivAsymDouble - pDivSymDiff;
    
       for(int i = 0; i < nsetState; i++){
         if(isSubstr(setState[i], "3")){
           if(isSubstr(setState[i], "1")){
               pcm[il][i].updateProbDiv(pDeathM1, pDiv, pDivSymM1, pDivAsym, pDivAsymSingle, pDivAsymDouble, pDivSymDiff, pDivAsymDiff);
           }else{
               pcm[il][i].updateProbDiv(pDeath, pDiv, pDivSym, pDivAsym, pDivAsymSingle, pDivAsymDouble, pDivSymDiff, pDivAsymDiff);
           }
         }
       }
}

void Simulator::SetStemMut(double pM1 , double pM2 , double pM3, double pM4, double pM5, double mfacMutRate){
    // the mutation rate effect is on M2
    for(int i = 0; i < nsetState; i++){
      if(isSubstr(setState[i], "2")){
        scm[i].updateProbMutation(pM1 * mfacMutRate , pM2 * mfacMutRate , pM3 * mfacMutRate , pM4 * mfacMutRate, pM5 * mfacMutRate);
      }else{
        scm[i].updateProbMutation(pM1 , pM2 , pM3, pM4, pM5);
      }
    }
    
}

void Simulator::SetProgMut(double pM1 , double pM2 , double pM3, double pM4, double pM5, double mfacMutRate, int il){
    for(int i = 0; i < nsetState; i++){
      if(isSubstr(setState[i], "2")){
        pcm[il][i].updateProbMutation(pM1 * mfacMutRate , pM2 * mfacMutRate , pM3 * mfacMutRate , pM4 * mfacMutRate, pM5 * mfacMutRate);
      }else{
        pcm[il][i].updateProbMutation(pM1 , pM2 , pM3, pM4, pM5);
      }
    }
}


void Simulator::getMChar(){
      //std::string MChar = "0"; // DO NOT re-initialize variable if it is a member !!! 
      if(dynamics_effect_option == "M1"){
        MChar = "1";
      }else if(dynamics_effect_option == "M1M3"){
        MChar = "3";
      }
}

void Simulator::StemDivToNext(unsigned int s, std::string thisType){
    for(int i = 0; i < nsetState; i++){    
        if(thisType == "stemM"){
          if(!isSubstr(setState[i], MChar)){continue;}
        }else if(thisType == "stemNonM"){
          if(isSubstr(setState[i], MChar)){continue;}
        }
        scm[i].getCellNumberByBehavior(s);
    }
    for(int i = 0; i < nsetState; i++){
        if(thisType == "stemM"){
          if(!isSubstr(setState[i], MChar)){continue;}
        }else if(thisType == "stemNonM"){
          if(isSubstr(setState[i], MChar)){continue;}
        }  
        scm[i].updateCellNumberNext();
    }    
}

void Simulator::StemNumUpdateFromPrev(std::string thisType, std::string prevType){
    for(int i = 0; i < nsetState; i++){
        if(thisType == "stemM"){
          //std::cout<<" MChar "<<MChar<<std::endl;
          if(!isSubstr(setState[i], MChar)){
            //std::cout<<" ! stemM skiped "<<std::endl;
            continue;}
        }else if(thisType == "stemNonM"){
          if(isSubstr(setState[i], MChar)){continue;}
        }
        double total = 0;
        for(int j = 0; j < nsetState; j++){
            if(prevType == "stemM"){
              if(!isSubstr(setState[j], MChar)){continue;}
            }else if(prevType == "stemNonM"){
              if(isSubstr(setState[j], MChar)){continue;}
            }
            total += scm[j].m_cellNumberNextM[i];
        }
        if(thisType == prevType){
          scm[i].updateCellNumber(total);
        }else{
          scm[i].addCellNumber(total); 
        }        
    }
}



void Simulator::ProgDivToNext(unsigned int s, int thisLayer, std::string thisType){ // +++ the activity that progcell of layer il divide one generation
    assert(thisLayer >= 0 && thisLayer <= number_max_layer_P);
     for(int i = 0; i < nsetState; i++){
        if(thisType == "progM"){
          if(!isSubstr(setState[i], MChar)){continue;}
        }else if(thisType == "progNonM"){
          if(isSubstr(setState[i], MChar)){continue;}
        }
        pcm[thisLayer][i].getCellNumberByBehavior(s);
     }
     for(int i = 0; i < nsetState; i++){
        if(thisType == "progM"){
          if(!isSubstr(setState[i], MChar)){continue;}
        }else if(thisType == "progNonM"){
          if(isSubstr(setState[i], MChar)){continue;}
        }
        pcm[thisLayer][i].updateCellNumberNext();
     }
}

void Simulator::ProgNumUpdateFromPrev(int thisLayer, int prevLayer, std::string thisType, std::string prevType){ // +++ update prog cell number from prev generation's division
    assert(thisLayer >= 0 && thisLayer <= number_max_layer_P);
    assert(prevLayer >= 0 && prevLayer <= number_max_layer_P);

    for(int i = 0; i < nsetState; i++){
        if(thisType == "progM"){
          if(!isSubstr(setState[i], MChar)){continue;}
        }else if(thisType == "progNonM"){
          if(isSubstr(setState[i], MChar)){continue;}
        }
        double total = 0;
        for(int j = 0; j < nsetState; j++){
          if((prevType == "progM") || (prevType == "stemM")){
            //std::cout<<" stemM !! "<<std::endl;
            if(!isSubstr(setState[j], MChar)){
              //std::cout<<" ! stemM skiped "<<std::endl;
              continue;}
          }else if((prevType == "progNonM") || (prevType == "stemNonM")){
            if(isSubstr(setState[j], MChar)){continue;}
          }

          if((prevType == "progM") || (prevType == "progNonM") || (prevType == "all")){
            if(thisLayer == prevLayer){
                total += pcm[prevLayer][j].m_cellNumberNextM[i];
            }else{
                total += pcm[prevLayer][j].m_cellNumberNextM_NextLayer[i];
            }
            
          }else if((prevType == "stemM") || (prevType == "stemNonM")){
            total += scm[j].m_cellNumberNextM_P[i];
          }        
        }

        if( (thisType == prevType) && (thisLayer == prevLayer) ){
          pcm[thisLayer][i].updateCellNumber(total);
        }else if( (thisType == "progM") && (prevType == "stemM") ){
          pcm[thisLayer][i].updateCellNumber(total);
        }else{
          pcm[thisLayer][i].addCellNumber(total);
        }                    
  }
}

    
void Simulator::updateStemCellNumberFromStem(unsigned int s){
    // update stem cell number for all types of stem cells from current generation to the next
     for(int i = 0; i < nsetState; i++){
         scm[i].getCellNumberByBehavior(s);
     }
     for(int i = 0; i < nsetState; i++){
         scm[i].updateCellNumberNext();
     }
     for(int i = 0; i < nsetState; i++){
         double total = 0;
         for(int j = 0; j < nsetState; j++){
             total += scm[j].m_cellNumberNextM[i];
         }
         scm[i].updateCellNumber(total);
     }

}

void Simulator::updateProgCellNumberFromProg(unsigned int s, int il){ // +++ modify layerwise update, easy to make mistakes
     assert(il >= 0 && il <= number_max_layer_P);
     for(int i = 0; i < nsetState; i++){
         pcm[il][i].getCellNumberByBehavior(s);
     }
     for(int i = 0; i < nsetState; i++){
         pcm[il][i].updateCellNumberNext();
     }
     for(int i = 0; i < nsetState; i++){
         double total_sameLayer = 0;
         double total_nextLayer = 0;
         for(int j = 0; j < nsetState; j++){
             total_sameLayer += pcm[il][j].m_cellNumberNextM[i];
             total_nextLayer += pcm[il][j].m_cellNumberNextM_NextLayer[i];
         }
         pcm[il][i].updateCellNumber(total_sameLayer);
         if(il < number_max_layer_P){          
            pcm[il+1][i].addCellNumber(total_nextLayer);
         }
     }
}



void Simulator::updateProgCellNumberFromStem(){
    for(int i = 0; i < nsetState; i++){
         double total = 0;
         for(int j = 0; j < nsetState; j++){
             total += scm[j].m_cellNumberNextM_P[i];
         }
         pcm[0][i].updateCellNumber(total); // layer must be 0
    }
}

// ++++++
//::::: paused ::::: continued :::::://
   
void Simulator::updateStemCellNumberFromStemNonM1(unsigned int s){
    // NonM1 stem cell to NonM1 stem cell
    // and NonM1 stem cell to M1 stem cell
     for(int i = 0; i < nsetState; i++){
         if(isSubstr(setState[i], "1")){continue;}
         scm[i].getCellNumberByBehavior(s);
     }
    
     for(int i = 0; i < nsetState; i++){
         if(isSubstr(setState[i], "1")){continue;}
         scm[i].updateCellNumberNext();
     }

     for(int i = 0; i < nsetState; i++){
         double total = 0;
         for(int j = 0; j < nsetState; j++){
             if(isSubstr(setState[j], "1")){continue;}
             total += scm[j].m_cellNumberNextM[i];
         }
         if(isSubstr(setState[i], "1")){
             scm[i].updateCellNumberFromNonM(total);
         }else{
             scm[i].updateCellNumber(total);
         }   
     }
}




//+++
void Simulator::updateProgCellNumberFromProgNonM1(unsigned int s, int il){
    // NonM1 stem cell to NonM1 stem cell
    // and NonM1 stem cell to M1 stem cell
    assert(il >= 0 && il <= number_max_layer_P);
    for(int i = 0; i < nsetState; i++){
         if(isSubstr(setState[i], "1")){continue;}
         pcm[il][i].getCellNumberByBehavior(s);
     }
    
     for(int i = 0; i < nsetState; i++){
         if(isSubstr(setState[i], "1")){continue;}
         pcm[il][i].updateCellNumberNext();
     }

     for(int i = 0; i < nsetState; i++){
         double total_sameLayer = 0;
         double total_nextLayer = 0;
         for(int j = 0; j < nsetState; j++){
             if(isSubstr(setState[j], "1")){continue;}
             total_sameLayer += pcm[il][j].m_cellNumberNextM[i];
             total_nextLayer += pcm[il][j].m_cellNumberNextM_NextLayer[i];
         }
         if(isSubstr(setState[i], "1")){
             pcm[il][i].updateCellNumberFromNonM(total_sameLayer);
             if(il < number_max_layer_P){
              pcm[il+1][i].updateCellNumberFromNonM(total_nextLayer);
             }             
         }else{
             pcm[il][i].updateCellNumber(total_sameLayer);
             if(il < number_max_layer_P){
              pcm[il+1][i].updateCellNumber(total_nextLayer);
             }             
         }
     }
}






// +++
void Simulator::updateStemCellNumberFromStemNonM3(unsigned int s){

    for(int i = 0; i < nsetState; i++){
         if(isSubstr(setState[i], "3")){continue;}
         scm[i].getCellNumberByBehavior(s);
     }
    
     for(int i = 0; i < nsetState; i++){
         if(isSubstr(setState[i], "3")){continue;}
         scm[i].updateCellNumberNext();
     }

     for(int i = 0; i < nsetState; i++){
         double total = 0;
         for(int j = 0; j < nsetState; j++){
             if(isSubstr(setState[j], "3")){continue;}
             total += scm[j].m_cellNumberNextM[i];
         }
         if(isSubstr(setState[i], "3")){
             scm[i].updateCellNumberFromNonM(total);
         }else{
             scm[i].updateCellNumber(total);
         } 
     }
}




void Simulator::updateProgCellNumberFromProgNonM3(unsigned int s, int il){

    for(int i = 0; i < nsetState; i++){
         if(isSubstr(setState[i], "3")){continue;}
         pcm[il][i].getCellNumberByBehavior(s);
     }
    
     for(int i = 0; i < nsetState; i++){
         if(isSubstr(setState[i], "3")){continue;}
         pcm[il][i].updateCellNumberNext();
     }

     for(int i = 0; i < nsetState; i++){
         double total_sameLayer = 0;
         double total_nextLayer = 0;
         for(int j = 0; j < nsetState; j++){
             if(isSubstr(setState[j], "3")){continue;}
             total_sameLayer += pcm[il][j].m_cellNumberNextM[i];
             total_nextLayer += pcm[il][j].m_cellNumberNextM_NextLayer[i];
         }
         if(isSubstr(setState[i], "3")){
             pcm[il][i].updateCellNumberFromNonM(total_sameLayer);
             if(il < number_max_layer_P){
              pcm[il+1][i].updateCellNumberFromNonM(total_nextLayer);
             }             
         }else{
             pcm[il][i].updateCellNumber(total_sameLayer);
             if(il < number_max_layer_P){
              pcm[il+1][i].updateCellNumber(total_nextLayer);
             }             
         }
     }
}


void Simulator::updateProgCellNumberFromStemNonM1(){

    for(int i = 0; i < nsetState; i++){
         double total = 0;
         for(int j = 0; j < nsetState; j++){
             if(isSubstr(setState[j], "1")){continue;}
             total += scm[j].m_cellNumberNextM_P[i];
         }
         if(isSubstr(setState[i], "1")){
             pcm[0][i].updateCellNumberFromNonM(total);
         }else{
             pcm[0][i].updateCellNumber(total);
         }
     }
}




void Simulator::updateProgCellNumberFromStemNonM3(){
    for(int i = 0; i < nsetState; i++){
         double total = 0;
         for(int j = 0; j < nsetState; j++){
             if(isSubstr(setState[j], "3")){continue;}
             total += scm[j].m_cellNumberNextM_P[i];
         }
         if(isSubstr(setState[i], "3")){
             pcm[0][i].updateCellNumberFromNonM(total);
         }else{
             pcm[0][i].updateCellNumber(total);
         }
     }
}



// :::::::::
void Simulator::updateStemCellM1NumberFromStemM1(unsigned int s){
    // stem M1 to stem M1
    for(int i = 0; i < nsetState; i++){
         if(isSubstr(setState[i], "1")){
             scm[i].getCellNumberByBehavior(s);
         }
     }
    
     for(int i = 0; i < nsetState; i++){
         if(isSubstr(setState[i], "1")){
             scm[i].updateCellNumberNext();
         }
     }

     for(int i = 0; i < nsetState; i++){
         if(!isSubstr(setState[i], "1")){continue;}
         double total = 0;
         for(int j = 0; j < nsetState; j++){
             if(!isSubstr(setState[j], "1")){continue;}
             total += scm[j].m_cellNumberNextM[i];
         }
         scm[i].updateCellNumber(total);
     }
}

void Simulator::updateStemCellM3NumberFromStemM3(unsigned int s){
    // stem M3 to stem M3
    for(int i = 0; i < nsetState; i++){
         if(isSubstr(setState[i], "3")){
             scm[i].getCellNumberByBehavior(s);
         }
     }
    
     for(int i = 0; i < nsetState; i++){
         if(isSubstr(setState[i], "3")){
             scm[i].updateCellNumberNext();
         }
     }

     for(int i = 0; i < nsetState; i++){
         if(!isSubstr(setState[i], "3")){continue;}
         double total = 0;
         for(int j = 0; j < nsetState; j++){
             if(!isSubstr(setState[j], "3")){continue;}
             total += scm[j].m_cellNumberNextM[i];
         }
         scm[i].updateCellNumber(total);
     }
}
    
void Simulator::updateProgCellM1NumberFromProgM1(unsigned int s, int il){
    // prog M1 to prog M1
    for(int i = 0; i < nsetState; i++){
         if(isSubstr(setState[i], "1")){
             pcm[il][i].getCellNumberByBehavior(s);
         }
     }
    
     for(int i = 0; i < nsetState; i++){
         if(isSubstr(setState[i], "1")){
             pcm[il][i].updateCellNumberNext();
         }
     }

     for(int i = 0; i < nsetState; i++){
         if(!isSubstr(setState[i], "1")){continue;}
         double total_sameLayer = 0;
         double total_nextLayer = 0;
         for(int j = 0; j < nsetState; j++){
             if(!isSubstr(setState[j], "1")){continue;}
             total_sameLayer += pcm[il][j].m_cellNumberNextM[i];
             total_nextLayer += pcm[il][j].m_cellNumberNextM_NextLayer[i];
         }
         pcm[il][i].updateCellNumber(total_sameLayer);
         if(il < number_max_layer_P){
          pcm[il+1][i].addCellNumber(total_nextLayer);
         }         
     }
}

void Simulator::updateProgCellM3NumberFromProgM3(unsigned int s, int il){
    // prog M3 to prog M3
    for(int i = 0; i < nsetState; i++){
         if(isSubstr(setState[i], "3")){
             pcm[il][i].getCellNumberByBehavior(s);
         }
     }
    
     for(int i = 0; i < nsetState; i++){
         if(isSubstr(setState[i], "3")){
             pcm[il][i].updateCellNumberNext();
         }
     }

     for(int i = 0; i < nsetState; i++){
         if(!isSubstr(setState[i], "3")){continue;}
         double total_sameLayer = 0;
         double total_nextLayer = 0;
         for(int j = 0; j < nsetState; j++){
             if(!isSubstr(setState[j], "3")){continue;}
             total_sameLayer += pcm[il][j].m_cellNumberNextM[i];
             total_nextLayer += pcm[il][j].m_cellNumberNextM_NextLayer[i];
         }
         pcm[il][i].updateCellNumber(total_sameLayer);
         if(il < number_max_layer_P){
          pcm[il+1][i].addCellNumber(total_nextLayer);
         }         
     }
}
    
void Simulator::updateProgCellM1NumberFromStemM1(){
     // stem M1 to prog M1
     for(int i = 0; i < nsetState; i++){
         if(!isSubstr(setState[i], "1")){continue;}
         double total = 0;
         for(int j = 0; j < nsetState; j++){
             if(!isSubstr(setState[j], "1")){continue;}
             total += scm[j].m_cellNumberNextM_P[i];
         }
         pcm[0][i].updateCellNumber(total);
     }
}

void Simulator::updateProgCellM3NumberFromStemM3(){
     // stem M3 to prog M3
     for(int i = 0; i < nsetState; i++){
         if(!isSubstr(setState[i], "3")){continue;}
         double total = 0;
         for(int j = 0; j < nsetState; j++){
             if(!isSubstr(setState[j], "3")){continue;}
             total += scm[j].m_cellNumberNextM_P[i];
         }
         pcm[0][i].updateCellNumber(total);
     }
}


double Simulator::getTotalStemCellNumber(std::string thisType){
     double total = 0;
     for(int i = 0; i < nsetState; i++){
        if(thisType == "stemM"){
          if(!isSubstr(setState[i], MChar)){continue;}
        }else if(thisType == "stemNonM"){
          if(isSubstr(setState[i], MChar)){continue;}
        }
        total += scm[i].getCellNumber();
     }
     return total;
}

double Simulator::getTotalProgCellNumber(std::string thisType){
    double total = 0;
     for(int i = 0; i < nsetState; i++){
        if(thisType == "progM"){
          if(!isSubstr(setState[i], MChar)){continue;}
        }else if(thisType == "progNonM"){
          if(isSubstr(setState[i], MChar)){continue;}
        }
        for(int il = 0; il <= number_max_layer_P; il++){
          total += pcm[il][i].getCellNumber();
        }         
     }
     return total;
}

double Simulator::getTotalProgCellNumber_byLayer(int il, std::string thisType){
    double total = 0;
     for(int i = 0; i < nsetState; i++){
        if(thisType == "progM"){
          if(!isSubstr(setState[i], MChar)){continue;}
        }else if(thisType == "progNonM"){
          if(isSubstr(setState[i], MChar)){continue;}
        }
        total += pcm[il][i].getCellNumber();         
     }
     return total;
}

double Simulator::getTotalTermCellNumber(){
    double total = 0;
     for(int i = 0; i < nsetState; i++){
         total += pcm[number_max_layer_P][i].getChildrenTermCellNumber();
     }
     return total;

}

double Simulator::getTotalStemCellM1Number(){
    double total = 0;
     for(int i = 0; i < nsetState; i++){
         if(!isSubstr(setState[i], "1")){continue;}
         total += scm[i].getCellNumber();
     }
     return total;
}

double Simulator::getTotalStemCellM3Number(){
    double total = 0;
     for(int i = 0; i < nsetState; i++){
         if(!isSubstr(setState[i], "3")){continue;}
         total += scm[i].getCellNumber();
     }
     return total;
}

double Simulator::getTotalProgCellM1Number(){
    double total = 0;
     for(int i = 0; i < nsetState; i++){
         if(!isSubstr(setState[i], "1")){continue;}
         for(int il = 0; il <= number_max_layer_P; il++){
            total += pcm[il][i].getCellNumber();
          }
     }
     return total;
}
double Simulator::getTotalProgCellM3Number(){
    double total = 0;
     for(int i = 0; i < nsetState; i++){
         if(!isSubstr(setState[i], "3")){continue;}
         for(int il = 0; il <= number_max_layer_P; il++){
            total += pcm[il][i].getCellNumber();
          }
     }
     return total;
}

// ::::: paused here::::: continued

double Simulator::getCancerStemCellNumber(){
    double total = 0;
     for(int i = 0; i < nsetState; i++){
         if(isSubstr(setState[i], "123")){
             total += scm[i].getCellNumber();
         }
     }
     return total;
}

double Simulator::getCancerProgCellNumber(){
    double total = 0;
     for(int i = 0; i < nsetState; i++){
         if(isSubstr(setState[i], "1234")){
          for(int il = 0; il <= number_max_layer_P; il++){
            total += pcm[il][i].getCellNumber();
          }
         }
     }
     return total;
}

double Simulator::getCancerTermCellNumber(){
    double total = 0;
     for(int i = 0; i < nsetState; i++){
         for(int j = 0; j < nsetState; j++){
             if(isSubstr(setState[j], "12345")){ // temp change from 12345 to 1234
                 total += pcm[number_max_layer_P][i].m_cellNumberNextM_T[j];
             }
         }
     }
     return total;
}

void Simulator::clearProgCell(std::string thisType){
  for(int i = 0; i < nsetState; i++){
    if(thisType == "progM"){
      if(!isSubstr(setState[i], MChar)){continue;}
      }else if(thisType == "progNonM"){
      if(isSubstr(setState[i], MChar)){continue;}
    }
    for(int il = 0; il <= number_max_layer_P; il++){
        pcm[il][i].clearCell();
    }
  }
}



void Simulator::cellNumberDebugger(){
    //std::cout<<" index for type M12 "<<setState2ind["12"]<<std::endl;
    std::cout<<" cell number for M23 "<<scm[setState2ind["5"]].getCellNumber()<<std::endl;
}



void Simulator::memFree(){
    std::cout<<"deleting dynamic allocated memory "<<std::endl;
    for(int i = 0; i < nsetState; i++){
      scm[i].memDelete();
    }

    for(int i = 0; i < nsetState; i++){
      for(int il = 0; il <= number_max_layer_P; il++){
        pcm[il][i].memDelete();
      }
    }

    delete[] scm;

    for (int il = 0; il <= number_max_layer_P; ++il){
      delete[] pcm[il];
    }
    delete[] pcm;
}














