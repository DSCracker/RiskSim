#ifndef LOOP_H
#define LOOP_H
#include "Simulator.h"
#include <time.h>
#include <string>

class Loop: public Params
{
protected:
    Simulator sim;
    
    int countCancer = 0;
    int countCancerSym = 0;
    int countCancerAsym = 0;
    int* countCancerGeneration = NULL;
    int* countStemNumberGeneration = NULL;
    int* countMCellGeneration = NULL;
    
    double* countStemNonMNumber = NULL; // 365 * 80, per natural day count of stem Non M Number
    double* countStemMNumber = NULL; // 365 * 80, per natural day count of stem M Number
    double* countStemNumber = NULL; // 365 * 80, per natural day count of stem Number
    double* countProgNonMNumber = NULL; // 365 * 80, per natural day count of prog Non M Number
    double* countProgMNumber = NULL; // 365 * 80, per natural day count of prog M Number 
    double* countProgNumber = NULL; // 365 * 80, per natural day count of prog Number 
    double* countTermNumber = NULL; // 365 * 80, per natural day count of term Number
    //double* countTempNumber = NULL; // 365 * 80, per natural day placeholder count

    unsigned int s0 = 0; // initial seed before each lifetime run
    unsigned int s = 0; 
    clock_t t;
    
    bool isFirstCancerCell = false;
    bool isFirstMCell = false;
    int genFirstCancerCell = 100000;
    int genFirstMCell = 100000;
        
    int iterStemGen = 0;
    int iterStemMGen = 0; // temp stemMGen used for cell number tracking in prog inside stemM
    
    int dbug_iroll = -1;
public:
    Loop();
    void initCountCellNumber(); // init per natural day cell numbers
    void computeCountCellNumber(double value, int startTp, int endTp, std::string cellType); 
    // update per natural day cell numbers each generation

    void initcountCancerGeneration();
    void initcountStemNumberGeneration();
    void initcountMCellGeneration();
    void initCell();
    void updateOneStemGeneration();
    
    void loopStemCellMLifeTime(int num_stem_div_left_sym, int num_stem_div_left_asym);
    void loopProgCellLifeTime();
    void loopProgCellMLifeTime();
    
    // +++
    //void evolProgToProg();//
    //void evolProgToTerm();//
    //void evalProgMToProgM();//
    //void evolProgMToTermM();//
    
    
    void loopOneTissueLifetime();
    void loopProgCellLifeTime_simple();
    bool isCancer();
    bool isMCell();
    
    void computecountCancer(int ngen);
    void computecountMCell(int ngen);
    void loopRolls();
    void writeResults();
    
    void resetFirstCancer(){
        isFirstCancerCell = false;
        genFirstCancerCell = 100000;
        }
    void resetFirstMCell(){
        isFirstMCell = false;
        genFirstMCell = 100000;
        }
    
    void resetiterStemGen(){
        iterStemGen = 0;
    }
    
};








#endif
