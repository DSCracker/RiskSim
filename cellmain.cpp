#include <iostream>
#include <string>
#include <random>
#include <fstream>
#include <time.h>
#include <math.h>

#include "Cell.h"
#include "stemCell.h"
#include "progCell.h"

#include "Simulator.h"
#include "Params.h"
#include "Loop.h"
#include "strdiff.h"

//g++ -O3 -std=c++11 -o run cellmain.cpp Params.cpp Cell.cpp 
//stemCell.cpp progCell.cpp Simulator.cpp Params.cpp Loop.cpp


int main(){
    // ::::::::::::::: formal run :::::::::::::::::: 

    Loop loop;
    loop.loopRolls();
    loop.writeResults();
    
    // ::::::::::::::: end of formal run :::::::::::


    

return 0;
}







