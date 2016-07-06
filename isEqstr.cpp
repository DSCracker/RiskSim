
#include <iostream>
#include <assert.h>
#include "isSubstr.h"
#include <string>
bool isEqstr(std::string s1, std::string s2){
    bool iseq = false;
    if((s1 == "0") || (s2 == "0")){
        if(s1 == s2){
            iseq = true;
        }else{
            iseq = false;
        }
    }else{
        if ((s1.length() == s2.length()) && (isSubstr(s1, s2))) {
            iseq = true;
        }else{
            iseq = false;
        }
    }
    return iseq;
}


