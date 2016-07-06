#include <iostream>
#include <assert.h>
#include <string>

bool isSubstr(std::string s, std::string ssub){
    if(ssub == "0"){
        return true;
    }
    bool issubstr = true;
    for(int i = 0; i < ssub.length(); i++){
        if(s.find(ssub[i]) == std::string::npos){
            issubstr = false;
        }
    }
    return issubstr;
}
