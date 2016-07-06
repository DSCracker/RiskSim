
#include <iostream>
#include <assert.h>
#include <memory.h>
#include <string>
std::string strdiff(std::string sbig, std::string ssub){
// return the the chars in s but not in ssub (normally ssub is a subset of set s;)
// the order of chars in string for this function does not matter
    if(ssub == "0" && sbig != "0"){
        return sbig;
    }
        
    assert(ssub.length() < sbig.length());
    //std::string ch(sbig.length() - ssub.length(),'0');
    char *ch = new char [sbig.length() - ssub.length() + 1];
    memset(ch, '0', sbig.length() - ssub.length() + 1);
    int counter = 0;
    for(int i = 0; i < sbig.length(); i++){
        if(ssub.find(sbig[i]) == std::string::npos){
        // not found
            ch[counter] = sbig[i];
            ++counter;
        }
    }
    std::string strdiff = ch;
    delete [] ch;
    return strdiff;
}
