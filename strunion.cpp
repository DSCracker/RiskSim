
#include <iostream>
#include <assert.h>
#include <memory.h>
#include <string>

std::string strunion(std::string s1, std::string s2){
    std::string stru = "";
    int ind = 0;
    std::string strset[5] = {"1","2","3","4","5"};
    
    if(s1 == "0"){
       stru = s2;
    }else if (s2 == "0"){
       stru = s1;
    }else{
       int freq[5] = {0}; // freqs for {1,2,3,4,5}
       for(int i = 0; i < s1.length(); i++){
           ind = s1[i] - '1';
           freq[ind] ++;
       }
       for(int i = 0; i < s2.length(); i++){
           ind = s2[i] - '1';
           freq[ind] ++;
       }
       for(int i = 0; i < 5; i++){
           if(freq[i] == 0){continue;}
           stru += strset[i];
       }
    }
return stru;
}


