//
// Created by msh on 22-4-22.
//

#ifndef BOOTSTRAPPING_MYUTILS_H
#define BOOTSTRAPPING_MYUTILS_H

#include <NTL/ZZX.h>
#include <helib/helib.h>

void mean_var_mag(const NTL::ZZX& poly, double& mean_dst, double& var_dst, double& mag_dst, int len= 0){
    int d = NTL::deg(poly);
    if(len == 0)
        len = d;
    double mean = 0, var = 0, mag = 0, tmp;
    for(int i = 0; i < d; i++){
        tmp = to_double(poly[i]);
        mag = std::max(mag, std::abs(tmp));
        mean += tmp;
        var += tmp * tmp;
    }
    mean /= len;
    var /= len;
    var -= mean * mean;
//    mean_dst = std::log2(std::abs(mean)); // may be negative
	mean_dst = mean; // NOTE: need to keep the sign!
    var_dst = std::log2(var);
    mag_dst = std::log2(mag);
}

void ZZX_mod(NTL::ZZX& poly, long q){
    for(auto &c : poly.rep){ // NOTE: nonzero constant poly has deg 0
        c = c % q;
    }

}


#endif //BOOTSTRAPPING_MYUTILS_H
