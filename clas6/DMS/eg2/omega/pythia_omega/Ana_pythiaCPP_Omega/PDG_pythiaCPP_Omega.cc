#include "PDG_pythiaCPP_Omega.h"

PDG_pythiaCPP_Omega::PDG_pythiaCPP_Omega(){
    
}

string PDG_pythiaCPP_Omega::Get_PDGname(int PDGcode){
    string ret;
    switch (PDGcode) {
        case 1: ret = "down"; break;
        case -1: ret = "anti-down"; break;
        case 2: ret = "up"; break;
        case -2: ret = "anti-up"; break;
        case 3: ret = "strange"; break;
        case -3: ret = "anti-strange"; break;
        case 11: ret = "e-"; break;
        case -11: ret = "e+"; break;
        case 21: ret = "gluon"; break;
        case 22: ret = "gamma"; break;
        case 111: ret = "pi0"; break;
        case 211: ret = "pi+"; break;
        case -211: ret = "pi-"; break;
        case 221: ret = "eta"; break;
        case 113: ret = "rho0"; break;
        case 213: ret = "rho+"; break;
        case -213: ret = "rho-"; break;
        case 130: ret = "KL0"; break;
        case 310: ret = "KS0"; break;
        case 311: ret = "K0"; break;
        case -311: ret = "anti-K0"; break;
        case -313: ret = "anti-K*(892)0"; break;
        case 313: ret = "K*(892)0"; break;
        case 323: ret = "K*(892)+"; break;
        case -323: ret = "K*(892)-"; break;
        case 321: ret = "K+"; break;
        case -321: ret = "K-"; break;
        case 331: ret = "eta'(958)"; break;
        case 333: ret = "phi (1020)"; break;
        case 2101: ret = "(ud)0"; break;
        case 2103: ret = "(ud)1"; break;
        case 2203: ret = "(uu)1"; break;
        case 2112: ret = "n"; break;
        case 2212: ret = "p"; break;
        case 1114: ret = "Delta-"; break;
        case 2114: ret = "Delta0"; break;
        case 2214: ret = "Delta+"; break;
        case 2224: ret = "Delta++"; break;
        case 3114: ret = "Sigma*-"; break;
        case 3122: ret = "Lambda"; break;
        case 3112: ret = "Sigma-"; break;
        case 3212: ret = "Sigma0"; break;
        case 3214: ret = "Sigma*0"; break;
        case 3222: ret = "Sigma+"; break;
        case 3224: ret = "Sigma(1385)"; break;
        default: ret = "unknown"; break;
    }
    return ret;
}

