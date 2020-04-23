#include <cassert>
#include <iostream>
#include <algorithm>
#include "LexiConDMIOHelpers.h"
#include "LexiConDMHelperFunctions.h"


using namespace std;


namespace LexiConDMIOHelpers{
  LexiConDMIOHelpers::IOBasisType getIOBasisFromString(std::string const& str){
    std::string str_lower;
    LexiConDMHelperFunctions::lowercase(str, str_lower);

    for (unsigned short ibase=0; ibase<static_cast<unsigned int const>(nIOBases); ibase++){
      IOBasisType btype = static_cast<IOBasisType>(ibase);
      std::string strtype = "";
      switch (btype){
      case bAmplitude_JHUGen:
        strtype = "amp_jhu";
        break;
      case bEFT_JHUGen:
        strtype = "eft_jhu";
        break;
      case bEFT_HiggsBasis:
        strtype = "eft_hbasis";
        break;
      default:
        cerr << "LexiConDMIOHelpers::getIOBasisFromString: Basis type " << btype << " is not implemented! Please revise the implementation." << endl;
        assert(0);
      }

      if (strtype==str_lower) return btype;
    }

    cerr << "LexiConDMIOHelpers::getIOBasisFromString: Basis type string " << str << " is not recognized." << endl;
    assert(0);
    return nIOBases; // Dummy return
  }

}
