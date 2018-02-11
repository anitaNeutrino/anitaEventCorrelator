#ifndef HC2_HKTREE
#define HC2_HKTREE

#include "TObject.h"
#include "TROOT.h"


//tree for the HC HK data, all that ANITAns care about probably
class HCHKTree{
 public:
  HCHKTree();
  double lat;
  double lon;
  double alt;
  double time;
  //ClassDef(HCHKTree, 0);
};
//ClassImp(HCHKTree)
inline HCHKTree::HCHKTree(){}

#endif
