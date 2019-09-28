
TH2F * N = 0;
TGraph * gref = 0;
TMultiGraph * direct = 0;
TMultiGraph * reflected = 0;


void setupTree(TTree * t) 
{
  t->SetBranchAddress("hN",&N); 
  t->SetBranchAddress("ref",&gref); 
  t->SetBranchAddress("gdirect",&direct); 
  t->SetBranchAddress("greflect",&reflected); 
  t->GetEntry(0); 
}
