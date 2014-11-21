void pedToVCF(vector <string> chrom, vector <int> pos, vector <string> geno){

vector <string> inVCF, FAMIDinVCF;
vector <int> colMOT, colFAT, oriMOT, oriFAT, upToMOT, upToFAT, FOUNDERinVCF, SEXinVCF, SEQUENCEinVCF, DISEASEinVCF;


string chr = "1";

//initialize values for sampling from VCF
while(inVCF.size() < FAMID.size()){ 
for(int i = 0; i < FAMID.size(); i++){
if(find(inVCF.begin(), inVCF.end(), ID[i]) == inVCF.end()){
if(FOUNDER[i]==1){ inVCF.push_back(ID[i]); FOUNDERinVCF.push_back(FOUNDER[i]); SEXinVCF.push_back(SEX[i]); SEQUENCEinVCF.push_back(SEQUENCE[i]); DISEASEinVCF.push_back(DISEASE[i]); FAMIDinVCF.push_back(FAMID[i]);
colMOT.push_back(pickRefSample()); colFAT.push_back(pickRefSample()); oriMOT.push_back(bernoulliTrial()); oriFAT.push_back(bernoulliTrial()); upToMOT.push_back(pickCrossoverPoint()); upToFAT.push_back(pickCrossoverPoint());}

if(find(inVCF.begin(), inVCF.end(), FATID[i]) != inVCF.end() && find(inVCF.begin(), inVCF.end(), MOTID[i]) != inVCF.end() ){ 
inVCF.push_back(ID[i]); FOUNDERinVCF.push_back(FOUNDER[i]); SEXinVCF.push_back(SEX[i]); SEQUENCEinVCF.push_back(SEQUENCE[i]); DISEASEinVCF.push_back(DISEASE[i]); FAMIDinVCF.push_back(FAMID[i]);
colMOT.push_back(myfindString(inVCF,MOTID[i])); colFAT.push_back(myfindString(inVCF,FATID[i]));
oriMOT.push_back(bernoulliTrial()); oriFAT.push_back(bernoulliTrial()); upToMOT.push_back(pickCrossoverPoint()); upToFAT.push_back(pickCrossoverPoint()); }
}}
}

cout << "pedigree reordering complete" << endl;


//construct fake VCF line-by-line

int iter;
//start by reading in ref VCF and creating first 9 fields of my VCF
for(int k=0; k<numVariants;k++){
string myVCFgeno;
iter = k;

//genotype each individual and add to my VCF line
for(int i=0;i<inVCF.size();i++){


if(FOUNDERinVCF[i]==1){
if(chr!= chrom[k]){ chr=chrom[k]; oriMOT[i]=bernoulliTrial(); oriFAT[i]=bernoulliTrial(); upToMOT[i]=pickCrossoverPoint(); upToFAT[i]=pickCrossoverPoint();}
if(upToMOT[i] < pos[k]){oriMOT[i] = (oriMOT[i]+1) % 2; upToMOT[i]=upToMOT[i]+pickCrossoverPoint(); }
if(upToFAT[i] < pos[k]){oriFAT[i] = (oriFAT[i]+1) % 2; upToFAT[i]=upToFAT[i]+pickCrossoverPoint(); }
if(chr!="X"){myVCFgeno += geno[k].substr(colMOT[i]*2+oriMOT[i],1); myVCFgeno += geno[k].substr(colFAT[i]*2+oriFAT[i],1);}
else if(SEXinVCF[i]==1){upToFAT[i]=500000000; myVCFgeno += geno[k].substr(colMOT[i]*2+oriMOT[i],1); myVCFgeno += geno[k].substr(colFAT[i]*2,1);}
else if(SEXinVCF[i]==0){myVCFgeno += geno[k].substr(colMOT[i]*2+oriMOT[i],1)+string(".");}
}



if(FOUNDERinVCF[i]==0){
string m = myVCFgeno.substr(colMOT[i]*2+oriMOT[i],1);
string f = myVCFgeno.substr(colFAT[i]*2+oriFAT[i],1);

if(chr!= chrom[k]){ chr=chrom[k]; oriMOT[i]=bernoulliTrial(); oriFAT[i]=bernoulliTrial(); upToMOT[i]=pickCrossoverPoint(); upToFAT[i]=pickCrossoverPoint();}
if(upToMOT[i] < pos[k]){oriMOT[i] = (oriMOT[i]+1) % 2; upToMOT[i]=upToMOT[i]+pickCrossoverPoint(); }
if(upToFAT[i] < pos[k]){oriFAT[i] = (oriFAT[i]+1) % 2; upToFAT[i]=upToFAT[i]+pickCrossoverPoint(); }
if(chr!="X"){myVCFgeno += m; myVCFgeno += f;}
else if(SEXinVCF[i]==1){upToFAT[i]=500000000; myVCFgeno += m; f = myVCFgeno.substr(colFAT[i]*2,1); myVCFgeno += f;}
else if(SEXinVCF[i]==0){myVCFgeno += m; myVCFgeno += ".";}
}

}

//define AC for this variant

int numSeq = 0;

double myMAF;

int myAC = 0;
for(int i=0; i<inVCF.size(); i++){ 
if(SEQUENCEinVCF[i]==1 && (myVCFgeno.substr(i,2) == ".1" || myVCFgeno.substr(i,2) == "1." || myVCFgeno.substr(i,2) == "01" || myVCFgeno.substr(i,2) == "10")){myAC++;}
if(SEQUENCEinVCF[i]==1 && myVCFgeno.substr(i,2) == "11"){myAC++; myAC++;}
}

//remove if monomorphic
if(myAC==0){continue;}


for(int i=0; i < myVCFgeno.size(); i++){if(myVCFgeno.substr(i,1) != "."){ numSeq++;}}


if(myAC*2 <= numSeq){myMAF = (double)myAC/(double)numSeq;}
else{myMAF = 1 - (double)myAC/(double)numSeq;}




applyFilters(myVCFgeno, inVCF, DISEASEinVCF, SEQUENCEinVCF, iter, myMAF); 

} 

}
