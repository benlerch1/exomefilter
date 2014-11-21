// exomefilter v1
//
// Simulate whole-exome sequencing study to identify mendelian disease loci using 1000 genomes reference data
//
// Author: Ben Lerch
// Last Update: Oct 2014

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <float.h>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <time.h>

using namespace std;




int main(int argc, char** argv){

readArguments(argc, argv);

ifstream g (genofile.c_str());
if(g.is_open()){
while( getline(g,line) ){
geno.push_back(line);
}}

ifstream c (chromfile.c_str());
if(c.is_open()){
while( getline(c,line) ){
chrom.push_back(line);
}}

ifstream p (posfile.c_str());
if(p.is_open()){
while( getline(p,line) ){
pos.push_back(boost::lexical_cast<int>(line));
}}

int counter = 1;
for(int i=0;i<numReps;i++){ 
pedToVCF(chrom, pos, geno);
for(int i=0; i<filters.size(); i++){ results[i]=results[i]+string(",")+boost::lexical_cast<string>(passTotal[i]);}
for(int i=0; i<passTotal.size(); i++){passTotal[i]=0;}
rinkisults=rinkisults+string(",")+boost::lexical_cast<string>(rinkipass);
rinkipass = 0;
cout << "Iteration " << counter << " complete." << endl;
counter++;
}

for(int i=0; i<filters.size();i++){ cout << filters[i].substr(1,string::npos) << "<-c(" << results[i].substr(1,string::npos) << ")" << endl; }

cout << "rinki<-c(" << rinkisults.substr(1,string::npos) << ")" << endl;

//string rfile = string("output/") + out + string("/plots.R");

//ofstream myfile;
//myfile.open (rfile.c_str());

//for(int i=0; i<filters.size();i++){ myfile << filters[i].substr(1,string::npos) << "<-c(" << results[i].substr(1,string::npos) << ");"; }

//for(int i=0; i<filters.size();i++){ myfile << "png(filename = 'output/" << out << "/" << filters[i].substr(1,string::npos) << ".png'); max<-max("<< filters[i].substr(1,string::npos)<<"); min<-min("
//    << filters[i].substr(1,string::npos) << ");mybreaks<-max/(max-min)*100; hist(" << filters[i].substr(1,string::npos) << 
//    ",breaks=mybreaks,col=rgb(2/255, 112/255, 200/255, 0.90),xlab='Number of Variants',ylab='Iterations');dev.off();";}


//myfile.close();

}
