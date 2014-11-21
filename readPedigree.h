void readPedigree(string pedigree){

string line;
ifstream pedigreeFile (pedigree.c_str());

if(pedigreeFile.is_open()){
while( getline(pedigreeFile,line) ){

typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
boost::char_separator<char> sep("\t");
tokenizer tokens(line, sep);
tokenizer::iterator tok_iter = tokens.begin();

FAMID.push_back(*tok_iter); ++tok_iter; ID.push_back(*tok_iter); ++tok_iter; FATID.push_back(*tok_iter); ++tok_iter; MOTID.push_back(*tok_iter); ++tok_iter;
SEX.push_back(boost::lexical_cast<int>(*tok_iter)); ++tok_iter; DISEASE.push_back(boost::lexical_cast<int>(*tok_iter)); ++tok_iter; SEQUENCE.push_back(boost::lexical_cast<int>(*tok_iter));
if(FATID.back()=="." && MOTID.back()=="."){ FOUNDER.push_back(1);}
else{FOUNDER.push_back(0);}

}
}
else{cout << "ERROR: Unable to open pedigree file." << endl;}

for(int i=0; i < FAMID.size(); i++){cout << FAMID[i] << "\t" << ID[i] << "\t" << FATID[i] << "\t" << MOTID[i] << "\t" << SEX[i] << "\t" << DISEASE[i] << "\t" << SEQUENCE[i] << endl;}

}
