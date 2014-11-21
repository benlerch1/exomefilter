void readArguments(int argc, char** argv){

if(argc == 1){

cout << "\n\nstudyNoise v1\nBen Lerch\n\n" << endl;

cout << "description: Compute distribution of number of variants remaining after filtering a vcf containing genotypes assigned to provided list of individuals." << endl;
cout << "Genotype assignments for gen 1 of pedigree are randomly drawn from 1000G phase 1, and children are assigned genotypes by drawing from parent genotypes, switching haplotype based on exp(100 Million).\n\n" << endl;

cout << "-rep [Int]: number of times to sample from ref vcf" << endl;
cout << "-numRefIndv [Int]: number of individuals in reference vcf" << endl;
cout << "-numVariants [Int]: number of variants in ref vcf" << endl;
cout << "-genoFile [file]: genotype fields from ref VCF" << endl;
cout << "-chromFile [file]: chrom column from ref VCF" << endl;
cout << "-posFile [file]: pos column from ref VCF" << endl;


cout << "PEDIGREE ARGUMENTS" << endl;
cout << "-unrelatedStudy [Int]: construct vcf of Int unrelated individuals" << endl;
cout << "-pedigree [file]: construct study based on pedigree file\n" << endl;

cout << "FILTERING ARGUMENTS" << endl;

cout << "Frequency\n" << endl;
cout << "-AlleleFrequency [Double] [file]: filter out variants that are not rare in 1000G\n" << endl;
cout << "-novelindbSNP [file]: filter out variants that are present in dbSNP\n" << endl;
cout << "-novelin1000G: filter out variants that are present in 1000 Genomes with 1000 Genomes Allele Count >= Allele Count of founders in VCF\n" << endl;


cout << "-nonSyn [file]: only keep nonsynonymous variants\n\n" << endl; 


cout << "Conservation\n" << endl;
cout << "-minGERP [Double] [file]: filter out variants that are not highly conserved based on GERP score\n" << endl;
cout << "-minphastCons [Double] [file]: filter out variants that are not highly conserved based on phastCons score\n" << endl;
cout << "-minphyloP [Double] [file]: filter out variants that are not highly conserved based on phyloP score\n" << endl;

cout << "Impact\n" << endl;
cout << "-PolyPhen2 [Double] [file]: filter out variants that are not predicted to be damaging based on PolyPhen2 score\n" << endl;
cout << "-SIFT [Double] [file]: filter out variants that are not predicted to be damaging based on SIFT score\n" << endl;
cout << "-CADD [Double] [file]: filter out variants that are not predicted to be damaging based on CADD score\n" << endl;
cout << "-MutationTasterScore [Double] [file]: filter out variants that are not predicted to be damaging based on MutationTaster score\n" << endl;
cout << "-MutationTasterPred [file]: filter out variants that are not predicted to be deleterious by MutationTaster (N or P)\n" << endl;
cout << "-MutationAssessorScore [Double] [file]: filter out variants that are not predicted to be damaging based on MutationAssessor score\n" << endl;
cout << "-MutationAssessorPred [String] [file]: filter out variants that are not predicted to be deleterious by MutationAssessor\n" << endl;
cout << "-LRT [Double] [file]: filter out variants that are not predicted to be damaging based on LRT score\n" << endl;
cout << "-LRTpred [Double] [file]: filter out variants that are not predicted to be damaging based on LRT score\n" << endl;

}

for(int i = 1; i < argc; i++){

if(string(argv[i])=="-rep"){numReps = boost::lexical_cast<int>(argv[i+1]);}
if(string(argv[i])=="numRefIndv"){numRefIndv = boost::lexical_cast<int>(argv[i+1]);}
if(string(argv[i])=="numVariants"){numVariants = boost::lexical_cast<int>(argv[i+1]);}
if(string(argv[i])=="-genoFile"){genofile = argv[i+1];}
if(string(argv[i])=="-chromFile"){chromfile = argv[i+1];}
if(string(argv[i])=="-posFile"){posfile = argv[i+1];}


if(string(argv[i])=="-unrelatedStudy"){ int numunrelatedStudy = boost::lexical_cast<int>(argv[i+1]); 
for(int i=0;i<numunrelatedStudy;i++){ FAMID.push_back(string("FAM")+boost::lexical_cast<string>(i)); ID.push_back(string("ID")+boost::lexical_cast<string>(i)); 
FATID.push_back("."); MOTID.push_back("."); SEX.push_back(bernoulliTrial()); DISEASE.push_back(bernoulliTrial()); SEQUENCE.push_back(1); FOUNDER.push_back(1);} 
} 

if(string(argv[i]) == "-pedigree"){ string pedfile = argv[i+1]; readPedigree(pedfile); }

if(string(argv[i]) == "-out" ){ out = argv[i+1];}

if(string(argv[i]) == "-Nonsynonymous"){filters.push_back(argv[i]);  passTotal.push_back(0); results.push_back(""); ifstream n (Annofile.c_str()); if(n.is_open()){ while( getline(n,line) ){ anno.push_back(line);}}}

if(string(argv[i]) == "-AlleleFrequency"){filters.push_back(argv[i]); maxAF = boost::lexical_cast<double>(argv[i+1]); passTotal.push_back(0); results.push_back("");
ifstream a (AFfile.c_str()); if(a.is_open()){ while( getline(a,line) ){ af.push_back(boost::lexical_cast<double>(line));}}}

if(string(argv[i]) == "-novelin1000G"){filters.push_back(argv[i]); AFfile = argv[i+1]; passTotal.push_back(0); results.push_back("");
ifstream a (AFfile.c_str()); if(a.is_open()){ while( getline(a,line) ){ af.push_back(boost::lexical_cast<double>(line));}} cout << "length of maf vector: " << af.size() << endl;}


if(string(argv[i]) == "-GERP"){filters.push_back(argv[i]); minGERP = boost::lexical_cast<double>(argv[i+1]); passTotal.push_back(0); results.push_back("");
ifstream g (GERPfile.c_str()); if(g.is_open()){ while( getline(g,line) ){ GERP.push_back(boost::lexical_cast<double>(line));}}}

if(string(argv[i]) == "-phastCons"){filters.push_back(argv[i]); minphastCons = boost::lexical_cast<double>(argv[i+1]); passTotal.push_back(0); results.push_back("");
ifstream g (phastConsfile.c_str()); if(g.is_open()){ while( getline(g,line) ){ phastCons.push_back(boost::lexical_cast<double>(line));}}}

if(string(argv[i]) == "-phyloP"){filters.push_back(argv[i]); minphyloP = boost::lexical_cast<double>(argv[i+1]); passTotal.push_back(0); results.push_back("");
ifstream g (phyloPfile.c_str()); if(g.is_open()){ while( getline(g,line) ){ phyloP.push_back(boost::lexical_cast<double>(line));}}}

if(string(argv[i]) == "-PolyPhen2"){filters.push_back(argv[i]); minPolyPhen2 = boost::lexical_cast<double>(argv[i+1]); string PolyPhen2file = argv[i+2]; passTotal.push_back(0); results.push_back("");
ifstream g (PolyPhen2file.c_str()); if(g.is_open()){ while( getline(g,line) ){ PolyPhen2.push_back(boost::lexical_cast<double>(line));}}}

if(string(argv[i]) == "-SIFT"){filters.push_back(argv[i]); maxSIFT = boost::lexical_cast<double>(argv[i+1]); string SIFTfile = argv[i+2]; passTotal.push_back(0); results.push_back("");

ifstream g (SIFTfile.c_str()); if(g.is_open()){ while( getline(g,line) ){ SIFT.push_back(boost::lexical_cast<double>(line));}}}


if(string(argv[i]) == "-CADD"){filters.push_back(argv[i]); minCADD = boost::lexical_cast<int>(argv[i+1]); passTotal.push_back(0); results.push_back("");
ifstream g (CADDfile.c_str()); if(g.is_open()){ while( getline(g,line) ){ CADD.push_back(boost::lexical_cast<int>(line));}}}

if(string(argv[i]) == "-MutationTasterScore"){filters.push_back(argv[i]); minMutationTasterScore = boost::lexical_cast<double>(argv[i+1]); passTotal.push_back(0); results.push_back("");
ifstream g (MutationTasterScorefile.c_str()); if(g.is_open()){ while( getline(g,line) ){ MutationTasterScore.push_back(boost::lexical_cast<double>(line));}}}

if(string(argv[i]) == "-MutationTasterPred"){filters.push_back(argv[i]); string MutationTasterPredfile = argv[i+1]; passTotal.push_back(0); results.push_back("");
ifstream g (MutationTasterPredfile.c_str()); if(g.is_open()){ while( getline(g,line) ){ MutationTasterPred.push_back(line);}}}

if(string(argv[i]) == "-MutationAssessorScore"){filters.push_back(argv[i]); minMutationAssessorScore = boost::lexical_cast<double>(argv[i+1]); passTotal.push_back(0); results.push_back("");
ifstream g (MutationAssessorScorefile.c_str()); if(g.is_open()){ while( getline(g,line) ){ MutationAssessorScore.push_back(boost::lexical_cast<double>(line));}}}

if(string(argv[i]) == "-MutationAssessorPred"){filters.push_back(argv[i]); string MutationAssessorPredfile = argv[i+1]; passTotal.push_back(0); results.push_back("");
ifstream g (MutationAssessorPredfile.c_str()); if(g.is_open()){ while( getline(g,line) ){ MutationAssessorPred.push_back(line);}}}

if(string(argv[i]) == "-LRT"){filters.push_back(argv[i]); minLRT = boost::lexical_cast<int>(argv[i+1]); passTotal.push_back(0); results.push_back("");
ifstream g (LRTfile.c_str()); if(g.is_open()){ while( getline(g,line) ){ LRT.push_back(boost::lexical_cast<int>(line));}}}

if(string(argv[i]) == "-LRTpred"){filters.push_back(argv[i]); string LRTpredfile = argv[i+1]; passTotal.push_back(0); results.push_back("");
ifstream g (LRTpredfile.c_str()); if(g.is_open()){ while( getline(g,line) ){ LRTpred.push_back(line);}}}


if(string(argv[i]) == "-segregatesDom"){filters.push_back(argv[i]); passTotal.push_back(0); results.push_back("");}
if(string(argv[i]) == "-segregatesRec"){filters.push_back(argv[i]); passTotal.push_back(0); results.push_back("");}}

}
