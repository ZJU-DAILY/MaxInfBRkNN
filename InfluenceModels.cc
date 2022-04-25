
#include "InfluenceModels.h"


InfluenceModels::InfluenceModels() {
	cout<<"InfluenceModels()!"<<endl;

}

InfluenceModels::~InfluenceModels() {


	if(mc!=NULL) delete mc;
	if(opt!=NULL) delete opt;
	if(imm!=NULL) {
		delete imm;
	}
	else if(opp!=NULL) delete opp;

	/*string pid = intToStr(unsigned(getpid()));
	string outfile = "temp/tmp_" + pid + ".txt";
	
	string command = string("rm -f ") + outfile ;
	system(command.c_str());*/
	
	//delete input;
}


// -c ../config_test.txt
AnyOption* InfluenceModels::readOptions(int argc, char* argv[] ) {
	// read the command line options
	//cout<<"new AnyOption"<<endl;
	AnyOption *opt = new AnyOption();
	//cout<<"after new AnyOption"<<endl;

	// ignore POSIX style options
	opt->noPOSIX();
	//cout<<"opt->noPOSIX"<<endl;
	opt->setVerbose(); /* print warnings about unknown options */
	//cout<<"opt->setVerbose"<<endl;
	opt->autoUsagePrint(true); /* print usage for bad options */
	//cout<<"opt->autoUsagePrint"<<endl;

	/* 3. SET THE USAGE/HELP   */
	//cout<<"SET THE USAGE/HELP"<<endl;
	opt->addUsage( "" );
	opt->addUsage( "Usage: " );
	opt->addUsage( "" );
	opt->addUsage( " -help                            Prints this help " );
	opt->addUsage( " -c <config_file>                 Specify config file " );
	opt->addUsage( " -outdir <output_dir>             Output directory ");
	opt->addUsage( " -propModel <propagation model>   IC or LT ");
	opt->addUsage( " -mcruns <MC runs>                Number of Monte Carlo runs ");
	opt->addUsage( " -budget <budget>                 Size of Seed Set ");
	opt->addUsage( " -celfPlus <0 or 1>               0 for CELF and 1 for CELF++ ");
	opt->addUsage( " -probGraphFile <0 or 1>          Input file for influence network ");
//	opt->addUsage( " -seedFileName <Seed Set File Name if th>   Number of Monte Carlo runs ");
	opt->addUsage( "" );

	/* 4. SET THE OPTION STRINGS/CHARACTERS */
	/* by default all  options  will be checked on the command line and from option/resource file */
	//cout<<"SET THE OPTION STRINGS"<<endl;

	opt->setOption("propModel");
	opt->setOption("outdir");
	opt->setOption("seedFileName");
	//cout<<"seedFileName"<<endl;


	opt->setOption("mcruns");
    opt->setOption("budget");
    opt->setOption("celfPlus");
    opt->setOption("startIt");

	/* for options that will be checked only on the command line not in option/resource file */

	//cout<<"for options that will be checked only on the command line not in option/resource file "<<endl;
	opt->setCommandFlag("help");
	opt->setCommandOption("c");

	/* for options that will be checked only from the option/resource file */
	//cout<<"for options that will be checked only from the option/resource file"<<endl;
	opt->setOption("phase");
	opt->setOption("probGraphFile");
	
	// for static or dynamic friendship graph setting

	// option for dataset size
	// 0 : movielens (very small)
	// 1 : yahoo movies (intermediate)
	// 2 : extremely large (flickr)

	//cout<<"opt->processCommandArgs"<<endl;

	/* go through the command line and get the options  */
	opt->processCommandArgs( argc, argv );

	/* 6. GET THE VALUES */
	if(opt->getFlag( "help" )) {
		opt->printUsage();
		delete opt;
		exit(0);
	}

	const char* configFile = opt->getValue("c");
	if (configFile == NULL) {
		cout << "Config file not mentioned" << endl;
		opt->printUsage();
		delete opt;
		exit(0);
	}

	//cout << "Config file : " << configFile << endl;

	opt->processFile(configFile);
	opt->processCommandArgs( argc, argv );

	phase = strToInt(opt->getValue("phase"));


	//training_dir = opt->getValue("training_dir");
	outdir = opt->getValue("outdir");

	//	cout << "doGenuineLeaders: " << doGenuineLeaders << endl;
	//	cout << "doTribeLeaders: " << doTribeLeaders << endl;
	//cout << "training directory: " << training_dir << endl;
	cout << "output directory: " << outdir << endl;
	cout << endl << endl;
	return opt;

}


int InfluenceModels:: getCurrentSize (int alg){
	if(alg==CELF)
		return mc->getAMSize();
	if(alg==IMM)
		return imm->social_graph->in_neighbors.size();
	if(alg==OPC)
		return opp->social_graph->in_neighbors.size();
}

void InfluenceModels:: seperateThenPMC (vector<int>& stores, BatchResults& results, int b, int R){
#ifdef SEED
	vector<int> seeds = pmc -> runIM(b,R);
	for (int i = 0; i < b; i++) {
		cout << i << "-th u =\t" << seeds[i] << endl;
	}

#else
	cout<<"R="<<R<<endl;
	vector<int> pois =   pmc -> minPOIbyPMC(stores,results,b,R);
	for (int i = 0; i < b; i++) {
		cout << i << "-th poi =\t" << pois[i] << endl;
	}

	//cout<<"R="<<R<<endl;
#endif

}


void InfluenceModels:: batchThenCelf (vector<int>& stores, BatchResults& results, int b){
	mc -> minePOI_CLEFPlus(stores, results, b);
	//mc -> mineSeedSetPlus();
	//mc -> MCSimulation4Seeds(10000);
}

void InfluenceModels:: batchThenIMM (vector<int>& stores, BatchResults& results, int b){
	imm ->checkReachable(stores, results);


	imm -> minePOIsIMM();

	//imm -> MCSimulation4Seeds(10000);
}



int  InfluenceModels:: batchThenOPIMC (vector<int>& stores, BatchResults& results, int b, ofstream& exp_resultOutput){

	opp ->checkReachable(stores, results);
#ifdef SEED
	opp -> InfluenceMaximize_by_opimc(b);
#else
	double runtime = 0; double opimc_time; double hop_time;
	clock_t start_time; clock_t end_time;
	start_time = clock();
	int iter_num =0; int rrset_num = 0; vector<int> poi_seeds;
	opp ->minePOIsOPIMC(b,iter_num,rrset_num,poi_seeds);
	end_time = clock();
	runtime = (double)(end_time- start_time);
	opimc_time += runtime;
	double selection_time = (double)(end_time- start_time)/ CLOCKS_PER_SEC * 1000;
	cout << "*minePOIsOPIMC runtime is: " << selection_time << "ms" << endl;
	exp_resultOutput<< "minePOIsOPIMC runtime is: " << selection_time << "ms, ";
	exp_resultOutput<< "iteration last " << iter_num  << "rounds, rr set number:" << rrset_num << endl;
	exp_resultOutput<<"seleted poi seed set: ";
	for(int poi_selected: poi_seeds){
		exp_resultOutput<<"p"<<poi_selected<<",";
	}
	//统计poi seed对应的潜在用户个数
	set<int> potential_users;
	for(int poi_selected: poi_seeds){
		exp_resultOutput<<"p"<<poi_selected<<",";
		for(ResultDetail rr: results[poi_selected]){
			int _u_id = rr.usr_id;
			potential_users.insert(_u_id);
		}
	}
	exp_resultOutput<<"influence set size="<<potential_users.size()<<endl;
	exp_resultOutput<<endl;


#endif

	return  rrset_num;
}


int InfluenceModels:: batchThenOPIMC_Hybrid(vector<int>& stores, BatchResults& results, int b,ofstream& exp_resultOutput){

	opp ->checkReachable(stores, results);
#ifdef SEED
	opp -> InfluenceMaximize_by_opimc(b);
#else
	double runtime = 0; double opimc_time; double hop_time;
	int iter_num =0; int rrset_num = 0; vector<int> poi_seeds;
	clock_t start_time; clock_t end_time;
	start_time = clock();
	opp ->minePOIsbyHybrid(b, iter_num,rrset_num,poi_seeds);
	end_time = clock();
	runtime = (double)(end_time- start_time);
	cout << "minePOIsbyHybrid runtime is: " << runtime / CLOCKS_PER_SEC * 1000 << "ms" << endl;
	exp_resultOutput<< "*minePOIsbyHybrid runtime is: " << runtime / CLOCKS_PER_SEC * 1000 << "ms" << endl;
	exp_resultOutput<< "iteration last " << iter_num  << "rounds, rr set number:" << rrset_num << endl;
	exp_resultOutput<<"seleted poi seed set: ";
	set<int> potential_users;
	for(int poi_selected: poi_seeds){
		exp_resultOutput<<"p"<<poi_selected<<",";
		for(ResultDetail rr: results[poi_selected]){
			int _u_id = rr.usr_id;
			potential_users.insert(_u_id);
		}
	}
	exp_resultOutput<<"influence set size="<<potential_users.size()<<endl;

#endif
	return rrset_num;

}

void InfluenceModels:: batchThenHopbased (vector<int>& stores, BatchResults& results, int b,ofstream& exp_resultOutput){

	opp ->checkReachable(stores, results);
#ifdef SEED
	opp -> InfluenceMaximize_by_opimc(b);
#else
	double runtime = 0; double opimc_time; double hop_time;
	clock_t start_time; clock_t end_time;

	start_time = clock();
	opp ->minePOIsHopBased(b);
	end_time = clock();
	runtime = (double)(end_time- start_time);
	hop_time += runtime;
	cout << "minePOIsHopBased runtime is: " << runtime / CLOCKS_PER_SEC * 1000 << "ms" << endl;
	exp_resultOutput<< "minePOIsHopBased runtime is: " << runtime / CLOCKS_PER_SEC * 1000 << "ms" << endl;
#endif

	//opp -> MCSimulation4Seeds(10000);
}




void InfluenceModels:: test_batchThenOPIMC (vector<int>& stores, BatchResults& results, int b){

	opp ->checkReachable(stores, results);
#ifdef SEED
	opp -> InfluenceMaximize_by_opimc(b);
#else
	double runtime = 0; double opimc_time; double hop_time;
	clock_t start_time; clock_t end_time;
	int iter_num =0; int rrset_num = 0; vector<int> poi_seeds;
	for(int i=0;i<5;i++){
		start_time = clock();
		opp ->minePOIsOPIMC(b, iter_num, rrset_num, poi_seeds);
		end_time = clock();
		runtime = (double)(end_time- start_time);
		opimc_time += runtime;
		cout << "minePOIsOPIMC runtime is: " << runtime / CLOCKS_PER_SEC * 1000 << "ms" << endl;
		//getchar();
		start_time = clock();
		opp ->minePOIsHopBased(b);
		end_time = clock();
		runtime = (double)(end_time- start_time);
		hop_time += runtime;
		cout << "minePOIsHopBased runtime is: " << runtime / CLOCKS_PER_SEC * 1000 << "ms" << endl;
	}
	cout << "minePOIsOPIMC average runtime is: " << (opimc_time /5.0)/ CLOCKS_PER_SEC * 1000 << "ms" << endl;
	cout << "minePOIsHopBased average runtime is: " << (hop_time /5.0)/ CLOCKS_PER_SEC * 1000 << "ms" << endl;
#endif

	//opp -> MCSimulation4Seeds(10000);
}







