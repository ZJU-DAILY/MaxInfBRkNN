#include <iomanip>
#include "MC.h"

namespace _MC {

MC::MC(AnyOption* opt1) {
	opt = opt1;

	//
	setSimulatorType(CELF);

	cout << "In testing phase" << endl;

    outdir = opt->getValue("outdir");	
	probGraphFile = opt->getValue("probGraphFile");
	//eta = strToInt(opt->getValue("eta"));
	graphType = DIRECTED;
    //flag_in = opt->getValue(""); 

	setModel(opt);

	AM = new HashTreeCube(89041);
	revAM = NULL;
	
	countIterations = strToInt(opt->getValue("mcruns"));
//	testingActionsFile = (const char*) opt->getValue("testingActionsFile");	

	cout << "User specified options: " << endl;
	cout << "model : " << m << " or " << model << endl;
	cout << "outdir : " << outdir << endl;
	cout << "probGraphFile : " << probGraphFile << endl;
	cout << "Number of iterations in MC : " << countIterations << endl;

	time (&startTime);
	srand ( startTime );
	time (&stime_mintime);

}

MC::~MC() {
    if(AM != NULL)
        delete AM;
    if (outFile.is_open()) {
        outFile.close();
    }
	cout << "total time taken : " << getTime() << endl;
}


void MC::setModel(AnyOption* opt) {
	m = opt->getValue("propModel");
	
	model = LT;

	if (m.compare("LT") == 0) {
		model = LT;
		
	} else if (m.compare("IC") == 0) {
		model = IC;
	}
}

	;

void MC::addSocialLink(int u1, int u2, double prob) {
	users.insert(u1);
	users.insert(u2);
	FriendsMap* neighbors = AM->find(u1);
	if (neighbors == NULL) {
		neighbors = new FriendsMap();
		neighbors->insert(pair<UID, float>(u2, prob));
		AM->insert(u1, neighbors);
	} else {

		neighbors->insert(pair<UID, float>(u2, prob));

	}
}



void MC::computeCov() {
	// read the seed set from the file
	// for each radius R, generate the table 
	// <S_alg, R, cov^R>
    curSeedSet.clear();
    vector<UID> seedsVec;
    seedsVec.clear();

	//int seeds = -1; 
	int seeds = 0; 

	//	str budget = intToStr(curSeedSet.size());
	string seedFileName = opt->getValue("seedFileName");
	string filename = seedFileName + "_cov.txt" ;
	ofstream outFile (filename.c_str(), ios::out);

	cout << "Seedfilename: " << seedFileName << "; outputfile: " << filename << endl;

	ifstream myfile (seedFileName.c_str(), ios::in);
	string delim = " \t";	
	if (myfile.is_open()) {
		while (! myfile.eof() )	{
			std::string line;
			getline (myfile,line);
			if (line.empty()) continue;
			seeds++;
			//if (seeds == 0) continue; // ignore the first line
			std::string::size_type pos = line.find_first_of(delim);
			int	prevpos = 0;

			UID u = 0;
			// get the user
			pos = line.find_first_of(delim, prevpos);
			if (pos == string::npos) 
				u = strToInt(line.substr(prevpos));
			else
				u = strToInt(line.substr(prevpos, pos-prevpos));

            seedsVec.push_back(u); // push into the vector
            cout << "Seed No. " << seedsVec.size() << ": " << u << endl;
			curSeedSet.insert(u);

		}
	}

	myfile.close();

	float cov = 0;
	if (model == IC) {
		cout << "Computing the cov under IC model" << endl;
		cov = ICCov(curSeedSet);
	} else if (model == LT) {
		cout << "Computing the cov under LT model" << endl;

        UserList S;
        S.clear();
        for (vector<UID>::iterator i = seedsVec.begin(); i != seedsVec.end(); ++i) {
            UID v = *i;
            S.insert(v);

            if (S.size() >= 51)
                break;
            // 
            if (S.size() == 1 || S.size() % 5 == 0) {
                float cov1 = LTCov(S);
                outFile << S.size() << " " << v << " " << cov1 << endl;
                cout << S.size() << " " << v << " " << cov1 << endl;
            }

        }

        //outFile <<
        //cout <<
		//cout << "(CovNew, Cov): " << covNew << ", " << cov << endl;
	}

	//cout << "Coverage achieved: " << cov << endl;

	//outFile << cov << endl;

	outFile.close();

}

/**************** NEW ************************************/

 float MC::mineSeedSet(int t_ub){
        cout << "In mineSeedSet" << endl;
        clear();
        //int budget = 50;
        int budget = strToInt(opt->getValue("budget"));

        multimap<float, UID> covQueue; // needed to implement CELF
        totalCov = 0;

        int countUsers = 0 ;

        // first pass, without CELF
        for (UserList::iterator i = users.begin(); i!=users.end(); ++i) {
            UID v = *i;
            UserList S;
            S.insert(v);

            float newCov = 0;

            if (model == IC) {
                //newCov = ICCov(S,t_ub);
                newCov = ICCov(S);
            } else if (model == LT) {
                newCov = LTCov(S);
            }

            covQueue.insert(pair<float, UID>(newCov, v));

            countUsers++;

            if (countUsers % 1000 == 0) {
                cout << "Number of users done in this iteration: " << countUsers << endl;
                cout << "Time, cov, v : " << getTime() << ", " << totalCov << ", " << *curSeedSet.begin() << endl;
            }

            if (totalCov < newCov) {
                totalCov = newCov;
                curSeedSet.clear();
                curSeedSet.insert(v);
            }
        }

        cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^ time: " << getTime() << endl;

        if (curSeedSet.size() != 1) {
            cout << "Some problem" << endl ;
        }

        float actualCov = totalCov;

        float prevActualCov = actualCov;

        writeInFile(*curSeedSet.begin(), totalCov, totalCov, t_ub, actualCov, actualCov, countUsers);
        cout << "Number of nodes examined in this iteration: " << countUsers << endl;
//	cout << "True Cov " << trueCov << endl;


        //cout << "debug stop" << endl;
        //exit(0);

        // remove the last element from covQueue
        // we have already picked the best node
        multimap<float, UID>::iterator i = covQueue.end();
        i--;
        covQueue.erase(i);

        // CELF
        countUsers = 0;
        UserList usersExamined;

//	while (totalCov < eta) {
        while (curSeedSet.size() < budget) {
            bool flag = false;

            multimap<float, UID>::iterator i = covQueue.end();
            i--;

            UID v = i->second;

            curSeedSet.insert(v);
            float newCov = 0;

            if (usersExamined.find(v) == usersExamined.end()) {
                if (model == IC) {
                    newCov = ICCov(curSeedSet);
                } else if (model == LT) {
                    newCov = LTCov(curSeedSet);
                }

                usersExamined.insert(v);
                countUsers++;

            } else {
                newCov = i->first + totalCov;
                flag = true;
            }


            i--;
            float oldCov_nextElement = i->first;
            UID nextElement = i->second;

            if (flag || newCov - totalCov >= oldCov_nextElement) {
                // pick v in the seed set


                writeInFile(v, newCov, newCov - totalCov, t_ub, actualCov, actualCov - prevActualCov, countUsers);
                cout << "Number of nodes examined in this iteration: " << usersExamined.size() << endl;
//			cout << "True Cov " << trueCov << endl;

                prevActualCov = actualCov;

                totalCov = newCov;
                // reset parameters
                countUsers = 1;
                usersExamined.clear();

                i++;
                covQueue.erase(i);
                continue;
            }

            // move the element to another place
            i++;
            covQueue.erase(i);
            covQueue.insert(pair<float, UID>(newCov - totalCov, v));

            // its not in seed set
            curSeedSet.erase(v);
        }

        float trueCov = totalCov;

        cout << "True Cov " << trueCov << endl;

        return totalCov;

}


 float MC::mineSeedSetPlus() {   //选seed的函数，Algorithm 1:CELF++
        cout << "In mineSeedSetPlus" <<endl;
        // cout << "Total # users in this dataset = " << users.size() << endl;

        clear();
//	int budget = 50;
        int budget = strToInt(opt->getValue("budget"));  //kegai

        totalCov = 0;
        Gains gainTable; // the five-column table we maintain (multimap<float, MGStruct*>)

        MGStruct* pBestMG = NULL; // points to the current best node, initially nothing
        UID lastSeedID = -1;
        int countUsers = 0;

        // Do 1st pass without CELF.  Here only the top-1 will be added to Seed Set
        for (UserList::iterator i = users.begin(); i!=users.end(); ++i) {
            curSeedSet.clear(); // seed set is empty for now

            MGStruct* pMG = new MGStruct;
            pMG->nodeID = *i; // set "1st column"

            if (model == IC) {
                // IC model with CELF++
                bool isBest; // compute MG(u|empty)
                isBest = ICCovPlus(pMG, pBestMG);

                if (isBest == true) {
                    pBestMG = pMG;
                } //endif

            }  else if (model == LT) {
                // LT model with CELF++
                bool isBest = false;
                isBest = LTCovPlus(pMG, pBestMG, curSeedSet);

                if (isBest == true) {
                    pBestMG = pMG;
                }//endif
            } //endif

            pMG->flag = 1; // set "5th column" -- only one seed, so flag is only 1
            gainTable.insert(pair<float, MGStruct*>(pMG->gain, pMG));  //将节点影响力加入 表
            countUsers++;

            if (countUsers % 50 == 0) {
                cout << "Number of users done in this iteration: " << countUsers << endl;
                // cout << "Time, cov, v : " << getTime() << ", " << pBestMG->gain << ", " << pBestMG->nodeID << endl;
                Gains::iterator iter = gainTable.end();
                iter--;
                MGStruct* pTop = (MGStruct *) (iter->second);
                if (pTop != NULL) {
                    cout << "Time, cov, v :" << getTime() << ", " <<  pTop->gain << ", " << pTop->nodeID << endl;
                }
            }

        }

        // pick the 1st seed
        Gains::iterator iter = gainTable.end();
        iter--;
        MGStruct* pTop = (MGStruct *) (iter->second);
        if (pTop == NULL) {
            cout << "Some problem: pTop is NULL pointer" << endl;
        } else {
            lastSeedID = pTop->nodeID; // record the last selected seed

            curSeedSet.insert(pTop->nodeID); // the top-1 becomes a seed
            totalCov += pTop->gain; // update the current seeds' coverage

            writeInFile(pTop->nodeID, totalCov, totalCov, 0, 0, 0, countUsers);

            gainTable.erase(iter); // remove that seed from queue（从列表中删去已选点）
        }//第一个节点选择完毕


        // Optimized CELF, from 2nd pass onwards
        int totalUserNum = users.size();
        pBestMG = NULL;
        covBestNode.clear();
        MGStruct * pCurr = NULL;
        MGStruct * pNext = NULL;
        countUsers = 0; // # of users need recomputation!
        //剩余k-1个seed的选择
        while (curSeedSet.size() < budget && totalCov < totalUserNum) {

            Gains::iterator iter = gainTable.end(); //当前float最大的出列
            iter--;
            pCurr = (MGStruct *) (iter->second); // pointing to the current node

            if (pCurr == NULL) {
                cout << "Some problem: pCurr is NULL pointer" << endl;
            }

            // cout << "Scanning node: " << pCurr->nodeID << " with MG = " << pCurr->gain << " ..... And MG_next = " << pCurr->gain_next << "....." << endl;

            // if it has been visited again, pick it as seed   //当前真正实际最优node（line 7)
            if ( pCurr->flag == (curSeedSet.size() + 1) ) {	 // Alg 1 line 7-9
                //cout << "Node " << pCurr->nodeID << " seen twice....." << endl;

                lastSeedID = pCurr->nodeID;
                curSeedSet.insert(pCurr->nodeID);
                totalCov += pCurr->gain;

                writeInFile(pCurr->nodeID, totalCov, pCurr->gain, 0, 0, 0, countUsers);
                // cout << "Seeds picked at No.1 place" << endl;

                gainTable.erase(iter);
                pBestMG = NULL;
                covBestNode.clear();
                countUsers = 0;

                continue;

            } else {  // it has NOT been vistied yet, in this iteration of seeking new seed
                pCurr->flag = curSeedSet.size() + 1;  //之前没在S情况下算marginal,现在要S纳入
                // countUsers++;

                iter--;
                pNext = (MGStruct *) (iter->second); // pointing to the next

                if (pCurr->v_best == lastSeedID && pCurr->gain_next >= pNext->gain) { //line 10,之前预计算有效，且依然为当前最佳 （CELF++优化的地方）
                    // pick it as seed

                    lastSeedID = pCurr->nodeID;

                    curSeedSet.insert(pCurr->nodeID); // update the seed set
                    totalCov += pCurr->gain_next; // update current seeds' coverage

                    writeInFile(pCurr->nodeID, totalCov, pCurr->gain_next, 0, 0, 0, countUsers);
                    //cout << "Seeds picked at No.2 place" << endl;

                    iter++;
                    gainTable.erase(iter); // delete new seed from queue

                    pBestMG = NULL;
                    covBestNode.clear();
                    countUsers = 0;

                    continue;

                } else {
                    // recompute MG, or, just resort if 3rd column is the last picked seed

                    // IF v_best is still the last seed, just move the previouly 4th column value to the 2nd column
                    // To avoid re-computation
                    if (pCurr->v_best == lastSeedID) {  //line 10

                        pCurr->gain = pCurr->gain_next;
                        pCurr->v_best = 0;
                        pCurr->gain_next = 0;
                        pCurr->flag = curSeedSet.size() + 1;


                        if (pBestMG == NULL) {
                            pBestMG = pCurr;
                        } else {
                            if (pCurr->gain > pBestMG->gain) {
                                pBestMG = pCurr;
                            }
                        }

                        iter++;
                        gainTable.erase(iter);
                        gainTable.insert(pair<float, MGStruct *>(pCurr->gain, pCurr));

                    } else {
                        bool isBest; // compute Cov(S+u) and Cov(S+u+prevBest)

                        if (model == IC) {  //利用MC计算 S+u 的 influence
                            isBest = ICCovPlus(pCurr, pBestMG); //moni

                        } else if (model == LT) {
                            isBest = LTCovPlus(pCurr, pBestMG, curSeedSet);
                        } // endif

                        pCurr->flag = curSeedSet.size() + 1;
                        countUsers++;

                        if(pCurr->gain >= pNext->gain) { // pick it as seed
                            lastSeedID = pCurr->nodeID; // record it to be the last added seed

                            curSeedSet.insert(pCurr->nodeID); // add it to seed set
                            totalCov += pCurr->gain; // update current seeds' coverage

                            writeInFile(pCurr->nodeID, totalCov, pCurr->gain, 0, 0, 0, countUsers);
                            // cout << "Seeds picked at No.4 place" << endl;
                            iter++;
                            gainTable.erase(iter); // remove new seed from queue

                            pBestMG = NULL;
                            covBestNode.clear();
                            countUsers = 0;

                            continue;
                        } else {
                            // cannot be seed for now, re-insert it into table (at another position)
                            if (isBest == true) {
                                pBestMG = pCurr;   //当前迭代发现的最优node更新
                            }
                            iter++;
                            gainTable.erase(iter); //擦除旧条目，更新列表
                            gainTable.insert(pair<float, MGStruct *>(pCurr->gain, pCurr)); //加入新条目 lin14

                            continue;
                        } //endif

                    }//endif
                } //endif

            }//endif

        } //endwhile

        return totalCov;

    }



float MC:: minePOI_CLEFPlus(vector<int>& stores, BatchResults& results, int b) {   //选seed的函数，Algorithm 1:CELF++
        cout << "In minePOI_CLEFPlus" <<endl;

        clear();
        //openOutputFiles();
        openExpResutlsFiles();

        int budget = strToInt(opt->getValue("budget"));  //kegai

        totalCov = 0;
        Gains gainTable; // the five-column table we maintain (multimap<float, MGStruct*>)

        MGStruct* pBestMG = NULL; // points to the current best node, initially nothing
        UID lastSeedID = -1;
        int countPOIs = 0;

        // Do 1st pass without CELF.  Here only the top-1 will be added to Seed Set
        for (int store: stores) {
            curSeedSet.clear(); // seed set is empty for now
            if(results[store].size()>0){
                countPOIs++;
                MGStruct* pMG = new MGStruct;
                pMG->nodeID = store; // set "1st column"

                if (model == IC) {
                    // IC model with CELF++
                    bool isBest; // compute MG(u|empty)
                    isBest = ICPOICovPlus(pMG, results, pBestMG);

                    if (isBest == true) {
                        pBestMG = pMG;
                    } //endif

                }

                pMG->flag = 1; // set "5th column" -- only one seed, so flag is only 1
                gainTable.insert(pair<float, MGStruct*>(pMG->gain, pMG));  //将节点影响力加入 表


                if (countPOIs % 5 == 0) {
                    cout << "Number of poi done in this iteration: " << countPOIs << endl;
                    // cout << "Time, cov, store : " << getTime() << ", " << pBestMG->gain << ", " << pBestMG->nodeID << endl;
                    Gains::iterator iter = gainTable.end();
                    iter--;
                    MGStruct* pTop = (MGStruct *) (iter->second);
                    if (pTop != NULL) {
                        cout << "Time, cov, v :" << getTime() << ", " <<  pTop->gain << ", " << pTop->nodeID << endl;
                    }
                }

            }
            else continue;

        }

        // pick the 1st seed
        Gains::iterator iter = gainTable.end();
        iter--;
        MGStruct* pTop = (MGStruct *) (iter->second);
        if (pTop == NULL) {
            cout << "Some problem: pTop is NULL pointer" << endl;
        } else {
            lastSeedID = pTop->nodeID; // record the last selected seed

            curPOISet.insert(pTop->nodeID); // the top-1 becomes a seed
            for(ResultDetail u: results[pTop->nodeID]){
                curSeedSet.insert(u.usr_id);
            }


            totalCov += pTop->gain; // update the current seeds' coverage

            writePOIChooseToFile(pTop->nodeID, totalCov, totalCov, 0, 0, 0, countPOIs);

            gainTable.erase(iter); // remove that seed from queue（从列表中删去已选点）
        }//第一个POI选择完毕


        // Using CELF++, from 2nd pass onwards
        int totalUserNum = users.size();
        pBestMG = NULL;
        covBestNode.clear();
        MGStruct * pCurr = NULL;
        MGStruct * pNext = NULL;
        countPOIs = 0; // # of users need recomputation!
        //进行剩余k-1个POI的选择
        while (curPOISet.size() < budget && totalCov < totalUserNum) {

            Gains::iterator iter = gainTable.end(); //当前float最大的出列
            iter--;
            pCurr = (MGStruct *) (iter->second); // pointing to the current node

            if (pCurr == NULL) {
                cout << "Some problem: pCurr is NULL pointer" << endl;
            }

            // cout << "Scanning node: " << pCurr->nodeID << " with MG = " << pCurr->gain << " ..... And MG_next = " << pCurr->gain_next << "....." << endl;

            // if it has been visited again, pick it as seed   //当前真正实际最优node（line 7)
            if ( pCurr->flag == (curPOISet.size() + 1) ) {	 // Alg 1 line 7-9
                //cout << "Node " << pCurr->nodeID << " seen twice....." << endl;

                lastSeedID = pCurr->nodeID;
                curPOISet.insert(pCurr->nodeID);
                for(ResultDetail u: results[pCurr->nodeID]){
                    curSeedSet.insert(u.usr_id);
                }
                totalCov += pCurr->gain;

                writePOIChooseToFile(pCurr->nodeID, totalCov, pCurr->gain, 0, 0, 0, countPOIs);
                // cout << "Seeds picked at No.1 place" << endl;

                gainTable.erase(iter);
                pBestMG = NULL;
                covBestNode.clear();
                countPOIs = 0;

                continue;

            } else {  // it has NOT been vistied yet, in this iteration of seeking new seed
                pCurr->flag = curPOISet.size() + 1;  //之前没在S情况下算marginal,现在要S纳入
                // countUsers++;

                iter--;
                pNext = (MGStruct *) (iter->second); // pointing to the next

                if (pCurr->v_best == lastSeedID && pCurr->gain_next >= pNext->gain) { //line 10,之前预计算有效，且依然为当前最佳 （CELF++优化的地方）
                    // pick it as a candidate poi

                    lastSeedID = pCurr->nodeID;

                    // update the poi set
                    curPOISet.insert(pCurr->nodeID);
                    for(ResultDetail u: results[pCurr->nodeID]){
                        curSeedSet.insert(u.usr_id);
                    }



                    totalCov += pCurr->gain_next; // update current seeds' coverage

                    writePOIChooseToFile(pCurr->nodeID, totalCov, pCurr->gain_next, 0, 0, 0, countPOIs);
                    //cout << "Seeds picked at No.2 place" << endl;

                    iter++;
                    gainTable.erase(iter); // delete new seed from queue

                    pBestMG = NULL;
                    covBestNode.clear();
                    countPOIs = 0;

                    continue;

                } else {
                    // recompute MG, or, just resort if 3rd column is the last picked seed

                    // IF v_best is still the last seed, just move the previouly 4th column value to the 2nd column
                    // To avoid re-computation
                    if (pCurr->v_best == lastSeedID) {  //line 10

                        pCurr->gain = pCurr->gain_next;
                        pCurr->v_best = 0;
                        pCurr->gain_next = 0;
                        pCurr->flag = curPOISet.size() + 1;


                        if (pBestMG == NULL) {
                            pBestMG = pCurr;
                        } else {
                            if (pCurr->gain > pBestMG->gain) {
                                pBestMG = pCurr;
                            }
                        }

                        iter++;
                        gainTable.erase(iter);
                        gainTable.insert(pair<float, MGStruct *>(pCurr->gain, pCurr));

                    } else {
                        bool isBest; // compute Cov(S+u) and Cov(S+u+prevBest)

                        if (model == IC) {  //利用MC计算 S+u 的 influence
                            isBest = ICPOICovPlus(pCurr, results, pBestMG); //moni

                        } else if (model == LT) {
                            //isBest = LTCovPlus(pCurr, pBestMG, curSeedSet);
                        } // endif

                        pCurr->flag = curPOISet.size() + 1;
                        countPOIs++;

                        if(pCurr->gain >= pNext->gain) { // pick it as seed
                            lastSeedID = pCurr->nodeID; // record it to be the last added seed

                            // add it to seed set
                            curPOISet.insert(pCurr->nodeID);
                            for(ResultDetail u: results[pCurr->nodeID]){
                                curSeedSet.insert(u.usr_id);
                            }


                            totalCov += pCurr->gain; // update current seeds' coverage

                            writePOIChooseToFile(pCurr->nodeID, totalCov, pCurr->gain, 0, 0, 0, countPOIs);
                            // cout << "Seeds picked at No.4 place" << endl;
                            iter++;
                            gainTable.erase(iter); // remove new seed from queue

                            pBestMG = NULL;
                            covBestNode.clear();
                            countPOIs = 0;

                            continue;
                        } else {
                            // cannot be seed for now, re-insert it into table (at another position)
                            if (isBest == true) {
                                pBestMG = pCurr;   //当前迭代发现的最优node更新
                            }
                            iter++;
                            gainTable.erase(iter); //擦除旧条目，更新列表
                            gainTable.insert(pair<float, MGStruct *>(pCurr->gain, pCurr)); //加入新条目 lin14

                            continue;
                        } //endif

                    }//endif
                } //endif

            }//endif

        } //endwhile
        cout<<"selected poi set:";
        printPOISelections();

        return totalCov;

    }


/************ New for www2011 poster **********/
/* pMG: the MG structue of node curV.  curV's ID is included */ //IC模型下的MC模拟
bool MC::ICCovPlus(MGStruct* pMG, MGStruct* pBestMG) {

	startIt = strToInt(opt->getValue("startIt"));
	//cout << "startIt : " << startIt << endl;
	
	if(pMG == NULL) {
		cout << "Some problem 3" << endl;
	}
	
	bool ret = false; // indicates if the current node is better than the previous best
	
	float cov = 0; // compute cov(S + currentNode)
	float covX = 0; // compute cov(S + currentNode + prevBestNode)
	
	int seedSetSize = (int) curSeedSet.size();
	
	// for each MC run:
	for (int b = 0; b < countIterations; b++) {
		queue<UID> Q;
		UserList activeNodes;
		
		// initialize Q and activeNodes using the Seed Set //xg
		for(UserList::iterator i = curSeedSet.begin(); i != curSeedSet.end(); i++) {
			UID v = *i;
			Q.push(v);
			activeNodes.insert(v);
		}
		Q.push(pMG->nodeID); // also put the node being examined into queue //xg
		activeNodes.insert(pMG->nodeID);   //source node 一定是激活的
		
		while (Q.empty() == false) {
			UID v = Q.front();
			
			FriendsMap* neighbors = AM->find(v);
			if (neighbors != NULL) {
				for (FriendsMap::iterator j = neighbors->begin(); j != neighbors->end(); ++j) {
					UID u = j->first;
				
					if (activeNodes.find(u) == activeNodes.end()) {  //邻居节点尚未激活，则激活之
						float toss = ((float)(rand() % 1001))/(float)1000;
						float p = j->second;

						if (p >= toss) {
							activeNodes.insert(u);  //此次激活成功，jihuo
							Q.push(u);
						}
					}			
				}
			}
			
			Q.pop();  //队首节点出列，失去活性
		}//一次MC模拟结束
		
		cov += (float)activeNodes.size()/countIterations;		
		//这是什么？
		// if x has not been activated, and this iteration allows computing 4th column
		if(pBestMG != NULL && seedSetSize >= (startIt - 1)  && activeNodes.find(pBestMG->nodeID) == activeNodes.end()) {
	
			Q.push(pBestMG->nodeID); // Q starts with x only
			activeNodes.insert(pBestMG->nodeID);			

			while (Q.empty() == false) {
				UID v = Q.front();

				FriendsMap* neighbors = AM->find(v);
				if (neighbors != NULL) {

					// for each inactive neighbour of v
					for (FriendsMap::iterator j = neighbors->begin(); j != neighbors->end(); ++j) {
						UID u = j->first;

						if (activeNodes.find(u) == activeNodes.end()) {
							float toss = ((float)(rand() % 1001))/(float)1000;
							float p = j->second;

							if (p >= toss) {
								activeNodes.insert(u);
								Q.push(u);
							}
						}			
					}
				}

				Q.pop();
			}
		
		} 
		
		covX += (float)activeNodes.size()/countIterations;

	} //endfor -- each MC run
	
	// After entire MC process
//	cov = cov / countIterations;
//	covX = covX / countIterations;
	
	pMG->gain = cov - totalCov; // set "2nd column", MG(u|S)...
	
	if (pBestMG == NULL) {
		// for first run of each iteration, Best Node is not specified
		pMG->v_best = pMG->nodeID;
		pMG->gain_next = 0;
	
		if (seedSetSize >= (startIt - 1)) {
			return true;
		} else {  
			return false;  // when it's not allowed to compute 4th column
		}
		
	} else { // if Best Node has been specified

		if (pMG->gain > pBestMG->gain) { 
			// this node is better than the previous best
			pMG->v_best = pMG->nodeID; // set "3rd column" to be itself
			pMG->gain_next = 0; // set "4th column" to be zero
			ret = true;
		} else { 
			// this node is NO better than the previous best
			pMG->v_best = pBestMG->nodeID;
			pMG->gain_next = covX - totalCov - pBestMG->gain; // MG(u|S+x)
			ret = false;
		}
	}
	
	return ret;
}

bool MC::ICPOICovPlus(MGStruct* pMG,  BatchResults& results, MGStruct* pBestMG) {

        startIt = strToInt(opt->getValue("startIt"));
        //cout << "startIt : " << startIt << endl;

        if(pMG == NULL) {
            cout << "Some problem 3" << endl;
        }

        bool ret = false; // indicates if the current node is better than the previous best

        float cov = 0; // compute cov(S + currentNode)
        float covX = 0; // compute cov(S + currentNode + prevBestNode)

        int seedSetSize = (int) curSeedSet.size();

        // for each MC run:
        for (int b = 0; b < countIterations; b++) {
            queue<UID> Q;
            UserList activeNodes;

            // initialize Q and activeNodes using the Seed Set //xg
            for(UserList::iterator i = curSeedSet.begin(); i != curSeedSet.end(); i++) {
                UID v = *i;
                Q.push(v);
                activeNodes.insert(v);
            }
            //shuangse
            for(ResultDetail u: results[pMG->nodeID]){
                Q.push(u.usr_id); // also put the node being examined into queue //xg
                activeNodes.insert(u.usr_id);   //perspective users are source node 一定是激活的
            }

            while (Q.empty() == false) {
                UID v = Q.front();

                FriendsMap* neighbors = AM->find(v);
                if (neighbors != NULL) {
                    for (FriendsMap::iterator j = neighbors->begin(); j != neighbors->end(); ++j) {
                        UID u = j->first;

                        if (activeNodes.find(u) == activeNodes.end()) {  //邻居节点尚未激活，则激活之
                            float toss = ((float)(rand() % 1001))/(float)1000;
                            float p = j->second;

                            if (p >= toss) {
                                activeNodes.insert(u);  //此次激活成功，jihuo
                                Q.push(u);
                            }
                        }
                    }
                }

                Q.pop();  //队首节点出列，失去活性
            }//一次MC模拟结束

            cov += (float)activeNodes.size()/countIterations;
            //这是什么？
            // if x has not been activated, and this iteration allows computing 4th column
            if(pBestMG != NULL && seedSetSize >= (startIt - 1)  && activeNodes.find(pBestMG->nodeID) == activeNodes.end()) {

                Q.push(pBestMG->nodeID); // Q starts with x only
                activeNodes.insert(pBestMG->nodeID);

                while (Q.empty() == false) {
                    UID v = Q.front();

                    FriendsMap* neighbors = AM->find(v);
                    if (neighbors != NULL) {

                        // for each inactive neighbour of v
                        for (FriendsMap::iterator j = neighbors->begin(); j != neighbors->end(); ++j) {
                            UID u = j->first;

                            if (activeNodes.find(u) == activeNodes.end()) {
                                float toss = ((float)(rand() % 1001))/(float)1000;
                                float p = j->second;

                                if (p >= toss) {
                                    activeNodes.insert(u);
                                    Q.push(u);
                                }
                            }
                        }
                    }

                    Q.pop();
                }

            }

            covX += (float)activeNodes.size()/countIterations;

        } //endfor -- each MC run

        // After entire MC process
//	cov = cov / countIterations;
//	covX = covX / countIterations;

        pMG->gain = cov - totalCov; // set "2nd column", MG(u|S)...

        if (pBestMG == NULL) {
            // for first run of each iteration, Best Node is not specified
            pMG->v_best = pMG->nodeID;
            pMG->gain_next = 0;

            if (seedSetSize >= (startIt - 1)) {
                return true;
            } else {
                return false;  // when it's not allowed to compute 4th column
            }

        } else { // if Best Node has been specified

            if (pMG->gain > pBestMG->gain) {
                // this node is better than the previous best
                pMG->v_best = pMG->nodeID; // set "3rd column" to be itself
                pMG->gain_next = 0; // set "4th column" to be zero
                ret = true;
            } else {
                // this node is NO better than the previous best
                pMG->v_best = pBestMG->nodeID;
                pMG->gain_next = covX - totalCov - pBestMG->gain; // MG(u|S+x)
                ret = false;
            }
        }

        return ret;
    }

double MC:: MCSimulation4Seeds(int r) {

        cout<<"MCSimulation4Seeds, for "<<r<<"round..."<<endl;
        double cov = 0; // compute cov(S + currentNode)

        // for each MC run:
        for (int b = 1; b < r+1; b++) {
            queue<UID> Q;
            UserList activeNodes;

            // initialize Q and activeNodes using the Seed Set //xg
            for(UserList::iterator i = curSeedSet.begin(); i != curSeedSet.end(); i++) {
                UID v = *i;
                Q.push(v);
                activeNodes.insert(v);
            }


            while (Q.empty() == false) {
                UID v = Q.front();

                FriendsMap* neighbors = AM->find(v);
                if (neighbors != NULL) {
                    for (FriendsMap::iterator j = neighbors->begin(); j != neighbors->end(); ++j) {
                        UID u = j->first;

                        if (activeNodes.find(u) == activeNodes.end()) {  //邻居节点尚未激活，则激活之
                            float toss = ((float)(rand() % 1001))/(float)1000;
                            float p = j->second;

                            if (p >= toss) {
                                activeNodes.insert(u);  //此次激活成功，jihuo
                                Q.push(u);
                            }
                        }
                    }
                }

                Q.pop();  //队首节点出列，失去活性
            }//一次MC模拟结束

            cov += activeNodes.size();
            if(b%2000==0)
                cout<<"round"<<b<<",activated="<<activeNodes.size()<<"current average ="<< (cov / b) << endl;
                //getchar();

        }
        cov = cov /r;
        cout<<"the number of attracted users = "<< curSeedSet.size()<<" 具体为："<<endl;
        printRealSeeds();
        cout<<"their total social influence = "<<cov<<endl;

        return cov;
    }




float MC::ICCov(UserList& S) {
	float cov = 0;
				
	for (int b = 0; b < countIterations; ++b) {
		// Q is the queue in the depth/breadth first search
		queue<UID> Q;
		// activeNodes is the set of nodes that are activated in the current
		// run of Monte Carlo simulation
		UserList activeNodes;

		// S is the seed set
		// for each seed node v is S, 
		// add it to activeNodes
		// add it to Q as well
		for (UserList::iterator i=S.begin(); i!=S.end(); ++i) {
			UID v = *i;
			Q.push(v);
			activeNodes.insert(v);
		}

		while(Q.empty() == false) {
			UID v = Q.front(); 
			
			// AM is adjacency matrix
			FriendsMap* neighbors = AM->find(v);
			if (neighbors != NULL) {
				for (FriendsMap::iterator j = neighbors->begin(); j!=neighbors->end(); ++j) {
					UID u = j->first;

					if (activeNodes.find(u) == activeNodes.end()) {
						float toss = ((float)(rand() % 1001))/(float)1000;
						float p = j->second;

						if (p >= toss) {
							activeNodes.insert(u);
							Q.push(u);
						}
					}

				}
			}
			
			Q.pop(); 

		}

		cov += (float)activeNodes.size()/countIterations;

	}

	// compute two things: cov(S) and cov(S+x)
//	cov = cov/countIterations;
	return cov;
	

}

void MC::clear() {
	covBestNode.clear();
	curSeedSet.clear();
	seedSetNeighbors.clear();
	totalCov = 0;
}



float MC::LTCov(UserList& S) {
//    cout << "In LTCov" << endl;
    float tol = 0.00001;
    float cov = 0;

//    map<UID, float> ppIn; // what is the prob with which the node is covered
   
    for (int b = 0; b < countIterations; ++b) {
        /* initialize random seed: */

        float cov1 = 0;
        // T is the set of nodes that are to be processed
        queue<UID> T;
        // Q is the set of nodes that have been seen until now
        // Thus, Q is a superset of T
        // Q contains the nodes that are processed and to be processed
        map<UID, NodeParams> Q;

//        cov += S.size();
        // cov1 == coverage in one run .. that is, number of nodes reachable
        cov1 += S.size();
       
        for (UserList::iterator i=S.begin(); i!=S.end(); ++i) {
            UID v = *i;
 //           ppIn[v] += 1;

            FriendsMap* neighbors = AM->find(v);
            if (neighbors != NULL) {
                for (FriendsMap::iterator j = neighbors->begin(); j!=neighbors->end(); ++j) {
                    UID u = j->first;

                    if (S.find(u) == S.end()) {

                        if (Q.find(u) == Q.end()) {
                            // if the node u has not been seen before
                            // create a new NodeParams
                            // create a random threhsold
                            NodeParams& np = Q[u];
                            np.active = false;
                            np.inWeight = j->second;

                            /* generate secret number: */
                            np.threshold = ((float)(rand() % 1001))/(float)1000;
                            T.push(u);
                        } else {
                            NodeParams& np = Q[u];
                            np.inWeight += j->second;                   
                        } //endif
                    } //endif
                } //endfor
            } //endif
        } //endfor

        while (!T.empty()) {
            UID u = T.front();

//            cout << "T.size " << T.size() << endl;

            NodeParams& np = Q.find(u)->second;
            if (np.active == false && np.inWeight >= np.threshold + tol) {
//                ppIn[u] += 1;
                np.active = true;
//                cov++;
                cov1++;

                // add u's neighbors to T
                FriendsMap* neighbors = AM->find(u);
                if (neighbors != NULL) {
                    // for each neighbor w of u
                    for (FriendsMap::iterator k = neighbors->begin(); k!=neighbors->end(); ++k) {
                        UID w = k->first;
                        // is w is in S, no need to do anything
                        if (S.find(w) != S.end()) continue;

                        // if w is not in S, locate it in Q
                        map<UID, NodeParams>::iterator it = Q.find(w);

                        if (it == Q.end()) {
                            // if it is not in Q, then
                            NodeParams& np_w = Q[w];
                            np_w.threshold = ((float)(rand() % 1001))/(float)1000;
//                            np_w.threshold = (float)rand()/RAND_MAX;
                            np_w.active = false;
                            np_w.inWeight = k->second;
                            T.push(w);
                        } else {
                            // if w is in Q, then
                            NodeParams& np_w = it->second;
                            if (np_w.active == false) {
                                T.push(w);
                                np_w.inWeight += k->second;

                                if (np_w.inWeight - 1 > tol) {
                                    cout << "Something wrong, the inweight for a node is > 1. (w, inweight) = " << w << ", " << np_w.inWeight - 1<< endl;
                                }
                            }
                        }
                    }
                }
            }

            // deletes the first element
            T.pop();

        } //endwhile
//        cout << "Coverage in this iteration: " << cov1 << endl;
		cov += cov1/countIterations;
    }

	/*
    // coverage from ppIn
    float cov1 = 0;


    for (map<UID, float>::iterator i = ppIn.begin(); i!=ppIn.end(); ++i) {
        cov1 += i->second;
        cout << "ppIN: " << i->first << " " << i->second/countIterations << endl;
    }
    
    cov1 = cov1/countIterations;
	*/
//    double retCov = (double) cov/countIterations;

//    cout << "(cov, retCov, cov1) = " << cov << ", " << retCov << ", " << cov1 << endl;

    return cov;
   

}


/* pMG: the MG structue of node curV.  curV's ID is included */
bool MC::LTCovPlus(MGStruct *pMG, MGStruct *pBestMG, UserList &S) {
// S is the current seed set, without the current candidate node or its previous best node

    startIt = strToInt(opt->getValue("startIt"));

    int seedSetSize = S.size();

    if(pMG == NULL) {
		cout << "Some problem 3" << endl;
	}

    UID x = pMG->nodeID; // current candidate seed
    S.insert(x); // temporarily insert x into seed set (should remove it later in this function)

    bool ret = false; // indicates if the current node is better than the previous best
    double cov = 0; // compute cov(S + currentNode)
	double covPlus = 0; // compute cov(S + currentNode + prevBestNode)
    float tol = 0.00001;

    for (int b = 0; b < countIterations; ++b) {
        
        queue<UID> T; // T is the set of nodes that are to be processed
        // Q is the set of nodes that have been seen until now
        // Thus, Q is a superset of T
        // Q contains the nodes that are processed and to be processed
        map<UID, NodeParams> Q;

//        cov += S.size(); // including S+x themselves
//        covPlus += S.size();
       	double cov1 = S.size(); // cov1 == coverage in one run .. that is, number of nodes reachable
       	double covPlus1 = S.size(); // cov1 == coverage in one run .. that is, number of nodes reachable
        
        // add neighorhood of S (current seeds, plus x) to queue
        for (UserList::iterator i = S.begin(); i != S.end(); ++i) {
            UID v = *i;
            
            FriendsMap *neighbors = AM->find(v);
            if (neighbors != NULL) {
                for (FriendsMap::iterator j = neighbors->begin(); j!=neighbors->end(); ++j) {
                    UID u = j->first;
                    
                    if (S.find(u) == S.end()) {
                        if (Q.find(u) == Q.end()) {
                            // if the node u has not been seen before
                            // create a new NodeParams
                            // create a random threhsold
                            NodeParams& np = Q[u];
                            np.active = false;
                            np.inWeight = j->second;


                            /* generate secret number: */
                            np.threshold = ((float)(rand() % 1001))/(float)1000;
                            T.push(u);
                        } else {
                            NodeParams& np = Q[u];
                            np.inWeight += j->second;                   
                        } //endif
                    } //endif
                } //endfor
            } //endif
  
        } //endfor
        

        bool is_best_covered = false;

        while ( !T.empty()) {
            UID u = T.front();
//            cout << "T.size " << T.size() << endl;
            NodeParams& np = Q.find(u)->second;
            if (np.active == false && np.inWeight >= np.threshold + tol) {
                np.active = true;

				if (pBestMG != NULL && u == pBestMG->nodeID) {
					is_best_covered = true;
				}

                cov1++;
                covPlus1++;

                // add u's neighbors to T
                FriendsMap* neighbors = AM->find(u);
                if (neighbors != NULL) {
                    // for each neighbor w of u
                    for (FriendsMap::iterator k = neighbors->begin(); k!=neighbors->end(); ++k) {
                        UID w = k->first;
                        // is w is in S, no need to do anything
                        if (S.find(w) != S.end()) continue;

                        // if w is not in S, locate it in Q
                        map<UID, NodeParams>::iterator it = Q.find(w);

                        if (it == Q.end()) {
                            // if it is not in Q, then
                            NodeParams& np_w = Q[w];
                            np_w.threshold = ((float)(rand() % 1001))/(float)1000;
                            np_w.active = false;
                            np_w.inWeight = k->second;
                            T.push(w);

							
                        } else {
                            // if w is in Q, then
                            NodeParams& np_w = it->second;
                            if (np_w.active == false) {
                                T.push(w);
                                np_w.inWeight += k->second;


                                if (np_w.inWeight - 1 > tol) {
                                    cout << "Something wrong, the inweight for a node is > 1. (w, inweight) = " << w << ", " << np_w.inWeight - 1<< endl;
                                } //endif
                            } //endfor
                        } //endif
                    } //endfor
                } //endif
            } ///////////////////////////////////////////// 
            // deletes the first element
            T.pop();

        } //endwhile

        // CELF++ come here!!!!
        // if prev_best is not covered yet..
        if(pBestMG != NULL && seedSetSize >= (startIt-1)  && is_best_covered == false) {
            UID bestID = pBestMG->nodeID;
            covPlus1++;

            FriendsMap *neighbors = AM->find(bestID);
            if (neighbors != NULL) {
                // add all neighbors of bestNode to T and Q.
                for (FriendsMap::iterator j = neighbors->begin(); j != neighbors->end(); ++j) {
                    UID u = j->first;
                    if (S.find(u) != S.end()) continue;

                    if (Q.find(u) == Q.end()) {
                        // if the node u has not been seen before
                        // create a new NodeParams
                        // create a random threhsold
                        NodeParams& np = Q[u];
                        np.active = false;
                        np.inWeight = j->second;

                        /* generate secret number: */
                        np.threshold = ((float)(rand() % 1001))/(float)1000;
                        T.push(u);
                    } else {
                        NodeParams& np = Q[u];
                        if (np.active == false) {
                            np.inWeight += j->second;
                            T.push(u);
                        }//endif
                    }//endif 
                } //endfor

            } //endif
        
            // propagate through bestNode's neighborhood
             
			while ( !T.empty()) {
				UID u = T.front();
				NodeParams& np = Q.find(u)->second;
				if (np.active == false && np.inWeight >= np.threshold + tol) {
					np.active = true;

					covPlus1++;

					// add u's neighbors to T
					FriendsMap* neighbors = AM->find(u);
					if (neighbors != NULL) {
						// for each neighbor w of u
						for (FriendsMap::iterator k = neighbors->begin(); k!=neighbors->end(); ++k) {
							UID w = k->first;
							// is w is in S, no need to do anything
							if (w == bestID || S.find(w) != S.end()) continue;

							// if w is not in S, locate it in Q
							map<UID, NodeParams>::iterator it = Q.find(w);

							if (it == Q.end()) {
								// if it is not in Q, then
								NodeParams& np_w = Q[w];
								np_w.threshold = ((float)(rand() % 1001))/(float)1000;
								np_w.active = false;
								np_w.inWeight = k->second;
								T.push(w);
							} else {
								// if w is in Q, then
								NodeParams& np_w = it->second;
								if (np_w.active == false) {
									T.push(w);
									np_w.inWeight += k->second;
									if (np_w.inWeight - 1 > tol) {
										cout << "Something 2 wrong, the inweight for a node is > 1. (w, inweight) = " << w << ", " << np_w.inWeight - 1<< endl;
									} //endif
								} //endif
							} //endif
						} //endfor
					} //endif
				} //endif
				// deletes the first element
				T.pop();

			} //endwhile
        } //endif
   
		cov += (double) 1.0 * cov1 / countIterations;
		covPlus += (double) 1.0 * covPlus1 / countIterations;

    } //endfor -- each MC run

    // After the whole MC process
	S.erase(x); // remove x from current seed set now

	
//    double cov1 = (double)cov / countIterations;
//    double covPlus1 = (double) covPlus / countIterations;

    pMG->gain = cov - totalCov;

    if (pBestMG == NULL) {
        // for first run of each iteration, Best Node is not specified
        pMG->v_best = pMG->nodeID;
        pMG->gain_next = 0;

        if (seedSetSize >= (startIt-1)) {
            return true;
        } else {
            return false;
        } //endif

    } else {
        // if Best Node has been specified before
        if (pMG->gain > pBestMG->gain) {
            // this node is better than its prev. best
            pMG->v_best = pMG->nodeID;
            pMG->gain_next = 0;
            ret = true;
        } else {
            // this node is NO better than its prev. best
            pMG->v_best = pBestMG->nodeID;
            pMG->gain_next = covPlus - totalCov - pBestMG->gain; // MG(u|S+x)
            ret = false;
        } //endif

    } //endif

    return ret;

} // end of function



void MC::printVector(vector<UID>& vec, float pp) {
	cout << "AMIT " << pp << " " ;
	for (vector<UID>::iterator i=vec.begin(); i!=vec.end(); ++i) {

		cout << *i << " ";
	}

	cout << endl;

}


void MC::writeInFile(UID v, float cov, float marginal_gain, int curTimeStep, float actualCov, float actualMG, int countUsers) {
	cout << endl << endl << "Picked a seed node: " << v << ", total: " << curSeedSet.size() << endl;
	outFile << v << " " << cov << " " << marginal_gain << " " << curTimeStep << " " << getCurrentMemoryUsage() << " " << getTime() <<  " " << getTime_cur() << " " << actualCov << " " << actualMG << " " << countUsers << endl;
	cout << v << " " << cov << " " << marginal_gain << " " << curTimeStep << " " << getCurrentMemoryUsage() << " " << getTime() << " " << getTime_cur() <<  " " << actualCov << " " << actualMG << " " << countUsers << endl;
	cout << endl << endl;
}

void MC::writePOIChooseToFile(int p, float cov, float marginal_gain, int curTimeStep, float actualCov, float actualMG, int countPOIs) {
        cout << endl << endl << "Select a candidate poi: " << p << ", total: " << curPOISet.size() << endl;
        outFile << p << " " << cov << " " << marginal_gain << " " << curTimeStep << " " << getCurrentMemoryUsage() << " " << getTime() <<  " " << getTime_cur() << " " << actualCov << " " << actualMG << " " << countPOIs << endl;
        cout << p << " " << cov << " " << marginal_gain << " " << curTimeStep << " " << getCurrentMemoryUsage() << " " << getTime() << " " << getTime_cur() <<  " " << actualCov << " " << actualMG << " " << countPOIs << endl;
        cout << endl << endl;
}



void MC::openOutputFiles() {

	if (outFile.is_open()) {
		outFile.close();   //qudiao
	}
	
	string algorithm = "Greedy";
	cout << "problem : " << "MaxInfBRGSTkNN" << endl;
	string filename;
	
	int celfPlus = strToInt(opt->getValue("celfPlus"));	
	if (celfPlus == 1) {  //sfxz
		filename = outdir + "/" +"MaxInfBRGSTQ_Greedy_CELF_PLUS.txt";
	} else {
        filename = outdir + "/" +"MaxInfBRGSTQ_Greedy_CELF.txt";
	}

	outFile.open (filename.c_str());

	if (outFile.is_open() == false) {
		cout << "Can't open file " << filename << " for writing" << endl;
		exit(1);
	}
}

void MC::openExpResutlsFiles() {

        if (outFile.is_open()) {
            outFile.close();   //qudiao
        }

        string algorithm = "Greedy";
        cout << "solve problem : " << "MaxInfBRGSTkNN" << endl;
        string filename;

        int celfPlus = strToInt(opt->getValue("celfPlus"));
        if (celfPlus == 1) {  //sfxz
            filename = outdir + "/" +"MaxInfBRGSTQ_Greedy_CELF_PLUS.txt";
        } else {
            filename = outdir + "/" +"MaxInfBRGSTQ_Greedy_CELF.txt";
        }

        outFile.open (filename.c_str());

        if (outFile.is_open() == false) {
            cout << "Can't open file " << filename << " for writing" << endl;
            exit(1);
        }
    }


PropModels MC::getModel() {
	return model;
}

float MC::getTime_cur() const {
	time_t curTime;
	time(&curTime);

	float min = ((float)(curTime - stime_mintime))/60;
	return min;
}

float MC::getTime() const {
	time_t curTime;
	time(&curTime);

	float min = ((float)(curTime - startTime))/60;
	return min;
}


void MC::readInputData(float alpha) {
	cout << "in readInputData for model " << model << " with alpha " << alpha << endl;

	unsigned int edges = -1;

	string probGraphFile = opt->getValue("probGraphFile");
	cout << "Reading file " << probGraphFile << endl;  //dutu
	ifstream myfile (probGraphFile.c_str(), ios::in);
	string delim = " \t";
	if (myfile.is_open()) {
		while (! myfile.eof() )	{
			std::string line;
			getline (myfile,line);
			if (line.empty()) continue;
			edges++;
			if (edges == 0) continue; // ignore the first line
			std::string::size_type pos = line.find_first_of(delim);
			int	prevpos = 0;

			// get first user
			string str = line.substr(prevpos, pos-prevpos);
			UID u1 = strToInt(str);

			// get the second user
			prevpos = line.find_first_not_of(delim, pos);
			pos = line.find_first_of(delim, prevpos);
			UID u2 = strToInt(line.substr(prevpos, pos-prevpos));

			if (u1 == u2) continue;

			// get the parameter
			float parameter1 = 0;
			prevpos = line.find_first_not_of(delim, pos);
			pos = line.find_first_of(delim, prevpos);
			if (pos == string::npos)
				parameter1 = strToFloat(line.substr(prevpos));
			else
				parameter1 = strToFloat(line.substr(prevpos, pos-prevpos));

			if (parameter1 == 0) continue;

			parameter1 = parameter1 ;
			users.insert(u1);
			users.insert(u2);


			if (edges % 1000000 == 0) {
				cout << "(node1, node2, weight,  AM size till now, edges till now, mem) = " << u1 << ", " << u2 << ", " << parameter1 << ", " << AM->size() << ", " << edges << ", " << getCurrentMemoryUsage() << endl;
			}


			FriendsMap* neighbors = AM->find(u1);
			if (neighbors == NULL) {
				neighbors = new FriendsMap();
				neighbors->insert(pair<UID, float>(u2, parameter1));
				AM->insert(u1, neighbors);
			} else {
				//FriendsMap::iterator it = neighbors->find(u2);
				//if (it == neighbors->end()) {
					neighbors->insert(pair<UID, float>(u2, parameter1));
				//} else {
				//	cout << "WARNING: Edge redundant between users " << u1 << " and " << u2 << endl;
				//}
			}

            if (revAM != NULL) {
                multimap<float, UID> *revNeighbors = revAM->find(u1);
                if (revNeighbors == NULL) {
                    revNeighbors = new multimap<float, UID>();
                    revNeighbors->insert(std::make_pair(parameter1, u2));
                    revAM->insert(u1, revNeighbors);
                } else {
                    revNeighbors->insert(std::make_pair(parameter1, u2));
                }
            }

            // also add the edges u2->u1 but done allocate Edge class to them
			// .. it is just to find friends efficiently
			if (graphType == UNDIRECTED) {
				neighbors = AM->find(u2);
				if (neighbors == NULL) {
					neighbors = new FriendsMap();
					neighbors->insert(pair<UID, float>(u1, parameter1 ));
					AM->insert(u2, neighbors);
				} else {
					FriendsMap::iterator it = neighbors->find(u1);
					if (it == neighbors->end()) {
						neighbors->insert(pair<UID, float>(u1, parameter1 ));
					} else {
						//	cout << "WARNING: Edge redundant between users " << u1 << " and " << u2 << endl;
					}
				}
			}
		}
		myfile.close();
	} else {
		cout << "Can't open friendship graph file " << probGraphFile << endl;
	}

    this->numEdges = edges;

    cout << "BuildAdjacencyMatFromFile done" << endl;
	cout << "Size of friendship graph hashtree is : " << AM->size() << endl;
    //cout << "Size of NEW friendship graph hashtree is : " << revAM->size() << endl;
	cout << "Number of users are: " << users.size() << endl;
	cout << "Number of edges in the friendship graph are: " << edges << endl;
}


}



