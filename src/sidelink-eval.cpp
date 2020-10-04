//============================================================================
// Name        : sidelink-eval.cpp
// Author      : Angelo
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <stdlib.h>
#include <stdio.h>
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <algorithm>    // std::find
#include <vector>       // std::vector
#include <list>       // std::list
#include <stack>
#include <map>       // std::list
#include <cstdlib>
#include <ctime>
#include <random>
#include <chrono>

#include <boost/algorithm/string.hpp>

#include "MyCoord.h"

using namespace std;

class InputParser{
public:
	InputParser (int &argc, char **argv){
		for (int i=1; i < argc; ++i)
			this->tokens.push_back(std::string(argv[i]));
	}
	const std::string& getCmdOption(const std::string &option) const{
		std::vector<std::string>::const_iterator itr;
		itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
		if (itr != this->tokens.end() && ++itr != this->tokens.end()){
			return *itr;
		}
		static const std::string empty_string("");
		return empty_string;
	}
	bool cmdOptionExists(const std::string &option) const{
		return std::find(this->tokens.begin(), this->tokens.end(), option)
		!= this->tokens.end();
	}
private:
	std::vector <std::string> tokens;
};

class Arc {
public:
	int id;
	int nodeStart;
	int nodeEnd;
	int txTime;
	int channel;
	int genPoI;
	int genTime;
};

class Stat {
public:
	Stat(){
		probRcv = -1;
		probRcv_nolimit = -1;
		delay = 0;
		isRcv = 0;
		nHop = 0;
	};
public:
	double probRcv;
	double probRcv_nolimit;
	double delay;
	double isRcv;
	double nHop;
};

bool compare_packets (const Arc *first, const Arc *second) {
	return (first->txTime < second->txTime);
}

int main(int argc, char **argv) {

	map<int, MyCoord> uavPos;

	map<pair<int, int>, list<Arc *>> arcMap;
	map<int, list<Arc *>> arcMap_tx;

	map<pair<int, int>, Stat *> risMap;

	string fin_pos = string("input_pos.txt");
	string fin = string("input.txt");
	string fout = string("ris.txt");

	double distMaxUAV = 1000.0;
	double distMaxBS = 1000.0;

	InputParser input(argc, argv);

	const std::string &in_string = input.getCmdOption("-fin");
	if (!in_string.empty()) {
		fin = in_string;
	}
	const std::string &inpos_string = input.getCmdOption("-fpos");
	if (!inpos_string.empty()) {
		fin_pos = inpos_string;
	}

	ifstream infile_pos;
	infile_pos.open (fin_pos, std::ifstream::in);
	if (infile_pos.is_open()){
		bool continue_read = true;
		std::string line;

		cout << "FILE " << fin_pos << endl;

		while ( (std::getline(infile_pos, line)) && (continue_read) ) {
			std::string delimiter_field = ";";
			std::string delimiter_eq = ":";

		    cout << "Line: " << line << endl;

		    vector<string> strs;
		    boost::split(strs, line, boost::is_any_of(";"));

		    int uavIdx = -1;

		    for (auto& var : strs) {
		    	cout << var << endl;

		    	vector<string> strs_var;
		    	boost::split(strs_var, var, boost::is_any_of(":"));

		    	for (auto& el : strs_var) {
		    		cout << el << endl;
		    	}

		    	if (strs_var.size() == 2) {
		    		if (strs_var[0].compare("U") == 0) {
		    			uavIdx = stoi(strs_var[1]);
		    			uavPos[uavIdx] = MyCoord::ZERO;
		    		}
		    		else if (strs_var[0].compare("x") == 0) {
		    			uavPos[uavIdx].x = stod(strs_var[1]);
		    		}
		    		else if (strs_var[0].compare("y") == 0) {
		    			uavPos[uavIdx].y = stod(strs_var[1]);
		    		}
		    		else if (strs_var[0].compare("POI") == 0) {
		    			break;
		    		}
		    		else if (strs_var[0].compare("BS") == 0) {
		    			break;
		    		}
		    		else if (strs_var[0].compare("DR") == 0) {

		    		}
		    		else if (strs_var[0].compare("DB") == 0) {
		    			distMaxBS = stod(strs_var[1]);
		    		}
		    		else if (strs_var[0].compare("DM") == 0) {
		    			distMaxUAV = stod(strs_var[1]);

		    			continue_read = false;
		    		}
		    	}
		    }

		}

		infile_pos.close();
	}

	for (auto& u : uavPos) {
		cout << "UAV" << u.first << " at pos: " << u.second << endl;
	}

	ifstream infile;
	infile.open (fin, std::ifstream::in);
	if (infile.is_open()){
		int ididx = 0;
		std::string line;

		while (std::getline(infile, line)) {
			std::string delimiter_field = ";";
			std::string delimiter_eq = ":";

		    cout << "Line: " << line << endl;

		    vector<string> strs;
		    boost::split(strs, line, boost::is_any_of(";"));

		    Arc *a = new Arc();
		    a->id = ididx++;

		    int colIdx = 0;
		    for (auto& var : strs) {
		    	cout << var << endl;

		    	vector<string> strs_var;
		    	boost::split(strs_var, var, boost::is_any_of(":"));

		    	for (auto& el : strs_var) {
		    		cout << el << endl;
		    	}

		    	if (strs_var.size() == 2) {
		    		if (strs_var[0].compare("U") == 0) {
		    			if (colIdx == 0) {
		    				a->nodeStart = stoi(strs_var[1]);
		    			}
		    			else {
		    				a->nodeEnd = stoi(strs_var[1]);
		    			}
		    		}
		    		else if (strs_var[0].compare("B") == 0) {
		    			a->nodeEnd = stoi(strs_var[1]);
		    		}
		    		else if (strs_var[0].compare("TS") == 0) {
		    			a->txTime = stoi(strs_var[1]);
		    		}
		    		else if (strs_var[0].compare("CH") == 0) {
		    			a->channel = stoi(strs_var[1]);
		    		}
		    		else if (strs_var[0].compare("POI") == 0) {
		    			a->genPoI = stoi(strs_var[1]);
		    		}
		    		else if (strs_var[0].compare("TX") == 0) {
		    			a->genTime = stoi(strs_var[1]);
		    		}
		    	}

		    	++colIdx;
		    }

		    //if (arcMap.count(a->gen) == 0) {
		    //	arcMap[a->gen] = std::list<Arc *>();
		    //}
		    arcMap[make_pair(a->genPoI, a->genTime)].push_back(a);
		    arcMap_tx[a->txTime].push_back(a);
		}

		infile.close();
	}

	//sort and print and check
	for (auto& m : arcMap) {
		int start = -1;
		m.second.sort(compare_packets);
		cout << "For packet generated at PoI " << m.first.first << " at time " << m.first.second << " we have: " << endl;
		for (auto& l : m.second) {
			cout << "   From " << l->nodeStart << " to " << l->nodeEnd << " using channel " << l->channel << " at time slot " << l->txTime << endl;

			if ( ((start >= 0) && (start != l->nodeStart)) || (l->nodeStart == l->nodeEnd)) {
				cerr << "Error in calculating the routing" << endl;
				cerr << "   From " << l->nodeStart << " to " << l->nodeEnd << " using channel " << l->channel << " at time slot " << l->txTime << endl;
			}
			start = l->nodeEnd;
		}
		cout << endl;

		risMap[m.first] = new Stat();
	}

	// Calculate Statistics for each packet
	for (auto& m : arcMap) {
		double totProb = 1;
		double totProb_noLimit = 1;
		double totDelay = 0;
		bool rcvBS = false;

		for (auto& l : m.second) {
			double linkProb = 1;
			double linkProb_noLimit = 1;

			if (l->nodeEnd == 0) {
				rcvBS = true;
				totDelay = l->txTime - l->genTime;
			}
			else if (arcMap_tx[l->txTime].size() > 1){
				linkProb = 0.5;	//TODO
				linkProb = calculate
				linkProb_noLimit = 0.5;	//TODO
			}

			totProb = totProb * linkProb;
			totProb_noLimit = totProb_noLimit * linkProb_noLimit;
		}

		risMap[m.first]->delay = totDelay;
		risMap[m.first]->probRcv = totProb;
		risMap[m.first]->probRcv_nolimit = totProb_noLimit;
		risMap[m.first]->nHop = m.second.size();
		risMap[m.first]->isRcv = (rcvBS ? 1 : 0);

		cout << "For packet generated at PoI " << m.first.first << " at time " << m.first.second << " we have:"
				<< " Probability: " << risMap[m.first]->probRcv
				<< " Probability(NoLimits): " << risMap[m.first]->probRcv_nolimit
				<< " Delay: " << risMap[m.first]->delay
				<< " nHops: " << risMap[m.first]->nHop
				<< " RCV?: " << risMap[m.first]->isRcv
				<< endl;

	}

	// Calculate Statistics for each packet
	double countEl = 0;
	double sumProbability = 0;
	double sumProbability_nolimit = 0;
	double sumDelay = 0;
	double sumNHops = 0;
	double sumRcv = 0;

	for (auto& m : risMap) {
		sumProbability += m.second->probRcv;
		sumProbability_nolimit += m.second->probRcv_nolimit;
		sumDelay += m.second->delay;
		sumNHops += m.second->nHop;
		sumRcv += m.second->isRcv;

		++countEl;
	}

	cout << "Total statistics are:"
			<< " Probability: " << (sumProbability / countEl)
			<< " Probability(NoLimits): " << (sumProbability_nolimit / countEl)
			<< " Delay: " << (sumDelay / countEl)
			<< " nHops: " << (sumNHops / countEl)
			<< " RCV?: " << (sumRcv / countEl)
			<< endl;

	//cout << endl << "End" << endl; // prints !!!Hello World!!!
	return EXIT_SUCCESS;
}
