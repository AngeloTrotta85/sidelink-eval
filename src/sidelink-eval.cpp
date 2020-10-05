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
#include <complex>      // std::complex, std::real

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

	map<int, double> rcvRSS_fromStart;
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

double linear2dBm(double x) {
	return (10.0 * log10(x * 1000.0));
}

double dbm2linear(double x) {
	return pow(10.0, (x-30.0)/10.0);
	//return (10.0 * log10(x * 1000.0));
}

double getRealNormal (double mean, double stdev, default_random_engine &generator_rand) {
	std::normal_distribution<double> n_distribution(mean, stdev);
	return n_distribution(generator_rand);
}

double getUAVChannelLoss (double freq, double tx_height, double rx_height, double dist) {
	//double C = 0;
	//double temp = rx_height * (1.1 * log10(freq) - 0.7) - (1.56 * log10(freq) - 0.8);
	//double sigma = 8; // Standard Deviation for Shadowing Effect

	double path_loss;
	path_loss = 41.1 + 20.9 * log10(dist);
	//double path_loss = 41.1 + 41.8 * log10(dist);
	//path_loss = 41.1 + 41.8 * log10(dist);

	//path_loss = (20.0 * log10(dist)) + (20.0 * log10(freq)) - 27.55;	//free space path loss
	//double alpha = 2.05;
	//path_loss = ((20.0 * log10(1.0)) + (20.0 * log10(freq)) - 27.55) + (10.0 * alpha * log10(dist));	//log distance path loss (free-space at d_0=1)

	/*if (randomness) {
		//double path_loss = 46.3 + 33.9 * log10(freq) - 13.82 * log10(tx_height) - temp + log10(dist/1000.0)*(44.9 - 6.55 * log10(tx_height))+C;
		double channel_loss = -path_loss + (-1 * sigma * RandomGenerator::getInstance().getRealNormal(0, 1));

		return channel_loss;
	}
	else {*/
		return -path_loss;
	//}
}

double rss(MyCoord src, MyCoord dst) {

	double distance = src.distance(dst);

	double UAV_TX_pt_db = 23; //24;
	//double UAV_TX_pt_watt = 0.2512;
	//double GAIN_ag_watt = pow(10.0, 0.6);  //% Transmiter Antenna Gain (6 dB)
	double GAIN_ag_db = 6;
	double freq = 3410;
	double tx_height = 30;
	double rx_height = 30;
	//double distance = u_src->actual_coord.distance(u_dst->actual_coord);
	double loss_dB = getUAVChannelLoss(freq,tx_height,rx_height,distance);

	double receivedPower_db = UAV_TX_pt_db + GAIN_ag_db + loss_dB;
	double receivedPower = pow(10.0, receivedPower_db / 10.0) / 1000.0;
	//double receivedPower = UAV_TX_pt * GAIN_ag * pow(10.0, loss_dB/10.0) ;//(10.^((loss_dB)/10)); // received power including path loss,shadowing

	return receivedPower;
}

double rss_with_fading(MyCoord src, MyCoord dst, default_random_engine &generator_rand) {

	double receivedPower = rss(src, dst);

	double fading_variance = 1.59; // Fading model is Complex Gaussian Random Variable
	double chan_real = sqrt(fading_variance/2) * getRealNormal(0, 1, generator_rand);
	double chan_complex = getRealNormal(0, 1, generator_rand);
	std::complex<double> chan (chan_real, chan_complex);
	//chan = sqrt(fading_variance/2) * RandomGenerator::getInstance().getRealNormal(0, 1) + 1i*RandomGenerator::getInstance().getRealNormal(0, 1);
	double chan_value = std::abs(chan); // fading loss

	receivedPower *= chan_value;


	return receivedPower;
}

double getProb_sigmoid(double sinr) {
	return (1.0 / (1.0 + exp((10.0 - sinr)*0.25)));
}

double getProb_linear(double sinr) {
	if (sinr <= -5.0) {
		return 0.0;
	}
	else if (sinr >= 25.0) {
		return 1.0;
	}
	else {
		return ((sinr+5.0)/30.0);
	}
}

double calculateProbability(default_random_engine &generator_rand, Arc *actTx, list<Arc *> &uavTxList, map<int, MyCoord> &uavPosMap, double dMaxUAV) {
	double pRis = 1;

	// Calculate Noise Parameters
	double temperature = 290; // Kelvin
	double k = 1.3806488 * pow(10.0, -23.0); // Boltzman Constant
	double bw = 9*1e6; // Effective Bandwidth of channel (9 MHz)
	double ue_noise_figure = 7 ; // 7 dB noise figure is considered
	double noise = linear2dBm(k * temperature * bw);
	double total_noise_dBm = ue_noise_figure + noise;
	double total_noise = pow(10.0, total_noise_dBm/10.0) / 1000.0;

	double rcvPow = actTx->rcvRSS_fromStart[actTx->nodeEnd]; //rss_with_fading(uavPosMap[actTx->nodeStart], uavPosMap[actTx->nodeEnd], generator_rand);

	cout << "      Dist: " << uavPosMap[actTx->nodeStart].distance(uavPosMap[actTx->nodeEnd]) << "; RcvPow: " << rcvPow << "; Noise: " << total_noise;

	double sumInterf = 0;
	cout << "; [";
	for (auto& interfLink : uavTxList) {
		if (interfLink->id != actTx->id) {
			double interActLink = 0;
			if (	(interfLink->channel == actTx->channel) &&
					(uavPosMap[interfLink->nodeStart].distance(uavPosMap[actTx->nodeEnd]) <= dMaxUAV)
			){
				interActLink = interfLink->rcvRSS_fromStart[actTx->nodeEnd];//rss_with_fading(uavPosMap[interfLink->nodeStart], uavPosMap[actTx->nodeEnd], generator_rand);
				sumInterf += interActLink;
				cout << "!";
			}
			cout << "U:" << interfLink->nodeStart << "=" << interActLink << ";d:" << uavPosMap[interfLink->nodeStart].distance(uavPosMap[actTx->nodeEnd]) << "|";
		}
	}
	cout << "]";

	double sinr = linear2dBm(rcvPow / (total_noise + sumInterf));

	cout << "; TotInt: " << sumInterf << "; SINR: " << sinr;

	pRis = getProb_linear(sinr);

	cout << endl;

	return pRis;
}

int main(int argc, char **argv) {

	default_random_engine generator_rand = std::default_random_engine();

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

	const std::string &seedUser = input.getCmdOption("-seed");
	if (!seedUser.empty()) {
		int seedR = atoi(seedUser.c_str());
		generator_rand.seed(seedR);
	}
	else {
		unsigned seedR = std::chrono::system_clock::now().time_since_epoch().count();
		generator_rand.seed(seedR);
	}
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

	const std::string &dm_string = input.getCmdOption("-dm");
	if (!dm_string.empty()) {
		distMaxUAV = stod(dm_string);
	}
	const std::string &db_string = input.getCmdOption("-db");
	if (!db_string.empty()) {
		distMaxBS = stod(db_string);
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

		    //cout << "      Rcv Signal -> ";
		    for (auto& intUAV : uavPos) {
		    	if (intUAV.first != a->nodeStart) {
		    		double rss_rcv = rss_with_fading(uavPos[a->nodeStart], intUAV.second, generator_rand);
		    		//a->rcvRSS_fromStart[a->nodeEnd] = rss_rcv;
		    		a->rcvRSS_fromStart[intUAV.first] = rss_rcv;

		    		//cout << "U:" << intUAV.first << "=" << rss_rcv << " ";
		    	}
		    }
		    //cout << endl;

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

			cout << "      Rcv Signal -> ";
			for (auto& dstUAV : l->rcvRSS_fromStart) {
				cout << "U:" << dstUAV.first << "=" << dstUAV.second << " ";
			}
			/*for (auto& intUAV : uavPos) {
				if (intUAV.first != l->nodeStart) {
					double rss_rcv = rss_with_fading(uavPos[l->nodeStart], intUAV.second, generator_rand);
					l->rcvRSS_fromStart[l->nodeEnd] = rss_rcv;

					cout << "U:" << intUAV.first << "=" << rss_rcv << " ";
				}
			}*/
			cout << endl;
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

		cout << "For packet generated at PoI " << m.first.first << " at time " << m.first.second << " we have: " << endl;

		for (auto& l : m.second) {
			double linkProb = 1;
			double linkProb_noLimit = 1;

			cout << "   From " << l->nodeStart << " to " << l->nodeEnd << " using channel " << l->channel << " at time slot " << l->txTime << endl;

			if (l->nodeEnd == 0) {
				totDelay = l->txTime - l->genTime;
				if (uavPos[l->nodeStart].length() > distMaxBS) {
					linkProb = 0;
					rcvBS = false;
				}
				else {
					rcvBS = true;
				}
			}
			else if (arcMap_tx[l->txTime].size() > 0){
				//linkProb = 0.5;	//TODO
				linkProb = calculateProbability(generator_rand, l, arcMap_tx[l->txTime], uavPos, distMaxUAV);
				//linkProb_noLimit = 0.5;	//TODO
				linkProb_noLimit = calculateProbability(generator_rand, l, arcMap_tx[l->txTime], uavPos, numeric_limits<double>::max());
			}

			cout << "   Final lProb: " << linkProb << " lProbNL: " << linkProb_noLimit << endl;

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
				<< endl << endl;

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
