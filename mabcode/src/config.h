/*
 * config.h
 *
 *  Created on: Dec 28, 2018
 *      Author: jmshinju
 */

#ifndef CONFIG_H_
#define CONFIG_H_

#ifdef WIN32
#define FILESEP "\\"
#else
#define FILESEP "/"
#endif

#include <iostream>

static bool verbose = 0;

const std::string QUERY = "query";
const std::string GENERATE_ARMS = "generate-arms";

const std::string EGDE = "egde";
const std::string DELTA_E = "delta_e";
const std::string TOPK_DELTA_E = "topk_delta_e";
const std::string D_DELTA_E = "d_delta_e";
const std::string TOPK_D_DELTA_E = "topk_d_delta_e";

const std::string UNIFORM = "uniform";
const std::string SEGMENT = "segment";
const std::string NORMAL = "normal";
class Config {
public:
	int num_arms = 0;
	double epsilon = 0.01;
	double delta = 0.05;
	int k = 1;
	double c = 8.0;
	double Q = 0.0;
	int R = 2;

	std::string action = ""; // query/generate index, etc..

	std::string arm_dist = UNIFORM;

	std::string prefix;

	std::string get_arm_file() {
		if (arm_dist == UNIFORM || arm_dist == NORMAL)
			return this->prefix + "arms_" + std::to_string(num_arms) + "_" + arm_dist + ".txt";
		else
			return this->prefix + "arms_" + std::to_string(num_arms) + "_" + arm_dist + "_" + std::to_string(k) + ".txt";
	}

	std::string algo;
	int repeat = 1;
	std::string result_dir;

//	int adjust_scale = 0;//bool 0 or 1

	double scale_factor_ME = 1;
	double scale_factor_US = 1;

	double segment_high = 0.5;
	double segment_low = 0.49;

	void display() {
		std::cerr << "\naction " << action << "\nalgo " << algo << "num_arms " << num_arms << "\nepsilon " << epsilon
				<< "\ndelta " << delta << "\nk " << k << "\nc " << c << "\nQ " << Q << "\nR " << R << "\narm_dist "
				<< arm_dist << "\nprefix " << prefix << "\nrepeat " << repeat << std::endl;
	}
};

extern Config config;

#endif /* CONFIG_H_ */
