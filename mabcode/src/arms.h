/*
 * arms.h
 *
 *  Created on: Dec 28, 2018
 *      Author: jmshinju
 */

#ifndef ARMS_H_
#define ARMS_H_

#include <vector>
#include <tuple>

#include "config.h"
#include "utils.h"

std::vector<std::tuple<int, double, double> > load_arms() {
	std::vector<std::tuple<int, double, double> > arms;

	std::string filename = config.get_arm_file();
	std::cerr << filename << std::endl;
	if (!file_exists_test(filename)) {
		std::cerr << "file does not exist, please generate files first" << std::endl;
		exit(0);
	}
	std::ifstream armfile(filename);
	int arm_id = 0;
	double real_avg;
	double emprical_avg = 0.0;
	while (armfile >> real_avg) {
		arms.push_back(std::make_tuple(arm_id, real_avg, emprical_avg));
		arm_id++;
	}
//	std::cerr << "loaded " << arms.size() << " arms" << std::endl;
	return arms;
}

#endif /* ARMS_H_ */
