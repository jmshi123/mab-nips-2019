/*
 * gen.h
 *
 *  Created on: Dec 28, 2018
 *      Author: jmshinju
 */

#ifndef GEN_H_
#define GEN_H_

void generate_uniform_arms() {
	std::string filename = config.get_arm_file();
	std::cout << filename << std::endl;
	std::ofstream armfile(filename);

	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(0.0, 1.0); //[0,1)
	for (int n = 0; n < config.num_arms; ++n) {
		double arm_avg = dis(gen);
		armfile << arm_avg << std::endl;
	}
}

void generate_normal_arms() {
	std::string filename = config.get_arm_file();
	std::cout << filename << std::endl;
	std::ofstream armfile(filename);

	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0.5, 0.2);
	int cnt = 0;
	double max = 0;
	while (cnt < config.num_arms) {
		double number = distribution(generator);
		if ((number >= 0.0) && (number <= 1.0)) {
			cnt++;
			armfile << number << std::endl;
			if (number > max) max = number;
		}
	}

	std::cout << max << std::endl;
}

std::vector<std::tuple<int, double, double> > generate_uniform_arms_online() {
	std::vector<std::tuple<int, double, double> > arms;

	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(0.0, 1.0); //[0,1)
	double empirical_avg = 0.0;
	for (int arm_id = 0; arm_id < config.num_arms; ++arm_id) {
		double real_avg = dis(gen);
		arms.push_back(std::make_tuple(arm_id, real_avg, empirical_avg));
	}
	return arms;
}

void generate_segment_arms() {
	std::string filename = config.get_arm_file();
	std::cout << filename << std::endl;
	std::ofstream armfile(filename);

	for (int n = 0; n < config.k; ++n) {
		armfile << config.segment_high << std::endl;
	}

	for (int n = config.k; n < config.num_arms; ++n) {
		armfile << config.segment_low << std::endl;
	}
}

#endif /* GEN_H_ */
