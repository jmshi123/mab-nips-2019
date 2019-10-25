/*
 * Performance.h
 *
 *  Created on: Dec 29, 2018
 *      Author: jmshinju
 */

#ifndef PERFORMANCE_H_
#define PERFORMANCE_H_

#include <vector>
class Performance {

public:
	double total_pulls = 0;
	double iterations = 0;
	double itertaions_old_before_next_repeat = 0;
	double best_arm_accuracy = 0;
	double runtime = 0;
	double regret = 0;

	int round_max = 0;
	double candidate_size_per_iteration[1000] = { 0.0 };
	double cost_per_itertaion[1000] = { 0.0 };
	double cost_accumulate[1000] = { 0.0 };

	std::pair<double, double> correct_error_cnt_per_itertaion[1000];
	std::pair<double, double> correct_error_cnt_per_itertaion_finished[1000];
	std::pair<double, double> correct_error_cnt_per_itertaion_finished_accumulate[1000];

	void init_accuracy_per_iteration() {
		for (int i = 0; i < 1000; ++i) {
			correct_error_cnt_per_itertaion[i] = std::make_pair(0.0, 0.0);
			correct_error_cnt_per_itertaion_finished[i] = std::make_pair(0.0, 0.0);
			correct_error_cnt_per_itertaion_finished_accumulate[i] = std::make_pair(0.0, 0.0);
		}
	}

	void display() {
		std::cout << "avg perf: rounds: " << iterations << ", pulls: " << (long long)total_pulls << ", accuracy: " << best_arm_accuracy << ", time: " << runtime << std::endl;
//		std::cout << iterations << " " << (long long) total_pulls << " " << best_arm_accuracy << " " << runtime
//				<< std::endl;

//		for (int r = 0; r < 1000; ++r) {
//			if (candidate_size_per_iteration[r] > 0) {
//				std::cout << candidate_size_per_iteration[r] << " ";
//			}
//		}
//		std::cout << std::endl;

//		for (int r = 0; r <= round_max; ++r) {
//			std::cout << cost_per_itertaion[r] / (double) config.repeat << " ";
//		}
//		std::cout << std::endl;

//		for (int r = 0; r <= round_max; ++r) {
//			std::cout << (long long)cost_accumulate[r] << " ";
//		}
//		std::cout << std::endl;

//		for (int r = 0; r <= round_max; ++r) {
//			std::cout << correct_error_cnt_per_itertaion[r].first << " ";
//		}
//		std::cout << std::endl;
//		for (int r = 0; r <= round_max; ++r) {
//			std::cout << correct_error_cnt_per_itertaion[r].second << " ";
//		}
//		std::cout << std::endl;

//		for (int r = 0; r <= round_max; ++r) {
//			std::cout << correct_error_cnt_per_itertaion_finished[r].first << " ";
//		}
//		std::cout << std::endl;
//		for (int r = 0; r <= round_max; ++r) {
//			std::cout << correct_error_cnt_per_itertaion_finished[r].second << " ";
//		}
//		std::cout << std::endl;

//		for (int r = 0; r <= round_max; ++r) {
//			std::cout << correct_error_cnt_per_itertaion_finished_accumulate[r].first << " ";
//		}
//		std::cout << std::endl;
//		for (int r = 0; r <= round_max; ++r) {
//			std::cout << correct_error_cnt_per_itertaion_finished_accumulate[r].second << " ";
//		}
//		std::cout << std::endl;

//		for (int r = 0; r <= round_max; ++r) {
//			double correct_not_finished = correct_error_cnt_per_itertaion[r].first;
//			double total_not_finished = correct_error_cnt_per_itertaion[r].first
//					+ correct_error_cnt_per_itertaion[r].second;
//			double correct_finished_accumulate = correct_error_cnt_per_itertaion_finished_accumulate[r].first;
//			double total_finished = correct_error_cnt_per_itertaion_finished_accumulate[r].first
//					+ correct_error_cnt_per_itertaion_finished_accumulate[r].second;
//			double accuracy = 0;
//			if (total_not_finished + total_finished > 0)
//				accuracy = (correct_not_finished + correct_finished_accumulate) / (total_not_finished + total_finished);
//			std::cout << accuracy << " ";
//		}
//		std::cout << std::endl;
	}

};

#endif /* PERFORMANCE_H_ */
