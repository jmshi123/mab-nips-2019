/*
 * MEAS.h
 *
 *  Created on: Dec 28, 2018
 *      Author: jmshinju
 */

#ifndef UTILMAB_H_
#define UTILMAB_H_
#include "Performance.h"

bool evaluateAccuracyPerItr(std::vector<std::tuple<int, double, double> >& result_ground_truth,
		std::vector<std::tuple<int, double, double> >& result, Performance& performance, int itr, bool isFinished) {
	bool accurate = true;
	for (int i = 0; i < config.k; ++i) {
		if (std::abs(std::get < 1 > (result_ground_truth[i]) - std::get < 1 > (result[i])) > config.epsilon) {
			accurate = false; //any top k with large error, the whole process fail.
			break;
		}
	}

	if (isFinished) {
		if (accurate) {
			performance.correct_error_cnt_per_itertaion_finished[itr].first += 1;
		} else {
			performance.correct_error_cnt_per_itertaion_finished[itr].second += 1;
		}
	} else {
		if (accurate) {
			performance.correct_error_cnt_per_itertaion[itr].first += 1;
		} else {
			performance.correct_error_cnt_per_itertaion[itr].second += 1;
		}
	}

	if (itr > performance.round_max) {
		performance.round_max = itr;
	}
	return accurate;
}

double getEmpiricalAvgOfArm(double arm_real_avg, long long num_samples, std::uniform_real_distribution<> &dist,
		std::mt19937 &rand_gen) {
	//bernoulli distribution
	double sum = 0.0;
	for (long long i = 0; i < num_samples; ++i) {
		sum += dist(rand_gen) < arm_real_avg ? 1.0 : 0.0;
	}
	double emprical_avg = sum / (double) num_samples;
	return emprical_avg;
}

//long long computeNumPulls_1stItr_MEAS_without_scale(double epsilon, double delta, double mu, int k, int num_arms) {
//	double epsilon_l = epsilon / 16.0;
//	double delta_l = delta / 8.0;
//	double mu_l = mu;
//	long long num_samples_per_arm = (12 / epsilon_l / epsilon_l) * (1 / mu_l) * log(6 * k / delta_l);
//	return num_samples_per_arm * num_arms;
//}

long long computeNumPulls_1stItr_TopkDE(double delta, double epsilon, double num_arms, int k,
		double scale_factor_1_TopkDE) {
	double c = 8;
	double Q = c / epsilon / epsilon;
	double beta_l = 1;
	double delta_l = delta / 4.0;

	long long num_samples_per_arm_origin = beta_l * Q * log(((double) k) / delta_l);
	long long num_samples_per_arm = num_samples_per_arm_origin / scale_factor_1_TopkDE;
	return num_samples_per_arm * num_arms;
}

//std::vector<std::tuple<int, double, double>> MedianElimination(const std::vector<std::tuple<int, double, double>> &arms,
//		double epsilon, double delta, double mu, int k, std::uniform_real_distribution<> &dist, std::mt19937 &rand_gen,
//		Performance &performance, double scale_factor, std::vector<std::tuple<int, double, double>> &ground_truth, long long *total_pulls_ME) {
//
//	std::vector<std::tuple<int, double, double>> S_l(arms.begin(), arms.end());
//	double epsilon_l = epsilon / 16.0;
//	double delta_l = delta / 8.0;
//	double mu_l = mu;
//
//	while (S_l.size() > 4 * k) {
//		long long num_samples_per_arm = (12 / epsilon_l / epsilon_l) * (1 / mu_l) * log(6 * k / delta_l) / scale_factor;
//
//		std::vector<std::tuple<int, double, double>> S_empirical;
//
//		std::cerr << "ME pulling " << num_samples_per_arm << " times for each of " << S_l.size() << " arms, total: "
//				<< num_samples_per_arm * S_l.size() << std::endl;
//		for (int i = 0; i < S_l.size(); ++i) {
//			int arm_id = std::get < 0 > (S_l[i]);
//			double real_avg = std::get < 1 > (S_l[i]);
//			double emprical_avg = getEmpiricalAvgOfArm(real_avg, num_samples_per_arm, dist, rand_gen);
//			S_empirical.push_back(std::make_tuple(arm_id, real_avg, emprical_avg));
//		}
//		//sort arms in desc order based on emperical averages
//		std::sort(S_empirical.begin(), S_empirical.end(), sortDescbyThird);
//
//		performance.total_pulls += num_samples_per_arm * S_l.size();
//		performance.cost_per_itertaion[(int) (performance.iterations + 1 - performance.itertaions_old_before_next_repeat)] +=
//				num_samples_per_arm * S_l.size();
//		total_pulls_ME[0] += num_samples_per_arm * S_l.size();
//
//		int next_size = ceil(((double) S_l.size()) / 2.0);
//		S_l.clear();
//		S_l.assign(S_empirical.begin(), S_empirical.begin() + next_size); // copy [0, next_size) to S_l for next iteration.
//		std::cerr << "next_size " << S_l.size() << std::endl;
//
//		mu_l *= (1 - epsilon_l);
//		epsilon_l *= 0.75;
//		delta_l *= 0.5;
//
//		performance.iterations++;
//		performance.candidate_size_per_iteration[(int) (performance.iterations
//				- performance.itertaions_old_before_next_repeat)] += S_l.size();
//
//		std::vector<std::tuple<int, double, double>> S_topk_per_itr(S_l.begin(), S_l.end());
//		std::sort(S_topk_per_itr.begin(), S_topk_per_itr.end(), sortDescbyThird);
//		std::vector<std::tuple<int, double, double>> S_topk_per_itr2(S_topk_per_itr.begin(),
//				S_topk_per_itr.begin() + k);
//		evaluateAccuracyPerItr(ground_truth, S_topk_per_itr2, performance,
//				(int) (performance.iterations - performance.itertaions_old_before_next_repeat), false);
//
//	}
//
//	return S_l;
//}

std::vector<std::tuple<int, double, double>> UniformSampling(const std::vector<std::tuple<int, double, double>> &S,
		double epsilon, double delta, double mu_s, int k, std::uniform_real_distribution<> &dist,
		std::mt19937 &rand_gen, Performance &performance, double scale_factor,
		std::vector<std::tuple<int, double, double>> &ground_truth, long long *total_pulls_ME) {

	int ssize = S.size();
	long long num_samples_per_arm;
//	if (config.adjust_scale) {
//		num_samples_per_arm = total_pulls_ME[0] / S.size() / scale_factor;
//	} else {
		num_samples_per_arm = (96.0 / epsilon / epsilon) * (1.0 / mu_s) * log(4.0 * ssize / delta) / scale_factor;
//	}
	std::vector<std::tuple<int, double, double>> S_empirical;

//	std::cerr << "US pulling " << num_samples_per_arm << " times for each of " << S.size() << " arms, total: "
//			<< num_samples_per_arm * S.size() << std::endl;
	for (int i = 0; i < S.size(); ++i) {
		int arm_id = std::get < 0 > (S[i]);
		double real_avg = std::get < 1 > (S[i]);
		double emprical_avg = getEmpiricalAvgOfArm(real_avg, num_samples_per_arm, dist, rand_gen);
		S_empirical.push_back(std::make_tuple(arm_id, real_avg, emprical_avg));
	}

	performance.total_pulls += num_samples_per_arm * S.size();
	performance.iterations++;
	performance.cost_per_itertaion[(int) (performance.iterations - performance.itertaions_old_before_next_repeat)] +=
			num_samples_per_arm * S.size();

	//sort arms in desc order based on emperical averages
	std::sort(S_empirical.begin(), S_empirical.end(), sortDescbyThird);

	std::vector<std::tuple<int, double, double>> result;
	result.assign(S_empirical.begin(), S_empirical.begin() + k); // copy [0, next_size) to S_l for next iteration.
	performance.candidate_size_per_iteration[(int) (performance.iterations
			- performance.itertaions_old_before_next_repeat)] += result.size();

	std::vector<std::tuple<int, double, double>> S_topk_per_itr(result.begin(), result.end());
	evaluateAccuracyPerItr(ground_truth, S_topk_per_itr, performance,
			(int) (performance.iterations - performance.itertaions_old_before_next_repeat), true);
	return result;
}

//std::vector<std::tuple<int, double, double>> MEAS_algo(const std::vector<std::tuple<int, double, double>> &arms,
//		double epsilon, double delta, int k, std::uniform_real_distribution<> &dist, std::mt19937 &rand_gen,
//		Performance &performance, double scale_factor_ME, double scale_factor_US,
//		std::vector<std::tuple<int, double, double>> &ground_truth, int repeati) {
//	double mu = 1.0;
//	if (config.arm_dist == UNIFORM) {
//		mu = 1.0;
//	} else if (config.arm_dist == SEGMENT) {
//		mu = 0.5;
//	} else if (config.arm_dist == NORMAL) {
//		mu = 0.997539;
//	}
//	std::cerr << "mu " << mu << std::endl;
//
//	double scale_factor_adjusted_ME = 1;
//
//	if (config.adjust_scale) {
//		long long num_pulls_1st_itr_ME = computeNumPulls_1stItr_MEAS_without_scale(epsilon, delta, mu, k, arms.size());
//		std::cerr << "scale_factor_ME for TopkDE = " << config.scale_factor_ME << std::endl;
//		long long num_pulls_1st_itr_TopkDE = computeNumPulls_1stItr_TopkDE(delta, epsilon, arms.size(), k,
//				config.scale_factor_ME);
//
//		scale_factor_adjusted_ME = (double) num_pulls_1st_itr_ME / (double) num_pulls_1st_itr_TopkDE;
//		std::cerr << "adjusted scale factor ME: " << scale_factor_adjusted_ME << std::endl;
//	}
//
//	std::cerr << "final scale factors: " << scale_factor_adjusted_ME << " " << scale_factor_US << std::endl;
//	if (repeati == 0)
//		std::cout << "final scale factors: " << scale_factor_adjusted_ME << " " << scale_factor_US << std::endl;
//
//	// scale_factor_ME is not in use actually
//	long long *total_pulls_ME = new long long [1];
//	total_pulls_ME[0] = 0;
//	std::vector<std::tuple<int, double, double>> S_ME = MedianElimination(arms, epsilon, delta, mu, k, dist, rand_gen,
//			performance, scale_factor_adjusted_ME, ground_truth, total_pulls_ME);
//
//	double mu_s = mu * (1 - epsilon / 2.0);
//	std::vector<std::tuple<int, double, double>> S_US = UniformSampling(S_ME, epsilon, delta, mu_s, k, dist, rand_gen,
//			performance, scale_factor_US, ground_truth, total_pulls_ME);
//	return S_US;
//}

std::vector<std::tuple<int, double, double>> GetGroundTruth(const std::vector<std::tuple<int, double, double>> &arms,
		int k) {
	std::vector<std::tuple<int, double, double>> S(arms.begin(), arms.end());
	//sort arms in desc order based on emperical averages
	std::sort(S.begin(), S.end(), sortDescbySecond);

	std::vector<std::tuple<int, double, double>> result;
	result.assign(S.begin(), S.begin() + k); // copy [0, next_size) to S_l for next iteration.
	return result;
}

#endif /* UTILMAB_H_ */
