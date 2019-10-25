/*
 * ExponentialGap.h
 *
 *  Created on: May 17, 2019
 *      Author: jmshinju
 */

#ifndef EXPONENTIALGAP_H_
#define EXPONENTIALGAP_H_
#include "Performance.h"


std::tuple<int, double, double> ExponentialGapDeltaE(const std::vector<std::tuple<int, double, double>> &arms,
		double delta, std::uniform_real_distribution<> &dist, std::mt19937 &rand_gen, Performance &performance,
		std::vector<std::tuple<int, double, double>> &ground_truth, int repeati) {

	std::vector<std::tuple<int, double, double>> S_l(arms.begin(), arms.end());
	int l = 1;

	while (S_l.size() > 1) {
		double epsilon_l = pow(2.0, l * (-1.0)) / 4.0;
		double delta_l = delta / (50.0 * l * l * l);

		long long num_samples_per_arm = (2.0 / epsilon_l / epsilon_l) * log(2.0 / delta_l);

		std::vector<std::tuple<int, double, double>> S_empirical;

//		std::cerr << "EG pulling " << num_samples_per_arm << " times for each of " << S_l.size() << " arms, total: "
//				<< num_samples_per_arm * S_l.size() << std::endl;
		for (int i = 0; i < S_l.size(); ++i) {
			int arm_id = std::get < 0 > (S_l[i]);
			double real_avg = std::get < 1 > (S_l[i]);
			double emprical_avg = getEmpiricalAvgOfArm(real_avg, num_samples_per_arm, dist, rand_gen);
			S_empirical.push_back(std::make_tuple(arm_id, real_avg, emprical_avg));
		}

		double Q = config.c / epsilon_l / epsilon_l;
		std::tuple<int, double, double> bestArmByDE =
				DeltaEliminateBestArm(S_l, epsilon_l / 2.0, delta_l, Q, dist, rand_gen, performance, config.scale_factor_ME, config.scale_factor_US, ground_truth, repeati);

		double empiricalRewardStar = std::get < 2 > (bestArmByDE);
		//sort arms in desc order based on emperical averages
//		std::sort(S_empirical.begin(), S_empirical.end(), sortDescbyThird);

		performance.total_pulls += num_samples_per_arm * S_l.size();
		performance.cost_per_itertaion[(int) (performance.iterations + 1 - performance.itertaions_old_before_next_repeat)] +=
				num_samples_per_arm * S_l.size();

		S_l.clear();
		for (int i = 0; i < S_empirical.size(); ++i) {
			int arm_id = std::get < 0 > (S_empirical[i]);
			double real_avg = std::get < 1 > (S_empirical[i]);
			double emprical_avg = std::get < 2 > (S_empirical[i]);
			if (emprical_avg >= empiricalRewardStar - epsilon_l) {
				S_l.push_back(std::make_tuple(arm_id, real_avg, emprical_avg));
			}
		}
		performance.iterations++;
		performance.candidate_size_per_iteration[(int) (performance.iterations
				- performance.itertaions_old_before_next_repeat)] += S_l.size();

		std::vector<std::tuple<int, double, double>> S_topk_per_itr(S_l.begin(), S_l.end());
		std::sort(S_topk_per_itr.begin(), S_topk_per_itr.end(), sortDescbyThird);
		std::vector<std::tuple<int, double, double>> S_topk_per_itr2(S_topk_per_itr.begin(), S_topk_per_itr.begin()
				+ 1);
		evaluateAccuracyPerItr(ground_truth, S_topk_per_itr2, performance, (int) (performance.iterations
				- performance.itertaions_old_before_next_repeat), false);
		l++;
	}
	return S_l[0];
}

#endif /* EXPONENTIALGAP_H_ */
