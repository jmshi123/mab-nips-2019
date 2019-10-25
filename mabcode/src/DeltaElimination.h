/*
 * DeltaElimination.h
 *
 *  Created on: Dec 29, 2018
 *      Author: jmshinju
 */

#ifndef DELTAELIMINATION_H_
#define DELTAELIMINATION_H_

#include <set>

int countNonzeroEmpiricalValue(std::vector<std::tuple<int, double, double>> Ka) {
	int count = 0;
	for (auto & element : Ka) {
		if (std::get < 2 > (element) > 0) count++;
	}
	return count;
}

std::vector<std::tuple<int, double, double>> getNonzeroEmpiricalValueElements(
		std::vector<std::tuple<int, double, double>> Ka) {
	std::vector<std::tuple<int, double, double>> result;
	for (auto & element : Ka) {
		if (std::get < 2 > (element) > 0) {
			result.push_back(
					std::make_tuple(std::get < 0 > (element), std::get < 1 > (element), std::get < 2 > (element)));
		}
	}
	return result;
}

std::tuple<int, double, double> DeltaEliminateBestArm(const std::vector<std::tuple<int, double, double>> &arms,
		double epsilon, double delta, double Q, std::uniform_real_distribution<> &dist, std::mt19937 &rand_gen,
		Performance &performance, double scale_factor_1, double scale_factor_2,
		std::vector<std::tuple<int, double, double>> &ground_truth, int repeati) {

	std::vector<std::tuple<int, double, double>> S_l(arms.begin(), arms.end());
	std::vector<std::tuple<int, double, double>> Ka(arms.begin(), arms.end());
	double ka_max_empirical_avg = 0.0;
	double ka_max_arm_id = -1;

	double beta_l = 1;
	double delta_l = delta / 4.0;

	double l = 1;

	while (S_l.size() > 0) {
		long long num_samples_per_arm_origin = beta_l * Q * log(1.0 / delta_l);
		long long num_samples_per_arm = num_samples_per_arm_origin / scale_factor_1;
		std::vector<std::tuple<int, double, double>> S_empirical;
//		std::cerr << "  DE------" << l << "-------\n";
//		std::cerr << "pulling " << num_samples_per_arm << " times for each of " << S_l.size() << " arms, total: "
//				<< num_samples_per_arm * S_l.size() << std::endl;
		long long total_pulls_step1 = num_samples_per_arm * S_l.size();

		for (int i = 0; i < S_l.size(); ++i) {
			int arm_id = std::get < 0 > (S_l[i]);
			double real_avg = std::get < 1 > (S_l[i]);
			double emprical_avg = getEmpiricalAvgOfArm(real_avg, num_samples_per_arm, dist, rand_gen);
			S_empirical.push_back(std::make_tuple(arm_id, real_avg, emprical_avg));
		}
		performance.total_pulls += num_samples_per_arm * S_l.size();
		performance.iterations++;
		performance.cost_per_itertaion[(int) (performance.iterations - performance.itertaions_old_before_next_repeat)] +=
				num_samples_per_arm * S_l.size();

		performance.candidate_size_per_iteration[(int) (performance.iterations
				- performance.itertaions_old_before_next_repeat)] += S_l.size() + countNonzeroEmpiricalValue(Ka);

		std::vector<std::tuple<int, double, double>> S_topk_per_itr(S_empirical.begin(), S_empirical.end());
		std::vector<std::tuple<int, double, double>> ka_tmp = getNonzeroEmpiricalValueElements(Ka);
		S_topk_per_itr.insert(S_topk_per_itr.end(), ka_tmp.begin(), ka_tmp.end());
		std::sort(S_topk_per_itr.begin(), S_topk_per_itr.end(), sortDescbyThird);
		std::vector<std::tuple<int, double, double>> S_topk_per_itr2(S_topk_per_itr.begin(),
				S_topk_per_itr.begin() + 1);
		evaluateAccuracyPerItr(ground_truth, S_topk_per_itr2, performance,
				(int) (performance.iterations - performance.itertaions_old_before_next_repeat), false);

		std::sort(S_empirical.begin(), S_empirical.end(), sortDescbyThird);
		int top_rank = ceil(pow(delta_l, beta_l) * ((double) S_empirical.size()) / 2.0);
		int arm_j = floor(dist(rand_gen) * top_rank); // random choose arm_j from [0, top_rank)

		int arm_id = std::get < 0 > (S_empirical[arm_j]);
		double real_avg_j = std::get < 1 > (S_empirical[arm_j]);
		double empirical_avg_j_old = std::get < 2 > (S_empirical[arm_j]);

//		if (config.adjust_scale) {
//			num_samples_per_arm = total_pulls_step1 / scale_factor_2;
//		} else {
			num_samples_per_arm = num_samples_per_arm_origin / scale_factor_2;
//		}
		double emprical_avg_j = getEmpiricalAvgOfArm(real_avg_j, num_samples_per_arm, dist, rand_gen);
		std::get < 2 > (S_empirical[arm_j]) = emprical_avg_j;
//		std::cerr << "randomly choose arm " << arm_id << "(idx " << arm_j << ") < " << top_rank
//				<< ", update its empirical avg " << empirical_avg_j_old << " to " << std::get < 2
//				> (S_empirical[arm_j]) << "\n  by pulling it " << num_samples_per_arm << " times\n";
		performance.total_pulls += num_samples_per_arm;
		performance.iterations++;
		performance.cost_per_itertaion[(int) (performance.iterations - performance.itertaions_old_before_next_repeat)] +=
				num_samples_per_arm;

		S_topk_per_itr.clear();
		S_topk_per_itr.insert(S_topk_per_itr.end(), S_empirical.begin(), S_empirical.end());
		ka_tmp = getNonzeroEmpiricalValueElements(Ka);
		S_topk_per_itr.insert(S_topk_per_itr.end(), ka_tmp.begin(), ka_tmp.end());
		std::sort(S_topk_per_itr.begin(), S_topk_per_itr.end(), sortDescbyThird);
		S_topk_per_itr2.clear();
		S_topk_per_itr2.insert(S_topk_per_itr2.end(), S_topk_per_itr.begin(), S_topk_per_itr.begin() + 1);

		std::get < 2 > (Ka[arm_id]) = std::get < 2 > (S_empirical[arm_j]);
		double ka_chosen_empirical_avg = std::get < 2 > (Ka[arm_id]);
		if (ka_chosen_empirical_avg > ka_max_empirical_avg) {
			ka_max_empirical_avg = ka_chosen_empirical_avg; // always get the max empirical avg in Ka vector.
			ka_max_arm_id = arm_id;
		}

		//update S_l for next iteration
		double size_S_l_old = S_l.size();
		S_l.clear();
		double threshold = ka_max_empirical_avg + 0.5 * epsilon;
		for (auto &element : S_empirical) {
			if (std::get < 2 > (element) >= threshold) {
				S_l.push_back(std::make_tuple(std::get < 0 > (element), std::get < 1 > (element), 0.0));
			}
		}
//		std::cerr << "next size " << S_l.size() << ", old size " << size_S_l_old << std::endl;
		if (S_l.size() == 0) {
			performance.candidate_size_per_iteration[(int) (performance.iterations
					- performance.itertaions_old_before_next_repeat)] += 1; // after last iteration, there is only one candidate, i.e., result
			evaluateAccuracyPerItr(ground_truth, S_topk_per_itr2, performance,
					(int) (performance.iterations - performance.itertaions_old_before_next_repeat), true);
		} else {
			performance.candidate_size_per_iteration[(int) (performance.iterations
					- performance.itertaions_old_before_next_repeat)] += S_l.size() + countNonzeroEmpiricalValue(Ka);
			evaluateAccuracyPerItr(ground_truth, S_topk_per_itr2, performance,
					(int) (performance.iterations - performance.itertaions_old_before_next_repeat), false);
		}

		//update parameters
//		if (config.adjust_scale) {
//			beta_l *= size_S_l_old / (double) S_l.size();
//		} else {
			if (S_l.size() <= 2 * delta * size_S_l_old) {
				beta_l *= size_S_l_old / 2.0 / (double) S_l.size();
			} else {
				beta_l *= size_S_l_old / (double) S_l.size();
			}
//		}

		delta_l = delta / 2.0 / pow(2.0, l);
		l = l + 1;
	}

	return Ka[ka_max_arm_id];
}

std::vector<std::tuple<int, double, double>> TopkDeltaEliminate(
		const std::vector<std::tuple<int, double, double>> &arms, double epsilon, double delta, double Q, int k,
		std::uniform_real_distribution<> &dist, std::mt19937 &rand_gen, Performance &performance, double scale_factor_1,
		double scale_factor_2, std::vector<std::tuple<int, double, double>> &ground_truth, int repeati) {

	std::vector<std::tuple<int, double, double>> S_l(arms.begin(), arms.end());
	std::vector<std::tuple<int, double, double>> Ka(arms.begin(), arms.end());

	double beta_l = 1;
	double delta_l = delta / 4.0;

	double l = 1;

//	std::cerr << "k: " << k << ", Q: " << Q << ", delta " << delta << ", epsilon " << epsilon << ", scale_factor "
//			<< scale_factor_1 << " and " << scale_factor_2 << std::endl;

//	std::cerr << "final scale factors: " << scale_factor_1 << " " << scale_factor_2 << std::endl;
//	if (repeati == 0) std::cout << "final scale factors: " << scale_factor_1 << " " << scale_factor_2 << std::endl;

	while (S_l.size() > 0) {
//		std::cerr << "beta_l: " << beta_l << ", delta_l " << delta_l << std::endl;
		long long num_samples_per_arm_origin = beta_l * Q * log(((double) k) / delta_l);
		long long num_samples_per_arm = num_samples_per_arm_origin / scale_factor_1;
		std::vector<std::tuple<int, double, double>> S_empirical;
//		std::cerr << "------" << l << "-------\n";
//		std::cerr << "1 pulling " << num_samples_per_arm << " times for each of " << S_l.size() << " arms, total: "
//				<< num_samples_per_arm * S_l.size() << std::endl;
		long long total_pulls_step1 = num_samples_per_arm * S_l.size();

		for (int i = 0; i < S_l.size(); ++i) {
			int arm_id = std::get < 0 > (S_l[i]);
			double real_avg = std::get < 1 > (S_l[i]);
			double emprical_avg = getEmpiricalAvgOfArm(real_avg, num_samples_per_arm, dist, rand_gen);
			S_empirical.push_back(std::make_tuple(arm_id, real_avg, emprical_avg));
		}
		performance.total_pulls += num_samples_per_arm * S_l.size();
		performance.iterations++;
		performance.cost_per_itertaion[(int) (performance.iterations - performance.itertaions_old_before_next_repeat)] +=
				num_samples_per_arm * S_l.size();
		performance.candidate_size_per_iteration[(int) (performance.iterations
				- performance.itertaions_old_before_next_repeat)] += S_l.size() + countNonzeroEmpiricalValue(Ka);

		std::vector<std::tuple<int, double, double>> S_topk_per_itr(S_empirical.begin(), S_empirical.end());
		std::vector<std::tuple<int, double, double>> ka_tmp = getNonzeroEmpiricalValueElements(Ka);
		S_topk_per_itr.insert(S_topk_per_itr.end(), ka_tmp.begin(), ka_tmp.end());
		std::sort(S_topk_per_itr.begin(), S_topk_per_itr.end(), sortDescbyThird);
		std::vector<std::tuple<int, double, double>> S_topk_per_itr2(S_topk_per_itr.begin(),
				S_topk_per_itr.begin() + k);
		evaluateAccuracyPerItr(ground_truth, S_topk_per_itr2, performance,
				(int) (performance.iterations - performance.itertaions_old_before_next_repeat), false);

		int k_min = k < S_l.size() ? k : S_l.size();
		std::sort(S_empirical.begin(), S_empirical.end(), sortDescbyThird);

		int top_rank = ceil(pow(delta_l / (double) k, beta_l) * ((double) S_empirical.size()) / 2.0) + k_min - 1;
		std::vector<std::tuple<int, double, double>> S_k_min;
		std::set<int> arm_ids_S_k_min;
		for (int i = 0; i < k_min; ++i) {
			int arm_j = floor(dist(rand_gen) * top_rank); // random choose arm_j from [0, top_rank)

			int arm_id = std::get < 0 > (S_empirical[arm_j]);
			double real_avg_j = std::get < 1 > (S_empirical[arm_j]);
			double empirical_avg_j_old = std::get < 2 > (S_empirical[arm_j]);
			S_k_min.push_back(std::make_tuple(arm_id, real_avg_j, empirical_avg_j_old));
			arm_ids_S_k_min.insert(arm_id);

			S_empirical.erase(S_empirical.begin() + arm_j); //remove selected j
			top_rank--;
		}

			num_samples_per_arm = num_samples_per_arm_origin / scale_factor_2;
		std::vector<std::tuple<int, double, double>> S_k_min_empirical;
//		std::cerr << "2 pulling " << num_samples_per_arm << " times for each of " << S_k_min.size() << " arms, total: "
//				<< num_samples_per_arm * S_k_min.size() << std::endl;
		for (int i = 0; i < S_k_min.size(); ++i) {
			int arm_id = std::get < 0 > (S_k_min[i]);
			double real_avg = std::get < 1 > (S_k_min[i]);
			double emprical_avg = getEmpiricalAvgOfArm(real_avg, num_samples_per_arm, dist, rand_gen);
			S_k_min_empirical.push_back(std::make_tuple(arm_id, real_avg, emprical_avg));
		}
		performance.total_pulls += num_samples_per_arm * S_k_min.size();
		performance.iterations++;
		performance.cost_per_itertaion[(int) (performance.iterations - performance.itertaions_old_before_next_repeat)] +=
				num_samples_per_arm * S_k_min.size();

		//update Ka's empirical avg of the arms in S_k_min_empirical
//		std::sort(S_k_min_empirical.begin(), S_k_min_empirical.end(), sortDescbyThird);
		for (int i = 0; i < S_k_min_empirical.size(); ++i) {
			int arm_id = std::get < 0 > (S_k_min_empirical[i]);
			double emprical_avg = std::get < 2 > (S_k_min_empirical[i]);
			std::get < 2 > (Ka[arm_id]) = emprical_avg;
//			std::cout << arm_id << ": " << std::get < 2 > (Ka[arm_id]) << ",";
		}
//		std::cout << std::endl;

		S_topk_per_itr.clear();
		S_topk_per_itr.insert(S_topk_per_itr.end(), S_empirical.begin(), S_empirical.end());
		ka_tmp = getNonzeroEmpiricalValueElements(Ka);
		S_topk_per_itr.insert(S_topk_per_itr.end(), ka_tmp.begin(), ka_tmp.end());
		std::sort(S_topk_per_itr.begin(), S_topk_per_itr.end(), sortDescbyThird);
		S_topk_per_itr2.clear();
		S_topk_per_itr2.insert(S_topk_per_itr2.end(), S_topk_per_itr.begin(), S_topk_per_itr.begin() + k);

		//find the k-th empirical value in Ka using a temp Ka
		std::vector<std::tuple<int, double, double>> Ka_tmp(Ka.begin(), Ka.end());
		std::sort(Ka_tmp.begin(), Ka_tmp.end(), sortDescbyThird);


		std::tuple<int, double, double> kth_arm_in_Ka = Ka_tmp[k - 1];
		double kth_empirical_value = std::get < 2 > (kth_arm_in_Ka);
//		std::cout << "kth_empirical_value " << kth_empirical_value << std::endl;

		//update S_l for next iteration
		double size_S_l_old = S_l.size();
		S_l.clear();
		double threshold = kth_empirical_value + 0.5 * epsilon;
		for (auto &element : S_empirical) {
			if (std::get < 2 > (element) >= threshold
					&& arm_ids_S_k_min.find(std::get < 0 > (element)) == arm_ids_S_k_min.end()) {
				S_l.push_back(std::make_tuple(std::get < 0 > (element), std::get < 1 > (element), 0.0));
			}
		}
//		std::cerr << "next size " << S_l.size() << ", old size " << size_S_l_old << std::endl;
		if (S_l.size() == 0) {
			performance.candidate_size_per_iteration[(int) (performance.iterations
					- performance.itertaions_old_before_next_repeat)] += k; // after last iteration, there is only one candidate, i.e., result
			evaluateAccuracyPerItr(ground_truth, S_topk_per_itr2, performance,
					(int) (performance.iterations - performance.itertaions_old_before_next_repeat), true);
		} else {
			performance.candidate_size_per_iteration[(int) (performance.iterations
					- performance.itertaions_old_before_next_repeat)] += S_l.size() + countNonzeroEmpiricalValue(Ka);
			evaluateAccuracyPerItr(ground_truth, S_topk_per_itr2, performance,
					(int) (performance.iterations - performance.itertaions_old_before_next_repeat), false);
		}

		//update parameter


			if (S_l.size() <= 2 * delta * size_S_l_old / (double) k) {
				beta_l *= size_S_l_old / 2.0 / (double) S_l.size();
			} else {
				beta_l *= size_S_l_old / (double) S_l.size();
			}
//		}

		delta_l = delta / 2.0 / pow(2.0, l);
		l = l + 1;
	}

	//now sort Ka and get top k.
	std::sort(Ka.begin(), Ka.end(), sortDescbyThird);
	std::vector<std::tuple<int, double, double>> result;
	result.assign(Ka.begin(), Ka.begin() + k);
	return result;
}

inline double computeBeta_l(int num_arms, int R, double delta, int k) {
	//always assuem num_arms >= 1/delta
	double ilogn = num_arms;
	for (int i = 0; i < R; ++i) {
		ilogn = log(ilogn);
	}

	double beta_l = 1;
	double ratio = ilogn / log((double) k / delta);
	if (ratio >= 1) {
		beta_l += ratio;
	}

	return beta_l;
}

std::vector<std::tuple<int, double, double>> UniformSampling_DE(std::vector<std::tuple<int, double, double>> &S_prime,
		const std::vector<std::tuple<int, double, double>> &S_R, double Q, int R, double epsilon, double delta,
		double beta_R, std::uniform_real_distribution<> &dist, std::mt19937 &rand_gen, Performance &performance,
		double scale_factor_ME, double scale_factor, int k, std::vector<std::tuple<int, double, double>> &ground_truth,
		long long *total_pulls_ME) {

//	std::cerr << "US1 Q: " << Q << ", beta_R " << beta_R << ", R " << R << ", delta " << delta << ", epsilon "
//			<< epsilon << std::endl;

	std::vector<std::tuple<int, double, double>> S_R_empirical;
	long long num_samples_per_arm = Q * beta_R * log(8.0 * k * pow(2.0, R) / delta) / scale_factor_ME;
//	std::cerr << "US1 pulling " << num_samples_per_arm << " times for each of " << S_R.size() << " arms in S_R, total: "
//			<< num_samples_per_arm * S_R.size() << std::endl;
	long long total_pulls_US1 = num_samples_per_arm * S_R.size();
	if (num_samples_per_arm > 0) {
		for (int i = 0; i < S_R.size(); ++i) {
			int arm_id = std::get < 0 > (S_R[i]);
			double real_avg = std::get < 1 > (S_R[i]);
			double emprical_avg = getEmpiricalAvgOfArm(real_avg, num_samples_per_arm, dist, rand_gen);
			S_R_empirical.push_back(std::make_tuple(arm_id, real_avg, emprical_avg));
		}
		std::sort(S_R_empirical.begin(), S_R_empirical.end(), sortDescbyThird);
		performance.total_pulls += num_samples_per_arm * S_R.size();
		performance.cost_per_itertaion[(int) (performance.iterations + 1 - performance.itertaions_old_before_next_repeat)] +=
				num_samples_per_arm * S_R.size();
		total_pulls_ME[0] += num_samples_per_arm * S_R.size();
	}

	int k_min = k < S_R_empirical.size() ? k : S_R_empirical.size();

	std::vector<std::tuple<int, double, double>> S_prime_empirical;

		num_samples_per_arm = Q * log(2.0 * S_prime.size() / delta) / scale_factor;
//	}
//	std::cerr << "US2 Q: " << Q << ", S_prime " << S_prime.size() << ", delta " << delta << std::endl;
//	std::cerr << "US2 pulling " << num_samples_per_arm << " times for each of " << S_prime.size()
//			<< " arms in S_prime, total: " << num_samples_per_arm * S_prime.size() << std::endl;
	for (int i = 0; i < S_prime.size(); ++i) {
		int arm_id = std::get < 0 > (S_prime[i]);
		double real_avg = std::get < 1 > (S_prime[i]);
		double emprical_avg = getEmpiricalAvgOfArm(real_avg, num_samples_per_arm, dist, rand_gen);
		S_prime_empirical.push_back(std::make_tuple(arm_id, real_avg, emprical_avg));
	}
	performance.total_pulls += num_samples_per_arm * S_prime.size();
	performance.iterations++;
	performance.cost_per_itertaion[(int) (performance.iterations - performance.itertaions_old_before_next_repeat)] +=
			num_samples_per_arm * S_prime.size();
	performance.candidate_size_per_iteration[(int) (performance.iterations
			- performance.itertaions_old_before_next_repeat)] += k;

	S_prime_empirical.insert(S_prime_empirical.end(), S_R_empirical.begin(), S_R_empirical.begin() + k_min);

	//sort arms in desc order based on emperical averages
	std::sort(S_prime_empirical.begin(), S_prime_empirical.end(), sortDescbyThird);

	std::vector<std::tuple<int, double, double>> result;
	result.assign(S_prime_empirical.begin(), S_prime_empirical.begin() + k); // copy [0, next_size) to S_l for next iteration.

	std::vector<std::tuple<int, double, double>> S_topk_per_itr(result.begin(), result.end());
	evaluateAccuracyPerItr(ground_truth, S_topk_per_itr, performance,
			(int) (performance.iterations - performance.itertaions_old_before_next_repeat), true);

	return result;

}

long long computeNumPulls_1stItr_BestDE_fix_without_scale(double epsilon, double delta, double num_arms, int R) {
	double c = 57;
	double Q = c / epsilon / epsilon;

	double delta_l = delta / 4.0;
	double beta_l = computeBeta_l(num_arms, R, delta, 1);

	long long num_samples_per_arm = beta_l * Q * log((double) 1 / delta_l);
	return num_arms * num_samples_per_arm;
}

std::tuple<int, double, double> DeltaEliminateBestArm_Deterministic(
		const std::vector<std::tuple<int, double, double>> &arms, double epsilon, double delta, double Q, int R,
		std::uniform_real_distribution<> &dist, std::mt19937 &rand_gen, Performance &performance,
		double scale_factor_DE, double scale_factor_US, std::vector<std::tuple<int, double, double>> &ground_truth,
		int repeati) {

	long long *total_pulls_ME = new long long[1];
	total_pulls_ME[0] = 0;

	double scale_factor_adjusted_ME = 1;


	int l = 1;
	double delta_l = delta / 4.0; //16.0;
	double beta_l = computeBeta_l(arms.size(), R, delta, 1);

	std::vector<std::tuple<int, double, double>> S_l(arms.begin(), arms.end());
	std::vector<std::tuple<int, double, double>> S_prime;

//	std::cerr << "DE Q: " << Q << ", R " << R << ", delta " << delta << ", epsilon " << epsilon << std::endl;
//	std::cerr << " beta_" << l << ": " << beta_l << " delta_" << l << ": " << delta_l << std::endl;

	for (l = 1; l <= R - 1; ++l) {
		long long num_samples_per_arm = beta_l * Q * log(1.0 / delta_l) / scale_factor_adjusted_ME;
		std::vector<std::tuple<int, double, double>> S_empirical;
//		std::cerr << "------" << l << "-------\n";
//		std::cerr << "pulling " << num_samples_per_arm << " times for each of " << S_l.size() << " arms, total: "
//				<< num_samples_per_arm * S_l.size() << std::endl;
		for (int i = 0; i < S_l.size(); ++i) {
			int arm_id = std::get < 0 > (S_l[i]);
			double real_avg = std::get < 1 > (S_l[i]);
			double emprical_avg = getEmpiricalAvgOfArm(real_avg, num_samples_per_arm, dist, rand_gen);
			S_empirical.push_back(std::make_tuple(arm_id, real_avg, emprical_avg));
		}
		performance.total_pulls += num_samples_per_arm * S_l.size();
		performance.iterations++;
		performance.cost_per_itertaion[(int) (performance.iterations - performance.itertaions_old_before_next_repeat)] +=
				num_samples_per_arm * S_l.size();
		total_pulls_ME[0] += num_samples_per_arm * S_l.size();

		std::vector<std::tuple<int, double, double>> S_topk_per_itr(S_empirical.begin(), S_empirical.end());
		std::sort(S_topk_per_itr.begin(), S_topk_per_itr.end(), sortDescbyThird);
		std::vector<std::tuple<int, double, double>> S_topk_per_itr2(S_topk_per_itr.begin(),
				S_topk_per_itr.begin() + 1);
		evaluateAccuracyPerItr(ground_truth, S_topk_per_itr2, performance,
				(int) (performance.iterations - performance.itertaions_old_before_next_repeat), false);

		std::sort(S_empirical.begin(), S_empirical.end(), sortDescbyThird);

		int top_rank = ceil(pow(delta_l, beta_l) * ((double) S_empirical.size()) / 2.0);
		int arm_j = floor(dist(rand_gen) * top_rank); // random choose arm_j from [0, top_rank)
		int arm_id_j = std::get < 0 > (S_empirical[arm_j]);
		double real_avg_j = std::get < 1 > (S_empirical[arm_j]);
		double empirical_avg_j = std::get < 2 > (S_empirical[arm_j]);

		S_prime.push_back(std::make_tuple(arm_id_j, real_avg_j, empirical_avg_j));
		S_empirical.erase(S_empirical.begin() + arm_j); //remove selected j

		int size_S_l_old = S_l.size();
		int next_size = ceil(2.0 * pow(delta_l, beta_l) * S_l.size()) - 1;
		if (scale_factor_adjusted_ME > 1) {
//			std::cout << delta_l << " " << beta_l << " " << S_l.size() << std::endl;
			int next_size = ceil(100.0 * pow(delta_l, beta_l) * S_l.size()) - 1;
		}
		S_l.clear();
		S_l.assign(S_empirical.begin(), S_empirical.begin() + next_size); // copy [0, next_size) to S_l for next iteration.
//		std::cerr << "next_size " << S_l.size() << std::endl;
		performance.candidate_size_per_iteration[(int) (performance.iterations
				- performance.itertaions_old_before_next_repeat)] += S_l.size() + S_prime.size();
		//update parameters
		beta_l *= size_S_l_old / 2.0 / (double) S_l.size();

//		std::cerr << "beta_" << l + 1 << ": " << beta_l << std::endl;

		delta_l = delta / 8.0 / pow(2.0, l);
	}

	std::vector<std::tuple<int, double, double>> top1 = UniformSampling_DE(S_prime, S_l, Q, R, epsilon, delta, beta_l,
			dist, rand_gen, performance, scale_factor_adjusted_ME, scale_factor_US, 1, ground_truth, total_pulls_ME);

	delete total_pulls_ME;
	return top1[0];
}

long long computeNumPulls_1stItr_TopkDE_fix_without_scale(double epsilon, double delta, double num_arms, int R, int k) {
	double c = 57;
	double Q = c / epsilon / epsilon;

	double delta_l = delta / 4.0;
	double beta_l = computeBeta_l(num_arms, R, delta, k);

	long long num_samples_per_arm = beta_l * Q * log((double) k / delta_l);
	return num_arms * num_samples_per_arm;
}

std::vector<std::tuple<int, double, double>> TopKDeltaEliminate_Deterministic(
		const std::vector<std::tuple<int, double, double>> &arms, double epsilon, double delta, double Q, int R, int k,
		std::uniform_real_distribution<> &dist, std::mt19937 &rand_gen, Performance &performance,
		double scale_factor_DE, double scale_factor_US, std::vector<std::tuple<int, double, double>> &ground_truth,
		int repeati) {

	long long *total_pulls_ME = new long long[1];
	total_pulls_ME[0] = 0;

	double scale_factor_adjusted_ME = 1;


	int l = 1;
	double delta_l = delta / 4.0;
	double beta_l = computeBeta_l(arms.size(), R, delta, k);

	std::vector<std::tuple<int, double, double>> S_l(arms.begin(), arms.end());
	std::vector<std::tuple<int, double, double>> S_prime;


	for (l = 1; l <= R - 1; ++l) {
		long long num_samples_per_arm = beta_l * Q * log((double) k / delta_l) / scale_factor_adjusted_ME;
		std::vector<std::tuple<int, double, double>> S_empirical;
//		std::cerr << "------" << l << "-------\n";
//		std::cerr << "pulling " << num_samples_per_arm << " times for each of " << S_l.size() << " arms, total: "
//				<< num_samples_per_arm * S_l.size() << std::endl;
		for (int i = 0; i < S_l.size(); ++i) {
			int arm_id = std::get < 0 > (S_l[i]);
			double real_avg = std::get < 1 > (S_l[i]);
			double emprical_avg = getEmpiricalAvgOfArm(real_avg, num_samples_per_arm, dist, rand_gen);
			S_empirical.push_back(std::make_tuple(arm_id, real_avg, emprical_avg));
		}
		performance.total_pulls += num_samples_per_arm * S_l.size();
		performance.iterations++;
		performance.cost_per_itertaion[(int) (performance.iterations - performance.itertaions_old_before_next_repeat)] +=
				num_samples_per_arm * S_l.size();
		total_pulls_ME[0] += num_samples_per_arm * S_l.size();;

		int k_min = k < S_l.size() ? k : S_l.size();

		std::sort(S_empirical.begin(), S_empirical.end(), sortDescbyThird);

		int top_rank = ceil(pow(delta_l / ((double) k), beta_l) * ((double) S_empirical.size()) / 2.0) + k_min - 1;

		//randomly choose k_min arms from top_rank of S_empirical, i.e., S_l, and put them into S_prime
		for (int i = 0; i < k_min; ++i) {
			int arm_j = floor(dist(rand_gen) * top_rank); // random choose arm_j from [0, top_rank)
			int arm_id_j = std::get < 0 > (S_empirical[arm_j]);
			double real_avg_j = std::get < 1 > (S_empirical[arm_j]);
			double empirical_avg_j = std::get < 2 > (S_empirical[arm_j]);

			S_prime.push_back(std::make_tuple(arm_id_j, real_avg_j, empirical_avg_j));
			S_empirical.erase(S_empirical.begin() + arm_j); //remove selected j

			top_rank--; //shrink since erased
		}

		int size_S_l_old = S_l.size();
		int next_size = ceil(2.0 * pow(delta_l / ((double) k), beta_l) * S_l.size()) + k_min - 1 - k_min;

		S_l.clear();
		if (next_size != 0) {
			S_l.assign(S_empirical.begin(), S_empirical.begin() + next_size); // copy [0, next_size) to S_l for next iteration.
		}
//		std::cerr << "next_size " << S_l.size() << std::endl;
		performance.candidate_size_per_iteration[(int) (performance.iterations
				- performance.itertaions_old_before_next_repeat)] += S_l.size() + S_prime.size();

		std::vector<std::tuple<int, double, double>> S_topk_per_itr(S_l.begin(), S_l.end());
		S_topk_per_itr.insert(S_topk_per_itr.end(), S_prime.begin(), S_prime.end());
		std::sort(S_topk_per_itr.begin(), S_topk_per_itr.end(), sortDescbyThird);
		if (S_prime.size() + S_l.size() <= k) {
			evaluateAccuracyPerItr(ground_truth, S_topk_per_itr, performance,
					(int) (performance.iterations - performance.itertaions_old_before_next_repeat), true);
		} else {
			evaluateAccuracyPerItr(ground_truth, S_topk_per_itr, performance,
					(int) (performance.iterations - performance.itertaions_old_before_next_repeat), false);
		}

		//update parameters
		if (next_size == 0) {
			beta_l = 0;
		} else {
			beta_l *= size_S_l_old / 2.0 / (double) S_l.size();
		}

//		std::cerr << "beta_" << l + 1 << ": " << beta_l << " actual next size " << S_l.size() << std::endl;

		delta_l = delta / 8.0 / pow(2.0, l);

		if (next_size == 0 || S_l.size() + S_prime.size() <= k) {
			break;
		}
	}

	if (S_prime.size() + S_l.size() <= k) { //candidate set <=k, just return
//		std::cerr << S_prime.size() << "+" << S_l.size() << "<=" << k << std::endl;
		S_prime.insert(S_prime.end(), S_l.begin(), S_l.end());
		std::sort(S_prime.begin(), S_prime.end(), sortDescbyThird);
		return S_prime;
	} else {
		return UniformSampling_DE(S_prime, S_l, Q, R, epsilon, delta, beta_l, dist, rand_gen, performance,
				scale_factor_adjusted_ME, scale_factor_US, k, ground_truth, total_pulls_ME);
	}
}

#endif /* DELTAELIMINATION_H_ */
