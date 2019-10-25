

#include <chrono>
#include <iostream>
#include <tuple>

#include "arms.h"
#include "config.h"
#include "utils.h"
#include "gen.h"
#include "UTILMAB.h"
#include "DeltaElimination.h"
#include "ExponentialGap.h"

using namespace std;

Config config;

void DisplayResult(std::vector<std::tuple<int, double, double> > result) {
	for (auto& kth : result) {
		std::cerr << tuple_to_string(kth) << ";";
	}
	std::cerr << std::endl;
}

bool evaluateAccuracy(std::vector<std::tuple<int, double, double> >& result_ground_truth,
		std::vector<std::tuple<int, double, double> >& result, Performance& performance) {
	bool accurate = true;
	for (int i = 0; i < config.k; ++i) {
		if (std::abs(std::get < 1 > (result_ground_truth[i]) - std::get < 1 > (result[i])) > config.epsilon) {
			accurate = false; //any top k with large error, the whole process fail.
			break;
		}
	}
	if (accurate) {
		performance.best_arm_accuracy += 1;
	} else {
		DisplayResult(result);
		DisplayResult(result_ground_truth);
	}
	return accurate;
}

void computeAvgPerf(Performance& performance) {
	performance.best_arm_accuracy /= (double) (config.repeat);
	performance.iterations /= (double) (config.repeat);
	performance.total_pulls /= (double) (config.repeat);
	performance.runtime /= (double) (config.repeat);

	for (int r = 0; r < 1000; ++r) {
		if (performance.candidate_size_per_iteration[r] > 0) {
			performance.candidate_size_per_iteration[r] /= (double) (config.repeat);
		}
	}

	performance.cost_accumulate[0] = performance.cost_per_itertaion[0];
	for (int r = 1; r < 1000; ++r) {
		performance.cost_accumulate[r] = performance.cost_accumulate[r - 1] + performance.cost_per_itertaion[r];
	}
	for (int r = 0; r < 1000; ++r) {
		performance.cost_accumulate[r] /= (double) (config.repeat);
	}

	performance.correct_error_cnt_per_itertaion_finished_accumulate[0].first =
			performance.correct_error_cnt_per_itertaion_finished[0].first;
	performance.correct_error_cnt_per_itertaion_finished_accumulate[0].second =
			performance.correct_error_cnt_per_itertaion_finished[0].second;
	for (int r = 1; r < 1000; ++r) {
		performance.correct_error_cnt_per_itertaion_finished_accumulate[r].first =
				performance.correct_error_cnt_per_itertaion_finished_accumulate[r - 1].first
						+ performance.correct_error_cnt_per_itertaion_finished[r].first;
		performance.correct_error_cnt_per_itertaion_finished_accumulate[r].second =
				performance.correct_error_cnt_per_itertaion_finished_accumulate[r - 1].second
						+ performance.correct_error_cnt_per_itertaion_finished[r].second;
	}


}


int main(int argc, char *argv[]) {
	program_start(argc, argv);

	for (int i = 0; i < argc; i++) {
		std::string help_str = ""
				"mab query --algo <algo> [options]\n"
				"mab generate-arms [options]\n"

				"algo: \n"
				"  delta_e (Delta Elimination)\n"
				"  topk_delta_e (Topk Delta Elimination)\n"
				"  d_delta_e (Delta Elimination with Round Limit)\n"
				"  topk_d_delta_e (Topk Delta Elimination with Round Limit)\n"
				"  egde (Exponential Gap Delta Elimination)\n"
				"options: \n"
				"  --prefix <prefix>\n"
				"  --num_arms <number of arms, 2000>\n"
				"  --epsilon <epsilon>\n"
				"  --delta <delta>\n"
				"  --k <top k>\n"
				"  --R <round limit>\n"
				"  --arm_dist <normal, uniform, segment>\n"
				"  --repeat <repeat times>\n"
				;

		if (std::string(argv[i]) == "--help") {
			std::cerr << help_str << std::endl;
			exit(0);
		}
	}

	config.action = argv[1];
	std::cout << "action: " << config.action << std::endl;

	for (int i = 0; i < argc; i++) {
		if (std::string(argv[i]) == "--prefix") {
			config.prefix = argv[i + 1];
		}
	}

	for (int i = 0; i < argc; i++) {
		std::string arg = argv[i];
		if (arg == "--algo") {
			config.algo = std::string(argv[i + 1]);
		} else if (arg == "--num_arms") {
			config.num_arms = atoi(argv[i + 1]);
			std::cout << "num_arms " << config.num_arms << std::endl;
		} else if (arg == "--epsilon") {
			config.epsilon = atof(argv[i + 1]);
			std::cout << "epsilon " << config.epsilon << std::endl;
		} else if (arg == "--delta") {
			config.delta = atof(argv[i + 1]);
			std::cout << "delta " << config.delta << std::endl;
		} else if (arg == "--result_dir") {
			config.result_dir = std::string(argv[i + 1]);
		} else if (arg == "--k") {
			config.k = atoi(argv[i + 1]);
			std::cout << "k " << config.k << std::endl;
		}
//		else if (arg == "--c") {
//			config.c = atof(argv[i + 1]);
//			std::cout << "c " << config.c << std::endl;
//		}
		else if (arg == "--R") {
			config.R = atoi(argv[i + 1]);
			std::cout << "R " << config.R << std::endl;
		} else if (arg == "--repeat") {
			config.repeat = atoi(argv[i + 1]);
			std::cout << "repeat " << config.repeat << std::endl;
		} else if (arg == "--arm_dist") {
			config.arm_dist = std::string(argv[i + 1]);
			std::cout << "arm_dist " << config.arm_dist << std::endl;
		} else if (arg == "--prefix") {
			// pass
		} else if (arg.substr(0, 2) == "--") {
			std::cout << "command not recognize " << arg << std::endl;
			exit(1);
		}
	}

	std::vector<std::string> possibleAlgo = { DELTA_E, TOPK_DELTA_E, D_DELTA_E, TOPK_D_DELTA_E, EGDE };

	if (config.action == QUERY) {
		auto f = find(possibleAlgo.begin(), possibleAlgo.end(), config.algo);
		if (f == possibleAlgo.end()) {
			std::cerr << "Wrong algo param: " << config.algo << std::endl;
			exit(1);
		}

		std::random_device rd;  //Will be used to obtain a seed for the random number engine
		std::mt19937 rand_gen(rd());  //Standard mersenne_twister_engine seeded with rd()
		std::uniform_real_distribution<> dist(0.0, 1.0);

		Performance performance;
		performance.init_accuracy_per_iteration();

		if (config.algo == DELTA_E) {
			config.c = 8;
			config.Q = config.c / config.epsilon / config.epsilon;

			for (int repi = 0; repi < config.repeat; ++repi) {
				std::cerr << "-----------repeat " << repi << "------------\n";
				std::vector<std::tuple<int, double, double> > arms = load_arms();
				std::vector<std::tuple<int, double, double> > result_ground_truth = GetGroundTruth(arms, config.k);

				//compute
				performance.itertaions_old_before_next_repeat = performance.iterations;
				performance.candidate_size_per_iteration[(int) (performance.iterations
						- performance.itertaions_old_before_next_repeat)] += arms.size();
				auto start = std::chrono::high_resolution_clock::now();
				std::tuple<int, double, double> best_arm =
						DeltaEliminateBestArm(arms, config.epsilon, config.delta, config.Q, dist, rand_gen, performance, config.scale_factor_ME, config.scale_factor_US, result_ground_truth, repi);
				auto stop = std::chrono::high_resolution_clock::now();
				auto duration = std::chrono::duration_cast < std::chrono::milliseconds > (stop - start);
				performance.runtime += duration.count();
				std::cerr << "runtime: " << duration.count() << "ms\n";

				//output performance metric
				std::vector<std::tuple<int, double, double> > result;
				result.push_back(best_arm);

				evaluateAccuracy(result_ground_truth, result, performance);
			}

			computeAvgPerf(performance);
			performance.display();
		} else if (config.algo == EGDE) {
			config.c = 8;
//			config.display();

			for (int repi = 0; repi < config.repeat; ++repi) {
				std::cerr << "-----------repeat " << repi << "------------\n";
				std::vector<std::tuple<int, double, double> > arms = load_arms();
				std::vector<std::tuple<int, double, double> > result_ground_truth = GetGroundTruth(arms, config.k);

				//compute
				performance.itertaions_old_before_next_repeat = performance.iterations;
				performance.candidate_size_per_iteration[(int) (performance.iterations
						- performance.itertaions_old_before_next_repeat)] += arms.size();
				auto start = std::chrono::high_resolution_clock::now();
				std::tuple<int, double, double> best_arm =
						ExponentialGapDeltaE(arms, config.delta, dist, rand_gen, performance, result_ground_truth, repi);
				auto stop = std::chrono::high_resolution_clock::now();
				auto duration = std::chrono::duration_cast < std::chrono::milliseconds > (stop - start);
				performance.runtime += duration.count();
				std::cerr << "runtime: " << duration.count() << "ms\n";

				//output performance metric
				std::vector<std::tuple<int, double, double> > result;
				result.push_back(best_arm);

				evaluateAccuracy(result_ground_truth, result, performance);
			}

			computeAvgPerf(performance);
			performance.display();
		} else if (config.algo == TOPK_DELTA_E) {
			config.c = 8;
			config.Q = config.c / config.epsilon / config.epsilon;

			for (int repi = 0; repi < config.repeat; ++repi) {
				std::cerr << "-----------repeat " << repi << "------------\n";
				//generate arms
				std::vector<std::tuple<int, double, double> > arms = load_arms();
				std::vector<std::tuple<int, double, double> > result_ground_truth = GetGroundTruth(arms, config.k);

				//compute
				performance.itertaions_old_before_next_repeat = performance.iterations;
				performance.candidate_size_per_iteration[(int) (performance.iterations
						- performance.itertaions_old_before_next_repeat)] += arms.size();
				auto start = std::chrono::high_resolution_clock::now();
				std::vector<std::tuple<int, double, double> > result =
						TopkDeltaEliminate(arms, config.epsilon, config.delta, config.Q, config.k, dist, rand_gen, performance, config.scale_factor_ME, config.scale_factor_US, result_ground_truth, repi);
				auto stop = std::chrono::high_resolution_clock::now();
				auto duration = std::chrono::duration_cast < std::chrono::milliseconds > (stop - start);
				performance.runtime += duration.count();
				std::cerr << "runtime: " << duration.count() << "ms\n";

				evaluateAccuracy(result_ground_truth, result, performance);
			}

			computeAvgPerf(performance);
			performance.display();
		} else if (config.algo == D_DELTA_E) {
			//config.c has to be 57
			config.c = 57;

			config.Q = config.c / config.epsilon / config.epsilon;

			for (int repi = 0; repi < config.repeat; ++repi) {
				std::cerr << "-----------repeat " << repi << "------------\n";
				//generate arms
				std::vector<std::tuple<int, double, double> > arms = load_arms();

				std::vector<std::tuple<int, double, double> > result_ground_truth = GetGroundTruth(arms, config.k);

				//compute
				performance.itertaions_old_before_next_repeat = performance.iterations;
				performance.candidate_size_per_iteration[(int) (performance.iterations
						- performance.itertaions_old_before_next_repeat)] += arms.size();
				auto start = std::chrono::high_resolution_clock::now();
				std::tuple<int, double, double> best_arm =
						DeltaEliminateBestArm_Deterministic(arms, config.epsilon, config.delta, config.Q, config.R, dist, rand_gen, performance, config.scale_factor_ME, config.scale_factor_US, result_ground_truth, repi);
				auto stop = std::chrono::high_resolution_clock::now();
				auto duration = std::chrono::duration_cast < std::chrono::milliseconds > (stop - start);
				performance.runtime += duration.count();
				std::cerr << "runtime: " << duration.count() << "ms\n";

				//output performance metric
				std::vector<std::tuple<int, double, double> > result;
				result.push_back(best_arm);

				evaluateAccuracy(result_ground_truth, result, performance);
			}
			computeAvgPerf(performance);
			performance.display();
		} else if (config.algo == TOPK_D_DELTA_E) {
			//config.c has to be 57
			config.c = 57;

			config.Q = config.c / config.epsilon / config.epsilon;
//			std::cerr << config.c / config.epsilon / config.epsilon << std::endl;
//			config.display();

			for (int repi = 0; repi < config.repeat; ++repi) {
				std::cerr << "-----------repeat " << repi << "------------\n";
				//generate arms
				std::vector<std::tuple<int, double, double> > arms = load_arms();

				std::vector<std::tuple<int, double, double> > result_ground_truth = GetGroundTruth(arms, config.k);

				//compute
				performance.itertaions_old_before_next_repeat = performance.iterations;
				performance.candidate_size_per_iteration[(int) (performance.iterations
						- performance.itertaions_old_before_next_repeat)] += arms.size();
				auto start = std::chrono::high_resolution_clock::now();
				std::vector<std::tuple<int, double, double> > result =
						TopKDeltaEliminate_Deterministic(arms, config.epsilon, config.delta, config.Q, config.R, config.k, dist, rand_gen, performance, config.scale_factor_ME, config.scale_factor_US, result_ground_truth, repi);
				auto stop = std::chrono::high_resolution_clock::now();
				auto duration = std::chrono::duration_cast < std::chrono::milliseconds > (stop - start);
				performance.runtime += duration.count();
				std::cerr << "runtime: " << duration.count() << "ms\n";

				evaluateAccuracy(result_ground_truth, result, performance);

//				DisplayResult(result);
//				DisplayResult(result_ground_truth);
			}

			computeAvgPerf(performance);
			performance.display();
		}

	} else if (config.action == GENERATE_ARMS) {
		if (config.arm_dist == UNIFORM) {
			generate_uniform_arms();
		} else if (config.arm_dist == SEGMENT) {
			generate_segment_arms();
		} else if (config.arm_dist == NORMAL) {
			generate_normal_arms();
		}
	}

	return 0;
}
