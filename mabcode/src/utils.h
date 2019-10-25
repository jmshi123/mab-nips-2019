/*
 * utils.h
 *
 *  Created on: Dec 28, 2018
 *      Author: jmshinju
 */

#ifndef UTILS_H_
#define UTILS_H_
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <string>
#include <algorithm>
#include <list>
#include <unordered_map>
#include <utility>

typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::microseconds microsec;
typedef std::chrono::milliseconds millisec;
typedef std::chrono::duration<float> floatsec;

static bool sortDescbyThird(const std::tuple<int, double, double>& a, const std::tuple<int, double, double>& b) {
	return (std::get < 2 > (a) > std::get < 2 > (b));
}

static bool sortDescbySecond(const std::tuple<int, double, double>& a, const std::tuple<int, double, double>& b) {
	return (std::get < 1 > (a) > std::get < 1 > (b));
}

template<typename T>
static void print_arr(T *arr, size_t size, std::string tag) {
	std::cout << tag << ": ";
	for (size_t idx = 0; idx < size; idx++) {
		std::cout << arr[idx] << ", ";
	}
	std::cout << std::endl;
}

static std::string tuple_to_string(const std::tuple<int, double, double> &a) {
	std::string result = "";
	result += std::to_string(std::get<0>(a)) + " " + std::to_string(std::get<1>(a));
	return result;
}

template<typename T>
static void print_arr(T *arr, size_t size, T threshold, std::string tag) {
	std::cout << tag << ": ";
	for (size_t idx = 0; idx < size; idx++) {
		if (arr[idx] >= threshold) std::cout << idx << "(" << arr[idx] << "), ";
	}
	std::cout << std::endl;
}

template<typename T>
static T sum_arr(T *arr, size_t size) {
	T sum = 0;
	for (size_t idx = 0; idx < size; idx++) {
		sum += arr[idx];
	}
	return sum;
}

void generate_ss_query(int num_vertices, int query_size) {
	std::string filename = "ssquery.txt";
	std::ofstream queryfile(filename);
	for (int i = 0; i < query_size; i++) {
		int v = rand() % num_vertices;
		queryfile << v << std::endl;
	}
}

bool file_exists_test(const std::string &name) {
	std::ifstream f(name.c_str());
	if (f.good()) {
		f.close();
		return true;
	} else {
		f.close();
		return false;
	}
}

template<typename IndexType>
void load_ss_query(std::vector<IndexType>& queries, std::string arm_dir) {
	std::string filename = arm_dir + "arms.txt";
	if (!file_exists_test(filename)) {
		std::cerr << "query file does not exist, please generate ss query files first" << std::endl;
		exit(0);
	}
	std::ifstream queryfile(filename);
	IndexType v;
	while (queryfile >> v) {
		queries.push_back(v);
	}
	std::cout << "loaded queries from " << filename << std::endl;
}

static std::string get_current_time_str() {
	time_t rawtime;
	struct tm *timeinfo;
	char buffer[80];

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", timeinfo);
	std::string str(buffer);

	return str;
}

const std::string Green = "\033[0;32m";
const std::string Reset = "\033[0m";
const std::string Red = "\033[0;31m";

static void program_start(int argc, char **argv) {

	std::cerr << Green << "--------------start------------" << get_current_time_str() << Reset << std::endl;
	std::string combine = "";
	for (int i = 1; i < argc; i++) {
		combine += argv[i];
		combine += " ";
	}
	std::cerr << Green << "args:" << combine << Reset << std::endl;
}

#endif /* UTILS_H_ */
