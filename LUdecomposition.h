#pragma once
#include <vector>
#include <string>

void createTest(std::vector<std::vector<double>>& A, std::vector<double>& b, std::string const& dir_prefix = "test");
void completeTest(int n, std::string const& dir_prefix = "test");
std::vector<double> readResult(int n, int size, std::string const& dir_prefix = "test");
