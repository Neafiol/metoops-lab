#include <iostream>
#include <fstream>
#include <vector>
#include <regex>
#include <filesystem>

namespace fs = std::filesystem;
using fs::directory_iterator;
using namespace std;

void UxSolution(vector<double>* x, int size, vector<double>* ia, vector<double>* al, vector<double>* au, vector<double>* di, vector<double>* y) {
	//Ux = y
	(*x)[size - 1] = (*y)[size - 1] / (*di)[size - 1];
	for (int k = size - 2; k >= 0; k--) {
		double sum = (*y)[k];
		for (int j = size - 1; j > k; j--) {
			int joffset = (*ia)[j + 1] - (*ia)[j] - j + k;
			if (joffset >= 0) {
				sum -= (*x)[j] * (*au)[(*ia)[j] - 1 + joffset];
			}
		}
		(*x)[k] = sum / (*di)[k];
	}
}

void LySolution(vector<double>* y, int size, vector<double>* ia, vector<double>* al, vector<double>* au, vector<double>* di, vector<double>* b) {
	//Ly = b
	(*y)[0] = (*b)[0];
	for (int k = 1; k < size; k++) {
		double sum = (*b)[k];
		for (int j = 0; j < k; j++) {
			int ioffset = (*ia)[k + 1] - (*ia)[k] - k + j;
			if (ioffset >= 0) {
				sum -= (*y)[j] * (*al)[(*ia)[k] - 1 + ioffset];
			}
		}
		(*y)[k] = sum;
	}
}

void solution(int size, vector<double>* ia, vector<double>* al, vector<double>* au, vector<double>* di, vector<double>* b, vector<double>* x) {
	vector<double> y(size);
	LySolution(&y, size, ia, al, au, di, b);
	UxSolution(x, size, ia, al, au, di, &y);
}

void LUprofile(int size, vector<double>* ia, vector<double>* al, vector<double>* au, vector<double>* di) {
	for (int k = 1; k < size; k++) {
		for (int i = k; i < size; i++) { // i - строка
			if ((*ia)[i + 1] - (*ia)[i] >= i - k + 1) {
				int ioffset = (*ia)[i + 1] - (*ia)[i] - i + k - 1;
				(*al)[(*ia)[i] + ioffset - 1] /= (*di)[k - 1];
				for (int j = k; j < size; j++) { // j - столбец
					if ((*ia)[j + 1] - (*ia)[j] >= j - k + 1) {
						int joffset = (*ia)[j + 1] - (*ia)[j] - j + k - 1;
						if (i == j) {
							(*di)[i] -= (*al)[(*ia)[i] + ioffset - 1] * (*au)[(*ia)[j] + joffset - 1];
						}
						else if (i < j) {
							(*au)[(*ia)[j] + joffset + i - k] -= (*al)[(*ia)[i] + ioffset - 1] * (*au)[(*ia)[j] + joffset - 1];
						}
						else {
							(*al)[(*ia)[i] + ioffset + j - k] -= (*al)[(*ia)[i] + ioffset - 1] * (*au)[(*ia)[j] + joffset - 1];
						}
					}
				}
			}
		}
	}
}

void printProfileForm(int size, vector<double>* ia, vector<double>* al, vector<double>* au, vector<double>* di) {
	cout << "\nia is: ";
	for (int i = 0; i < ia->size(); i++) {
		cout << (*ia)[i] << ' ';
	}
	cout << "\nal is: ";
	for (int i = 0; i < al->size(); i++) {
		cout << (*al)[i] << ' ';
	}
	cout << "\nau is: ";
	for (int i = 0; i < au->size(); i++) {
		cout << (*au)[i] << ' ';
	}
	cout << "\ndi is: ";
	for (int i = 0; i < di->size(); i++) {
		cout << (*di)[i] << ' ';
	}
	cout << '\n';
}

void completeTest(int n) {
	cout << "Test " << n << '\n';
	int size;
	vector<double> ia;
	vector<double> al;
	vector<double> au;
	vector<double> di;
	vector<double> b;
	cout << fs::current_path() << '\n';
	string fileIn = fs::current_path().string() + "\\test" + to_string(n) + "\\test.in";
	cout << fileIn << "\n";
	ifstream in(fileIn);
	if (!in.good()) {
		cout << "Test doesn't exist\n";
		return;
	}
	in >> size;
	for (int i = 0; i < size + 1; i++) {
		double v;
		in >> v;
		ia.push_back(v);
	}
	for (int i = 0; i < ia[size] - 1; i++) {
		double v;
		in >> v;
		al.push_back(v);
	}
	for (int i = 0; i < ia[size] - 1; i++) {
		double v;
		in >> v;
		au.push_back(v);
	}
	for (int i = 0; i < size; i++) {
		double v;
		in >> v;
		di.push_back(v);
	}
	for (int i = 0; i < size; i++) {
		double v;
		in >> v;
		b.push_back(v);
	}
	in.close();
	//printProfileForm(size, &ia, &al, &au, &di);
	LUprofile(size, &ia, &al, &au, &di);
	//printProfileForm(size, &ia, &al, &au, &di);
	vector<double> x(size);
	solution(size, &ia, &al, &au, &di, &b, &x);
	string fileOut = fs::current_path().string() + "\\test" + to_string(n) + "\\test.out";
	ofstream out(fileOut);
	for (int i = 0; i < x.size(); i++) {
		out << x[i] << '\n';
	}
	out.close();
}

void createTest(vector<vector<double>>* A, vector<double>* b) {
	const regex reg("test([0-9]*)");
	smatch match;
	int max = 0;
	for (const auto& file : directory_iterator(fs::current_path())) {
		string name(file.path().filename().string());
		regex_search(name, match, reg);
		if (match.size() != 0 && match[1] > max) max = stoi(match[1]);
	}
	fs::create_directory(fs::current_path().string() + "\\test" + to_string(max + 1));
	string fileName = fs::current_path().string() + "\\test" + to_string(max + 1) + "\\test.in";
	ofstream file(fileName);
	int size = A->size();
	vector<int> count;
	for (int i = 1; i < size; i++) {
		int zeros = 0;
		int shift = 0;
		while (i != shift && (*A)[i][shift] == 0 && (*A)[shift][i] == 0) {
			zeros++;
			shift++;
		}
		count.push_back(i - zeros);
	}
	vector<double> ia(size + 1);
	vector<double> di(size);
	vector<double> al, au;
	ia[0] = 1;
	ia[1] = 1;
	for (int i = 2; i < size + 1; i++) {
		ia[i] = ia[i - 1] + count[i - 2];
	}
	for (int i = 1; i < size; i++) {
		for (int j = i - count[i - 1]; j < i; j++) {
			al.push_back((*A)[i][j]);
		}
	}
	for (int i = 1; i < size; i++) {
		for (int j = i - count[i - 1]; j < i; j++) {
			au.push_back((*A)[j][i]);
		}
	}
	for (int i = 0; i < size; i++) {
		di[i] = (*A)[i][i];
	}
	file << size << '\n';
	for (int i = 0; i < ia.size(); i++) {
		file << ia[i] << ' ';
	}
	file << '\n';
	for (int i = 0; i < al.size(); i++) {
		file << al[i] << ' ';
	}
	file << '\n';
	for (int i = 0; i < au.size(); i++) {
		file << au[i] << ' ';
	}
	file << '\n';
	for (int i = 0; i < di.size(); i++) {
		file << di[i] << ' ';
	}
	file << '\n';
	for (int i = 0; i < b->size(); i++) {
		file << (*b)[i] << ' ';
	}
	file.close();
}

int main() {
	/*vector<vector<double>> A = { {3, 0, 13, 0, 0},
								   {0, 11, 1, 0, 0},
								   {2, 0, 25, 4, 0},
								   {0, 9, 5, 7, 0},
								   {0, 0, 0, 0, 19} };
	vector<double> b = { 0, 0, 0, 0, 0 };
	createTest(&A, &b);*/
	completeTest(1);
	return 0;
}