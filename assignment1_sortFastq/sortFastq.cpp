#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include <sys/time.h>

using namespace std;

struct Sample {
  string identifier;
  string sequence;
  string description;
  string q_score;
};

void read_data(vector<Sample> &sample_array, string file_path) {
  ifstream fin;
  string line;
  int i = 0;

  fin.open(file_path, ios::in);
  assert(!fin.fail());

  while(!fin.eof()) {
    sample_array.push_back(Sample());
    getline(fin, line);
    sample_array[i].identifier = line;
    getline(fin, line);
    sample_array[i].sequence = line;
    getline(fin, line);
    sample_array[i].description = line;
    getline(fin, line);
    sample_array[i].q_score = line;
    i++;
  }
  fin.close();
}

void swap(vector<Sample> &sample_array, int idx, int jdx) {
  Sample tmp = sample_array[idx];
  sample_array[idx] = sample_array[jdx];
  sample_array[jdx] = tmp; 
}

int partition(vector<Sample> &sample_array, int low, int high) {
  Sample pivot = sample_array[high];
  int j = low - 1;

  for (int i = low; i <= high; i++) {
    if (sample_array[i].sequence < pivot.sequence) {
      j++;
      swap(sample_array, j, i);
    }
  }
  swap(sample_array, j + 1, high);
  return (j + 1);
}

void quicksort(vector<Sample> &sample_array, int low, int high) {
  if (low < high) {
    int p = partition(sample_array, low, high);
    // Left side and then right side recursively
    quicksort(sample_array, low, p - 1);
    quicksort(sample_array, p + 1, high);
  }
}

double time_to_double(timeval *t)
{
    return (t->tv_sec + (t->tv_usec/1000000.0)) * 1000.0;
}

double time_diff(timeval *t1, timeval *t2)
{
    return time_to_double(t2) - time_to_double(t1);
}

int main(int argc, char **argv) {
  // Tests to check that the number of arguments and file type are adequate
  if (argc > 2) {
    cout << "sortFastq takes no more than one command-line argument" << endl;
    exit(-1);
  }
  string delimiter = ".fastq";
  string const file_path = argv[1];
  if (file_path.find(delimiter) == string::npos) {
    cout << "Invalid file type provided. File must be a fastq file." << endl;
    exit(-1);
  }

  vector<Sample> sample_array;

  timeval t1, t2;
  gettimeofday(&t1, NULL);
  read_data(sample_array, argv[1]);
  gettimeofday(&t2, NULL);
  cout << "read_data1" << ": " << time_diff(&t1, &t2) << "ms" << endl;

  timeval t3, t4;
  gettimeofday(&t3, NULL);
  quicksort(sample_array, 0, sample_array.size() - 1);
  gettimeofday(&t4, NULL);
  cout << "quicksort:" << time_diff(&t1, &t2) << "ms" << endl;
  
  return 0;
}





