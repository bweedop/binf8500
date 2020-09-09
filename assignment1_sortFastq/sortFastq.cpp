#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstdio>

using namespace std;

struct Sample {
  string identifier;
  string sequence;
  string description;
  string q_score;
};

void read_data(vector<Sample> &sample_array, vector<Sample *> &ptr_array, const char* file_path) {
  ifstream fin(file_path);
  string line;
 
  while(!fin.eof()) {
    Sample *ptr = new Sample;
    getline(fin, line);
    ptr->identifier = line;
    getline(fin, line);
    ptr->sequence = line;
    getline(fin, line);
    ptr->description = line;
    getline(fin, line);
    ptr->q_score = line;
    sample_array.push_back(*ptr);
    ptr_array.push_back(ptr);
  }
  fin.close();
}

void swap(vector<Sample *> &ptr_array, int idx, int jdx) {
  Sample* tmp = ptr_array[idx];
  ptr_array[idx] = ptr_array[jdx];
  ptr_array[jdx] = tmp;
}

int partition(vector<Sample *> &ptr_array, int low, int high) {
  int mid = (low + high) / 2;
  if (ptr_array[mid]->sequence < ptr_array[low]->sequence) {
    swap(ptr_array, low, mid);
  } else if (ptr_array[high]->sequence < ptr_array[low]->sequence) {
    swap(ptr_array, low, high);
  } else if (ptr_array[mid]->sequence < ptr_array[high]->sequence){
    swap(ptr_array, mid, high);
  }
  Sample* pivot = ptr_array[high];
  int j = low - 1;

  for (int i = low; i <= high; i++) {
    if (ptr_array[i]->sequence < pivot->sequence) {
      j++;
      swap(ptr_array, j, i);
    }
  }
  swap(ptr_array, j + 1, high);
  return (j + 1);
}

void quicksort(vector<Sample *> &ptr_array, int low, int high) {
  if (low < high) {
    int p = partition(ptr_array, low, high);
    // Left side and then right side recursively
    quicksort(ptr_array, low, p - 1);
    quicksort(ptr_array, p + 1, high);
  }
}

int main(int argc, const char **argv) {
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
  vector<Sample*> ptr_array;
  
  read_data(sample_array, ptr_array, argv[1]);
  
  quicksort(ptr_array, 0, sample_array.size() - 1);

  for (int i = 1; i < sample_array.size(); i++) {
    cout << ptr_array[i]->identifier << "\n";
    cout << ptr_array[i]->sequence << "\n";
    cout << ptr_array[i]->description << "\n";
    cout << ptr_array[i]->q_score << "\n";
  }

  return 0;
}
