#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>    //JM: added for atoi to work

using namespace std;

// kmers_permutations will be called recursively to get all permutations of base pairs loaded into the kmers vector. 
void kmers_permutations (char *str, char* data, int k, int index, int start, int end, vector<string> &kmers) {
  // If the index has reached the desired kmer length, load into vector
  if (index == k) {
    kmers.push_back(data);
    return;
  }
  // Loop through array of base pairs and concatenate the index value to data. Call function again with new parameters.
  for (int i = 0; i < end; i++) {
    data[index] = str[i];
    kmers_permutations(str, data, k, index + 1, i, end, kmers);
  }
  return;
}
  
void possible_kmers (int k, vector<string> &kmers) {
  char bp[] = {'A', 'C', 'G', 'T'}; // Array holding all base pairs
  int n = sizeof(bp)/sizeof(bp[0]); // Size of the array
  char data[k]; // temporary placeholder for the each of the kmers 
  data[k] = '\0'; // ending of sequence
  kmers_permutations(bp, data, k, 0, 0, n, kmers);
}

int main (int argc, char* argv[]) {

  int k; // length of the k-mers
  k = atoi(argv[1]);
  // User must specify k-mer length which is >= 2
  if (k < 2) {
    cout << "Error: length of desired k-mer must be >= 2" << endl;
    return 0;
  }

  string seq_id; // placeholder for sequence id
  string tmp; // temporary holder for line from file
  string seq; // placeholder for full sequence from file
  vector<string> kmers; // vector of all possible k-mers given a DNA seqeunce

  // Read file given by the user
  ifstream fasta_file (argv[2]);
  if (fasta_file.is_open()) {
    getline(fasta_file, seq_id);
    cout << seq_id << endl; // print sequence identifier
    while (getline(fasta_file, tmp)) {
      seq = seq + tmp; // concatenate sequence into a single string
    }
  } else {
    // Error if file cannot be opened
    cout << "Error: Could not read file" << endl;
  }

  possible_kmers(k, kmers); // Load all possible k-mers
  
  for (int i = 0; i < kmers.size(); i++) {
    int tmp = 0;
    for (int j = 0; j < seq.size()-(k-1); j++) {
      string tmp_kmer = seq.substr(j, k);
      
      if (kmers[i] == tmp_kmer) {
	tmp++;
      }
    }
    cout << kmers[i] << " " << tmp << endl;
    }
  
  return 0;
}
