#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <bits/stdc++.h>

using namespace std;

/*
Read FASTA file and add label and sequence to their respective vectors

@param filename Path to file to read
@param labels Vector where label will be stored
@param sequences Vector where sequence will be stored  
*/
void read_data (char *filename, vector<string> &labels, vector<string> &sequences) {
    // temporary string to read file line by line
    string s;
    string tmp_sequence = "";

    ifstream fin (filename);
    if (fin.is_open()) {
        // Get label and check if in FASTA format
        getline(fin, s);
        if (s[0] != '>') {
            throw domain_error("Error: file must be in FASTA format. Missing '>' before label");
            exit(2);
        }
        labels.push_back(s);
        // Done with label; Now get the sequence
        while(getline(fin, s)) {
            // Transform sequence to be all uppercase
            transform(s.begin(), s.end(), s.begin(), ::toupper);
            tmp_sequence += s;
        };
        // Now that full sequence is in tmp_sequence, add to sequences vector
        sequences.push_back(tmp_sequence);
        fin.close();
    } else {
        // Check if file can not be opened
        throw domain_error("Error: could not open data file...");
        exit(1);
    }
}

/*
Initialize grid with user-defined gap penalty

@param grid Container where grid is stored
@param gap_penalty User-defined gap penalty
*/
void initialize_grid (vector<vector <int> > &grid, int gap_penalty) {
    // Fill first row with incremental gap penalty values
    for (int i = 1; i < grid[0].size(); i++) {
        grid[0][i] = i*gap_penalty;
    }
    // Fill first column with gap penalty values 
    for (int j = 1; j < grid.size(); j++) {
        grid[j][0] = j*gap_penalty;
    }
}

/*
Find maximum value in grid

@param top_cell Integer value contained in cell above the current cell
@param side_cell Integer value from cell to side of the current cell
@param diag_cell Integer value from cell diagonal to the current cell
*/
int max_value (int top_cell, int side_cell, int diag_cell) {
    int ans = -9999999;
    int tmp[] = { top_cell, side_cell, diag_cell };
    for (int i = 0; i < 3; i++) {
        if (tmp[i] > ans) {
            ans = tmp[i];
        }
    }
    return ans;
}

/*
Assign match value at current cell

@param grid Container where grid will is stored
@param row Index of current row
@param col Index of current column
@param match_score User-defined match score
@param gap_penalty User-defined gap penalty
*/
void match (vector<vector <int> > &grid, int row, int col, int match_score, int gap_penalty) {
    // Calculate all possible values
    int side_cell = grid[row][col-1] + gap_penalty;
    int diag_cell = grid[row-1][col-1] + match_score;
    int top_cell = grid[row-1][col] + gap_penalty;
    // Assign the maximum value calculated to current cell 
    grid[row][col] = max_value(top_cell, side_cell, diag_cell);
}

/*
Assign mismatch value at current cell

@param grid Container where grid  is stored
@param row Index of current row
@param col Index of current column
@param mismatch_score User-defined mismatch score
@param gap_penalty User-defined gap penalty
*/
void mismatch (vector<vector <int> > &grid, int row, int col, int mismatch_score, int gap_penalty) {
    // Calculate all possible values
    int side_cell = grid[row][col-1] + gap_penalty;
    int diag_cell = grid[row-1][col-1] + mismatch_score;
    int top_cell = grid[row-1][col] + gap_penalty;
    // Assign the maximum value calculated to current cell 
    grid[row][col] = max_value(top_cell, side_cell, diag_cell);
}

/*
Fill initialized grid given matches and mismatches

@param grid Container where grid is stored
@param sequences Container where sequences are stored
@param current_row Index of current row
@param current_col Index of current column
@param match_score User-defined match score
@param mismatch_score User-defined mismatch score
@param gap_penalty User-defined gap penalty
*/
int fill_grid (vector<vector <int> > &grid, vector<string> sequences, int current_row, int current_col, int match_score, int mismatch_score, int gap_penalty) {
    // Check if all columns have been filled
    if (current_col == grid[0].size())
        return 1;
    // Check if all rows have been filled
    if (current_row == grid.size())
        return 0;
    
    // Check sequences for match
    if (sequences[0][current_row - 1] == sequences[1][current_col - 1]) {
        match(grid, current_row, current_col, match_score, gap_penalty);
    } else {
        mismatch(grid, current_row, current_col, mismatch_score, gap_penalty);
    }

    // Call fill_grid recursively until all columns 
    if (fill_grid(grid, sequences, current_row, current_col + 1, match_score, mismatch_score, gap_penalty) == 1) {
        fill_grid(grid, sequences, current_row + 1, 1, match_score, mismatch_score, gap_penalty);
    }
    return 0;
}

/*
Produce similarity score given current characters in sequences

@param s1 Character from sequence 1
@param s2 Character from sequence 2
@param match User-defined match score
@param mismatch User-defined mismatch score
*/
int similarity_score(char s1, char s2, int match, int mismatch) {
    if (s1 == s2) {
        return match;
    } else {
        return mismatch;
    }
}

/*
Traverse grid and produce aligned sequences along the way

@param grid Container where grid is stored
@param sequences Container where sequences are stored
@param current_row Index of current row
@param current_col Index of current column
@param aligned1 First aligned sequence 
@param aligned2 Second aligned sequence
@param gap_penalty User-defined gap penalty
@param match User-defined match score
@param mismatch User-defined mismatch score
*/
void grid_traversal(vector<vector <int> > &grid, vector<string> sequences, int current_row, int current_col, string &aligned1, string &aligned2, int gap_penalty, int match, int mismatch) {
    // Continue until reach left most cell
    while (current_row > 0 || current_col > 0) {
        // Get values in current cell and cells to side, top and diagonal to it
        int current_cell = grid[current_row][current_col];
        int top_cell = grid[current_row-1][current_col];
        int side_cell = grid[current_row][current_col-1];
        int diag_cell = grid[current_row-1][current_col-1];
        // Get similarity score for current characters
        int sim_score = similarity_score(sequences[0][current_row - 1], sequences[1][current_col - 1], match, mismatch);
        // Find match
        if (current_row > 0 && current_col > 0 && current_cell == (diag_cell +  sim_score)) {
            aligned1.insert(0, 1, sequences[0][current_row - 1]);
            aligned2.insert(0, 1, sequences[1][current_col - 1]);    
            current_row = current_row - 1;
            current_col = current_col - 1;
        // gap in sequence 2
        } else if (current_row > 0 && current_cell == top_cell + gap_penalty) {
            aligned1.insert(0, 1, sequences[0][current_row - 1]);
            aligned2.insert(0, 1, '-');
            current_row = current_row - 1;
        // gap in sequence 1
        } else {
            aligned1.insert(0, 1, '-');
            aligned2.insert(0, 1, sequences[1][current_col - 1]);
            current_col = current_col - 1;
        }
    }
}

int main (int argc, char **argv) {
    if (argc == 1) {
        cout << "Command-line argument must include the following (in order):\n" <<
                "\t1. Command to run program (e.g. './nw')\n" <<
                "\t2. Path to FASTA-formatted sequence (e.g. 'seq1.fasta')\n" <<
                "\t3. Path to FASTA-formatted sequence (e.g.'seq2.fasta')\n" <<
                "\t4. Gap penalty (e.g. '-1')\n" <<
                "\t5. Match score (e.g. '1')\n" <<
                "\t6. Mismatch score (e.g. '-1')\n";
        exit(-1);
    }
    // vector for labels
    vector<string> labels;
    // vector for sequences
    vector<string> seqs;

    // Sequences must be in separate files.
    // Read data into labels and sequences.
    read_data(argv[1], labels, seqs);
    read_data(argv[2], labels, seqs);

    int n = seqs[0].size() + 1; 
    int m = seqs[1].size() + 1;
  
    // Create a vector containing n vectors of size m.  
    vector<vector <int> > f_matrix(n, vector<int> (m, 0));

    /* 
    Initialize the first row and first column with numbers decreasing from 0 to negative value of the size of each sequence
    */
    initialize_grid(f_matrix, atoi(argv[3]));

    //Fill the matrix with values based on matches and mismatches
    fill_grid(f_matrix, seqs, 1, 1, atoi(argv[4]), atoi(argv[5]), atoi(argv[3]));

    // Traverse matrix from the bottom right to top left
    string aligned1;
    string aligned2;
    
    grid_traversal(f_matrix, seqs, f_matrix.size() - 1, f_matrix[0].size() - 1, aligned1, aligned2, atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));

    int idx = 0;
    while (idx < aligned1.size()) {

        ostringstream tmp;
        tmp << idx;
        
        cout << idx << ": ";
        for (int i = idx; i < aligned1.size(); i++) {
            cout << aligned1[i];
            if (i % 80 == 0 && i != 0) {
                cout << "\n";
                break;
            } else if (i == aligned1.size() - 1) {
                cout << "\n";
                break;
            }
        }
        
        cout << string(tmp.str().size() + 2, ' ');
        for (int j = idx; j < aligned1.size(); j++) {
            if (aligned1[j] == aligned2[j]) {
                cout << "|";
            } else {
                cout << " ";
            }
            if (j % 80 == 0 && j != 0) {
                cout << "\n";
                break;
            } else if (j == aligned1.size() - 1) {
                cout << "\n";
                break;
            }
        }
        cout << idx << ": ";
        for (int k = idx; k < aligned2.size(); k++) {
            cout << aligned2[k];
            if (k % 80 == 0 && k != 0) {
                idx = k + 1;
                cout << "\n";
                break;
            } else if (k == aligned2.size() - 1) {
                idx = k + 1;
                cout << "\n";
                break;
            }
        }
    }

    return 0;
}