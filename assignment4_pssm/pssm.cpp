#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <math.h>
#include <bits/stdc++.h>

using namespace std;

void scan_motifs(char *filepath, int motifs_info[])
{
    string s;

    ifstream fin(filepath);
    while (getline(fin, s))
    {
        if (s[0] != '>')
        {
            if (motifs_info[0] == 0) 
            {
                motifs_info[0] = s.size();
            }
            else if (motifs_info[0] != s.size())
            {
                cout << "Error: Motifs are not the same length\n";
                exit(-1);
            }
            motifs_info[1]++;
        }
    }
}

int scan_genome(char *filepath) {
    string s;
    int genome_size = 0;

    ifstream fin(filepath);
    getline(fin, s);
    if (s[0] != '>')
    {
        cout << "Error: Sequence file is not FASTA format\n";
        exit(-1);
    }
    while (getline(fin, s))
    {
        genome_size += s.size();
    }
    return genome_size;
}

void frequency_matrix(int motif[], int motif_size, vector<vector<float>> &f_matrix)
{

    for (int i = 0; i < motif_size; i++)
    {
        if (motif[i] > 3)
        {
            cout << "Ambiguous character in motif" << i << "\n";
            exit(-3);
        }
        else
        {
            f_matrix[motif[i]][i]++;
        }
    }
}

void convert_matrix(vector<vector<float>> &f_matrix)
{
    for (int col = 0; col < f_matrix[0].size(); col++)
    {
        float total = 0.0;
        for (int row = 0; row < f_matrix.size(); row++)
        {
            total += f_matrix[row][col];
        }
        for (int i = 0; i < f_matrix.size(); i++)
        {
            f_matrix[i][col] = f_matrix[i][col] / total;
        }
    }
}

void pssm(vector<vector<float>> &f_matrix, float gc_content, float at)
{
    // Background nucleotide probabilities that correspond to rows of frequency matrix
    vector<float> background_probs{at / 2, at / 2, gc_content / 2, gc_content / 2};

    for (int row = 0; row < f_matrix.size(); row++)
    {
        for (int col = 0; col < f_matrix[0].size(); col++)
        {
            f_matrix[row][col] = log2(f_matrix[row][col] / background_probs[row]);
        }
    }
}

float get_score(vector<vector<float>> f_matrix, int motif[], int motif_size, float ambiguous_content)
{
    if (motif_size != f_matrix[0].size())
    {
        cout << "Error: Sequence to be scored is not the same length as motif\n";
        exit(-4);
    }
    float score = 0.0;

    for (int i = 0; i < motif_size; i++)
    {
        if (motif[i] < 4) 
        {
            score += f_matrix[motif[i]][i];
        } 
        else 
        {
            score += log2(0.25);
        }
    }
    return score;
}

void segment_arr(int arr[], int start, int end, int genome[]){
    int idx = 0;
    for (int i = start; i < end; i++) 
    {
        arr[idx] = genome[i];
        idx++;
    }
}

void get_complement(int arr[], int motif_length, int genome_segment[]){
    int idx = 0;

    unordered_map<int, int> complement_map = {
        {0, 1},
        {1, 0},
        {2, 3},
        {3, 2}};

    for (int i = (motif_length - 1); i >= 0; i--) 
    {
        unordered_map<int, int>::const_iterator val = complement_map.find(genome_segment[i]);
        if (val != complement_map.end()) {
            arr[idx] = val->second;
        
        }
        else
        {
            arr[idx] = 4;
        }
        idx++;
    }
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        cout << "Position-Specific Scoring Matrix requires the following arguments: \n\n";
        cout << "\tProgram command-line call (e.g. './pssm')\n";
        cout << "\tPath to alignment file (e.g. 'FruR.fasta')\n";
        cout << "\tPath to sequence file (e.g. 'ecoK12-MG1655.fasta')\n";
        cout << "\tMinimum score cutoff (e.g. '13.25')\n";
    }
    char bp[] = {'A', 'T', 'C', 'G', 'N'};
    unordered_map<char, int> bp_map = {
        {'A', 0},
        {'T', 1},
        {'C', 2},
        {'G', 3}};
    
    // Get all motifs into vector
    int motifs_info[2] = { 0, 0 };
    scan_motifs(argv[1], motifs_info);
 
    int motifs[motifs_info[1]][motifs_info[0]];
    int motif_idx = 0;
    string s;

    ifstream motifs_in(argv[1]);
    while (getline(motifs_in, s))
    {
        if (s[0] != '>')
        {
            transform(s.begin(), s.end(), s.begin(), ::toupper);
            for (int c = 0; c < s.size(); c++)
            {
                unordered_map<char, int>::const_iterator val = bp_map.find(s[c]);
                if (val == bp_map.end())
                {
                    cout << "Error: Ambiguous character in motif" << "\n";
                    exit(-1);   
                }
                else
                {   
                    motifs[motif_idx][c] = val->second;
                }
            }
            motif_idx++;
        }
    }
    
    // Get genome size
    int genome_size = scan_genome(argv[2]);
    
    // Get genome and complement as array of ints
    int* genome = new int[genome_size];

    int genome_idx = 0;

    ifstream genome_in(argv[2]);
    getline(genome_in, s);
    if (s[0] != '>')
    {
        cout << "Error: Sequence file is not FASTA format\n";
        exit(-1);
    }
    while (getline(genome_in, s))
    {
        transform(s.begin(), s.end(), s.begin(), ::toupper);
        for (int c = 0; c < s.size(); c++) {
            unordered_map<char, int>::const_iterator genome_val = bp_map.find(s[c]);
            if (genome_val != bp_map.end())
            {
                genome[genome_idx] = genome_val->second;
            }
            else
            {
                // Dealing with ambiguous codes.
                genome[genome_idx] = 4;
            }
            genome_idx++;
        }
    }

    // Calculate GC content
    float gc = 0.0;
    float at = 0.0;
    for (int c = 0; c < genome_size; c++)
    {
        if (genome[c] == 2 || genome[c] == 3)
        {
            gc++;
        }
        else if (genome[c] == 0 || genome[c] == 1) 
        {
            at++;
        }
    }
    at = at / genome_size;
    gc = gc / genome_size;

 
    // Create a vector containing 4 vectors (ATCG) each being the size of the motif
    vector<vector<float>> f_matrix(4, vector<float>(motifs_info[0], 0));

    // Take each sequence and add to the frequency matrix
    for (int i = 0; i < motifs_info[1]; i++)
    {
        frequency_matrix(motifs[i], motifs_info[0], f_matrix);
    }

    // Print frequency matrix to output
    cout << "Frequency matrix:\n";
    for (int i = 0; i < f_matrix.size(); i++)
    {
        cout << bp[i] << ": ";
        for (int j = 0; j < f_matrix[0].size(); j++)
        {
            cout << f_matrix[i][j] << " ";
        }
        cout << "\n\n";
    }
 
    // Add pseudocounts to the frequency matrix
    for (int i = 0; i < f_matrix.size(); i++)
    {
        for (int j = 0; j < f_matrix[0].size(); j++)
        {
            f_matrix[i][j] += 0.25;
        }
    }

    // Convert f_matrix to probability matrix
    convert_matrix(f_matrix);
    //Print probability matrix
    cout << "Probability matrix:\n";
    cout << fixed << setprecision(3);
    for (int i = 0; i < f_matrix.size(); i++)
    {
        cout << bp[i] << ": ";
        for (int j = 0; j < f_matrix[0].size(); j++)
        {
            cout << f_matrix[i][j] << " ";
        }
        cout << "\n\n";
    }
       
    // Probability matrix to PSSM
    pssm(f_matrix, gc, at);
    cout << "PSSM:\n";
    for (int i = 0; i < f_matrix.size(); i++)
    {
        cout << bp[i] << ": ";
        for (int j = 0; j < f_matrix[0].size(); j++)
        {
            cout << f_matrix[i][j] << " ";
        }
        cout << "\n";
    }

    cout << "\n";

    //Training set scores
    float cutoff = 99999999.0;
    cout << "Training set motif scores: \n";
    for (int i = 0; i < motifs_info[1]; i++) {
        float tmp_score = get_score(f_matrix, motifs[i], motifs_info[0], 1 - (at + gc));
        for (int c = 0; c < motifs_info[0]; c++) {
            cout << bp[motifs[i][c]];
        }
        cout << "\t\t" << tmp_score << "\n";
        if (tmp_score < cutoff)
        {
            cutoff = tmp_score;
        }
    }
    //Add some space    
    cout << "\n";
    
    //Print header
    cout << "Matches with score " << cutoff << " or higher in " << argv[2] << " (length " << genome_size << " bp):\n\n";
    cout << "Start\tEnd\tStrand\tSequence\t\tScore\n";
    
    // Length of motifs and a placeholder for score
    int genome_segment[motifs_info[0]];
    int complement_segment[motifs_info[0]];
    float tmp_score;
    // Loop through entire genome and its complement
    for (int low = 0; low < genome_size - motifs_info[0]; low++) {
        
        // Get segments to test
        segment_arr(genome_segment, low, low + motifs_info[0], genome);
        get_complement(complement_segment, motifs_info[0], genome_segment);
        
        // Get score of genome segment
        tmp_score = get_score(f_matrix, genome_segment, motifs_info[0], 1 - (at + gc));
        if (tmp_score >= cutoff) {
            cout << low + 1 << "\t" << low + motifs_info[0] << "\t" << "+\t";
            for (int c = 0; c < motifs_info[0]; c++) 
            {
                cout << bp[genome_segment[c]];
            }
            cout << "\t" << tmp_score << "\n";
        }
        
        // Get score of complement sequence
        tmp_score = get_score(f_matrix, complement_segment, motifs_info[0], 1 - (at + gc));
        if (tmp_score >= cutoff) {
            cout << low + 1 << "\t" << low + motifs_info[0] << "\t" << "-\t";
            for (int c = 0; c < motifs_info[0]; c++) 
            {
                cout << bp[complement_segment[c]];
            }
            cout << "\t" << tmp_score << "\n";
        }
    }
    delete[] genome;
    return 0;
}