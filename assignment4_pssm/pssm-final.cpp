#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
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

void frequency_matrix(int motif[], int motif_size, float **f_matrix)
{
    for (int i = 0; i < motif_size; i++)
    {
        if (motif[i] > 3)
        {
            cout << "Error: Ambiguous character in motif" << i << "\n";
            exit(-3);
        }
        else
        {
            f_matrix[motif[i]][i]++;
        }
    }
}

void convert_matrix(float **f_matrix, int motif_size)
{
    for (int col = 0; col < motif_size; col++)
    {
        float total = 0.0;
        for (int row = 0; row < 4; row++)
        {
            total += f_matrix[row][col];
        }
        for (int i = 0; i < 4; i++)
        {
            f_matrix[i][col] = f_matrix[i][col] / total;
        }
    }
}

void pssm(float **f_matrix, int motif_size, float gc_content, float at)
{
    // Background nucleotide probabilities that correspond to rows of frequency matrix
    float background_probs[4] = {at / 2, at / 2, gc_content / 2, gc_content / 2};

    for (int row = 0; row < 4; row++)
    {
        for (int col = 0; col < motif_size; col++)
        {
            f_matrix[row][col] = log2(f_matrix[row][col] / background_probs[row]);
        }
    }
}

float get_score(float **f_matrix, int motif[], int motif_size, float ambiguous_content)
{
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
   int genome_idx = 0;

    for (int i = (motif_length - 1); i >= 0; i--) 
    {
        if (genome_segment[i] == 0)
        {
            arr[genome_idx] = 1;
        }
        else if (genome_segment[i] == 1)
        {
            arr[genome_idx] = 0;
        }
        else if (genome_segment[i] == 2)
        {
            arr[genome_idx] = 3;
        }
        else if (genome_segment[i] == 3)
        {
            arr[genome_idx] = 2;
        }
        else
        {
            arr[genome_idx] = 4;   
        }
        genome_idx++;
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
    }
    char bp[] = {'A', 'T', 'C', 'G', 'N'};
    
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
                if (s[c] == 'A')
                {
                    motifs[motif_idx][c] = 0;
                }
                else if (s[c] == 'T')
                {
                    motifs[motif_idx][c] = 1;
                }
                else if (s[c] == 'C')
                {
                    motifs[motif_idx][c] = 2;
                }
                else if (s[c] == 'G')
                {
                    motifs[motif_idx][c] = 3;
                }
                else
                {
                    cout << "Error: Ambiguous character in motif" << "\n";
                    exit(-1);   
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
            if (s[c] == 'A')
            {
                genome[genome_idx] = 0;
            }
            else if (s[c] == 'T')
            {
                genome[genome_idx] = 1;
            }
            else if (s[c] == 'C')
            {
                genome[genome_idx] = 2;
            }
            else if (s[c] == 'G')
            {
                genome[genome_idx] = 3;
            }
            else
            {
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
    float **f_matrix = new float*[4];
    for (int i = 0; i < 4; i++) {
        f_matrix[i] = new float[motifs_info[0]];
    }

    // Make sure that frequency matrix has all 0s
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < motifs_info[0]; j++)
        {
            f_matrix[i][j] = 0;
        }
    }

    // Take each sequence and add to the frequency matrix
    for (int i = 0; i < motifs_info[1]; i++)
    {
        frequency_matrix(motifs[i], motifs_info[0], f_matrix);
    }

    // Print frequency matrix to output
    cout << "Frequency matrix:\n";
    for (int i = 0; i < 4; i++)
    {
        cout << bp[i] << ": ";
        for (int j = 0; j < motifs_info[0]; j++)
        {
            cout << f_matrix[i][j] << " ";
        }
        cout << "\n\n";
    }
 
    // Add pseudocounts to the frequency matrix
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < motifs_info[0]; j++)
        {
            f_matrix[i][j] += 0.25;
        }
    }

    // Convert f_matrix to probability matrix
    convert_matrix(f_matrix, motifs_info[0]);
    //Print probability matrix
    cout << "Probability matrix:\n";
    cout << fixed << setprecision(3);
    for (int i = 0; i < 4; i++)
    {
        cout << bp[i] << ": ";
        for (int j = 0; j < motifs_info[0]; j++)
        {
            cout << f_matrix[i][j] << " ";
        }
        cout << "\n\n";
    }
       
    // Probability matrix to PSSM
    pssm(f_matrix, motifs_info[0], gc, at);
    cout << "PSSM:\n";
    for (int i = 0; i < 4; i++)
    {
        cout << bp[i] << ": ";
        for (int j = 0; j < motifs_info[0]; j++)
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
            cout << low << "\t" << low + motifs_info[0] - 1 << "\t" << "+\t";
            for (int c = 0; c < motifs_info[0]; c++) 
            {
                cout << bp[genome_segment[c]];
            }
            cout << "\t" << tmp_score << "\n";
        }
        
        // Get score of complement sequence
        tmp_score = get_score(f_matrix, complement_segment, motifs_info[0], 1 - (at + gc));
        if (tmp_score >= cutoff) {
            cout << low << "\t" << low + motifs_info[0] - 1 << "\t" << "-\t";
            for (int c = 0; c < motifs_info[0]; c++) 
            {
                cout << bp[complement_segment[c]];
            }
            cout << "\t" << tmp_score << "\n";
        }
    }
    
    // Free array holding genome
    delete[] genome;

    //Free each sub-array
    for(int i = 0; i < 4; i++) {
        delete[] f_matrix[i];
    }
    //Free the array of pointers
    delete[] f_matrix;
    
    return 0;
}