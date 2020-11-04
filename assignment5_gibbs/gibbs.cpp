#include <iostream>
#include <fstream>
#include <string>
#include <bits/stdc++.h>
#include <stdlib.h>
#include <time.h>


using namespace std;

int * scan_data(char *filepath)
{
    string s;
    static int seq_data[2] = { 0, 0 };
    int seq_length = 0;
    fstream fin(filepath);
    while(getline(fin, s))
    {
        if (s[0] == '>')
        {
            // If no length has been recorded, put it into the array to return
            if (seq_data[1] == 0) {
                seq_data[1] = seq_length;
            }
            // Otherwise, check current length against what has been recorded
            else if (seq_length != seq_data[1]) 
            {
                cout << "Error: Lengths of given sequences are not equal\n";
                exit(-1);
            }
            seq_length = 0;
            seq_data[0]++;
        }
        else
        {
            seq_length += s.size();            
        }
    }
    return seq_data;
}

void convert_sequence(string s, int arr[])
{
    transform(s.begin(), s.end(), s.begin(), ::toupper);
    for (int c = 0; c < s.size(); c++) {
        if (s[c] == 'A')
        {
            arr[c] = 0;
        }
        else if (s[c] == 'T')
        {
            arr[c] = 1;
        }
        else if (s[c] == 'C')
        {
            arr[c] = 2;
        }
        else if (s[c] == 'G')
        {
            arr[c] = 3;
        }
        else
        {
            arr[c] = 4;   
        }
    }
}

float gc_content(int sequence[], int sequence_size)
{
    // Calculate GC content
    float gc = 0.0;
    for (int c = 0; c < sequence_size; c++)
    {
        if (sequence[c] == 2 || sequence[c] == 3)
        {
            gc++;
        }
    }
    return gc / sequence_size;
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

void pssm(float **f_matrix, int motif_size, float gc_content)
{
    // Background nucleotide probabilities that correspond to rows of frequency matrix
    float background_probs[4] = {1 - (gc_content / 2), 1 - (gc_content / 2), gc_content / 2, gc_content / 2};

    for (int row = 0; row < 4; row++)
    {
        for (int col = 0; col < motif_size; col++)
        {
            f_matrix[row][col] = log2(f_matrix[row][col] / background_probs[row]);
        }
    }
}

float get_score(float **f_matrix, int motif[], int motif_size)
{
    float score = 0.0;

    for (int i = 0; i < motif_size; i++)
    {
        score += f_matrix[motif[i]][i];
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

int main(int argc, char **argv)
{
    // Check if all required arguments are provided...
    if (argc == 1)
    {
        cout << "Unsupervised motif finding/Gibbs Sampler\n";
        cout << "The following arguments are required to be specified by the user:\n";
        cout << "\t1. Function call (e.g. `./gibbs`)\n";
        cout << "\t2. Path to data file (e.g. `data/some-data-file.fasta`)\n";
        cout << "\t3. Length of motif (e.g. `10`)\n";
        return 0;
    }

    // initialize random seed:
    srand (time(NULL));
    // Assign variable to length of motifs
    int m = atoi(argv[2]);
    char bp[] = {'A', 'T', 'C', 'G', 'N'};
    // Get the number of sequences in the data file
    int *data_info;
    data_info = scan_data(argv[1]);
    // Create integer arrays to hold sequences
    int sequences[data_info[0]][data_info[1]];
    // Create an array for sequence labels
    int labels[data_info[0]];

    // Read file and put sequences into sequence array
    string line;
    string seq;
    int seq_idx = 0;
    fstream fin(argv[1]);
    while(getline(fin, line)) 
    {
        if (line[0] != '>') 
        {
            seq += line;
            if (seq.size() == data_info[1])
            {
                convert_sequence(seq, sequences[seq_idx]);
                seq_idx++;
                seq = "";
            }
        }
    }

    // Get random indexes where motifs begin
    int motif_locations[data_info[0]];
    int motifs[data_info[0]][m];
    for (int i = 0; i < data_info[0]; i++)
    {
        // Random index chosen from 0 to (length of sequence - length of motifs)
        motif_locations[i] = rand() % (data_info[1] - m);
        for (int j = motif_locations[i]; j < motif_locations[i] + m; j++)
        {
            motifs[i][j - motif_locations[i]] = sequences[i][j];
        }
    }

    float gc = gc_content(sequences[0], data_info[1]);
    
    // Create a vector containing 4 vectors (ATCG) each being the size of the motif
    float **f_matrix = new float*[4];
    for (int i = 0; i < 4; i++) {
        f_matrix[i] = new float[m];
    }

    // Make sure that frequency matrix has all 0s
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < m; j++)
        {
            f_matrix[i][j] = 0;
        }
    }
    
    // Take each sequence and add to the frequency matrix
    for (int i = 1; i < 6; i++)
    {
        frequency_matrix(motifs[i], m, f_matrix);
    }

    // Print frequency matrix to output
    cout << "Frequency matrix:\n";
    for (int i = 0; i < 4; i++)
    {
        cout << bp[i] << ": ";
        for (int j = 0; j < m; j++)
        {
            cout << f_matrix[i][j] << " ";
        }
        cout << "\n\n";
    }

    for (int i = 1; i < data_info[0]; i++)
    {
        frequency_matrix(motifs[i], m, f_matrix);
    }

    // Add pseudocounts to the frequency matrix
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < m; j++)
        {
            f_matrix[i][j] += 0.25;
        }
    }

    // Convert f_matrix to probability matrix
    convert_matrix(f_matrix, m);
    //Print probability matrix
    cout << "Probability matrix:\n";
    cout << fixed << setprecision(3);
    for (int i = 0; i < 4; i++)
    {
        cout << bp[i] << ": ";
        for (int j = 0; j < m; j++)
        {
            cout << f_matrix[i][j] << " ";
        }
        cout << "\n\n";
    }

    // Probability matrix to PSSM
    pssm(f_matrix, m, gc);
    cout << "PSSM:\n";
    for (int i = 0; i < 4; i++)
    {
        cout << bp[i] << ": ";
        for (int j = 0; j < m; j++)
        {
            cout << f_matrix[i][j] << " ";
        }
        cout << "\n";
    }

    cout << "\n";

    //Training set scores
    float cutoff = 99999999.0;
    cout << "Training set motif scores: \n";
    for (int i = 0; i < data_info[0]; i++) {
        float tmp_score = get_score(f_matrix, motifs[i], m);
        for (int c = 0; c < m; c++) {
            cout << bp[motifs[i][c]];
        }
        cout << "\t\t" << tmp_score << "\n";
        if (tmp_score < cutoff)
        {
            cutoff = tmp_score;
        }
    }

    //Print header
    cout << "Matches with score " << cutoff << " or higher";
    cout << "Start\tEnd\tStrand\tSequence\t\tScore\n";
    
    // Length of motifs and a placeholder for score
    int current_sequence[m];
    float tmp_score;
    // Loop through entire genome and its complement
    for (int low = 0; low < data_info[1] - m; low++) {
        
        // Get segments to test
        segment_arr(current_sequence, low, low + m, sequences[0]);
        
        // Get score of genome segment
        tmp_score = get_score(f_matrix, current_sequence, m);
        if (tmp_score >= 0) {
            cout << low << "\t" << low + m - 1 << "\t" << "+\t";
            for (int c = 0; c < m; c++) 
            {
                cout << bp[current_sequence[c]];
            }
            cout << "\t" << tmp_score << "\n";
        }
    }
    /*
    */
    //Free each sub-array
    for(int i = 0; i < 4; i++) {
        delete[] f_matrix[i];
    }
    //Free the array of pointers
    delete[] f_matrix;
    return 0;
}