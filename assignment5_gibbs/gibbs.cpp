#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <bits/stdc++.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// Max iterations to avoid infinite loops
#define MAX_ITER 2000
// Max iterations until plateau
#define ITER 400
//Max number of seeds
#define SEEDS 20
// Number of adjustments
#define NUM_ADJ 5

using namespace std;

int *scan_data(char *filepath)
{
    string s;
    static int seq_data[2] = {0, 0};
    int seq_length = 0;
    fstream fin(filepath);
    while (getline(fin, s))
    {
        if (s[0] == '>')
        {
            // If no length has been recorded, put it into the array to return
            if (seq_data[1] == 0)
            {
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
    for (int c = 0; c < s.size(); c++)
    {
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

void frequency_matrix(int **sequences, int motif_locations[], int data_info[], int motif_size, float **f_matrix, int exclude_idx)
{
    for (int i = 0; i < data_info[0]; i++)
    {
        if (i != exclude_idx)
        {
            for (int j = 0; j < motif_size; j++)
            {
                if (sequences[i][motif_locations[i] + j] > 3)
                {
                    cout << "Error: Ambiguous character in motif" << j << "\n";
                    exit(-3);
                }
                else
                {
                    f_matrix[sequences[i][motif_locations[i] + j]][j]++;
                }
            }
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
    float background_probs[4] = {(1 - gc_content) / 2, (1 - gc_content) / 2, gc_content / 2, gc_content / 2};

    for (int row = 0; row < 4; row++)
    {
        for (int col = 0; col < motif_size; col++)
        {
            f_matrix[row][col] = log2(f_matrix[row][col] / background_probs[row]);
        }
    }
}

float get_score(float **f_matrix, int sequence[], int low, int high)
{
    float score = 0.0;
    int matrix_idx = 0;
    for (int i = low; i < high; i++)
    {
        score += f_matrix[sequence[i]][matrix_idx];
        matrix_idx++;
    }
    return score;
}

double objective_fx(int **sequences, int motif_locations[], int data_info[], int m, float gc_values[])
{
    // Create a vector containing 4 vectors (ATCG) each being the size of the motif
    float **f_matrix = new float *[4];
    for (int i = 0; i < 4; i++)
    {
        f_matrix[i] = new float[m];
    }

    double total_score = 0.0;
    // Loop through sequences, forget sequences[i], add motif score to total
    for (int i = 0; i < data_info[0]; i++)
    {
        // Make sure frequency matrix is zeroed out
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < m; k++)
            {
                f_matrix[j][k] = 0;
            }
        }
        // Fill freq_matrix with sequences except for sequences[exclude_idx]
        frequency_matrix(sequences, motif_locations, data_info, m, f_matrix, i);
        // Add pseudocounts to the frequency matrix
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < m; k++)
            {
                f_matrix[j][k] += 0.25;
            }
        }
        // Convert f_matrix to probability matrix
        convert_matrix(f_matrix, m);
        // Probability matrix to PSSM
        pssm(f_matrix, m, gc_values[i]);
        // Get score at current motif_locations[i]
        total_score += get_score(f_matrix, sequences[i], 0, m);
    }

    // Cleaning
    //Free each sub-array
    for (int i = 0; i < 4; i++)
    {
        delete[] f_matrix[i];
    }
    //Free the array of pointers
    delete[] f_matrix;

    return total_score;
}

void run_pssm(int **sequences, int data_info[], int motif_locations[], int m, float gc_values[])
{
    // Create a vector containing 4 vectors (ATCG) each being the size of the motif
    float **f_matrix = new float *[4];
    for (int i = 0; i < 4; i++)
    {
        f_matrix[i] = new float[m];
    }

    for (int exclude_idx = 0; exclude_idx < data_info[0]; exclude_idx++)
    {
        // Make sure frequency matrix is zeroed out
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < m; j++)
            {
                f_matrix[i][j] = 0;
            }
        }
        // Fill freq_matrix with sequences except for sequences[exclude_idx]
        frequency_matrix(sequences, motif_locations, data_info, m, f_matrix, exclude_idx);

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
        // Probability matrix to PSSM
        pssm(f_matrix, m, gc_values[exclude_idx]);

        // Placeholder for scores; Use length of sequence as length of array
        float score_array[data_info[1]];
        int index_arr[data_info[1]];
        // Where to put score in score_array if > 0
        int score_idx = 0;
        // High score at index and placeholder for index
        float high_score = 0.0;
        int high_index = 0;
        // Score at index
        float tmp_score = 0.0;
        // Sum of scores > 0
        float total_score = 0.0;
        // Loop through entire genome
        for (int low = 0; low < (data_info[1] - m); low++)
        {
            // Get score at position `low`
            tmp_score = get_score(f_matrix, sequences[exclude_idx], low, low + m);
            // Eliminate negative scores...
            if (tmp_score > 0)
            {
                total_score += tmp_score;
                // store current total score in score array
                score_array[score_idx] = total_score;
                // store the index of the sequence that produced the score above
                index_arr[score_idx] = low;
                score_idx++;
            }
        }

        // Normalize scores
        for (int i = 0; i < score_idx; i++)
        {
            score_array[i] = score_array[i] / total_score;
        }

        int new_index = 0;
        //Generate random probability to select new index
        float random_prob = (double)rand() / (RAND_MAX);
        // Loop through probs until I get value lower than random_prob
        while (random_prob >= score_array[new_index])
        {
            new_index++;
        }

        // Assign new location for motif_location[exclude_idx]
        motif_locations[exclude_idx] = index_arr[new_index];
    }
    // Cleaning
    //Free each sub-array
    for (int i = 0; i < 4; i++)
    {
        delete[] f_matrix[i];
    }
    //Free the array of pointers
    delete[] f_matrix;
}

void adjuster(int **sequences, int motif_locations[], int data_info[], int &m, float gc_values[])
{
    double optimum = -999999.9;
    double current_score = 0.0;
    int tmp_locations[data_info[0]];
    int adjusted_locations[data_info[0]];
    int tmp_m = m;

    // Left side adjustment
    for (int i = -1; i < 2; i++)
    {
        for (int j = 0; j < data_info[0]; j++)
        {
            adjusted_locations[j] = motif_locations[j] + i;
            if ((adjusted_locations[j] + m - 1) > data_info[1])
            {
                break;
            }
        }       
        // Current score from objective_fx
        current_score = objective_fx(sequences, adjusted_locations, data_info, m, gc_values);
        if (current_score > optimum)
        {
            optimum = current_score;
            //tmp_m = tmp_m - i;
            for (int j = 0; j < data_info[0]; j++)
            {
                tmp_locations[j] = adjusted_locations[j];
            }
        }
    }
    
    // Right side adjustment
    int tmp = -1;
    bool oneAtEnd = false;
    while(tmp < 2 || oneAtEnd == false)
    {
        for (int j = 0; j < data_info[0]; j++)
        {
            if ((motif_locations[j] + m - 1) + tmp)
            {
                oneAtEnd = true;
            }
        }

        // Current score from objective_fx
        current_score = objective_fx(sequences, motif_locations, data_info, m + tmp, gc_values);

        if (current_score > optimum)
        {
            optimum = current_score;
            //tmp_m = m + tmp;
        }
        tmp++;
    }
    /*
    //Both sides
    tmp = -1;
    oneAtEnd = false;
    while(tmp < 2 || oneAtEnd == false)
    {
        // Don't need to do 0
        if (tmp != 0)
        {
            for (int j = 0; j < data_info[0]; j++)
            {
                adjusted_locations[j] = motif_locations[j] + tmp;
                if ((adjusted_locations[j] + m - 1) + (m + tmp))
                {
                    oneAtEnd = true;
                }
            }
        
            // Current score from objective_fx
            current_score = objective_fx(sequences, tmp_locations, data_info, m + tmp, gc_values);

            if (current_score > optimum)
            {
                optimum = current_score;
                //tmp_m = m - i;
                for (int j = 0; j < data_info[0]; j++)
                {
                    tmp_locations[j] = adjusted_locations[j];
                }
            }
        }
    }
    
    tmp = -1;
    oneAtEnd = false;
    while(tmp < 2 || oneAtEnd == false)
    {
        // Don't need to do 0
        if (tmp != 0)
        {
            for (int j = 0; j < data_info[0]; j++)
            {
                adjusted_locations[j] = motif_locations[j] + tmp;
                if ((adjusted_locations[j] + m - 1) + (m - tmp))
                {
                    oneAtEnd = true;
                }
            }
        
            // Current score from objective_fx
            current_score = objective_fx(sequences, tmp_locations, data_info, m - tmp, gc_values);

            if (current_score > optimum)
            {
                optimum = current_score;
                tmp_m = m - (2 * tmp);
                for (int j = 0; j < data_info[0]; j++)
                {
                    tmp_locations[j] = adjusted_locations[j];
                }
            }
        }
    }
    */
    //m = tmp_m;
    for (int j = 0; j < data_info[0]; j++)
    {
        motif_locations[j] = tmp_locations[j];
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
    srand(time(NULL));
    // Assign variable to length of motifs
    int m = atoi(argv[2]);
    char bp[] = {'A', 'T', 'C', 'G', 'N'};
    // Get the number of sequences in the data file
    int *data_info;
    data_info = scan_data(argv[1]);
    // Create integer arrays to hold sequences
    int **sequences = new int *[data_info[0]];
    for (int i = 0; i < data_info[0]; i++)
    {
        sequences[i] = new int[data_info[1]];
    }
    // Create an array for sequence labels
    vector<string> labels;

    // Read file and put sequences into sequence array
    string line;
    string seq;
    int seq_idx = 0;
    fstream fin(argv[1]);
    while (getline(fin, line))
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
        else
        {
            labels.push_back(line);
        }
    }

    // Get all GC content values for sequences
    float gc_values[data_info[0]];
    for (int i = 0; i < data_info[0]; i++)
    {
        gc_values[i] = gc_content(sequences[i], data_info[1]);
    }

    int num_seeds = 0;
    double optimum = -999999.0;
    int optimum_locations[data_info[0]];
    int motif_locations[data_info[0]];
    while (num_seeds < SEEDS)
    {
        for (int i = 0; i < data_info[0]; i++)
        {
            // Random index from 0 to (length of sequence - length of motifs)
            motif_locations[i] = rand() % (data_info[1] - m);
        }

        int iterations = 0;
        int max_iterations = 0;
        int until_adjust = 0;
        double seed_optimum = -9999999.0;
        double current_score = 0.0;
        while (iterations < ITER || max_iterations < 2000)
        {
            if (until_adjust > NUM_ADJ)
            {
                // Run adjustment step
                adjuster(sequences, motif_locations, data_info, m, gc_values);
            }
            // Run gibbs sampler for all sequences and reassign motif_locations
            run_pssm(sequences, data_info, motif_locations, m, gc_values);

            // Current score from objective_fx
            current_score = objective_fx(sequences, motif_locations, data_info, m, gc_values);
            // Check if current_score is higher than optimum
            if (current_score > seed_optimum)
            {
                seed_optimum = current_score;
                iterations = 0;
            }
            iterations++;
            until_adjust++;
            max_iterations++;
        }
        // Check if seed optimum is greater than total optimum
        if (seed_optimum > optimum)
        {
            // Set new optimum
            optimum = seed_optimum;
            // Set optimum motif locations
            for (int i = 0; i < data_info[0]; i++)
            {
                optimum_locations[i] = motif_locations[i];
            }
        }
        num_seeds++;
    }

    // Print output
    for (int i = 0; i < data_info[0]; i++)
    {
        for (int j = 0; j < m; j++)
        {
            cout << bp[sequences[i][optimum_locations[i] + j]];
        }
        cout << " " << labels[i] << " " << optimum_locations[i] << "\n";
    }
    cout << optimum << "\n";

    //Cleaning house
    //Free each sub-array
    for (int i = 0; i < data_info[0]; i++)
    {
        delete[] sequences[i];
    }
    //Free the array of pointers
    delete[] sequences;

    return 0;
}