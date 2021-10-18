# estimate_spatial_information
In this repository, we use two methods to correct the upward bias in the naive calculation of information content, seeking a bias-free estimation of the information carried by neuronal activity about a given encoded variable.


## Implementation
The repository includes two methods for correcting the bias in the naive calculation of information content for limited sample sizes:

**1. The scaled shuffle reduction (SSR) method** - estimating the ratio between the bias in the naive and shuffle information and subtracting from the naïve information the shuffle information scaled by the estimated bias ratio.

**2. The bounded asymptotic extrapolation (BAE) method** - fitting the function of how the naive information changes with sample size and extrapolating it to infinity.

The naive estimation of information content is based either on the Skaggs information index (SI) or on the mutual information (MI).

## Usage and documentation
Scripts are provided in the *Scripts* directory.
Sample data sets for both real and simulated data are provided in the *sample data* directory.

The full data sets for both real and simulated data are provided in the *Full data set* directory.

To perform the analysis on the real or simulated sample data, use the *demo_real_data.m* or *demo_simulated_data.m* script.  
Before running the script, change the *data_pathway* to the *sample data* directory on your computer.

## Inputs
1. spike_train - a matrix of size TxN, where T is the number of time bins and N is the number of neurons. Each element is the spike count of a given neuron in a given time bin.
2. stimulus_trace - a vector of size T, where each element is the stimulus value in a given time bin.

## Outputs

1. General parameters (dt, number of shuffles, subsampling repetitions, etc.)
2. Data statistics (average rates, average active time bins, fraction of significant cells, etc.)
3. Estimated information (SI and MI based on the SSR and BE methods)

## Additional comments
The bias correction methods can be modified to accomodate additional information theoretic measures.
The obtained results of the SSR and BAE methods are compared against shuffle reduction (SR) and asymptotic extrapolation (AE). 

