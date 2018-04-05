This code was developed under Julia v0.5 by Ali Hassan in collaboration with Robert Mieth and Deepjyoti Deka. 

The following packages must be installed:

  - Distributions
  - Ipopt
  - JuMP
  - JuMPChance
  - DataFrames
  - Requests

To run the code, execute main.jl, or include() it from a Julia prompt.

The input data is given as follows:

  - File Input.jl load the system data from the files generators.csv, lines.csv
    and nodes.csv. Additional system data is specified in lines 19-21. 
  - File main.jl contains the input data for TCL ensembles (MDP) in lines 23-27,
    43-57.

The Spatio-Temporal Dual Decompostion (ST-D2) algorithm is formulated in the file main.jl.
The algorithm solves model for TCl ensemble and CC-OPF by calling functions from the files
mdp.jl and opf.jl respectively. The lagrange multipliers are updated in the file main.jl.
 
