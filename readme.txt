The .m files are associated with figure numbers in "Using a Cell Flux Model to Investigate Patterns of
Phytoplankton Growth and Macromolecular Allocation Under Iron Limitation"

All codes use optimized values of afe and A_-^Fe_pho, which were determined using the attached python codes (opt_subplots.py and dependencies) and saved to the excel file "optimized_params"

specific growth rates, elemental stoichiometry, and functional/macromolecular allocation are calculated using the function CFM_Fe, which requires DSolver and Qc_essential_computation to run

Toolboxes/functions needed to generate figures:
tight_subplot needed to generate all figures (https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w)
m_map is needed to generate figure 5. (https://www-old.eoas.ubc.ca/~rich/map.html)
* All of those 'addpath(genpath(' lines at the start of each .m file are just directing MATLAB to look through my files and find those toolboxes/functions

Questions about code/output should be directed to mbernish@gmail.com
