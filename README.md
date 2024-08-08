# RandomizedGreedyMPS
Code associated with the publication "Randomized greedy magic point selection schemes for nonlinear model reduction" by R. Zimmermann/Kai Cheng

The code implements 
* a full-order model
* a POD-based reduced order model
* a POD-DEIM based reduced order model
for the nonlinear Schrödinger equation.

The theoretical background is given in Section 3.3. of the associated paper.

To execute the code, go to folder "Matlab/NonLinSchroedinger" and execute the script
">> script_nonlinSchroedinger_PODDEIM.m".



# Detailed instructions: 
Open the script "script_nonlinSchroedinger_PODDEIM.m" and read the topmost comment block.
 
 1) First run the script in FOM mode (set model =1):
    ">> script_nonlinSchroedinger_PODDEIM"
    This will execute the full-order model and in particular,
    the snapshots will be generated.
    The reduced order models rely on these snapshots.

 2) Change the model to POD (2) or DEIM (3) and set the other user
    parameters according to your preferences.
    The default set up for reproducing the results from the paper is:
    * tstep_scheme  = 2;
    * MPE_mode      = 1; %1:onlinefast MPE, 2:leverage score sampling
    * r             = 40;
    * rDEIM         = 40;
    * nr_points     = 2*rDEIM;
    * nr_randruns   = 100;
    
    Then rerun 
    ">> script_nonlinSchroedinger_PODDEIM"
    
   This reproduces the first column of the table "TableResults.txt".
   Be aware that the script will run the ROM 100 times to compute
   performance statistics. This may take a while.
