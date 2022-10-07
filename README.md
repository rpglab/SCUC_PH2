# SCUC_PH2
 Benefit Analysis of Hydrogen Networks to the Renewable Power Grids. SCUC - day-ahead scheduling - for a hybrid super energy grid consisting of power systems and hydrogen networks

This program implements a security-constrained unit commitment (SCUC) model for day-ahead scheduling of Super Energy GRIDs consisting of power systems and hydrogen networks. Three other benchmark SCUC models are also implemented for comparison.

The four SCUC models are as follows:
* Model 1: the traditional SCUC (T-SCUC) - Benchmark
* Model 2: line capacity relaxed SCUC (R-SCUC) - Benchmark
* Model 3: SCUC for power systems with local hydrogen energy hub (EH-SCUC) - Benchmark
* Model 4: SCUC for power systems coupled with hydrogen conversion and transmission system (H-SCUC) - The Proposed Model

The detailed formulations of these four SCUC models are explained in our paper listed below.

### Environment Setting:
* Python version: Python 3.8.13 - Other versions may also work but test is needed.
* Required packages: numpy (1.22.2), pyomo (6.4.1) and pypower (5.1.15). Other versions may also work but test is needed.
* A third party solver (e.g. glpk-4.65) is required to solve the problem with Pyomo.


### Main Program:
To run all SCUC models for different wind penetration levels: run python file 'run_all_SCUC.py'. You can also run selected SCUC model(s) by commenting out other models in the for loop.

### Input data:
* case24.mat: the IEEE 24-bus test case in matpower format 
* wind_hrly.txt: the hourly wind output power profile
* load_dl.txt: the hourly load profile

### Source Codes:
* power_mod.py: create a Class for matpower/pypower case
* formpyom_UC_simplified.py: function to form the pyomo data file for the SCUC
* UC_simplified_function.py: functions to build T-SCUC and R-SCUC model
* UC_simp_hdg.py: functions to build EH-SCUC and H-SCUC model
* UC_nohdg_cases.py: to run the T-SCUC and R-SCUC simulation for modified IEEE 24-bus case
* UC_hdg_cases.py: to run the EH-SCUC and H-SCUC simulation for modified IEEE 24-bus case
* run_all_SCUC.py: run all types of SCUC for different wind penetration level


### Test Power System Case:
The test case used for this work is modified from the following IEEE 24-bus test system:
* IEEE 24-bus system is one of the three areas of the IEEE 73-bus system. They are described in this reference: "The IEEE Reliability Test System-1996. A report prepared by the Reliability Test System Task Force of the Application of Probability Methods Subcommittee" and link is <a class="" target="_blank" href="https://ieeexplore.ieee.org/document/780914">here</a>. 



### Simulation Results:

The simulation results are (automatically) saved in folder 'UC_results'.
SubFolders are for different wind penetration levels.
For example, the folder for 30% wind penetration is 'wdlv_0.3'.
Each subfolder includes 24 files for all 4 SCUC models; each SCUC model has 6 result files explained as follows:

* 'Typename_SCUC_lmp.txt': results of locational marginal price ($/MWh).
			                     row is bus, column is time interval (24 hours).
* 'Typename_SCUC_Opcost.txt': results of total cost ($).
* 'Typename_SCUC_pf.txt': results of power flow (MW), 
								positive/negative indicates the direction of the power flow.
								row is line, column is time interval (24 hours).
* 'Typename_SCUC_pfpct.txt': results of pecentage of the power flow to the line capacity, 
								   positive/negative indicates the direction of the power flow
								   row is line, column is time interval (24 hours).
* 'Typename_SCUC_pgt.txt': results of generator active power output (MW)
								 row is generator, column is time interval (24 hours).
* 'Typename_SCUC_wcur.txt': results of wind curtailment (MW)
								  row is wind power plant, column is time interval (24 hours).

For example, total cost of EH-SCUC is in 'EH_SCUC_Opcost.txt'

Note that each subfolder needs to be created before runing the codes.

## Citation:
If you use these codes for your work, please cite the following paper:

Jin Lu and Xingpeng Li, “The Benefits of Hydrogen Energy Transmission and Conversion Systems to the Renewable Power Grids: Day-ahead Unit Commitment”, *54th North American Power Symposium*, Salt Lake City, UT, USA, Oct. 2022.

Paper website: <a class="off" href="/papers/JinLu-BnftAnlys-H2Grid-SCUC/"  target="_blank">https://rpglab.github.io/papers/JinLu-BnftAnlys-H2Grid-SCUC/</a>


## Contributions:
This program was developed by Jin Lu. Xingpeng Li supervised this work.


## Contact:
If you need any techinical support, please feel free to reach out to Jin Lu at jlu27@CougarNet.UH.EDU.

For collaboration, please contact Dr. Xingpeng Li at xli83@central.uh.edu.

Website: https://rpglab.github.io/


## License:
This work is licensed under the terms of the <a class="off" href="https://creativecommons.org/licenses/by/4.0/"  target="_blank">Creative Commons Attribution 4.0 (CC BY 4.0) license.</a>


## Disclaimer:
The author doesn’t make any warranty for the accuracy, completeness, or usefulness of any information disclosed and the author assumes no liability or responsibility for any errors or omissions for the information (data/code/results etc) disclosed.
