# This file includes functions to build and run the T-SCUC and R-SCUC simulation for modified IEEE 24-bus case
# Author: Jin Lu (University of Houston)

# Import
from UC_simplified_function import *

### Case 1: run the T-SCUC for the modified IEEE 24-bus case
def run_T_SCUC(wd_lv):  # wd_lv is the wind penetration level
    ## In IEEE 24-bus case, the generators which costs are near 0 is set to wind generators
    renew_gen_list = [15,25,26,27,28,29,30] # gen_num for wind gens
    # read test case
    case_name = 'case24.mat'
    case_inst = GenericModel(case_name)
    # adjust line capacity
    for l in range(case_inst.branchtotnum):
        case_inst.branch.rateA[l] = case_inst.branch.rateA[l] * 0.7 # adjustment: line capacity is 70% of the origin value
    # origin generators in the test case: fuel type set to 'traditional'
    case_inst.gen.fuel_type = []
    for g in range(case_inst.gentotnum):
        case_inst.gen.fuel_type.append('traditional')
    # wind power plants in the test case: fuel type set to 'Wind'
    for g_num in renew_gen_list:
        case_inst.gen.fuel_type[g_num-1] = 'Wind'
    # gen cost for original generators
    # the polynomial gencost data => linear gencost used in UC
    case_inst.gen.c0 = []
    case_inst.gen.c1 = []
    case_inst.gen.c_su = []
    for g in range(case_inst.gentotnum):
        case_inst.gen.c0.append(case_inst.gencost_array[g,6])   # c2*(p^2)+c1*p+c0: use c0 as no-load cost
        case_inst.gen.c1.append(case_inst.gencost_array[g,5])   # c2*(p^2)+c1*p+c0: use c1 (ignore the c2 term)
        case_inst.gen.c_su.append(case_inst.gencost.startup[g])
    # form the ramping rate data
    for g in range(case_inst.gentotnum):
        case_inst.gen.ramp_agc[g] = case_inst.gen.Pmax[g]/60  # ramping rate in MW/min, assume ramping rate is large

    ## calculate the wind penetration level and may scale the wind power based on wind penetration level
    load_dl = np.loadtxt('load_dl.txt') # the load profile is saved in 'load_dl.txt'
    load_tot = 0
    for t in range(24):
        for b in range(case_inst.bustotnum):
            load_tot += load_dl[b,t]
    pwind_tot = 0
    wind_hrly = np.loadtxt('wind_hrly.txt')  # the wind profile is saved in 'wind_hrly.txt'
    for t in range(24):
        for wg in range(wind_hrly.shape[0]):
            pwind_tot += wind_hrly[wg][t]
    wd_lv_origin = pwind_tot/load_tot   # origin penetration level is about 7%
    # scale the wind profile to match the input penetration level
    for t in range(24):
        for wg in range(wind_hrly.shape[0]):
            wind_hrly[wg][t] = wind_hrly[wg][t]*wd_lv/wd_lv_origin
    np.savetxt('wind_hrly_adj.txt', wind_hrly, fmt='%.2f')  # [wgen][h]

    ## first run of UC
    UC_case_r1 = build_UC_pyomo(case_inst,load_dl)
    solve_UC(UC_case_r1,'TSCUC_solution.dat')
    # write_UCresult_Texas(UC_case_r1, 'UCsolution_case2.dat')
    ## Second run of UC
    # fixed the u_g_t, solve again and obtain dual variable (electricity price)
    UC_case_r2 = build_UC_pyomo_Run2(case_inst, load_dl, UC_case_r1)
    solve_UC(UC_case_r2,"TSCUC_solution_r2.dat")
    # write_UCresult(UC_case_r2,'Nhdg')
    write_UCresult_wdlv(UC_case_r2,'T_SCUC',wd_lv)

### Case 2: run the R-SCUC for the modified IEEE 24-bus case
def run_R_SCUC(wd_lv):
    ## In IEEE 24-bus case, the generators which costs are near 0 is set to wind generators
    renew_gen_list = [15,25,26,27,28,29,30] # gen_num for renew gens
    # read test case
    case_name = 'case24.mat'
    case_inst = GenericModel(case_name)
    # adjust line capacity
    for l in range(case_inst.branchtotnum):
        case_inst.branch.rateA[l] = case_inst.branch.rateA[l] * 0.7
    # origin generators in the test case, fuel type set to 'traditional'
    case_inst.gen.fuel_type = []
    for g in range(case_inst.gentotnum):
        case_inst.gen.fuel_type.append('traditional')
    # wind power plants in the test case, fuel type set to 'Wind'
    for g_num in renew_gen_list:
        case_inst.gen.fuel_type[g_num-1] = 'Wind'
    # gen cost for original generators
    # the polynomial gencost data => linear gencost used in UC
    case_inst.gen.c0 = []
    case_inst.gen.c1 = []
    case_inst.gen.c_su = []
    for g in range(case_inst.gentotnum):
        case_inst.gen.c0.append(case_inst.gencost_array[g,6])   # c2*(p^2)+c1*p+c0: use c0 as no-load cost
        case_inst.gen.c1.append(case_inst.gencost_array[g,5])   # c2*(p^2)+c1*p+c0: use c1 (ignore the c2 term)
        case_inst.gen.c_su.append(case_inst.gencost.startup[g])
    # form the ramping rate data
    for g in range(case_inst.gentotnum):
        case_inst.gen.ramp_agc[g] = case_inst.gen.Pmax[g]/60  # ramping rate in MW/min, assume ramping rate is large
    ## calculate the wind penetration level and may scale the wind power based on the wind penetration level
    load_dl = np.loadtxt('load_dl.txt')
    load_tot = 0
    for t in range(24):
        for b in range(case_inst.bustotnum):
            load_tot += load_dl[b,t]
    pwind_tot = 0
    wind_hrly = np.loadtxt('wind_hrly.txt')  # read the renewable profiles created in case 2 [rgen][hr]
    for t in range(24):
        for wg in range(wind_hrly.shape[0]):
            pwind_tot += wind_hrly[wg][t]
    wd_lv_origin = pwind_tot/load_tot   # origin penetration level is about 7%
    # scale the wind profile to match the input penetration level
    for t in range(24):
        for wg in range(wind_hrly.shape[0]):
            wind_hrly[wg][t] = wind_hrly[wg][t]*wd_lv/wd_lv_origin
    np.savetxt('wind_hrly_adj.txt', wind_hrly, fmt='%.2f')  # [wgen][h]

    ## first run of UC
    UC_case_r1 = build_UC_PFrelax(case_inst,load_dl)
    solve_UC(UC_case_r1,'RSCUC_solution.dat')
    # write_UCresult_Texas(UC_case_r1, 'UCsolution_PFrelax.dat')
    ## Second run of UC
    # fixed the u_g_t, solve again and obtain dual variable (electricity price)
    UC_case_r2 = build_UC_PFrelax_Run2(case_inst, load_dl, UC_case_r1)
    solve_UC(UC_case_r2,"RSCUC_solution_r2.dat")
    # write_UCresult(UC_case_r2,'Nhdg_PFrelax')
    write_UCresult_wdlv(UC_case_r2,'R_SCUC',wd_lv)


### Example: run the R-SCUC case for 30% wind penetration level
if __name__=='__main__':
    run_R_SCUC(wd_lv=0.3)