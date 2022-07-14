# This file includes functions to build and run the EH-SCUC and H-SCUC model for modified IEEE 24-bus case
# Author: Jin Lu (University of Houston)

# Import
from UC_simplified_function import *
from UC_simp_hdg import *

### Case 3: run the EH-SCUC for the modified IEEE 24-bus case
### 2 electrolyzer& 2 fuel cells added
def run_EH_SCUC(wd_lv):   # wd_lv is wind penetration level
    ### data for electrolyzer & fuel cell
    elctro_bus = [14, 22]  # the buses that have wind power plants
    fc_bus = [14, 22]  # the buses are same as electrolyzer's buses
    elctro_Pmax = [100,300] # electrolyzer maximum power rating, match the renewable capacity at the location
    fc_Pmax = [100,300] # total fuel cell capacity is larger than electrolyzer total capacity, because the load at one location may need more renewable energy than other locations
    renew_gen_list = [15, 25, 26, 27, 28, 29, 30]  # gen_num for renew gens
    # read test case
    case_name = 'case24.mat'
    case_inst = GenericModel(case_name)
    # adjust line capacity
    for l in range(case_inst.branchtotnum):
        case_inst.branch.rateA[l] = case_inst.branch.rateA[l]*0.7
    ## In IEEE 24-bus case, the generators which costs are near 0 is set to wind generators
    # Assign generator fuel type
    case_inst.gen.fuel_type = []
    for g in range(case_inst.gentotnum):
        case_inst.gen.fuel_type.append('traditional')
    for g_num in renew_gen_list:
        case_inst.gen.fuel_type[g_num-1] = 'Wind'
    # gen cost for original generators
    # the polynomial gencost data => linear gencost used in UC
    case_inst.gen.c0 = []
    case_inst.gen.c1 = []
    case_inst.gen.c_su = []
    for g in range(case_inst.gentotnum):
        case_inst.gen.c0.append(case_inst.gencost_array[g, 6])  # c2*(p^2)+c1*p+c0: use c0 as no-load cost
        case_inst.gen.c1.append(case_inst.gencost_array[g, 5])  # c2*(p^2)+c1*p+c0: use c1 (ignore the c2 term)
        case_inst.gen.c_su.append(case_inst.gencost.startup[g])
    # form the ramping rate data
    for g in range(case_inst.gentotnum):
        case_inst.gen.ramp_agc[g] = case_inst.gen.Pmax[g] / 60  # ramping rate in MW/min, assume ramping rate is large
    ## calculate the wind penetration level and may scale the wind power based on wind penetration level
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
    UC_hdg_r1 = build_EHSCUC_pyomo(case_inst, load_dl, elctro_bus, elctro_Pmax, fc_bus, fc_Pmax)
    solve_UC(UC_hdg_r1, 'EHSCUC_solution.dat')
    # write_UCresult_Texas(UC_hdg_r1,'UCsolution_hdg_nHN.dat')
    ## Second run of UC
    # fixed the u_g_t, solve again and obtain dual variable (electricity price)
    UC_hdg_r2 = build_EHSCUC_pyomo_run2(case_inst, load_dl, elctro_bus, elctro_Pmax, fc_bus, fc_Pmax, UC_hdg_r1)
    solve_UC(UC_hdg_r2, "EHSCUC_solution_r2.dat")
    # write_UCresult(UC_hdg_r2,'Hdg_nHN_c1')  # save to 'UC_results' folder
    write_UCresult_wdlv(UC_hdg_r2,'EH_SCUC',wd_lv)

### Case 4: run the H-SCUC for the modified IEEE 24-bus case
### 2 electrolyzer & 2 fuel cells added, electrolyzer and fuel cell are not at the same location
def run_H_SCUC(wd_lv):   # wd_lv is wind penetration level
    ### data for electrolyzer & fuel cell
    elctro_bus = [14, 22]  # the buses that have renewable power plants
    fc_bus = [13, 15]  # the buses that peak load over 200MW   # fuel cell location is different
    # fc_bus = [9, 10]  # this line is used to test different fuel cell location
    elctro_Pmax = [100,300] # electrolyzer maximum power rating, match the renewable capacity at the location
    fc_Pmax = [100,300] # total fuel cell capacity is larger than electrolyzer total capacity, because the load at one location may need more renewable energy than other locations
    renew_gen_list = [15, 25, 26, 27, 28, 29, 30]  # gen_num for renew gens
    # read test case
    case_name = 'case24.mat'
    case_inst = GenericModel(case_name)
    # adjust line capacity
    for l in range(case_inst.branchtotnum):
        case_inst.branch.rateA[l] = case_inst.branch.rateA[l]*0.7
    ## In IEEE 24-bus case, the generators which costs are near 0 is set to wind generators
    # origin generators in the test case fuel type set to 'traditional'
    case_inst.gen.fuel_type = []
    for g in range(case_inst.gentotnum):
        case_inst.gen.fuel_type.append('traditional')
    for g_num in renew_gen_list:
        case_inst.gen.fuel_type[g_num-1] = 'Wind'
    # gen cost for original generators
    # the polynomial gencost data => linear gencost used in UC
    case_inst.gen.c0 = []
    case_inst.gen.c1 = []
    case_inst.gen.c_su = []
    for g in range(case_inst.gentotnum):
        case_inst.gen.c0.append(case_inst.gencost_array[g, 6])  # c2*(p^2)+c1*p+c0: use c0 as no-load cost
        case_inst.gen.c1.append(case_inst.gencost_array[g, 5])  # c2*(p^2)+c1*p+c0: use c1 (ignore the c2 term)
        case_inst.gen.c_su.append(case_inst.gencost.startup[g])
    # form the ramping rate data
    for g in range(case_inst.gentotnum):
        case_inst.gen.ramp_agc[g] = case_inst.gen.Pmax[g] / 60  # ramping rate in MW/min, assume ramping rate is large
    ## calculate the wind penetration level and may scale the wind power based on the wind penetration level
    load_dl = np.loadtxt('load_dl.txt')
    load_tot = 0
    for t in range(24):
        for b in range(case_inst.bustotnum):
            load_tot += load_dl[b,t]
    pwind_tot = 0
    wind_hrly = np.loadtxt('wind_hrly.txt')  # read the renewable profiles [rgen][hr]
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
    UC_hdg_r1 = build_HSCUC_pyomo(case_inst, load_dl, elctro_bus, elctro_Pmax, fc_bus, fc_Pmax)
    solve_UC(UC_hdg_r1, 'HSCUCsolution.dat')
    ## Second run of UC
    # fixed the u_g_t, solve again and obtain dual variable (electricity price)
    UC_hdg_r2 = build_HSCUC_pyomo_run2(case_inst, load_dl, elctro_bus, elctro_Pmax, fc_bus, fc_Pmax, UC_hdg_r1)
    solve_UC(UC_hdg_r2, "HSCUC_solution_r2.dat")
    write_UCresult_wdlv(UC_hdg_r2,'H_SCUC',wd_lv)

### Example: run the H-SCUC for 30% wind penetration level
if __name__=='__main__':
    run_H_SCUC(wd_lv=0.5)