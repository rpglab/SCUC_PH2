# This file includes functions to build EH-SCUC, H-SCUC models
# Electrolyzer & FC efficiency is set in the functions
# Author: Jin Lu (University of Houston)


# Import
from power_mod import *
from pyomo.environ import *
from formpyom_UC_simplified import *
import numpy as np

#### FC & Electrolyzer same location (No Hydrogen Network)
### function to build the pyomo instance based on the test case
def build_EHSCUC_pyomo(case_inst, load_dl, elctro_bus, elctro_Pmax, fc_bus, fc_Pmax):
    # wind profile
    wind_hrly = np.loadtxt('wind_hrly_adj.txt')  # read the renewable profiles created in case 2 [rgen][hr]
    ## list of non wind or solar generators
    reg_gen_list = []
    nonreg_gen_list = []
    nonreg_gen_bus = []
    wgen_bus = []
    for i in range(case_inst.gentotnum):
        if case_inst.gen.fuel_type[i] != 'Wind':
            if case_inst.gen.fuel_type[i] != 'Solar':
                gen_num = i + 1
                reg_gen_list.append(gen_num)
        if case_inst.gen.fuel_type[i] == 'Wind':
            gen_num = i + 1
            nonreg_gen_list.append(gen_num)
            nonreg_gen_bus.append(case_inst.gen.bus[i])
            wgen_bus.append(case_inst.gen.bus[i])
        if case_inst.gen.fuel_type[i] == 'Solar':
            gen_num = i + 1
            nonreg_gen_list.append(gen_num)
            nonreg_gen_bus.append(case_inst.gen.bus[i])

    ## generate data file
    time_totnum = 24  # total number of time period
    pyomodata_UC_simplified('formpyomo_UC_simplified.dat',case_inst, time_totnum, load_dl, reg_gen_list)

    ### Python set
    # 24 hour time set
    Time_24 = []
    for i in range(24):
        Time_24.append(i + 1)
    # elctrolyzer set
    Elctro = []
    for i in range(len(elctro_bus)):
        Elctro.append(i+1)
    # fuel cell set
    FC = []
    for i in range(len(fc_bus)):
        FC.append(i+1)
    # renew gen set
    WGen = []
    for i in range(wind_hrly.shape[0]):
        WGen.append(i+1)

    ### The pyomo abstract model for UC
    model = AbstractModel()

    # Set
    model.BUS = Set()
    model.LINE = Set()
    model.GEN = Set()
    model.TIME = Set()
    model.BUS_TIME = Set(dimen=2)

    ## Param
    # Bus param
    model.bus_num = Param(model.BUS)
    # Gen param
    model.gen_num = Param(model.GEN)
    model.gen_bus = Param(model.GEN)
    model.gen_cost_P = Param(model.GEN)
    model.gen_cost_NL = Param(model.GEN)
    model.gen_cost_SU = Param(model.GEN)
    model.gen_Pmin = Param(model.GEN)
    model.gen_Pmax = Param(model.GEN)
    model.gen_r10 = Param(model.GEN)
    model.gen_rr = Param(model.GEN)
    # Line param
    model.line_num = Param(model.LINE)
    model.line_x = Param(model.LINE)
    model.line_Pmax = Param(model.LINE)
    model.line_fbus = Param(model.LINE)
    model.line_tbus = Param(model.LINE)
    # Time param
    model.time_num = Param(model.TIME)
    # Bus_Time param
    model.bustime_num = Param(model.BUS_TIME)
    model.load_b_t = Param(model.BUS_TIME)

    # Variable
    # Gen_Time Var
    model.p_g_t = Var(model.GEN, model.TIME)
    model.u_g_t = Var(model.GEN, model.TIME, domain=Binary)
    model.v_g_t = Var(model.GEN, model.TIME, domain=Binary)  # v_g_t is simplified to non-binary variable
    model.r_g_t = Var(model.GEN, model.TIME)
    # Line_Time Var
    model.theta_k_t = Var(model.LINE, model.TIME)
    model.p_k_t = Var(model.LINE, model.TIME)
    # Bus_Time Var
    model.theta_b_t = Var(model.BUS, model.TIME)
    # Electrolyzer Var
    model.p_elct_t = Var(Elctro, model.TIME,domain=NonNegativeReals)
    # fuel cell Var
    model.p_fc_t = Var(FC, model.TIME,domain=NonNegativeReals)
    # hydrogen tank storage
    model.E_elct_t = Var(Elctro,model.TIME,domain=NonNegativeReals)
    # Renew Curtailment Var
    model.pcur_wind = Var(WGen,model.TIME,domain=NonNegativeReals)


    ## Objective function
    def objfunction(model):
        obj = sum(
            model.gen_cost_P[g] * model.p_g_t[g, t] + model.gen_cost_NL[g] * model.u_g_t[g, t] + model.gen_cost_SU[g] *
            model.v_g_t[g, t] for g in model.GEN for t in model.TIME)
        return obj
    model.object = Objective(rule=objfunction, sense=minimize)

    # Generator initial u_g_t status
    def gen_Uinit_f(model, g):
        return model.u_g_t[g, 0] == 0
    model.gen_Uinit_cons = Constraint(model.GEN, rule=gen_Uinit_f)

    ## Generator power and reserve constraints
    # # P_g_t minimum constraint
    def gen_Pmin_f(model, g, t):
        return model.gen_Pmin[g] * model.u_g_t[g, t] <= model.p_g_t[g, t]

    model.gen_Pmin_cons = Constraint(model.GEN, model.TIME, rule=gen_Pmin_f)

    # P_g_t maximum constraint:
    def gen_Pmax_f(model, g, t):
        return model.p_g_t[g, t] + model.r_g_t[g, t] <= model.gen_Pmax[g] * model.u_g_t[g, t]

    model.gen_Pmax_cons = Constraint(model.GEN, model.TIME, rule=gen_Pmax_f)

    # r_g_t ramping constraint 1:
    def reserve_rr1_f(model, g, t):
        return model.r_g_t[g, t] <= model.gen_r10[g] * model.u_g_t[g, t]

    model.reserve_rr1_cons = Constraint(model.GEN, model.TIME, rule=reserve_rr1_f)

    # r_g_t ramping constraint 2:
    def reserve_rr2_f(model, g, t):
        return model.r_g_t[g, t] >= 0

    model.reserve_rr2_cons = Constraint(model.GEN, model.TIME, rule=reserve_rr2_f)

    # total reserve constraint
    def reserve_tot_f(model, g, t):
        reserve_tot_left = sum(model.r_g_t[g_1, t] for g_1 in model.GEN)
        reserve_tot_right = model.p_g_t[g, t] + model.r_g_t[g, t]
        return reserve_tot_left >= reserve_tot_right

    model.reserve_tot_cons = Constraint(model.GEN, model.TIME, rule=reserve_tot_f)

    # Theta define constraint
    def theta_def_f(model, k, t):
        fbus_num = model.line_fbus[k]
        tbus_num = model.line_tbus[k]
        return model.theta_k_t[k, t] == model.theta_b_t[fbus_num, t] - model.theta_b_t[tbus_num, t]

    model.theta_def_cons = Constraint(model.LINE, model.TIME, rule=theta_def_f)

    # Power flow constraint
    def pf_theta_f(model, k, t):
        return model.p_k_t[k, t] == model.theta_k_t[k, t] / model.line_x[k]

    model.pf_theta_cons = Constraint(model.LINE, model.TIME, rule=pf_theta_f)

    # Power flow min constraint:
    def pf_min_f(model, k, t):
        return model.p_k_t[k, t] >= -1 * model.line_Pmax[k]  # relax line rating

    model.pf_min_cons = Constraint(model.LINE, model.TIME, rule=pf_min_f)

    # Power flow max constraint:
    def pf_max_f(model, k, t):
        return model.p_k_t[k, t] <= 1 * model.line_Pmax[k]  # relax line rating

    model.pf_max_cons = Constraint(model.LINE, model.TIME, rule=pf_max_f)

    # Nodal balance constraint  (*** changed for hydrogen integration)
    def nodal_balance_f(model, b, t):
        nodal_balance_left = sum(model.p_g_t[g, t] for g in model.GEN if model.gen_bus[g] == b)
        nodal_balance_left += sum(model.p_k_t[k, t] for k in model.LINE if model.line_tbus[k] == b)
        nodal_balance_left -= sum(model.p_k_t[k, t] for k in model.LINE if model.line_fbus[k] == b)
        nodal_balance_right = model.load_b_t[b, t]
        if t!=0:
            nodal_balance_right -= sum(wind_hrly[wgen_num-1,t-1] for wgen_num in WGen if wgen_bus[wgen_num-1] == b)
            nodal_balance_right += sum(model.pcur_wind[wgen_num,t] for wgen_num in WGen if wgen_bus[wgen_num-1] == b)
        # integrate the electrolyzers
        for elct in range(len(elctro_bus)):
            if elctro_bus[elct]==b:
                nodal_balance_right += model.p_elct_t[elct+1,t]
        # integrate the fuel cells
        for fc in range(len(fc_bus)):
            if fc_bus[fc]==b:
                nodal_balance_left += model.p_fc_t[fc+1,t]
        return nodal_balance_left == nodal_balance_right

    model.nodal_balance_cons = Constraint(model.BUS, model.TIME, rule=nodal_balance_f)

    # Generator ramping rate constraint 1
    # Assume normal/startup/shutdown ramping rates are the same
    def gen_rr1_f(model, g, t):
        t_previous = model.time_num[t] - 1
        return model.p_g_t[g, t] - model.p_g_t[g, t_previous] <= model.gen_rr[g]

    model.gen_rr1_cons = Constraint(model.GEN, Time_24, rule=gen_rr1_f)

    # Generator ramping rate constraint 2
    def gen_rr2_f(model, g, t):
        t_previous = model.time_num[t] - 1
        return model.p_g_t[g, t] - model.p_g_t[g, t_previous] >= -model.gen_rr[g]

    model.gen_rr2_cons = Constraint(model.GEN, Time_24, rule=gen_rr2_f)

    # Variable V constraint
    def var_v_f(model, g, t):
        if t == 0:
            return model.v_g_t[g, t] >= 0  # no v_g_t constraint for t=0
        else:
            return model.v_g_t[g, t] >= model.u_g_t[g, t] - model.u_g_t[g, t - 1]

    model.var_v_cons = Constraint(model.GEN, model.TIME, rule=var_v_f)

    ### Power curtailment constraints
    def pcur_w_f(model,wg, t):
        return model.pcur_wind[wg,t] <= wind_hrly[wg-1,t-1]
    model.pcur_w_cons = Constraint(WGen,Time_24,rule=pcur_w_f)

    # ### new constraints for hydro integration
    # # Hydrogen-Power constraint 1
    # def hdg_p_f(model,b,t):
    #     if b in elctro_bus: # if the bus has FC & Electrolyzer
    #         alpha = 0.8 # 50% efficiency (round trip)
    #         beta = 0.6
    #         hdg_p_left = sum(alpha*model.p_elct_t[elct,t] for elct in Elctro if elctro_bus[elct-1]==b)
    #         hdg_p_right = sum(model.p_fc_t[fc,t]/beta for fc in FC if fc_bus[fc-1]==b)
    #         return hdg_p_left == hdg_p_right
    #     else: return model.p_elct_t[1,1]>=0
    # model.hdg_p_cons1 = Constraint(model.BUS,model.TIME, rule=hdg_p_f)

    ### new constraints for hydro integration, assume only one elctrolyzer and one fuel cell on a bus
    # hydrogen tank storage constraint
    def HNS_E_f(model,elct, t):
        alpha = 0.8  # 50% efficiency (round trip)
        beta = 0.6
        if t >= 1:
            bus_num = elctro_bus[elct-1]
            for fc in FC:
                if fc_bus[fc-1] == bus_num:
                    fc_onElct = fc
            return model.E_elct_t[elct,t] == model.E_elct_t[elct,t - 1] + model.p_elct_t[elct,t]*alpha - model.p_fc_t[fc_onElct, t] / beta
        else:
            return model.E_elct_t[elct,t] == 0  # initial hour, the storage is 0
    model.HNS_E_cons1 = Constraint(Elctro,model.TIME, rule=HNS_E_f)

    # elctrolyzer capacity constraint
    def elct_Pmax_f(model,elct,t):
        return  model.p_elct_t[elct,t]<=elctro_Pmax[elct-1]
    model.elct_Pmax_cons = Constraint(Elctro,model.TIME, rule=elct_Pmax_f)

    # fuel cell capacity constraint
    def fc_Pmax_f(model, fc, t):
        return model.p_fc_t[fc, t] <= fc_Pmax[fc - 1]
    model.fc_Pmax_cons = Constraint(FC, model.TIME, rule=fc_Pmax_f)


    # load case data and create instance
    print('start creating the instance')
    case_pyomo = model.create_instance('formpyomo_UC_simplified.dat')
    # dual variable setting
    case_pyomo.dual = pyomo.environ.Suffix(direction=pyomo.environ.Suffix.IMPORT)
    print('finish creating the instance')
    # case_pyomo.pprint()
    return case_pyomo

def build_EHSCUC_pyomo_run2(case_inst, load_dl, elctro_bus, elctro_Pmax, fc_bus, fc_Pmax, UC_case):
    # wind profile
    wind_hrly = np.loadtxt('wind_hrly_adj.txt')  # read the renewable profiles created in case 2 [rgen][hr]
    ## list of non wind or solar generators
    reg_gen_list = []
    nonreg_gen_list = []
    nonreg_gen_bus = []
    wgen_bus = []
    for i in range(case_inst.gentotnum):
        if case_inst.gen.fuel_type[i] != 'Wind':
            if case_inst.gen.fuel_type[i] != 'Solar':
                gen_num = i + 1
                reg_gen_list.append(gen_num)
        if case_inst.gen.fuel_type[i] == 'Wind':
            gen_num = i + 1
            nonreg_gen_list.append(gen_num)
            nonreg_gen_bus.append(case_inst.gen.bus[i])
            wgen_bus.append(case_inst.gen.bus[i])
        if case_inst.gen.fuel_type[i] == 'Solar':
            gen_num = i + 1
            nonreg_gen_list.append(gen_num)
            nonreg_gen_bus.append(case_inst.gen.bus[i])

    ## generate data file
    time_totnum = 24  # total number of time period
    pyomodata_UC_simplified('formpyomo_UC_simplified.dat',case_inst, time_totnum, load_dl, reg_gen_list)

    ### Python set
    # 24 hour time set
    Time_24 = []
    for i in range(24):
        Time_24.append(i + 1)
    # elctrolyzer set
    Elctro = []
    for i in range(len(elctro_bus)):
        Elctro.append(i+1)
    # fuel cell set
    FC = []
    for i in range(len(fc_bus)):
        FC.append(i+1)
    # renew gen set
    WGen = []
    for i in range(wind_hrly.shape[0]):
        WGen.append(i+1)

    ### The pyomo abstract model for UC
    model = AbstractModel()

    # Set
    model.BUS = Set()
    model.LINE = Set()
    model.GEN = Set()
    model.TIME = Set()
    model.BUS_TIME = Set(dimen=2)

    ## Param
    # Bus param
    model.bus_num = Param(model.BUS)
    # Gen param
    model.gen_num = Param(model.GEN)
    model.gen_bus = Param(model.GEN)
    model.gen_cost_P = Param(model.GEN)
    model.gen_cost_NL = Param(model.GEN)
    model.gen_cost_SU = Param(model.GEN)
    model.gen_Pmin = Param(model.GEN)
    model.gen_Pmax = Param(model.GEN)
    model.gen_r10 = Param(model.GEN)
    model.gen_rr = Param(model.GEN)
    # Line param
    model.line_num = Param(model.LINE)
    model.line_x = Param(model.LINE)
    model.line_Pmax = Param(model.LINE)
    model.line_fbus = Param(model.LINE)
    model.line_tbus = Param(model.LINE)
    # Time param
    model.time_num = Param(model.TIME)
    # Bus_Time param
    model.bustime_num = Param(model.BUS_TIME)
    model.load_b_t = Param(model.BUS_TIME)

    # Variable
    # Gen_Time Var
    model.p_g_t = Var(model.GEN, model.TIME)
    # model.u_g_t = Var(model.GEN, model.TIME, domain=Binary)
    model.v_g_t = Var(model.GEN, model.TIME)  # v_g_t is simplified to non-binary variable
    model.r_g_t = Var(model.GEN, model.TIME)
    # Line_Time Var
    model.theta_k_t = Var(model.LINE, model.TIME)
    model.p_k_t = Var(model.LINE, model.TIME)
    # Bus_Time Var
    model.theta_b_t = Var(model.BUS, model.TIME)
    # Electrolyzer Var
    model.p_elct_t = Var(Elctro, model.TIME,domain=NonNegativeReals)
    # fuel cell Var
    model.p_fc_t = Var(FC, model.TIME,domain=NonNegativeReals)
    # hydrogen tank storage
    model.E_elct_t = Var(Elctro, model.TIME, domain=NonNegativeReals)
    # Renew Curtailment Var
    model.pcur_wind = Var(WGen, model.TIME, domain=NonNegativeReals)

    ## Objective function
    def objfunction(model):
        obj = sum(
            model.gen_cost_P[g] * model.p_g_t[g, t] + model.gen_cost_NL[g] * UC_case.u_g_t[g, t]() + model.gen_cost_SU[g] *
            model.v_g_t[g, t] for g in model.GEN for t in model.TIME)
        return obj

    model.object = Objective(rule=objfunction, sense=minimize)

    ## v_g_t constraint
    # v_g_t constraint 1
    def vgt_f_1(model, g, t):
        return model.v_g_t[g, t] >= 0

    model.vgt_cons_1 = Constraint(model.GEN, model.TIME, rule=vgt_f_1)

    # v_g_t constraint 2
    def vgt_f_2(model, g, t):
        return model.v_g_t[g, t] <= 1

    model.vgt_cons_2 = Constraint(model.GEN, model.TIME, rule=vgt_f_2)

    # ## u_g_t constraint
    # def ugt_f(model, g, t):
    #     return model.u_g_t[g, t] == 1
    #
    # model.ugt_cons = Constraint(model.GEN, model.TIME, rule=ugt_f)

    ## Generator power and reserve constraints
    # # P_g_t minimum constraint
    def gen_Pmin_f(model, g, t):
        return model.gen_Pmin[g] * UC_case.u_g_t[g, t]() <= model.p_g_t[g, t]

    model.gen_Pmin_cons = Constraint(model.GEN, model.TIME, rule=gen_Pmin_f)

    # P_g_t maximum constraint:
    def gen_Pmax_f(model, g, t):
        return model.p_g_t[g, t] + model.r_g_t[g, t] <= model.gen_Pmax[g] * UC_case.u_g_t[g, t]()

    model.gen_Pmax_cons = Constraint(model.GEN, model.TIME, rule=gen_Pmax_f)

    # r_g_t ramping constraint 1:
    def reserve_rr1_f(model, g, t):
        return model.r_g_t[g, t] <= model.gen_r10[g] * UC_case.u_g_t[g, t]()

    model.reserve_rr1_cons = Constraint(model.GEN, model.TIME, rule=reserve_rr1_f)

    # r_g_t ramping constraint 2:
    def reserve_rr2_f(model, g, t):
        return model.r_g_t[g, t] >= 0

    model.reserve_rr2_cons = Constraint(model.GEN, model.TIME, rule=reserve_rr2_f)

    # total reserve constraint
    def reserve_tot_f(model, g, t):
        reserve_tot_left = sum(model.r_g_t[g_1, t] for g_1 in model.GEN)
        reserve_tot_right = model.p_g_t[g, t] + model.r_g_t[g, t]
        return reserve_tot_left >= reserve_tot_right

    model.reserve_tot_cons = Constraint(model.GEN, model.TIME, rule=reserve_tot_f)

    # Theta define constraint
    def theta_def_f(model, k, t):
        fbus_num = model.line_fbus[k]
        tbus_num = model.line_tbus[k]
        return model.theta_k_t[k, t] == model.theta_b_t[fbus_num, t] - model.theta_b_t[tbus_num, t]

    model.theta_def_cons = Constraint(model.LINE, model.TIME, rule=theta_def_f)

    # Power flow constraint
    def pf_theta_f(model, k, t):
        return model.p_k_t[k, t] == model.theta_k_t[k, t] / model.line_x[k]

    model.pf_theta_cons = Constraint(model.LINE, model.TIME, rule=pf_theta_f)

    # Power flow min constraint:
    def pf_min_f(model, k, t):
        return model.p_k_t[k, t] >= -1 * model.line_Pmax[k]  # relax line rating

    model.pf_min_cons = Constraint(model.LINE, model.TIME, rule=pf_min_f)

    # Power flow max constraint:
    def pf_max_f(model, k, t):
        return model.p_k_t[k, t] <= 1 * model.line_Pmax[k]  # relax line rating

    model.pf_max_cons = Constraint(model.LINE, model.TIME, rule=pf_max_f)

    # Nodal balance constraint  (*** changed for hydrogen integration)
    def nodal_balance_f(model, b, t):
        nodal_balance_left = sum(model.p_g_t[g, t] for g in model.GEN if model.gen_bus[g] == b)
        nodal_balance_left += sum(model.p_k_t[k, t] for k in model.LINE if model.line_tbus[k] == b)
        nodal_balance_left -= sum(model.p_k_t[k, t] for k in model.LINE if model.line_fbus[k] == b)
        nodal_balance_right = model.load_b_t[b, t]
        if t!=0:
            nodal_balance_right -= sum(wind_hrly[wgen_num - 1, t - 1] for wgen_num in WGen if wgen_bus[wgen_num - 1] == b)
            nodal_balance_right += sum(model.pcur_wind[wgen_num, t] for wgen_num in WGen if wgen_bus[wgen_num - 1] == b)
        # integrate the electrolyzers
        for elct in range(len(elctro_bus)):
            if elctro_bus[elct]==b:
                nodal_balance_right += model.p_elct_t[elct+1,t]
        # integrate the fuel cells
        for fc in range(len(fc_bus)):
            if fc_bus[fc]==b:
                nodal_balance_left += model.p_fc_t[fc+1,t]
        return nodal_balance_left == nodal_balance_right

    model.nodal_balance_cons = Constraint(model.BUS, model.TIME, rule=nodal_balance_f)

    # Generator ramping rate constraint 1
    # Assume normal/startup/shutdown ramping rates are the same
    def gen_rr1_f(model, g, t):
        t_previous = model.time_num[t] - 1
        return model.p_g_t[g, t] - model.p_g_t[g, t_previous] <= model.gen_rr[g]

    model.gen_rr1_cons = Constraint(model.GEN, Time_24, rule=gen_rr1_f)

    # Generator ramping rate constraint 2
    def gen_rr2_f(model, g, t):
        t_previous = model.time_num[t] - 1
        return model.p_g_t[g, t] - model.p_g_t[g, t_previous] >= -model.gen_rr[g]

    model.gen_rr2_cons = Constraint(model.GEN, Time_24, rule=gen_rr2_f)

    # Variable V constraint
    def var_v_f(model, g, t):
        if t == 0:
            return model.v_g_t[g, t] >= 0  # no v_g_t constraint for t=0
        else:
            return model.v_g_t[g, t] >= UC_case.u_g_t[g, t]() - UC_case.u_g_t[g, t - 1]()

    model.var_v_cons = Constraint(model.GEN, model.TIME, rule=var_v_f)

    ### Power curtailment constraints
    def pcur_w_f(model, wg, t):
        return model.pcur_wind[wg, t] <= wind_hrly[wg - 1, t - 1]

    model.pcur_w_cons = Constraint(WGen, Time_24, rule=pcur_w_f)

    ### new constraints for hydro integration, assume only one elctrolyzer and one fuel cell on a bus
    # hydrogen tank storage constraint
    def HNS_E_f(model, elct, t):
        alpha = 0.8  # 50% efficiency (round trip)
        beta = 0.6
        if t >= 1:
            bus_num = elctro_bus[elct - 1]
            for fc in FC:
                if fc_bus[fc - 1] == bus_num:
                    fc_onElct = fc
            return model.E_elct_t[elct, t] == model.E_elct_t[elct, t - 1] + model.p_elct_t[elct, t] * alpha - \
                   model.p_fc_t[fc_onElct, t] / beta
        else:
            return model.E_elct_t[elct, t] == 0  # initial hour, the storage is 0

    model.HNS_E_cons1 = Constraint(Elctro, model.TIME, rule=HNS_E_f)

    # elctrolyzer capacity constraint
    def elct_Pmax_f(model,elct,t):
        return  model.p_elct_t[elct,t]<=elctro_Pmax[elct-1]
    model.elct_Pmax_cons = Constraint(Elctro,model.TIME, rule=elct_Pmax_f)

    # fuel cell capacity constraint
    def fc_Pmax_f(model, fc, t):
        return model.p_fc_t[fc, t] <= fc_Pmax[fc - 1]

    model.fc_Pmax_cons = Constraint(FC, model.TIME, rule=fc_Pmax_f)


    # load case data and create instance
    print('start creating the instance')
    case_pyomo = model.create_instance('formpyomo_UC_simplified.dat')
    # dual variable setting
    case_pyomo.dual = pyomo.environ.Suffix(direction=pyomo.environ.Suffix.IMPORT)
    print('finish creating the instance')
    # case_pyomo.pprint()
    return case_pyomo

#### FC & Electrolyzer different location (Hydrogen Network, Network has inifinite storage capacity)
### function to build the pyomo instance based on the test case
def build_HSCUC_pyomo(case_inst, load_dl, elctro_bus, elctro_Pmax, fc_bus, fc_Pmax):
    # wind profile
    wind_hrly = np.loadtxt('wind_hrly_adj.txt')  # read the renewable profiles created in case 2 [rgen][hr]
    ## list of non wind or solar generators
    reg_gen_list = []
    nonreg_gen_list = []
    nonreg_gen_bus = []
    wgen_bus = []
    for i in range(case_inst.gentotnum):
        if case_inst.gen.fuel_type[i] != 'Wind':
            if case_inst.gen.fuel_type[i] != 'Solar':
                gen_num = i + 1
                reg_gen_list.append(gen_num)
        if case_inst.gen.fuel_type[i] == 'Wind':
            gen_num = i + 1
            nonreg_gen_list.append(gen_num)
            nonreg_gen_bus.append(case_inst.gen.bus[i])
            wgen_bus.append(case_inst.gen.bus[i])
        if case_inst.gen.fuel_type[i] == 'Solar':
            gen_num = i + 1
            nonreg_gen_list.append(gen_num)
            nonreg_gen_bus.append(case_inst.gen.bus[i])

    ## generate data file
    time_totnum = 24  # total number of time period
    pyomodata_UC_simplified('formpyomo_UC_simplified.dat',case_inst, time_totnum, load_dl, reg_gen_list)

    ### Python set
    # 24 hour time set
    Time_24 = []
    for i in range(24):
        Time_24.append(i + 1)
    # elctrolyzer set
    Elctro = []
    for i in range(len(elctro_bus)):
        Elctro.append(i + 1)
    # fuel cell set
    FC = []
    for i in range(len(fc_bus)):
        FC.append(i + 1)
    # renew gen set
    WGen = []
    for i in range(wind_hrly.shape[0]):
        WGen.append(i + 1)

    ### The pyomo abstract model for UC
    model = AbstractModel()

    # Set
    model.BUS = Set()
    model.LINE = Set()
    model.GEN = Set()
    model.TIME = Set()
    model.BUS_TIME = Set(dimen=2)

    ## Param
    # Bus param
    model.bus_num = Param(model.BUS)
    # Gen param
    model.gen_num = Param(model.GEN)
    model.gen_bus = Param(model.GEN)
    model.gen_cost_P = Param(model.GEN)
    model.gen_cost_NL = Param(model.GEN)
    model.gen_cost_SU = Param(model.GEN)
    model.gen_Pmin = Param(model.GEN)
    model.gen_Pmax = Param(model.GEN)
    model.gen_r10 = Param(model.GEN)
    model.gen_rr = Param(model.GEN)
    # Line param
    model.line_num = Param(model.LINE)
    model.line_x = Param(model.LINE)
    model.line_Pmax = Param(model.LINE)
    model.line_fbus = Param(model.LINE)
    model.line_tbus = Param(model.LINE)
    # Time param
    model.time_num = Param(model.TIME)
    # Bus_Time param
    model.bustime_num = Param(model.BUS_TIME)
    model.load_b_t = Param(model.BUS_TIME)

    # Variable
    # Gen_Time Var
    model.p_g_t = Var(model.GEN, model.TIME)
    model.u_g_t = Var(model.GEN, model.TIME, domain=Binary)
    model.v_g_t = Var(model.GEN, model.TIME, domain=Binary)  # v_g_t is simplified to non-binary variable
    model.r_g_t = Var(model.GEN, model.TIME)
    # Line_Time Var
    model.theta_k_t = Var(model.LINE, model.TIME)
    model.p_k_t = Var(model.LINE, model.TIME)
    # Bus_Time Var
    model.theta_b_t = Var(model.BUS, model.TIME)
    # Electrolyzer Var
    model.p_elct_t = Var(Elctro, model.TIME, domain=NonNegativeReals)
    # fuel cell Var
    model.p_fc_t = Var(FC, model.TIME, domain=NonNegativeReals)
    # hydrogen network storage
    model.p_HNS_t = Var(model.TIME, domain=NonNegativeReals)
    # Renew Curtailment Var
    model.pcur_wind = Var(WGen, model.TIME, domain=NonNegativeReals)

    ## Objective function
    def objfunction(model):
        obj = sum(
            model.gen_cost_P[g] * model.p_g_t[g, t] + model.gen_cost_NL[g] * model.u_g_t[g, t] + model.gen_cost_SU[g] *
            model.v_g_t[g, t] for g in model.GEN for t in model.TIME)
        return obj

    model.object = Objective(rule=objfunction, sense=minimize)

    # Generator initial u_g_t status
    def gen_Uinit_f(model, g):
        return model.u_g_t[g, 0] == 0

    model.gen_Uinit_cons = Constraint(model.GEN, rule=gen_Uinit_f)

    ## Generator power and reserve constraints
    # # P_g_t minimum constraint
    def gen_Pmin_f(model, g, t):
        return model.gen_Pmin[g] * model.u_g_t[g, t] <= model.p_g_t[g, t]

    model.gen_Pmin_cons = Constraint(model.GEN, model.TIME, rule=gen_Pmin_f)

    # P_g_t maximum constraint:
    def gen_Pmax_f(model, g, t):
        return model.p_g_t[g, t] + model.r_g_t[g, t] <= model.gen_Pmax[g] * model.u_g_t[g, t]

    model.gen_Pmax_cons = Constraint(model.GEN, model.TIME, rule=gen_Pmax_f)

    # r_g_t ramping constraint 1:
    def reserve_rr1_f(model, g, t):
        return model.r_g_t[g, t] <= model.gen_r10[g] * model.u_g_t[g, t]

    model.reserve_rr1_cons = Constraint(model.GEN, model.TIME, rule=reserve_rr1_f)

    # r_g_t ramping constraint 2:
    def reserve_rr2_f(model, g, t):
        return model.r_g_t[g, t] >= 0

    model.reserve_rr2_cons = Constraint(model.GEN, model.TIME, rule=reserve_rr2_f)

    # total reserve constraint
    def reserve_tot_f(model, g, t):
        reserve_tot_left = sum(model.r_g_t[g_1, t] for g_1 in model.GEN)
        reserve_tot_right = model.p_g_t[g, t] + model.r_g_t[g, t]
        return reserve_tot_left >= reserve_tot_right

    model.reserve_tot_cons = Constraint(model.GEN, model.TIME, rule=reserve_tot_f)

    # Theta define constraint
    def theta_def_f(model, k, t):
        fbus_num = model.line_fbus[k]
        tbus_num = model.line_tbus[k]
        return model.theta_k_t[k, t] == model.theta_b_t[fbus_num, t] - model.theta_b_t[tbus_num, t]

    model.theta_def_cons = Constraint(model.LINE, model.TIME, rule=theta_def_f)

    # Power flow constraint
    def pf_theta_f(model, k, t):
        return model.p_k_t[k, t] == model.theta_k_t[k, t] / model.line_x[k]

    model.pf_theta_cons = Constraint(model.LINE, model.TIME, rule=pf_theta_f)

    # Power flow min constraint:
    def pf_min_f(model, k, t):
        return model.p_k_t[k, t] >= -1 * model.line_Pmax[k]  # relax line rating

    model.pf_min_cons = Constraint(model.LINE, model.TIME, rule=pf_min_f)

    # Power flow max constraint:
    def pf_max_f(model, k, t):
        return model.p_k_t[k, t] <= 1 * model.line_Pmax[k]  # relax line rating

    model.pf_max_cons = Constraint(model.LINE, model.TIME, rule=pf_max_f)

    # Nodal balance constraint  (*** changed for hydrogen integration)
    def nodal_balance_f(model, b, t):
        nodal_balance_left = sum(model.p_g_t[g, t] for g in model.GEN if model.gen_bus[g] == b)
        nodal_balance_left += sum(model.p_k_t[k, t] for k in model.LINE if model.line_tbus[k] == b)
        nodal_balance_left -= sum(model.p_k_t[k, t] for k in model.LINE if model.line_fbus[k] == b)
        nodal_balance_right = model.load_b_t[b, t]
        if t != 0:
            nodal_balance_right -= sum(
                wind_hrly[wgen_num - 1, t - 1] for wgen_num in WGen if wgen_bus[wgen_num - 1] == b)
            nodal_balance_right += sum(model.pcur_wind[wgen_num, t] for wgen_num in WGen if wgen_bus[wgen_num - 1] == b)
        # integrate the electrolyzers
        for elct in range(len(elctro_bus)):
            if elctro_bus[elct] == b:
                nodal_balance_right += model.p_elct_t[elct + 1, t]
        # integrate the fuel cells
        for fc in range(len(fc_bus)):
            if fc_bus[fc] == b:
                nodal_balance_left += model.p_fc_t[fc + 1, t]
        return nodal_balance_left == nodal_balance_right

    model.nodal_balance_cons = Constraint(model.BUS, model.TIME, rule=nodal_balance_f)

    # Generator ramping rate constraint 1
    # Assume normal/startup/shutdown ramping rates are the same
    def gen_rr1_f(model, g, t):
        t_previous = model.time_num[t] - 1
        return model.p_g_t[g, t] - model.p_g_t[g, t_previous] <= model.gen_rr[g]

    model.gen_rr1_cons = Constraint(model.GEN, Time_24, rule=gen_rr1_f)

    # Generator ramping rate constraint 2
    def gen_rr2_f(model, g, t):
        t_previous = model.time_num[t] - 1
        return model.p_g_t[g, t] - model.p_g_t[g, t_previous] >= -model.gen_rr[g]

    model.gen_rr2_cons = Constraint(model.GEN, Time_24, rule=gen_rr2_f)

    # Variable V constraint
    def var_v_f(model, g, t):
        if t == 0:
            return model.v_g_t[g, t] >= 0  # no v_g_t constraint for t=0
        else:
            return model.v_g_t[g, t] >= model.u_g_t[g, t] - model.u_g_t[g, t - 1]

    model.var_v_cons = Constraint(model.GEN, model.TIME, rule=var_v_f)

    ### Power curtailment constraints
    def pcur_w_f(model, wg, t):
        return model.pcur_wind[wg, t] <= wind_hrly[wg - 1, t - 1]

    model.pcur_w_cons = Constraint(WGen, Time_24, rule=pcur_w_f)

    ### new constraints for hydro integration
    # hydrogen network capacity constraint
    def HNS_Pmax_f(model, t):
        alpha = 0.8  # 50% efficiency (round trip)
        beta = 0.6
        if t >= 1:
            return model.p_HNS_t[t] == model.p_HNS_t[t - 1] + sum(
                alpha * model.p_elct_t[elct, t] for elct in Elctro) - sum(model.p_fc_t[fc, t] / beta for fc in FC)
        else:
            return model.p_HNS_t[t] == 0  # initial hour, the storage is 0

    model.HNS_Pmax_cons1 = Constraint(model.TIME, rule=HNS_Pmax_f)

    # elctrolyzer capacity constraint
    def elct_Pmax_f(model, elct, t):
        return model.p_elct_t[elct, t] <= elctro_Pmax[elct - 1]

    model.elct_Pmax_cons = Constraint(Elctro, model.TIME, rule=elct_Pmax_f)

    # fuel cell capacity constraint
    def fc_Pmax_f(model, fc, t):
        return model.p_fc_t[fc, t] <= fc_Pmax[fc - 1]

    model.fc_Pmax_cons = Constraint(FC, model.TIME, rule=fc_Pmax_f)


    # load case data and create instance
    print('start creating the instance')
    case_pyomo = model.create_instance('formpyomo_UC_simplified.dat')
    # dual variable setting
    case_pyomo.dual = pyomo.environ.Suffix(direction=pyomo.environ.Suffix.IMPORT)
    print('finish creating the instance')
    # case_pyomo.pprint()
    return case_pyomo

def build_HSCUC_pyomo_run2(case_inst, load_dl, elctro_bus, elctro_Pmax, fc_bus, fc_Pmax, UC_case):
    # wind profile
    wind_hrly = np.loadtxt('wind_hrly_adj.txt')  # read the renewable profiles created in case 2 [rgen][hr]
    ## list of non wind or solar generators
    reg_gen_list = []
    nonreg_gen_list = []
    nonreg_gen_bus = []
    wgen_bus = []
    for i in range(case_inst.gentotnum):
        if case_inst.gen.fuel_type[i] != 'Wind':
            if case_inst.gen.fuel_type[i] != 'Solar':
                gen_num = i + 1
                reg_gen_list.append(gen_num)
        if case_inst.gen.fuel_type[i] == 'Wind':
            gen_num = i + 1
            nonreg_gen_list.append(gen_num)
            nonreg_gen_bus.append(case_inst.gen.bus[i])
            wgen_bus.append(case_inst.gen.bus[i])
        if case_inst.gen.fuel_type[i] == 'Solar':
            gen_num = i + 1
            nonreg_gen_list.append(gen_num)
            nonreg_gen_bus.append(case_inst.gen.bus[i])

    ## generate data file
    time_totnum = 24  # total number of time period
    pyomodata_UC_simplified('formpyomo_UC_simplified.dat',case_inst, time_totnum, load_dl, reg_gen_list)

    ### Python set
    # 24 hour time set
    Time_24 = []
    for i in range(24):
        Time_24.append(i + 1)
    # elctrolyzer set
    Elctro = []
    for i in range(len(elctro_bus)):
        Elctro.append(i + 1)
    # fuel cell set
    FC = []
    for i in range(len(fc_bus)):
        FC.append(i + 1)
    # renew gen set
    WGen = []
    for i in range(wind_hrly.shape[0]):
        WGen.append(i + 1)

    ### The pyomo abstract model for UC
    model = AbstractModel()

    # Set
    model.BUS = Set()
    model.LINE = Set()
    model.GEN = Set()
    model.TIME = Set()
    model.BUS_TIME = Set(dimen=2)

    ## Param
    # Bus param
    model.bus_num = Param(model.BUS)
    # Gen param
    model.gen_num = Param(model.GEN)
    model.gen_bus = Param(model.GEN)
    model.gen_cost_P = Param(model.GEN)
    model.gen_cost_NL = Param(model.GEN)
    model.gen_cost_SU = Param(model.GEN)
    model.gen_Pmin = Param(model.GEN)
    model.gen_Pmax = Param(model.GEN)
    model.gen_r10 = Param(model.GEN)
    model.gen_rr = Param(model.GEN)
    # Line param
    model.line_num = Param(model.LINE)
    model.line_x = Param(model.LINE)
    model.line_Pmax = Param(model.LINE)
    model.line_fbus = Param(model.LINE)
    model.line_tbus = Param(model.LINE)
    # Time param
    model.time_num = Param(model.TIME)
    # Bus_Time param
    model.bustime_num = Param(model.BUS_TIME)
    model.load_b_t = Param(model.BUS_TIME)

    # Variable
    # Gen_Time Var
    model.p_g_t = Var(model.GEN, model.TIME)
    # model.u_g_t = Var(model.GEN, model.TIME, domain=Binary)
    model.v_g_t = Var(model.GEN, model.TIME)  # v_g_t is simplified to non-binary variable
    model.r_g_t = Var(model.GEN, model.TIME)
    # Line_Time Var
    model.theta_k_t = Var(model.LINE, model.TIME)
    model.p_k_t = Var(model.LINE, model.TIME)
    # Bus_Time Var
    model.theta_b_t = Var(model.BUS, model.TIME)
    # Electrolyzer Var
    model.p_elct_t = Var(Elctro, model.TIME, domain=NonNegativeReals)
    # fuel cell Var
    model.p_fc_t = Var(FC, model.TIME, domain=NonNegativeReals)
    # hydrogen network storage
    model.p_HNS_t = Var(model.TIME, domain=NonNegativeReals)
    # Renew Curtailment Var
    model.pcur_wind = Var(WGen, model.TIME, domain=NonNegativeReals)

    ## Objective function
    def objfunction(model):
        obj = sum(
            model.gen_cost_P[g] * model.p_g_t[g, t] + model.gen_cost_NL[g] * UC_case.u_g_t[g, t]() + model.gen_cost_SU[
                g] *
            model.v_g_t[g, t] for g in model.GEN for t in model.TIME)
        return obj

    model.object = Objective(rule=objfunction, sense=minimize)

    ## v_g_t constraint
    # v_g_t constraint 1
    def vgt_f_1(model, g, t):
        return model.v_g_t[g, t] >= 0

    model.vgt_cons_1 = Constraint(model.GEN, model.TIME, rule=vgt_f_1)

    # v_g_t constraint 2
    def vgt_f_2(model, g, t):
        return model.v_g_t[g, t] <= 1

    model.vgt_cons_2 = Constraint(model.GEN, model.TIME, rule=vgt_f_2)

    # ## u_g_t constraint
    # def ugt_f(model, g, t):
    #     return model.u_g_t[g, t] == 1
    #
    # model.ugt_cons = Constraint(model.GEN, model.TIME, rule=ugt_f)

    ## Generator power and reserve constraints
    # # P_g_t minimum constraint
    def gen_Pmin_f(model, g, t):
        return model.gen_Pmin[g] * UC_case.u_g_t[g, t]() <= model.p_g_t[g, t]

    model.gen_Pmin_cons = Constraint(model.GEN, model.TIME, rule=gen_Pmin_f)

    # P_g_t maximum constraint:
    def gen_Pmax_f(model, g, t):
        return model.p_g_t[g, t] + model.r_g_t[g, t] <= model.gen_Pmax[g] * UC_case.u_g_t[g, t]()

    model.gen_Pmax_cons = Constraint(model.GEN, model.TIME, rule=gen_Pmax_f)

    # r_g_t ramping constraint 1:
    def reserve_rr1_f(model, g, t):
        return model.r_g_t[g, t] <= model.gen_r10[g] * UC_case.u_g_t[g, t]()

    model.reserve_rr1_cons = Constraint(model.GEN, model.TIME, rule=reserve_rr1_f)

    # r_g_t ramping constraint 2:
    def reserve_rr2_f(model, g, t):
        return model.r_g_t[g, t] >= 0

    model.reserve_rr2_cons = Constraint(model.GEN, model.TIME, rule=reserve_rr2_f)

    # total reserve constraint
    def reserve_tot_f(model, g, t):
        reserve_tot_left = sum(model.r_g_t[g_1, t] for g_1 in model.GEN)
        reserve_tot_right = model.p_g_t[g, t] + model.r_g_t[g, t]
        return reserve_tot_left >= reserve_tot_right

    model.reserve_tot_cons = Constraint(model.GEN, model.TIME, rule=reserve_tot_f)

    # Theta define constraint
    def theta_def_f(model, k, t):
        fbus_num = model.line_fbus[k]
        tbus_num = model.line_tbus[k]
        return model.theta_k_t[k, t] == model.theta_b_t[fbus_num, t] - model.theta_b_t[tbus_num, t]

    model.theta_def_cons = Constraint(model.LINE, model.TIME, rule=theta_def_f)

    # Power flow constraint
    def pf_theta_f(model, k, t):
        return model.p_k_t[k, t] == model.theta_k_t[k, t] / model.line_x[k]

    model.pf_theta_cons = Constraint(model.LINE, model.TIME, rule=pf_theta_f)

    # Power flow min constraint:
    def pf_min_f(model, k, t):
        return model.p_k_t[k, t] >= -1 * model.line_Pmax[k]  # relax line rating

    model.pf_min_cons = Constraint(model.LINE, model.TIME, rule=pf_min_f)

    # Power flow max constraint:
    def pf_max_f(model, k, t):
        return model.p_k_t[k, t] <= 1 * model.line_Pmax[k]  # relax line rating

    model.pf_max_cons = Constraint(model.LINE, model.TIME, rule=pf_max_f)

    # Nodal balance constraint  (*** changed for hydrogen integration)
    def nodal_balance_f(model, b, t):
        nodal_balance_left = sum(model.p_g_t[g, t] for g in model.GEN if model.gen_bus[g] == b)
        nodal_balance_left += sum(model.p_k_t[k, t] for k in model.LINE if model.line_tbus[k] == b)
        nodal_balance_left -= sum(model.p_k_t[k, t] for k in model.LINE if model.line_fbus[k] == b)
        nodal_balance_right = model.load_b_t[b, t]
        if t != 0:
            nodal_balance_right -= sum(
                wind_hrly[wgen_num - 1, t - 1] for wgen_num in WGen if wgen_bus[wgen_num - 1] == b)
            nodal_balance_right += sum(model.pcur_wind[wgen_num, t] for wgen_num in WGen if wgen_bus[wgen_num - 1] == b)
        # integrate the electrolyzers
        for elct in range(len(elctro_bus)):
            if elctro_bus[elct] == b:
                nodal_balance_right += model.p_elct_t[elct + 1, t]
        # integrate the fuel cells
        for fc in range(len(fc_bus)):
            if fc_bus[fc] == b:
                nodal_balance_left += model.p_fc_t[fc + 1, t]
        return nodal_balance_left == nodal_balance_right

    model.nodal_balance_cons = Constraint(model.BUS, model.TIME, rule=nodal_balance_f)

    # Generator ramping rate constraint 1
    # Assume normal/startup/shutdown ramping rates are the same
    def gen_rr1_f(model, g, t):
        t_previous = model.time_num[t] - 1
        return model.p_g_t[g, t] - model.p_g_t[g, t_previous] <= model.gen_rr[g]

    model.gen_rr1_cons = Constraint(model.GEN, Time_24, rule=gen_rr1_f)

    # Generator ramping rate constraint 2
    def gen_rr2_f(model, g, t):
        t_previous = model.time_num[t] - 1
        return model.p_g_t[g, t] - model.p_g_t[g, t_previous] >= -model.gen_rr[g]

    model.gen_rr2_cons = Constraint(model.GEN, Time_24, rule=gen_rr2_f)

    # Variable V constraint
    def var_v_f(model, g, t):
        if t == 0:
            return model.v_g_t[g, t] >= 0  # no v_g_t constraint for t=0
        else:
            return model.v_g_t[g, t] >= UC_case.u_g_t[g, t]() - UC_case.u_g_t[g, t - 1]()

    model.var_v_cons = Constraint(model.GEN, model.TIME, rule=var_v_f)

    ### Power curtailment constraints
    def pcur_w_f(model, wg, t):
        return model.pcur_wind[wg, t] <= wind_hrly[wg - 1, t - 1]

    model.pcur_w_cons = Constraint(WGen, Time_24, rule=pcur_w_f)

    ### new constraints for hydro integration
    # Hydrogen network storage
    def HNS_Pmax_f(model, t):
        alpha = 0.8  # 50% efficiency (round trip)
        beta = 0.6
        if t >= 1:
            return model.p_HNS_t[t] == model.p_HNS_t[t - 1] + sum(
                alpha * model.p_elct_t[elct, t] for elct in Elctro) - sum(model.p_fc_t[fc, t] / beta for fc in FC)
        else:
            return model.p_HNS_t[t] == 0  # initial hour, the storage is 0

    model.HNS_Pmax_cons1 = Constraint(model.TIME, rule=HNS_Pmax_f)

    # elctrolyzer capacity constraint
    def elct_Pmax_f(model, elct, t):
        return model.p_elct_t[elct, t] <= elctro_Pmax[elct - 1]

    model.elct_Pmax_cons = Constraint(Elctro, model.TIME, rule=elct_Pmax_f)

    # fuel cell capacity constraint
    def fc_Pmax_f(model, fc, t):
        return model.p_fc_t[fc, t] <= fc_Pmax[fc - 1]

    model.fc_Pmax_cons = Constraint(FC, model.TIME, rule=fc_Pmax_f)



    # load case data and create instance
    print('start creating the instance')
    case_pyomo = model.create_instance('formpyomo_UC_simplified.dat')
    # dual variable setting
    case_pyomo.dual = pyomo.environ.Suffix(direction=pyomo.environ.Suffix.IMPORT)
    print('finish creating the instance')
    # case_pyomo.pprint()
    return case_pyomo
