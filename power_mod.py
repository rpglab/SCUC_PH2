# Create a Class for matpower/pypower case
'''
GenericModel(power module), Bus, Gen, Branch, Gencost
input is index(start from 0), output is num(start from 1)
GenericModel_add add the power flow for branch
'''
# Author: Jin Lu (University of Houston)
## link: https://rpglab.github.io/resources/

# Import
from pypower import loadcase


# classic power system case class
class GenericModel:
    def __init__(self, py_case):
        self.case_dict = loadcase.loadcase(py_case)
        self.baseMVA = self.case_dict['baseMVA']
        self.bus_array = self.case_dict['bus']
        self.gen_array = self.case_dict['gen']
        self.branch_array = self.case_dict['branch']
        if 'gencost' in self.case_dict:
            self.gencost_array = self.case_dict['gencost']  # the gen cost array of the matpower case
            self.gencost = self.loadgencost(self.gencost_array) # the gen cost data: type, startup, shutdown
        self.bus = self.loadbus(self.bus_array)
        self.gen = self.loadgen(self.gen_array)
        self.branch = self.loadbranch(self.branch_array)
        self.bustotnum = self.bus_array.shape[0]             #total bus number
        self.gentotnum = self.gen_array.shape[0]  # total bus number
        self.branchtotnum = self.branch_array.shape[0]  # total bus number

    # loadbus function can load the bus array data to the class Bus
    def loadbus(self,bus_data):
        bus_num = bus_data[:, 0].astype(int)
        bus_type = bus_data[:, 1].astype(int)
        bus_Pd = bus_data[:, 2]
        bus_Qd = bus_data[:, 3]
        bus_Gs = bus_data[:, 4]
        bus_Bs = bus_data[:, 5]
        bus_area = bus_data[:, 6]
        bus_Vm = bus_data[:, 7]
        bus_Va = bus_data[:, 8]
        bus_baseKV = bus_data[:, 9]
        bus_zone = bus_data[:, 10]
        bus_Vmax = bus_data[:, 11]
        bus_Vmin = bus_data[:, 12]
        bus_loaded = Bus(bus_num, bus_type, bus_Pd, bus_Qd, bus_Gs, bus_Bs, bus_area, bus_Vm, bus_Va, bus_baseKV,
                         bus_zone, bus_Vmax, bus_Vmin)
        return bus_loaded

    def loadgen(self,gen_data):
        gen_bus = gen_data[:, 0].astype(int)
        gen_Pg = gen_data[:, 1]
        gen_Qg = gen_data[:, 2]
        gen_Qmax = gen_data[:, 3]
        gen_Qmin = gen_data[:, 4]
        gen_Vg = gen_data[:, 5]
        gen_mBase = gen_data[:, 6]
        gen_status = gen_data[:, 7]
        gen_Pmax = gen_data[:, 8]
        gen_Pmin = gen_data[:, 9]
        gen_Pc1 = gen_data[:, 10]
        gen_Pc2 = gen_data[:, 11]
        gen_Qc1min = gen_data[:, 12]
        gen_Qc1max = gen_data[:, 13]
        gen_Qc2min = gen_data[:, 14]
        gen_Qc2max = gen_data[:, 15]
        gen_ramp_agc = gen_data[:, 16]
        gen_ramp_10 = gen_data[:, 17]
        gen_ramp_30 = gen_data[:, 18]
        gen_ramp_q = gen_data[:, 19]
        gen_apf = gen_data[:, 20]
        gen_loaded = Gen(gen_bus, gen_Pg, gen_Qg, gen_Qmax, gen_Qmin, gen_Vg, gen_mBase, gen_status, gen_Pmax, gen_Pmin,
                         gen_Pc1, gen_Pc2, gen_Qc1min, gen_Qc1max, gen_Qc2min, gen_Qc2max, gen_ramp_agc, gen_ramp_10,
                         gen_ramp_30, gen_ramp_q, gen_apf)
        return gen_loaded

    def loadbranch(self,branch_data):
        branch_fbus = branch_data[:, 0].astype(int)
        branch_tbus = branch_data[:, 1].astype(int)
        branch_r = branch_data[:, 2]
        branch_x = branch_data[:, 3]
        branch_b = branch_data[:, 4]
        branch_rateA = branch_data[:, 5]
        branch_rateB = branch_data[:, 6]
        branch_rateC = branch_data[:, 7]
        branch_ratio = branch_data[:, 8]
        branch_angle = branch_data[:, 9]
        branch_status = branch_data[:, 10]
        branch_angmin = branch_data[:, 11]
        branch_angmax = branch_data[:, 12]
        branch_loaded = Branch(branch_fbus, branch_tbus, branch_r, branch_x, branch_b, branch_rateA, branch_rateB,
                               branch_rateC, branch_ratio, branch_angle, branch_status, branch_angmin, branch_angmax)
        return branch_loaded

    def loadgencost(self,gencost_data):
        gencost_type = gencost_data[:, 0]
        gencost_startup = gencost_data[:, 1]
        gencost_shutdown = gencost_data[:, 2]
        gencost_loaded = Gencost(gencost_type, gencost_startup, gencost_shutdown)
        return gencost_loaded

class Bus:
    def __init__(self, num, type, Pd, Qd, Gs, Bs, area, Vm, Va, baseKV, zone, Vmax, Vmin):
        self.num = num
        self.type = type
        self.Pd = Pd
        self.Qd = Qd
        self.Gs = Gs
        self.Bs = Bs
        self.area = area
        self.Vm = Vm
        self.Va = Va
        self.baseKV = baseKV
        self.zone = zone
        self.Vmax = Vmax
        self.Vmin = Vmin

class Gen:
    def __init__(self,bus, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin, Pc1, Pc2, Qc1min, Qc1max, Qc2min, Qc2max, ramp_agc, ramp_10, ramp_30, ramp_q, apf):
        self.bus = bus
        self.Pg = Pg
        self.Qg = Qg
        self.Qmax = Qmax
        self.Qmin = Qmin
        self.Vg = Vg
        self.mBase = mBase
        self.status = status
        self.Pmax = Pmax
        self.Pmin = Pmin
        self.Pc1 = Pc1
        self.Pc2 = Pc2
        self.Qc1min = Qc1min
        self.Qc1max = Qc1max
        self.Qc2min = Qc2min
        self.Qc2max = Qc2max
        self.ramp_agc = ramp_agc
        self.ramp_10 = ramp_10
        self.ramp_30 = ramp_30
        self.ramp_q = ramp_q
        self.apf = apf

class Branch:
    def __init__(self, fbus, tbus, r, x, b, rateA, rateB, rateC, ratio, angle, status, angmin, angmax):
        self.fbus = fbus
        self.tbus = tbus
        self.r = r
        self.x = x
        self.b = b
        self.rateA = rateA
        self.rateB = rateB
        self.rateC = rateC
        self.ratio = ratio
        self.angle = angle
        self.status = status
        self.angmin = angmin
        self.angmax = angmax

class Gencost:
    def __init__(self, type, startup, shutdown):
        self.type = type
        self.startup = startup
        self.shutdown = shutdown

### class add the power flow data for branch
class GenericModel_add(GenericModel):
    def loadbranch(self,branch_data):
        branch_fbus = branch_data[:, 0].astype(int)
        branch_tbus = branch_data[:, 1].astype(int)
        branch_r = branch_data[:, 2]
        branch_x = branch_data[:, 3]
        branch_b = branch_data[:, 4]
        branch_rateA = branch_data[:, 5]
        branch_rateB = branch_data[:, 6]
        branch_rateC = branch_data[:, 7]
        branch_ratio = branch_data[:, 8]
        branch_angle = branch_data[:, 9]
        branch_status = branch_data[:, 10]
        branch_angmin = branch_data[:, 11]
        branch_angmax = branch_data[:, 12]
        branch_Pf = branch_data[:, 13]
        branch_Qf = branch_data[:, 14]
        branch_Pt = branch_data[:, 15]
        branch_Qt = branch_data[:, 16]
        branch_loaded = Branch_add(branch_fbus, branch_tbus, branch_r, branch_x, branch_b, branch_rateA, branch_rateB,
                               branch_rateC, branch_ratio, branch_angle, branch_status, branch_angmin, branch_angmax,branch_Pf,branch_Qf,branch_Pt,branch_Qt)
        return branch_loaded

class Branch_add(Branch):
    def __init__(self, fbus, tbus, r, x, b, rateA, rateB, rateC, ratio, angle, status, angmin, angmax,Pf,Qf,Pt,Qt):
        self.fbus = fbus
        self.tbus = tbus
        self.r = r
        self.x = x
        self.b = b
        self.rateA = rateA
        self.rateB = rateB
        self.rateC = rateC
        self.ratio = ratio
        self.angle = angle
        self.status = status
        self.angmin = angmin
        self.angmax = angmax
        self.Pf = Pf
        self.Qf = Qf
        self.Pt = Pt
        self.Qt = Qt
