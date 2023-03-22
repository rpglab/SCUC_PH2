# This code is used to form the pyomo data file for the Unit Commitment
# Jin Lu, University of Houston
## link: https://rpglab.github.io/resources/

def pyomodata_UC_simplified(file_name,GenericModel,time_totnum,load_b_t,reg_gen_list):    # load_b_t[i][j] where i,j start from 0
    bus_str = 'param: BUS : bus_num := \n'
    gen_str = 'param: GEN: gen_num gen_bus gen_cost_P gen_cost_NL gen_cost_SU gen_Pmin gen_Pmax gen_r10 gen_rr := \n'
    line_str = 'param: LINE: line_num line_x line_Pmax line_fbus line_tbus := \n'
    time_str = 'param: TIME: time_num := \n'
    bus_time_str = 'param: BUS_TIME: bustime_num load_b_t:= \n'
    semicolon_str = ';\n'
    gen_totnum = GenericModel.gentotnum
    bus_totnum = GenericModel.bustotnum
    line_totnum = GenericModel.branchtotnum
    bustime_totnum = bus_totnum*time_totnum
    f = open(file_name, 'w')
    # Write bus data
    f.write(bus_str)
    i = 0
    while(i<bus_totnum):
        bus_num = i + 1     # In pyomo, start from 1
        busdata_i = str(bus_num)+' '+str(bus_num) + '\n'
        f.write(busdata_i)
        i+=1
    f.write(semicolon_str)
    # Write gen data
    f.write(gen_str)
    for i_num in reg_gen_list:
        #i = 10*origin_i + 1 # reduce the generators total size
        i = i_num - 1
        gendata_i = str(i_num)+' '+str(i_num)+' '+str(GenericModel.gen.bus[i])+' '+str(GenericModel.gen.c1[i])+' '
        gendata_i += str(GenericModel.gen.c0[i])+' '+str(GenericModel.gen.c_su[i])+' '+str(0)+' '    # Pmin is set to zero
        # ramp_10 is MW/10min, ramp_agc is MW/min, gen_rr is MW/hr
        gendata_i += str(GenericModel.gen.Pmax[i])+' '+str(GenericModel.gen.ramp_agc[i]*10)+' '+str(GenericModel.gen.ramp_agc[i]*60)+'\n '
        f.write(gendata_i)
    f.write(semicolon_str)
    # Write line data
    f.write(line_str)
    i = 0
    while (i < line_totnum):
        line_num = i + 1
        linedata_i = str(line_num)+' '+str(line_num) + ' ' + str(GenericModel.branch.x[i]) + ' ' + str(GenericModel.branch.rateA[i]) + ' '  # rateA is MVA
        linedata_i += str(GenericModel.branch.fbus[i]) + ' ' + str(GenericModel.branch.tbus[i]) + '\n'
        f.write(linedata_i)
        i += 1
    f.write(semicolon_str)
    # Write time data
    f.write(time_str)
    i = 0
    while (i < time_totnum+1):  # 24 hour and initial hour
        time_num = i
        timedata_i = str(time_num)+' '+str(time_num) + '\n'
        f.write(timedata_i)
        i += 1
    f.write(semicolon_str)
    # Write Bus_Time data
    f.write(bus_time_str)
    i = 0
    k = 0
    while (i < bus_totnum):
        bus_num = i + 1
        j = 0
        while (j<time_totnum+1):    # 24 hour and initial hour
            time_num = j
            bustime_num = k + 1
            # initial hour load data
            if j==0:
                bustime_data_i = str(bus_num)+' '+str(time_num)+' '+str(bustime_num)+' '+str(0)+'\n'
            else:
                bustime_data_i = str(bus_num)+' '+str(time_num)+' '+str(bustime_num)+' '+str(load_b_t[i][j-1])+'\n'
            f.write(bustime_data_i)
            j += 1
            k += 1
        i += 1
    f.write(semicolon_str)

    f.close()
