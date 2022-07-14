# function to run all types of SCUC (T-SCUC, R-SCUC, EH-SCUC, H-SCUC) under different wind penetration levels
# Author: Jin Lu (University of Houston)
#    Website: https://rpglab.github.io/resources/SCUC-H2-Python/


# Import
from UC_nohdg_cases import run_T_SCUC,run_R_SCUC
from UC_hdg_cases import run_EH_SCUC,run_H_SCUC

# Run UC & Save the results for different wind penetration level
# The wind penetration level is 10%, 20%, 30%, 40% and 50%
def run_SCUC_wdlv():
    step = 0.1         # step size for wind power penetration level is 'step' (e.g. 10%)
    lv_init = 0.3      # initial wind power penetration level (first step, e.g., 10%)
    lv_totnum = 1      # Run 'lv_totnum' (e.g. 5) scenarios of different wind power penetration levels
    for s in range(lv_totnum):
        wd_lv = lv_init + s * step
        wd_lv = round(wd_lv, 2)
        run_T_SCUC(wd_lv)        # traditional SCUC model
        run_R_SCUC(wd_lv)        # relaxed SCUC model - ignoring network/line capacity limit
        run_EH_SCUC(wd_lv)       # SCUC considering LOCAL "P2H - Hydrogen Storage - H2P"
        run_H_SCUC(wd_lv)        # SCUC considering a hydrogen network where P2H, Hydrogen Storage,
                                      # and H2P may not be in the same location.


# Example: run all types of SCUC (T-SCUC, R-SCUC, EH-SCUC, and H-SCUC) under different wind penetration levels
run_SCUC_wdlv()


