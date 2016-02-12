#
# calculations for returning beer recipies
# 6Feb2016
# Ashlee J
# References
# http://homebrewmanual.com/home-brewing-calculations/
###########
## clear workspace
rm(list = ls())
ls()

########
# specify directory
# read in file with desired parameters
# Need:
# 1. the final batch volumne
# 2. the preferred gravity
# 3. the original gravity
# 4. the malts used and percentage of each
# 5. the IBU target

################################################################
# Calculation parameters

###############
# Gravity
# Gravity points = (or.grav -1)*1000*volume

###############
# Extract potential/

###############
# Water volume

# AAUs

# IBUs

# mashing


# Multi-rest Infustion Equation
# Palmer bok equation
G.lb = 13.5
G.kg <- G.lb*0.453592
T1.degC = 56
Wm.L = 15
T2.degC = 65

T.strike.degC = 100
Wa.lbL = (T2.degC-T1.degC)*(0.2*G.lb+Wm.L)/(T.strike.degC-T2.degC)



