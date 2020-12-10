from pyrad.lbl.tips import TotalPartitionFunction
from pyrad.lbl.hitran import Hitran, Voigt
from pyrad.optics.gas import Gas

for x in ["H2O", "CO2", "O3"]:
    Hitran(x, Voigt()).create_database("hitran.sqlite")
    TotalPartitionFunction(x).create_database("tips-2017.sqlite")
gas = Gas("H2O", hitran_database="hitran.sqlite", tips_database="tips-2017.sqlite")
gas = Gas("CO2", hitran_database="hitran.sqlite", tips_database="tips-2017.sqlite")
gas = Gas("O3", hitran_database="hitran.sqlite", tips_database="tips-2017.sqlite")
