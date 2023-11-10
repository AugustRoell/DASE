"""
Data on the relevant material properties of a number of materials
"""


# Sapphire (alpha Aluminium Oxide, Al2O3)
Sapphire = {
    "name": r"alpha - Al_2 O_3",
    "f_0": 6.29e+9,
    "v_l": 1.120e+4,
    "v_t_eff": 9.208e+3,
    "kappa": 0.0,
    "rho": 3980.0,
    "L": 1.0e-3,
}

# Quartz (alpha Silicon Oxide)
Quartz = {
    "name": r"alpha - Si O_2",
    "f_0": 12.645e+9,
    "v_l": 6.338e+3,
    "v_t_eff": 7.782e+3,
    "kappa": 0.0,
    "rho": 2648.5,
    "L": 1.0e-3,
}

# Calcium Fluoride
CaF2 = {
    "name": r"Ca F_2",
    "f_0": 13.303e+9,
    "v_l": 7.204e+3,
    "v_t_eff": 5.134e+3,
    "kappa": 0.0,
    "rho": 3179.0,
    "L": 1.0e-3,
}


available_materials = {
    "Sapphire": Sapphire,
    "Quartz": Quartz,
    "CaF2": CaF2,
}