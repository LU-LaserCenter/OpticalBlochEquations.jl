using Parameters

"""
J1: quantum number of ground state total electron angular momentum
J2: quantum number of excited state total electron angular momentum
nucI: quantum number of total nuclear angular momentum
exFmax: maximum quantum number of excited state total atomic angular momentum
exFmin: minimum quantum number of excited state total atomic angular momentum
grFmax: maximum quantum number of ground state total atomic angular momentum
grFmin: minimum quantum number of ground state total atomic angular momentum
S1: ground state electron spin
S2: excited state electron spin
L1: quantum number of ground state orbital angular momentum
L2: quantum number of excited state orbital angular momentum
gI: nuclear Lande factor
A_gr: ground state magnetic dipole constant
A_ex: excited state magnetic dipole constant
B_ex: excited state electric quadrupole constant
Efs: energy difference between the ground and excited level of the induced hyperfine structure transition 
"""
@with_kw struct cezijsD1
    J1::Real = 1/2
    J2::Real = 1/2
    nucI::Real = 7/2
    exFmax::Real = nucI+J2
    exFmin::Real = nucI-J2
    grFmax::Real = nucI+J1
    grFmin::Real = nucI-J1
    S1::Real = 1/2
    S2::Real = 1/2
    L1::Real = J1-S1
    L2::Real = J2+S2 

    gI::Float64 = -0.00039885395 

    #Sīkstruktūras konstantes
    A_gr::Float64 = 2298.1579425  

    A_ex::Float64 = 291.9201  
    B_ex::Float64 = 0.0 

    #Sīkstruktūras enerģijas starpība
    Efs::Float64 = 335116048.807 

end

@with_kw struct cezijsD2
    J1::Real = 1/2
    J2::Real = 3/2
    nucI::Real = 7/2
    exFmax::Real = nucI+J2
    exFmin::Real = nucI-J2
    grFmax::Real = nucI+J1
    grFmin::Real = nucI-J1
    S1::Real = 1/2
    S2::Real = 1/2
    L1::Real = J1-S1
    L2::Real = J2+S2 

    gI::Float64 = -0.00039885395

    #Sīkstruktūras konstantes
    A_gr::Float64 = 2298.1579425 

    A_ex::Float64 = 50.28827  
    B_ex::Float64 = -0.4934

    #Sīkstruktūras enerģijas starpība
    Efs::Float64 = 351725718.50

end

