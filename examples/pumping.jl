using PhysicalConstants: CODATA2018 as PCs
ee = PCs.e.val
kb = PCs.k_B.val
mp = PCs.m_p.val
kb_ee = kb / ee
S_pump = 40
C_pump = 30
r_pump = 18000
R_omp = 4.0
R_pump = 5.0
λ_plasma = 0.002
p_pump = 0.004
I_pump = 7e19 * p_pump * (S_pump + C_pump)
I_max = r_pump * 0.004 * 7e19
Γ_0 = I_max / (2 * pi * R_omp * λ_plasma) #I_pump/
n_up = 1e20
T_up = 200
P_up = n_up * ee * T_up
4 * mp * Γ_0^2 * ee * 200 / P_up^2
Q = Γ_0 * ee * T_up
sqrt(ee * T_up / mp)

10^20

## simple scaling for pumping power


