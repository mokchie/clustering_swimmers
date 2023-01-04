from math import *
l_ref = rc = 1.6
eta_ref = 6.25
d0 = 6.25
rho_ref = 6.25

m_ref = rho_ref*rc**2
k_ref = 1/rc
t_ref = rc**2*rho_ref/eta_ref
v_ref = rc/t_ref
f_ref = rc**3*rho_ref/t_ref**2
p_ref = rc**2*rho_ref/t_ref**2
zeta_ref = rc**2*rho_ref/t_ref**2
zeta_theta_ref = rc**4*rho_ref/t_ref**2
kappa_ref = rc**5*rho_ref/t_ref**2
I_ref = rho_ref*rc**4
E_ref = f_ref*rc
omega_ref = 1/t_ref

def ndm_mass(x):
    return x/m_ref
def ndm_length(x):
    return x/rc
def ndm_k(x):
    return x/k_ref
def ndm_v(x):
    return x/v_ref
def ndm_omega(x):
    return x/omega_ref
def ndm_rho(x):
    return x/rho_ref
def ndm_time(x):
    return x/t_ref
def ndm_force(x):
    return x/f_ref
def ndm_grad_force(x):
    return x/(f_ref/rc)
def ndm_eta(x):
    return x/eta_ref
def ndm_pressure(x):
    return x/p_ref
def ndm_zeta(x):
    return x/zeta_ref
def ndm_zt(x):
    return x/zeta_theta_ref
def ndm_kappa(x):
    return x/kappa_ref
def ndm_I(x):
    return x/I_ref
def ndm_E(x):
    return x/E_ref
def ndm_angle(x):
    return x/180



