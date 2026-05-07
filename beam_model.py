# import packages
import os
import numpy as np



# constants

# Units: dimensionless

# Units: m/s**2
g = 9.80665 # standard gravity


# Units: J/Kg

LS0 = 2.836e6  # latent heat of sublimation (-100<=T<= 0-deg-C)
LV0 = 2.5009e6 # latent heat of vaporization at 0-deg-C
LF0 = LS0-LV0  # latent heat of fusion at 0-deg-C (lv0-ls)
CI = 2.097e3   # heat capacity of ice 
CL = 4.219e3   # heat capacity of liquid water (above freezing)

RD= 287.04 # gas constant of dry air
RV = 461.50 # gas constant of water vapor
CVD = 719. # heat capacity at constant volume for dry air
CPD = 1005.7 # heat capacity at constant pressure for dry air
CVV = 1410. # heat capacity at constant volume of water vapor
CPV = 1870. # heat capacity at constant pressure of water vapor
CPVCMCL = CL-CPV # cpvmcl seems to be a common notation for this value
EPS = RD/RV

# Model constants
T0 = 273.16   # triple point temperature (K)
p0 = 1000e2

# User specified input
DT_HOM = 40 # T difference between triple point and homogenous ice T - determine linear partition of mixed-phase between ice and supercooled water in lam(T) function
T_END = 220 # controls exopnential decay of condensate given by exp(beta/(T0-Te_end)*dT)
# theta_e functions

def lam(T,dT_hom = DT_HOM):
    """Ice fraction of condensate 
    
    Units: dimensionless
        
    Args:
        T: air temperature in K.
        
    Returns:
        f (fraction of ice)
    """    
    T0 = 273.16
    T = np.asarray(T, dtype=float)
    f = np.zeros(T.shape)
    f[T<=T0-dT_hom] = 1
    f[np.logical_and(T<T0, T>T0 - dT_hom)] =  -(T[np.logical_and(T<T0, T>T0 - dT_hom)]-T0)/dT_hom
    return f

def lf(T):
    """Latent heat of fusion

    Units: J/kg
    Args:
        T (double): temperature
    
    """
    return ls(T) - lv(T)

def lv(T,lv0=LV0,cl=CL,cpv=CPV):
    """Latent heat of vaporization

    Units: J/kg
    Args:
        T (double): temperature
    
    """
    return lv0-(cl-cpv)*(T-T0)

def ls(T,ls0=LS0,ci=CI,cpv=CPV):
    """Latent heat of sublimation

    Units: J/kg
    Args:
        T (double): temperature
    
    """
    return ls0-(ci-cpv)*(T-T0)


def eval_q(e, p, qt=0,eps=EPS):
    """Calculate specific humidity.
    
    Units: kg/kg (unitless).
    
    Args:
        e (double): vapor pressure in Pa.
        p (double): air pressure in Pa.
        EPS (constant): ratio of gas constant of dry air to gas constant of water vapor.
        
    Returns:
        double: specific humidity.
    """
    # q = (EPS * e) / (p + ((EPS-1.)*e))

#     if nargin==2 #assuming condensate-free
#         q = EPS * e / ( p-e*(1-EPS) );
#     else #given qt (suitable for plume computations)
    q = eps* e * (1-qt) / (p-e)
    #q = EPS * e / ( p-e*(1-EPS) )
    return q
    
def eval_es(T,lv0 = LV0, cl= CL, cpv=CPV):
    """Calculate saturation vapor pressure with respect to liquid water.
    
    Units: Pa.
    
    Using 2 formulae for different temperature ranges (separated by -80-deg-C).
    For T>-80-deg-C, the difference between the 2 formulae is within 0.5#.
    
    Args:
        T (double): air temperature in K.
        
    Returns:
        double: saturation vapor pressure of water (liquid).
    """
#     Tc = T-273.15
    
#     if Tc < -80.: # Is < a typo? Using Bolton
#         return eval_es_bolton(T)
#     else:
#         return 0.6105851e+03 + Tc *( 
#                0.4440316e+02 + Tc *(
#                0.1430341e+01 + Tc *( 
#                0.2641412e-01 + Tc *( 
#                0.2995057e-03 + Tc *( 
#                0.2031998e-05 + Tc *(
#                0.6936113e-08 + Tc *(
#                0.2564861e-11 - Tc * 0.3704404e-13)))))))
    es0 = 611.2
    lv1 = lv0 + (cl-cpv)*T0
    lv2 = cl-cpv
    es = 1e2 * np.exp( 53.67957 - 6743.769/T - 4.8451*np.log(T) );

    return es#es0*np.exp(lv1/rv*(T-T0)/T0/T - lv2/rv*np.log(T/T0))


def eval_esi(T,ls0=LS0,ci=CI,cpv=CPV):
    """Calculate ice saturation vapor pressure.
    
    Units: Pa.
    
    Accuracy within 0.14# for -80<=T<=0-deg-C.
    
    Args:
        T (double): air temperature in K.
        
    Returns:
        double: ice saturation vapor pressure following Eq. (4.4.15) in Emanuel (1994).
    """    
    es0 = 611.2
    ls1 = ls0 + (ci-cpv)*T0
    ls2 = ci-cpv
    return 100. * np.exp( 23.33086 - 6111.72784/T + 0.15215*np.log(T) )



def eval_qs(T,p,qt=0):
    """Effective qsat by partial pressures of ice and liquid
    Args:
        T: temperature in K
        p: pressure in Pa

    Returns:
        qsat (kg/kg)
    """
    es = eval_es(T)
    esi = eval_esi(T)
    # qs_l = eval_q(eval_es2(T),p)
    # qs_i = eval_q(eval_esi2(T),p)
    return (1-lam(T))*eval_q(es,p,qt) + lam(T)*eval_q(esi,p,qt)

def eval_theta_e(T,qv,ql,qi,p,rd=RD,rv=RV,cpd=CPD,cl=CL,ci=CI):
    """  Calculate equivalent potential temperature.

      Units: K.

      Following Eq. (2.67) in Stevens & Siebesma (2020)
       which, with ice correction, is equivalent to Eq. (4.5.11) in Emanuel (1994).

      Args:
        T (double): air temperature in K.
        qv (double): specific humidity of water vapor in kg/kg.
        ql (double): specific liquid water content in kg/kg.
        qi (double): specific liquid water content in kg/kg.
        p (double): air pressure in Pa.

      Returns:
        double: equivalent potential temperature.
        #}"""
    T0 = 273.16
    p0 = 1000e2
    qt = qv+ql+qi
    Re = (1-qt)*rd
    eps = rd/rv
    # e = p*qv/(eps + qv*(1-eps) )
    qc = ql + qi
    e = p*qv/(qv+eps*(1-qt))
    # e = p*qv/((1-eps)*qv + eps*(1-ql-qi))
    es = (1-lam(T))*eval_es(T) + lam(T)*eval_esi(T)
    cpl = cpd + (cl - cpd)*qt
    R = Re + qv*rv;
    omega_e = (R/Re)**(Re/cpl)*(es/e)**(qv*rv/cpl)
    chi = Re/cpl
    gammi = -(cl-ci)*(qi)/cpl

    theta_e = T * (p0/p)**chi*omega_e* (T/T0)**gammi* np.exp(qv*lv(T)/T/cpl - lf(T)*qi/T0/cpl)
    return theta_e

def q_scheme(T,p,qt):
    """Partition total water between phases given T and p
    
    Units: kg/kg
        
    Args:
        T: air temperature in K.
        p: pressure in Pa 
        
    Returns:
        qv: vapor (kg/kg)
        ql: liquid water (kg/kg)
        qi: ice (kg/kg)
    """    
    qs = eval_qs(T,p,qt)
    if qs<0:
        qs = 0
        # raise Exception("q_sat is below 0")
    qc = qt - qs
    if qc>0:
        qv = qs
        qi = lam(T)*qc
        ql = (1-lam(T))*qc
    else:
        qv = qt
        ql = 0
        qi = 0
    qv = qv
    qt1 = qv + qi + ql
    return qv,ql,qi


def inv_T_s(s0, p, qt, T_guess=280, abs_tol=1e-5, max_iter=10000):
    """  Invert T from specific entropy and determine water partition

      Units: K.

      Args:
        s0 (double): entropy in J/Kg/K
        p (double): air pressure in Pa.
        qt (double): total water in kg/kg
        T_guess: initial T guess for Newton scheme
        abs_tol: tolerance of convergence for Newton scheme
        max_iter: max number of iterations for scheme to converge
      Returns:
        T (double): temperature in K
        qv: vapor (kg/kg)
        ql: liquid water (kg/kg)
        qi: ice (kg/kg)
        #}"""
    T = T_guess
    for n in range(max_iter):
        if np.isnan(T):
            T = T_guess
            qv, ql, qi = q_scheme(T, p, qt)

            break
        
        qv, ql, qi = q_scheme(T, p, qt)
        s_1 = comp_s(T, p, qv, ql, qi)

        # finite-diff derivative 
        dT = 1e-3
        s_hi = comp_s(T + dT, p, *q_scheme(T + dT, p, qt))
        s_lo = comp_s(T - dT, p, *q_scheme(T - dT, p, qt))
        dT_s = (s_hi - s_lo) / (2 * dT)

        # Newton step
        err = (s_1 - s0) / dT_s
        T_new = T - err

        if abs(err) < abs_tol:
            break
        T = T_new
        
    return T, qv, ql, qi


def comp_s(T,p,qv,ql,qi,rv = RV,rd=RD,cl=CL,ci=CI,Tr=T0):
    """  Composite specfic entropy
    s = (1-qt)*sd + qt*sl + qv*(sv - sl)  - qi*(sl - si) + qv*(sv - ssv)
    from eqn 2.64 in Stevens & Siebesma (2020)
      Units: J/kg/K.

      Args:
        T (double): temperature in K
        p (double): air pressure in Pa.
        qv (double): vapor in kg/kg
        ql (double): liquid in kg/kg
        qi (double): ice in kg/kg
      Returns:
        specific entropy (J/kg/K)
        #}"""
    # Tr = T0
    p0 = 1000e2
    eps = rd/rv

    qt = qv + ql + qi
    e = p*qv/(qv+eps*(1-qt))
    es = (1-lam(T))*eval_es(T) + lam(T)*eval_esi(T)
    re = (1-qt)*rd
    r = re + rv*qv
    sd = cpd*np.log(T/Tr) - rd*np.log(p/p0*re/r)
    svsl = lv(T)/T
    svssv = rv*np.log(es/e)
    slsi = (cl-ci)*np.log(T/T0) + lf(T)/T0
    sl = cl*np.log(T/Tr)
    s = (1-qt)*sd + qt*sl + qv*svsl  - qi*slsi + qv*svssv
    return s


def beam_model_calc(T_init,p,qt,T_env,q_env,mix_coeff,gam=1,abs_tol = 1e-3,loss_type = 'fractional',cap=0,beta = np.log(10e3)/(T0 - T_END)):
    """Entraining plume calculation, assuming pressure decreasing with index.
    
    Input: 
    T_init: initial temperature at p[0] in K
    p: pressure (Pa)
    qt: total moisture (qt)
    T_env: vertical profile of environmental temperature in K
    q_env: vertical profile of environmental moisture in K
    mix_coeff: mixing coefficients 
    gam: fractional condensate retention rate (between 0 and 1)
    abs_tol: tolerance for convergence in inv_T_s solver
    loss_type: options are 'fractional', 'cap', or'decay'
            1. 'fractional': fractional condenate retention by gam
            2. 'cap': condensate above cap is removed
            3. 'decay': fractional condenate retention by gam below the freezing level, above freezing level, total condensate exponentially decays                  according to exp(beta*dT).
                - beta nominally gives 1/1000 of freezing level condensate at T0 - T_end
    Output: ep: plume theta_e (K)
    Tp: plume T (K)
    qv: plume moisture
    ql: plume liquid water
    qi: plume ice
    """
    Tp = np.zeros(p.shape)*np.nan
    ep = np.zeros(p.shape)*np.nan
    qv = np.zeros(p.shape)
    ql = np.zeros(p.shape)
    qi = np.zeros(p.shape)
    s = np.zeros(p.shape)
    # env measures
    s_env = comp_s(T_env,p,q_env,0,0)
    # initialize quantities
    Tp[0] = T_init
    qv[0],ql[0],qi[0] = q_scheme(Tp[0],p[0],qt)

    s[0] = comp_s(Tp[0],p[0],qv[0],ql[0],qi[0])
    
    cpl = cpd + (cl - cpd)*qt
    ep[0] = T0*np.exp(s[0]/cpl)
    if gam==0:
        loss_type = 'fractional'
            
    for a,i in tqdm(enumerate(p[:-1])):
        j = p[a+1]
        dp = (j - i)
        s[a] = mix_coeff[a]*(s_env[a] - s[a])*dp + s[a] 
        qt = mix_coeff[a]*(q_env[a] - qt)*dp + qt
        
        Tp[a+1],qv_pre,ql_pre,qi_pre = inv_T_s(s[a],p[a+1],qt,T_guess = Tp[a],abs_tol=1e-6)
    
        dqc = 0

        if loss_type=='fractional':
            dql = (1-gam)*(ql_pre - ql[a])
            dqi = (1-gam)*(qi_pre - qi[a])
            dqc = dql + dqi
            ql[a+1] = ql[a] + gam*(ql_pre - ql[a])
            qi[a+1] = qi[a] + gam*(qi_pre - qi[a])

        elif loss_type=='decay':
            if Tp[a+1]>T0:
                dql = (1-gam)*(ql_pre - ql[a])
                dqi = (1-gam)*(qi_pre - qi[a])
                dqc = dql + dqi
                ql[a+1] = ql[a] + gam*(ql_pre - ql[a])
                qi[a+1] = qi[a] + gam*(qi_pre - qi[a])
            else:
                # orginal
                # dT = Tp[a+1] - Tp[a]
                # ql[a+1]  = (1+beta*dT)*ql_pre
                # qi[a+1] = (1+beta*dT)*qi_pre
                # dql = abs(ql[a+1] - ql_pre)
                # dqi = abs(qi[a+1] - qi_pre)
                # dqc = dql + dqi
                dT = Tp[a+1] - Tp[a]          # negative when cooling upward
                dT_cool = min(dT, 0.0)
                fac = np.exp(beta*dT_cool)       # beta > 0 gives fac < 1 for dT<0
                ql[a+1] = ql_pre * fac
                qi[a+1] = qi_pre * fac
                dql = ql_pre - ql[a+1]
                dqi = qi_pre - qi[a+1]
                dqc = dql + dqi
                ########################### FIX ###################
                # cold-phase: exponential decay with temperature drop
                # dT = Tp[a+1] - Tp[a]
        
                # # Only allow fallout when cooling (dT <= 0). If warming, don't "anti-fallout".
                # dT_cool = min(dT, 0.0)
        
                # # Retention factor (0<fac<=1). Using exp makes it step-size consistent.
                # fac = np.exp(beta * dT_cool)
        
                # # Retain a fraction of the *equilibrated* condensate at this level
                # ql_ret = fac * ql_pre
                # qi_ret = fac * qi_pre
        
                # # Removed condensate (guaranteed >=0)
                # dql = ql_pre - ql_ret
                # dqi = qi_pre - qi_ret
                # dqc = dql + dqi
        
                # # Store retained condensate
                # ql[a+1] = ql_ret
                # qi[a+1] = qi_ret
    
        elif loss_type=='cap':
            qc = ql[a+1] + qi[a+1]
            if qc>cap:
                qc = cap
                dqc = ql[a+1] + qi[a+1] - cap
                ql[a+1] = cap*(1-lam(Tp[a+1]))
                qi[a+1] = cap*(lam(Tp[a+1]))


        
        s_pre = comp_s(Tp[a+1],p[a+1],qv_pre,ql_pre,qi_pre)

        #entropies
        sl = cl * np.log(Tp[a+1]/Tr)
        slsi = (cl-ci)*np.log(Tp[a+1]/T0) + lf(Tp[a+1])/T0#Tp[a+1]
        
        # s that is removed
        ds = sl*(dql+dqi) - slsi*dqi 

        s_new = (s_pre - ds)/(1-dqc)

        # boosted moisture quantities
        qv[a+1] = qv_pre/(1-dqc)
        ql[a+1] = ql[a+1]/(1-dqc)
        qi[a+1] = qi[a+1]/(1-dqc)
        qt = qv[a+1] + ql[a+1] + qi[a+1]
        
        cpl_new = cpd + (cl-cpd)*qt
        
        Tp[a+1],qv[a+1],ql[a+1],qi[a+1] = inv_T_s(s_new,p[a+1],qt,T_guess = Tp[a+1],abs_tol=1e-5)

        ep[a+1] = T0*np.exp(s_new/cpl_new)
        s[a+1] = s_new

    return ep,Tp,qv,ql,qi