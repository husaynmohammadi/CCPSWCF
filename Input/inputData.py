"C-PSW/CF Section: Shafaei PP=136"
import sys
# from functions.ClassComposite import compo
exec(open("MAIN.py").readlines()[18]) # It SHOULD read and execute exec(open(f"Input/units{'US'}.py").read())
# exec(open("../Input/unitsSI.py").read()) # It SHOULD read and execute exec(open("Input/units    .py").read())
# exec(open("MAIN.py").readlines()[20]) # It SHOULD read and execute exec(open("Input/materialParameters.py").read())

#=============================================================================
#    Design Inputs
#=============================================================================
# Seismic Coefficients
Cd      = 5.5
Ie      = 1.
R       = 8
Omega0  = 2.5
Rho     = 1.
SDC     = "Dmax"  # "Dmax", "Dmin"

    

# Composite wall resistance factor
Fi_v = Fi_b = Fi_c = Fi_t = 0.9



#=============================================================================
#    Material Properties
#=============================================================================
# Material
Es      = 200       *GPa
Gs      = 77.2      *GPa
Fy      = 50 *ksi
Fu      = 65 *ksi
Ry      = 1.1

fpc     = 6 *ksi
Ec      = 4500 *ksi
Gc      = 1800 *ksi
Rc      = 1.3

linearity = 1

#_______________________IMK+Pinching_Hinge_Properties________________________#
# K0          = 12 *EIeff /L **3 Should be given in the FuncModel.py
C_K0        = 70
My_Plus     = 1250 *kN*m  *1.27         # Pick Strength
My_Neg      = -1 *My_Plus
as_Plus     = as_Neg      = 0.0002      # Slope of 2nd line
Res_Pos     = Res_Neg     = 0.2         # Residual Strength
theta_p_Plus= theta_p_Neg = 0.011       # Top width (2nd line width)
theta_pc_Plus=theta_pc_Neg= 0.04        # 3rd line width (or slope)
FprPos      = FprNeg      = A_pinch     = 0.3
Lamda_S     = Lamda_C     = Lamda_A     = Lamda_K     = 1.0
c_S         = c_C         = c_A         = c_K         = 1
theta_u_Plus= theta_u_Neg = 0.3 
D_Plus      = D_Neg       = 0.2

#_______________________Steel02_ShearingHinge_Properties_____________________#
slope = -0.0004
R0,cR1,cR2= 18.5, 0.91, 0.06
a1=a3= 0.00001
a2=a4= 0.7
c_ktrans = 0.516858918
b1=1/c_ktrans * slope
c_krot   = 100
#=============================================================================
#    Elements
#=============================================================================
Hw          = 2800 *mm
L_CB        = 2700 *mm
L_eff       = L_CB +Hw
# L_eff       = 5150 *mm
# L_CB        = L_eff -Hw
# L_CB        = 6.5 *m -2 *Hw
tw          = 12 *mm
RhoW        = 0.04
bf          = 1 /RhoW *tw *2
tc          = bf -2 *tw
tf          = tw
# tf          = 10 *mm
t_sc        = bf
# RhoW        = 2 *tw /t_sc

btie        = 240 *mm # Vertical Spacing
Stie        = 250 *mm # Horizontal Spacing
dtie        = 25 *mm
lsr         = btie /tw

typeSect    = "Composite" # Composite, I_Shaped
ductility   = "high"  # moderate, high
t_pwCB      = 10 *mm
t_pfCB      = 10 *mm
H_CB        = 420 *mm
# tc_CB       = tc
# bf_CB       = bf -2 *tw
# tc_CB       = 700 *mm
# bf_CB       = tc_CB +2 *t_pwCB
bf_CB       = 600 *mm
tc_CB       = bf_CB -2 *t_pwCB

Lw      = Hw
b_cCB   = bf_CB -2 *t_pwCB
t_cCB   = bf_CB -2 *t_pwCB
h_CB    = H_CB
h_cCB   = h_CB -2 *t_pfCB # Clear height of the web plate

b           = 114*mm
NfibeY      = 5
NIP         = 7

Section = {
    'wall': { # C-PSW/CF Wall Section
        #tags       = [tagSec, tagMatStFlange, tagMatStWeb, tagMatCtUnconf, tagMatCtConf]
        'tags'      : [1,      1,              2,           3,              4           ],
        #propStPart = [B,      H,         Es,      Fy,      Fu,      eps_sh, eps_ult, nu,   alpha, beta, gamma, Cf,  a1,  limit] 
        'propWeb'   : [tw,     (Hw-2*tf), Es,      Fy,      Fu,      0.007,  0.2,    0.28, 0.65,  1.0,  1.0,   0.5, 4.3, 0.01],
        'propFlange': [bf,     tf,        Es,      Fy,      Fu,      0.007,  0.2,    0.28, 0.65,  1.0,  1.0,   0.5, 4.3, 0.01],
        #propCore   = [tc,     fpc,       wc,      lamConf, lamUnconf]
        'propCore'  : [tc,     fpc,       0.2*mm,  0.05,     0.25    ]
    },
    'beam': { # Composite Beam Section
        #tags       = [tagSec, tagMatStFlange, tagMatStWeb, tagMatCtUnconf, tagMatCtConf]
        'tags'      : [2,      5,              6,           7,              8           ],
        #propStPart = [B,      H,         Es,      Fy,      Fu,      eps_sh, eps_ult, nu,   alpha, beta, gamma, Cf,  a1,  limit] 
        'propWeb'   : [t_pwCB, h_cCB,     Es,      Fy,      Fu,      0.007,  0.12,    0.28, 0.65,  1.0,  1.0,   0.5, 4.3, 0.01],
        'propFlange': [bf_CB,  t_pfCB,    Es,      Fy,      Fu,      0.007,  0.12,    0.28, 0.65,  1.0,  1.0,   0.5, 4.3, 0.01],
        #propCore   = [tc,     fpc,       wc,      lamConf, lamUnconf]
        'propCore'  : [tc_CB,  fpc,       0.2*mm,  0.05,     0.25     ]
    },
    }

#=============================================================================
#    Frame Data:
#=============================================================================
n_story         = 8
H_typical       = 4.25      *m#14        *ft
H_first         = 5.2       *m #17        *ft
LDR_CB          = L_CB /H_CB; print(f"LDR_CB = {LDR_CB:.3f}")
L_Bay           = Hw + L_CB #(Hw+2*tf) + L_CB
H_story_List    = [H_first, *((n_story-1)*[H_typical])]       # [Hstory1, *((numStories-1)*[HstoryTypical])]
L_Bay_List      = 2*[L_Bay]#, 5.*m, 5.*m, 5.*m]        # [*LBays]

# L               = L_CB /1
L               = H_first

# Building Geometry
Lf              = 200   *ft
Wf              = 120   *ft
h_typ           = H_typical
h_1             = H_first
A_SFRS          = Lf *Wf /2
#=============================================================================
#    Loading
#=============================================================================
#   Cantilever Loads
Pno             = 12776238.963599999 *N
# or:
ALR             = 0.02  # Axial Load Ratio
Py              = ALR * Pno
#   Frame Loads
load            = {}
DL_Floor        = 120   *psf #90 *psf
DL_PWalls       = 0 #25 *psf
LL_Floor        = 50    *psf
LL_Roof         = 0 #20 *psf

##  Tributary Loading
L_Bay_y         = (42.5 +35) *ft
L_Bay_x         = (42.5 +30) /2 *ft
A_Tributary     = 0.5*L_Bay_y * L_Bay_x
DL_Tributary    = A_Tributary * DL_Floor
LL_Tributary    = A_Tributary * LL_Floor
load["wall"]    = 1.0*DL_Tributary + 0.25*LL_Tributary # This is to calculate the effective weight of the building for earthquake force 
load["wallG"]   = 1.2*DL_Tributary + 1.60*LL_Tributary # This is for gravity analysis of the structure
# LoadG           = 72
# load["wall"]    = LoadG * kip

##  Loading the Leaning Columns
n_Bay_x         = len(L_Bay_List)
# A_SFRS          = (1.5 * L_Bay_y) * ((n_Bay_x+1) * L_Bay_x)
A_Leaning       = A_SFRS - A_Tributary*n_Bay_x
L_PWall         = L_Bay_y + ((n_Bay_x+1) * L_Bay_x) - n_Bay_x*Hw
DL_Leaning      = A_Leaning * DL_Floor + L_PWall*H_typical * DL_PWalls
LL_Leaning      = A_Leaning * LL_Floor
load["leaningColumn"] = 1.0*DL_Leaning + 0.25*LL_Leaning # This is to calculate the effective weight of the building for earthquake force 
load["leaningColumnG"]= 1.2*DL_Leaning + 1.60*LL_Leaning  # This is for gravity analysis of the structure
# load["leaningColumn"] = 0 * kip


if SDC == "Dmax":
    S_MS    = 1.5
    S_M1    = 0.9
else:
    S_MS    = 0.75
    S_M1    = 0.3
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                   Step 1: Input Data
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# Floor Loads
DL      = DL_Floor +DL_PWalls
LL      = LL_Floor +LL_Roof

We      = (1.0 *DL  +0.25 *LL) *A_SFRS *n_story



# Coupling Beam Properties
if typeSect == "Composite":
    As_CB   = 2 *t_pwCB *h_cCB + 2 *t_pfCB *bf_CB
    Ac_CB   = tc_CB *h_cCB
    Asw_CB  = 2 *h_CB *t_pwCB
    # Is_CB   = 2 *(1/12 *t_pwCB *h_cCB **3 + 1/12 *bf_CB *t_pfCB **3 +
    #               t_pfCB *bf_CB *(h_cCB /2 +t_pfCB /2) **2)
    # Ic_CB   = 1/12 *tc_CB *h_CB **3
    
    # EAuncrCB= 0.8 *(Es *As_CB +Ec *Ac_CB)
    # C3      = min(0.9, 0.45 +3 *(As_CB /bf_CB/ h_CB))
    # EIeff_CB= 0.64 *(Es *Is_CB +C3 *Ec *Ic_CB)
    # GAv_CB  = Gs *Asw_CB + Gc *Ac_CB
    
elif typeSect == "I_Shaped":
    As_CB   = t_pwCB *h_cCB +2 *t_pfCB *bf_CB
    Asw_CB  = 2 *h_CB *t_pwCB
    Z_CB    = (bf_CB * t_pfCB) *(h_cCB +t_pfCB) + (t_pwCB *h_cCB /2) *(h_cCB /2)
    Ix      = 1/12 *(bf_CB *H_CB **3) -1/12 *((bf_CB -t_pwCB) *h_cCB **3)
    Iy      = 1/12 *(2 *t_pfCB) *bf_CB **3 +1/12 *h_cCB *t_pwCB **3
    print(f"Iy/Ix = {Iy/Ix:.3f}")
    if Iy/Ix <= 0.67:
        print("Iy/Ix should be greater than 0.67!!!"); sys.exit()
    
# Planar SpeedCore Wall Properties
As      = 2 *(tw *(Lw -2 *tf) + tf *bf)
Ac      =     tc *(Lw -2 *tf)
Asw     = 2 *Lw *tw
# Is      = (2 *(1/12 *tw *(Lw -2 *tf) **3) + 
#            2 *(1/12 *bf *tf **3 + bf *tf *(Lw /2 +tf /2) **2))
# Ic      = 1/12 * tc * (Lw -2 *tw) **3

# EAeff   = Es *As  +0.45 *Ec *Ac
# EIeff   = Es *Is  +0.35 *Ec *Ic
# GAveff  = Gs *Asw +      Gc *Ac

#______________________________________________________________________________
#                   Check with Thresholds
"""'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''"""

# 1.    Check Minimum and Maximum Area of Steel
"""'''''''''''''''''''''''''''''''''''''''''"""
rho_min  = 0.01
rho_max  = 0.1

# 1.1.  Coupling Beams
if typeSect == "Composite":
    As_CBmin = rho_min *bf_CB *h_CB
    As_CBmax = rho_max *bf_CB *h_CB
    if As_CB < As_CBmin: 
        print(f"Area of Steel is {abs(As_CBmin-As_CB)/As_CBmin*100:.1f}% less than minimum for Coupling Beams!!!")
    elif As_CB > As_CBmax:
        print(f"Area of Steel is {abs(As_CBmax-As_CB)/As_CBmax*100:.1f}% greater than maximum for Coupling Beams!!!")

# 1.2.  SpeedCore Walls
Asmin = rho_min *t_sc *Lw
Asmax = rho_max *t_sc *Lw
if As < Asmin: 
    print(f"Area of Steel is {abs(Asmin-As)/Asmin*100:.1f}% less than minimum for Composite Walls!!!")
elif As > Asmax:
    print(f"Area of Steel is {abs(Asmax-As)/Asmax*100:.1f}% greater than maximum for Composite Walls!!!")


# 2. Check Plate Slenderness Ratios
"""'''''''''''''''''''''''''''''"""
# 1.1.  Coupling Beams
if typeSect == "Composite":
    lambdaP_fCB     = 2.37 *(Es /Ry /Fy) **0.5
    lambdaP_wCB     = 2.66 *(Es /Ry /Fy) **0.5
    
    lambda_fCB      = b_cCB /t_pfCB
    lambda_wCB      = h_cCB /t_pwCB
        
elif typeSect == "I_Shaped":
    # if      ductility == "high":
    #     lambdaP_fCB     = 0.38 *(Es /Fy) **0.5
    #     lambdaP_wCB     = 3.76 *(Es /Fy) **0.5
    # elif    ductility == "moderate":
    #     lambdaP_fCB     = 1.00 *(Es /Fy) **0.5
    #     lambdaP_wCB     = 5.70 *(Es /Fy) **0.5
        
    lambdaP_fCB     = 1.00 *(Es /Fy) **0.5
    lambdaP_wCB     = 3.76 *(Es /Fy) **0.5
        
    lambda_fCB      = ((bf_CB -t_pwCB) /2) /t_pfCB
    lambda_wCB      = h_cCB /t_pwCB

R_lambda_fCB    = lambda_fCB/lambdaP_fCB
R_lambda_wCB    = lambda_wCB/lambdaP_wCB
if R_lambda_fCB > 1.0: 
    print(f"Coupling Beam Flange Plate Slenderness OVERRATED!!!\n===>>>Lam_fCB/LamP_fCB = {R_lambda_fCB:.2f}")
if R_lambda_wCB > 1.0: 
    print(   f"Coupling Beam Web Plate Slenderness OVERRATED!!!\n===>>>Lam_wCB/LamP_wCB = {R_lambda_wCB:.2f}")

# 1.2.  SpeedCore Walls
lambda_btie     = btie/tw
lambdaP_btie    = 1.05 *(Es /Ry /Fy) **0.5
R_lambda_btie   = lambda_btie/lambdaP_btie

lambdaP_btieTop = 1.2 *(Es /Fy) **0.5
R_lambda_btieTop= lambda_btie/lambdaP_btieTop

if R_lambda_btie > 1.0:    
    print(f"Wall Faceplate Slenderness in flexural yielding zone not CHECKED!!!\n===>>>Lam_fW/LamP_fW = {R_lambda_btie:.2f}") 
if R_lambda_btieTop > 1.0: 
    print(f"Wall Faceplate Slenderness above flexural yielding zone not CHECKED!!!\n===>>>Lam_fW/LamP_fW = {R_lambda_btieTop:.2f}") 

Alpha           = 1.7 *(t_sc /tw -2) *(tw /dtie) **4
lambdaP_tbs     = (Es /(2 *Alpha +1)) **0.5
lambda_tbs      = Stie /tw
R_lambda_tbs    = lambda_tbs/lambdaP_tbs
if R_lambda_tbs > 1.0: 
    print(f"Wall Tie bar spacing not CHECKED!!!\n===>>>Lam_tbs/LamP_tbs = {R_lambda_tbs:.2f}") 


# 3. Check Shear Criticality of Coupling Beams
"""''''''''''''''''''''''''''''''''''''''''"""
def calc_C(fpc, Fy, tc, tw, tf, h, P=0, Ry=1, Rc=1):
    if   P < 0:
        print("Section is under compression.")
    elif P > 0:
        print("Section is under tension.")    
    C = ((0.85 *Rc  *fpc *tc *tf +2 *tw *Ry *Fy *h -P)/
         (0.85 *Rc  *fpc *tc     +4 *tw *Ry *Fy))
    return C
def calc_Mp(fpc, Fy, tc, tw, tf, h, bf, P=0, Ry=1, Rc=1):
    y       = calc_C(fpc, Fy, tc, tw, tf, h, P, Ry, Rc)
    C1 = T1 = bf *tf             *Ry *Fy
    C2      = 2 *tw *(y -tf)     *Ry *Fy
    C3      = tc    *(y -tf)     *Rc *(0.85 *fpc)
    T2      = 2 *tw *(h -y - tf) *Ry *Fy

    Mp = ( C1 *( y -tf     /2) 
          +C2 *((y -tf)    /2) 
          +C3 *((y -tf)    /2) 
          +T1 *( h -y -tf  /2) 
          +T2 *((h -y -tf) /2) 
          -P *(h /2 -y))
    return Mp
    
if typeSect == "Composite":
    C_CBexp = calc_C(fpc, Fy, tc_CB, t_pwCB, t_pfCB, h_CB, Ry=1.1, Rc=1.3)
    print(f"C_CBexp = {C_CBexp *1000:.1f} mm")
    M_exp   = calc_Mp(fpc, Fy, tc_CB, t_pwCB, t_pfCB, h_CB, bf_CB, Ry=1.1, Rc=1.3)
    V_exp   = 0.6 *Ry *Fy *Asw_CB +0.06 *Ac_CB *(Rc *fpc /MPa) **0.5 #!!! Take care NOT to put fpc in other units except MPa
    
elif typeSect == "I_Shaped":
    Ry      = 1.1
    M_exp   =      (1.25 *Ry) *(Fy *Z_CB)
    V_exp   = 0.6 *(1.25 *Ry) *(Fy *Asw_CB)

print(f"M_exp = {M_exp /1000:.1f} kN.m")
print(f"V_exp = {V_exp /1000:.1f} kN")

# Check Flexure-Criticality Condition
if      V_exp *L_CB /M_exp >= 2.4: 
    print(f"The Beams are Flexure-Critical, since V_exp *L_CB /M_exp is = {V_exp *L_CB /M_exp:.2f}>2.4")
    if typeSect == "I_Shaped": print("Use a Composite Shape for Flexure-Critical Coupling Beam then try again!!!"); sys.exit(); 
    shearCriticality=False
elif    V_exp *L_CB /M_exp <= 1.6: 
    print(f"The Beams are Shear-Critical, since V_exp *L_CB /M_exp is = {V_exp *L_CB /M_exp:.2f}<1.6")
    if typeSect == "Composite": print("Use a Composite Shape for Flexure-Critical Coupling Beam then try again!!!"); sys.exit();
    shearCriticality=True
else:
    print(f"The Beams are Flexure-Shear-Critical, since 1.6 < V_exp *L_CB /M_exp is = {V_exp *L_CB /M_exp:.2f} <2.4")
    print("Not acceptable! It should be either Shear-Critical or Flexure-Critical. \nTry Again!!!")
    sys.exit()
    shearCriticality=False



