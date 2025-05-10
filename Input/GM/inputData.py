"C-PSW/CF Section: Shafaei PP=136"
# import sys
# from functions.ClassComposite import compo
exec(open("MAIN.py").readlines()[18]) # It SHOULD read and execute exec(open(f"Input/units{'US'}.py").read())
# exec(open("../Input/unitsSI.py").read()) # It SHOULD read and execute exec(open("Input/units    .py").read())
# exec(open("MAIN.py").readlines()[20]) # It SHOULD read and execute exec(open("Input/materialParameters.py").read())

#=============================================================================
#    Elements
#=============================================================================
#       Element Length
Hw          = 980   *mm
H_CB        = 300   *mm
bf          = 200   *mm
RhoW        = 0.052
tw          = RhoW*bf/2
tf          = tw
tc          = bf - 2*tw
lsr         = 24.
b           = 114*mm
NfibeY      = 10

fpc         = 45
Fy          = 422
Fu          = 1.12*Fy

Section = {
    'wall': { # C-PSW/CF Wall Section
        #tags       = [tagSec, tagMatStFlange, tagMatStWeb, tagMatCtUnconf, tagMatCtConf]
        'tags'      : [1,      1,              2,           3,              4           ],
        #propStPart = [B,      H,         Es,      Fy,      Fu,      eps_sh, eps_ult, nu,   alpha, beta, gamma, Cf,  a1,  limit] 
        'propWeb'   : [tw,     Hw,        200*GPa, Fy *MPa, Fu *MPa, 0.007,  0.12,    0.28, 0.65,  1.0,  1.0,   0.5, 4.3, 0.01],
        'propFlange': [bf,     tf,        200*GPa, Fy *MPa, Fu *MPa, 0.007,  0.12,    0.28, 0.65,  1.0,  1.0,   0.5, 4.3, 0.01],
        #propCore   = [tc,     fpc,       wc,      lamConf, lamUnconf]
        'propCore'  : [tc,     fpc*MPa,   0.2*mm,  0.05,     0.25    ]
    },
    'beam': { # Composite Beam Section
        #tags       = [tagSec, tagMatStFlange, tagMatStWeb, tagMatCtUnconf, tagMatCtConf]
        'tags'      : [2,      5,              6,           7,              8           ],
        #propStPart = [B,      H,         Es,      Fy,      Fu,      eps_sh, eps_ult, nu,   alpha, beta, gamma, Cf,  a1,  limit] 
        'propWeb'   : [tw,     H_CB,      200*GPa, Fy *MPa, Fu *MPa, 0.007,  0.12,    0.28, 0.65,  1.0,  1.0,   0.5, 4.3, 0.01],
        'propFlange': [bf,     tf,        200*GPa, Fy *MPa, Fu *MPa, 0.007,  0.12,    0.28, 0.65,  1.0,  1.0,   0.5, 4.3, 0.01],
        #propCore   = [tc,     fpc,       wc,      lamConf, lamUnconf]
        'propCore'  : [tc,     fpc*MPa,   0.2*mm,  0.05,     0.25     ]
    },
    }

#=============================================================================
#    Frame Data:
#=============================================================================
n_story         = 2
H_typical       = 3.    *m
H_first         = 4.    *m
LDR_CB          = 2
L_CB            = LDR_CB * H_CB
L_Bay           = (Hw+2*tf) + L_CB
H_story_List    = [H_first, *((n_story-1)*[H_typical])]       # [Hstory1, *((numStories-1)*[HstoryTypical])]
L_Bay_List      = 2*[L_Bay]#, 5.*m, 5.*m, 5.*m]        # [*LBays]

L               = H_typical

#=============================================================================
#    Loading
#=============================================================================
#   Cantilever Loads
Pno             = 12776238.963599999 *N
# or:
ALR             = 0.02  # Axial Load Ratio
Py              = ALR * Pno
#   Frame Loads
load={}
DL_Floor        = 90 *psf
DL_PWalls       = 25 *psf
LL_Floor        = 50 *psf
LL_Roof         = 20 *psf

##  Tributary Loading
L_Bay_y         = 4 *m
L_Bay_x         = L_Bay
A_Tributary     = 0.5*L_Bay_y * L_Bay_x
DL_Tributary    = A_Tributary * DL_Floor
LL_Tributary    = A_Tributary * LL_Floor
load["wall"]    = 1.0*DL_Tributary + 0.25*LL_Tributary
# LoadG           = 0
# load["wall"]    = LoadG * kip

##  Loading the Leaning Columns
n_Bay_x         = len(L_Bay_List)
A_SFRS          = (1.5 * L_Bay_y) * ((n_Bay_x+1) * L_Bay_x)
A_Leaning       = A_SFRS - A_Tributary*n_Bay_x
L_PWall         = L_Bay_y + ((n_Bay_x+1) * L_Bay_x) - n_Bay_x*Hw
DL_Leaning      = A_Leaning * DL_Floor + L_PWall*H_typical * DL_PWalls
LL_Leaning      = A_Leaning * LL_Floor
load["leaningColumn"] = 1.0*DL_Leaning + 0.25*LL_Leaning
# load["leaningColumn"] = 0 * kip






















