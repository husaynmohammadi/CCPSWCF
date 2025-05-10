exec(open("MAIN.py").readlines()[18]) # It SHOULD read and execute exec(open("Input/units    .py").read())
# exec(open("../Input/unitsSI.py").read()) # It SHOULD read and execute exec(open("Input/units    .py").read()) 
import sys
import numpy                       as np
import openseespy.opensees         as ops
import matplotlib.pyplot           as plt
import opsvis                      as opv
nho_Ct=0.15; nho_St=0.3
class steel:
    def __init__(self, Es, Fy, Fu, eps_sh, eps_ult, nu, 
                 alpha, beta, gamma, Cf, a1, limit):
        self.Es     = Es
        self.Gs     = Es /2 /(1 + nho_St)
        self.Esh    = Es/30
        self.Fy     = Fy 
        self.Fu     = Fu 
        self.epsy   = Fy/Es
        self.eps_sh = eps_sh
        self.eps_ult= eps_ult
        self.nu     = nu
        self.alpha  = alpha
        self.beta   = beta
        self.gamma  = gamma
        self.Cf     = Cf
        self.a1     = a1
        self.limit  = limit
        self.A      = 1

class concrete:
    def __init__(self, fpc,  wc, lam):
        self.fpc        = fpc
        self.fpcu       = 1
        self.wc         = wc
        self.lam        = lam
        self.Ec         = (21.5e3*MPa * 1.0 * (abs(fpc/MPa)/10)**(1/3))            # Tangent Modulus of Elasticity with fpc in MPa  ==> CEB-FIB-2010 5.1.7.2 (Selected by Masoumeh Asgharpoor)
        self.Gc         = self.Ec /2 /(1 +nho_Ct)
        self.epsc0      = 2*fpc/self.Ec                                            # Crack width (Average of 0.15 to  0.25 mm)
        self.epscU      = 1
        self.Gf         = (73 * (fpc/MPa)**0.18)                                   # Fracture Energy CEB-FIB-2010 section 5.1.5.2
        self.fts        = 1.0*(0.3 * (fpc/MPa - 8)**(2/3))*MPa                     # CEB-FIB-2010 Eq. (5.1-5)
        self.Ets        = abs(1/(1-(2*self.Ec*self.Gf)/(wc*self.fts**2)))*self.Ec
        self.A          = 1
        
NfibeZ = 1
class compo:
    def __init__(self, typeEle, tagSec, tagMatStFlange, tagMatStWeb, tagMatCtUnconf, tagMatCtConf, P, lsr, b, NfibeY,
                 tw, Hw, Esw, Fyw, Fuw, eps_shw, eps_ultw, nuw, alphaw, betaw, gammaw, Cfw, a1w, limitw, 
                 Bf, tf, Esf, Fyf, Fuf, eps_shf, eps_ultf, nuf, alphaf, betaf, gammaf, Cff, a1f, limitf,
                 tc, fpc, wc, lamConf, lamUnconf, linearity
                 ):
        
        self.tagSec             = tagSec
        self.tagMatStFlange     = tagMatStFlange 
        self.tagMatStWeb        = tagMatStWeb    
        self.tagMatCtUnconf     = tagMatCtUnconf 
        self.tagMatCtConf       = tagMatCtConf   
        
        self.St_flange  = steel(Esw, Fyw, Fuw, eps_shw, eps_ultw, nuw, alphaw, betaw, gammaw, Cfw, a1w, limitw)
        self.St_web     = steel(Esf, Fyf, Fuf, eps_shf, eps_ultf, nuf, alphaf, betaf, gammaf, Cff, a1f, limitf)
        self.Ct_unconf  = concrete(fpc,  wc, lamUnconf)
        self.Ct_conf    = concrete(fpc,  wc, lamConf) # This is only to create Ct_conf
        
        self.d          = Hw + 2*tf
        self.Hc2        = Hc2 = (tc + 2*tw)/2           # Height of   confined concrete i.e. Region 2 in documentation # Masoumeh Asgharpoor
        self.Hc1        = Hc1 = Hw - 2*Hc2   # Height of unconfined concrete i.e. Region 1 in documentation
        if Hc1 < 0:
            print(f"Hw is {Hw}, but it cannot be less than {tc + 2*tw}"); #sys.exit()
        self.St_web.A   = 2*(tw*Hw)
        self.St_flange.A= 2*(tf*Bf)
        self.Ct_unconf.A= self.Hc1 * tc
        self.Ct_conf.A  = self.Hc2 * tc * 2
        self.St_A       = self.St_flange.A + self.St_web.A
        self.St_Asw     = self.d * 2*tw
        self.Ct_A       = self.Ct_unconf.A + self.Ct_conf.A
        self.Ag         = self.St_A + self.Ct_A
        self.P          = P
        self.Pno        = 0.85*(self.Ct_A*abs(fpc)) + (self.St_web.A*abs(Fyw) + self.St_flange.A*abs(Fyf))
        self.ALR        = ALR = P/self.Pno
        self.lsr        = lsr
        self.b          = b
        self.NfibeY     = NfibeY
        self.tw         = tw
        self.Hw         = Hw
        self.Bf         = Bf
        self.tf         = tf
        self.tc         = tc
        self.Es         = Es        = (self.St_flange.Es        + self.St_web.Es        )/2
        self.Gs         = Gs        = (self.St_flange.Gs        + self.St_web.Gs        )/2
        self.Esh        = Esh       = (self.St_flange.Esh       + self.St_web.Esh       )/2
        self.Fy         = Fy        = (self.St_flange.Fy        + self.St_web.Fy        )/2
        self.Fu         = Fu        = (self.St_flange.Fu        + self.St_web.Fu        )/2
        self.eps_sh     = eps_sh    = (self.St_flange.eps_sh    + self.St_web.eps_sh    )/2
        self.eps_ult    = eps_ult   = (self.St_flange.eps_ult   + self.St_web.eps_ult   )/2
        self.nu         = nu        = (self.St_flange.nu        + self.St_web.nu        )/2
        self.alpha      = alpha     = (self.St_flange.alpha     + self.St_web.alpha     )/2
        self.beta       = beta      = (self.St_flange.beta      + self.St_web.beta      )/2
        self.gamma      = gamma     = (self.St_flange.gamma     + self.St_web.gamma     )/2
        self.Cf         = Cf        = (self.St_flange.Cf        + self.St_web.Cf        )/2
        self.a1         = a1        = (self.St_flange.a1        + self.St_web.a1        )/2
        self.limit      = limit     = (self.St_flange.limit     + self.St_web.limit     )/2
        
        epsc0, fpcu         =           self.Ct_unconf.epsc0, self.Ct_unconf.fpcu
        self.rhos           = rhos    = (2*tw)/(tc + 2*tw) # Percentage of Steel
        self.R              = R       = b/tw * (12*(1-nu**2)/(4*np.pi**2))**0.5 * (Fy/Es)**0.5 # Masoumeh Asgarpoor
        self.frp            = frp     = (-6.5*R*((fpc/MPa)**1.46/(Fy/MPa)) + 0.12*(fpc/MPa)**1.03)*MPa # Masoumeh Asgarpoor
        self.fpcc           = fpcc    = ((fpc/MPa) + 4*(frp/MPa)*(1+0.8*ALR)**3.5)*MPa # Masoumeh Asgarpoor
        self.Ct_conf                  = concrete(fpcc, wc, lamConf)
        self.Ct_conf.Ec     = Ecc     = (21.5e3*MPa * 1.0 * (abs(fpcc/MPa)/10)**(1/3)) # Tangent Modulus of Elasticity with fpc in MPa  ==> CEB-FIB-2010 5.1.7.2 (Selected by Masoumeh Asgharpoor)
        self.Ct_conf.epsc0  = epscc0  = 2*fpcc/Ecc
        self.Ct_unconf.fpcu = fpcu    = (0.15 *(fpc/MPa)** 1.5 * (Fy/MPa)**0.01 * rhos**0.16 * (1-ALR)**(-0.5) * lsr**(-0.025))*MPa # Masoumeh Asgarpoor 
        self.Ct_unconf.epscU= epscU   = 0.157   *epsc0 **(-0.67) * (fpcu/MPa) **0.23 *(Fy/MPa)**(-1.4) * rhos**0.06 * (1-ALR)**(-0.17) * lsr**0.12 #Masoumeh Asgarpoor
        self.Ct_conf.fpcu   = fpccu   = (2.63 *(fpcc/MPa)**0.65 * (Fy/MPa)**0.001 * rhos**(-0.04) * (1-ALR)**(0.4) * lsr**(0.070))*MPa # Masoumeh Asgarpoor 
        self.Ct_conf.epscU  = epsccU  = 1.12e-6 *epscc0**(-0.33) * (fpccu/MPa)**0.83 *(Fy/MPa)**(0.77) * rhos**0.07 * (1-ALR)**( 0.15) * lsr**0.16 #Masoumeh Asgarpoor 
        self.r              = r       = (15*lsr**0.37 * (Fy/fpc)**0.26 * eps_sh**0.3 * eps_ult**0.85)**-1 # Masoumeh Asgarpoor 
        # self.r              = r       = 0.9
        self.Cd             = Cd      = (16 * (Fy/Fu)**0.78 * eps_ult * ((fpcu/MPa)**0.5/(fpc/MPa))**0.58 * r**0.37) # Masoumeh Asgarpoor 
        # self.Cd             = Cd      = 1
        
        
        # Section Strength
        self.yp= yp = (2*Fy*tw*(Hw+2*tf) + 0.85*fpc*tc*tf)/(4*Fy*tw + 0.85*fpc*tc)
        C_1=T_1 = Fy*Bf*tf; C_2 = 2*Fy*tw*(yp-tf); C_3 = 0.85*fpc*tc*(yp-tf)
        T_2 = 2*Fy*tw*(Hw+tf-yp)
        self.Mp = Mp = C_1*(yp-tf/2) + (C_2+C_3)*(yp-tf)/2 + T_1*(Hw+3/2*tf - yp) + T_2*(Hw+tf-yp)/2
        self.Vp = Vp = 0.6*Fy*(2*tw*(Hw+2*tf))
        self.eMax = 2.6*Mp/Vp; self.eMin = 1.6*Mp/Vp; self.eAve = 2.0*Mp/Vp
        
        # Effective Stiffnesses

        self.EIeff_webs     = self.St_web.Es    * (1/12 * (2*tw) * Hw**3)
        self.EIeff_flanges  = self.St_flange.Es * (1/12 * Bf * (self.d**3 - Hw**3))
        self.EIeff_St       = self.EIeff_webs + self.EIeff_flanges
        
        self.EIeff_unconf   = self.Ct_unconf.Ec * (1/12 * tc * self.Hc1**3)
        self.EIeff_conf     = self.Ct_conf.  Ec * (1/12 * tc * Hw**3 - 1/12 * tc * self.Hc1**3)
        self.EIeff_Ct       = self.EIeff_unconf + self.EIeff_conf
        
        self.C3             = min(0.45+3*(self.St_A/self.Ag), 0.9)
        if typeEle == "wall":
            self.EIeff      = self.EIeff_St + 0.35*self.EIeff_Ct
        else:
            self.EIeff      = 0.64*(self.EIeff_St + self.C3*self.EIeff_Ct)
            
        
        self.EAeff_webs     = self.St_web.Es    * self.St_web.A
        self.EAeff_flanges  = self.St_flange.Es * self.St_flange.A
        self.EAeff_St       = self.EAeff_webs   + self.EAeff_flanges
        
        self.EAeff_unconf   = self.Ct_unconf.Ec * self.Ct_unconf.A
        self.EAeff_conf     = self.Ct_conf.Ec   * self.Ct_conf.A
        self.EAeff_Ct       = self.EAeff_unconf + self.EAeff_conf
        
        self.EAeff          = self.EAeff_St + 0.45*self.EAeff_Ct
        
        self.GAveff         = Gs *self.St_Asw +((self.EAeff_Ct/self.Ct_A) /(2*(1 +nho_Ct))) *self.Ct_A
        
        # Define Materials
        if linearity == False:
            self.tagMatStFlange = tagMatStFlange
            ops.uniaxialMaterial('ReinforcingSteel', tagMatStFlange, Fyf, Fuf, Esf, Esh, eps_shf, eps_ultf, 
                                 '-GABuck', lsr, beta, r, gamma, 
                                 '-CMFatigue', Cf, alpha, Cd, 
                                 '-IsoHard', a1, limit
                                 )
            self.tagMatStWeb    = tagMatStWeb
            ops.uniaxialMaterial('ReinforcingSteel', tagMatStWeb,    Fyw, Fuw, Esw, Esh, eps_shw, eps_ultw, 
                                 '-GABuck', lsr, beta, r, gamma, 
                                 '-CMFatigue', Cf, alpha, Cd, 
                                 '-IsoHard', a1, limit
                                 )
            ops.uniaxialMaterial('Concrete02', tagMatCtUnconf, fpc,  epsc0,  fpcu,  
                                 self.Ct_unconf.epscU, lamUnconf, self.Ct_unconf.fts, self.Ct_unconf.Ets)
            ops.uniaxialMaterial('Concrete02', tagMatCtConf,   fpcc, epscc0, fpccu, 
                                 self.Ct_conf.epscU,   lamConf,   self.Ct_conf.fts,   self.Ct_conf.Ets)
        else: 
            self.tagMatStFlange = tagMatStFlange
            ops.uniaxialMaterial('Elastic', tagMatStFlange, Esf)
            self.tagMatStWeb    = tagMatStWeb
            ops.uniaxialMaterial('Elastic', tagMatStWeb,    Esw)
            ops.uniaxialMaterial('Elastic', tagMatCtUnconf, self.Ct_unconf.Ec)
            ops.uniaxialMaterial('Elastic', tagMatCtConf,   self.Ct_conf.Ec)
        
        
                    
        
    def defineSection(self, plot_section=False):
        
        # Section Geometry
        ##  Steel Parts
        ##1    Bottom   Flange
        crdsJ1 = [- self.Hw/2           ,  self.Bf/2            ]
        crdsI1 = [-(self.Hw/2 + self.tf), -self.Bf/2            ]
        ##2    Top      Flange
        crdsJ2 = [ (self.Hw/2 + self.tf),  self.Bf/2            ]
        crdsI2 = [  self.Hw/2           , -self.Bf/2            ]
        if self.Hw > 2 *self.Hc2:
            ##3    Bottom   BZ Web
            crdsJ3 = [-self.Hc1/2           ,  self.tw              ]
            crdsI3 = [-self.Hw/2            , -self.tw              ]
            ##4    Top      BZ Web
            crdsJ4 = [ self.Hw/2            ,  self.tw              ]
            crdsI4 = [ self.Hc1/2           , -self.tw              ]
            ##5    Mid Web
            crdsJ5 = [ self.Hc1/2           ,  self.tw              ]
            crdsI5 = [-self.Hc1/2           , -self.tw              ]
            ##  Concrete Parts
            ##6    Bottom   BZ Concrete (Confined)
            crdsJ6 = [-self.Hc1/2           ,  self.tc/2            ]
            crdsI6 = [-self.Hw/2            , -self.tc/2            ]
            ##7    Top      BZ Concrete (Confined)
            crdsJ7 = [ self.Hw/2            ,  self.tc/2            ]
            crdsI7 = [ self.Hc1/2           , -self.tc/2            ]
            ##8    Mid Concrete (Unconfined)
            crdsJ8 = [ self.Hc1/2           ,  self.tc/2            ]
            crdsI8 = [-self.Hc1/2           , -self.tc/2            ]
        else:
            ##5    Web
            crdsJ5 = [ self.Hw/2            ,  self.tw              ]
            crdsI5 = [-self.Hw/2            , -self.tw              ]
            ##8    Concrete (confined)
            crdsJ8 = [ self.Hw/2            ,  self.tc/2            ]
            crdsI8 = [-self.Hw/2            , -self.tc/2            ]
        
        # Define Fiber Section
        GJ          = 1e8
        N_CZone     = 2
        N_UZone     = 3
        ops.section('Fiber', self.tagSec, '-GJ', GJ)
        ops.patch('rect', self.tagMatStFlange,  self.NfibeY, 1, *crdsI1, *crdsJ1) #1 Bot Flange
        ops.patch('rect', self.tagMatStFlange,  self.NfibeY, 1, *crdsI2, *crdsJ2) #2 Top Flange
        if self.Hw > 2 *self.Hc2:
            ops.patch('rect', self.tagMatStWeb,     N_CZone, 1, *crdsI3, *crdsJ3) #3 Bot BZ Web
            ops.patch('rect', self.tagMatStWeb,     N_CZone, 1, *crdsI4, *crdsJ4) #4 Top BZ Web
            ops.patch('rect', self.tagMatStWeb,     N_UZone, 1, *crdsI5, *crdsJ5) #5 Mid Web
            ops.patch('rect', self.tagMatCtConf,    N_CZone, 1, *crdsI6, *crdsJ6) #6 Bot BZ Concrete (Confined)
            ops.patch('rect', self.tagMatCtConf,    N_CZone, 1, *crdsI7, *crdsJ7) #7 Bot BZ Concrete (Confined)
            ops.patch('rect', self.tagMatCtUnconf,  N_UZone, 1, *crdsI8, *crdsJ8) #8 Mid Concrete (Unconfined)
        else:
            ops.patch('rect', self.tagMatStWeb,     N_UZone, 1, *crdsI5, *crdsJ5) #5 Mid Web
            ops.patch('rect', self.tagMatCtConf,    N_UZone, 1, *crdsI8, *crdsJ8) #8 Mid Concrete (Confined)
        
        if plot_section == True: # True , False
            # this part of the code is just to sent out a varibale for plotting the fiber section
            if self.Hw > 2 *self.Hc2:
                fib_sec = [['section', 'Fiber', self.tagSec, '-GJ', GJ],
        
                         ['patch', 'rect', self.tagMatStFlange, self.NfibeY, 1, *crdsI1, *crdsJ1], #1 Bot Flange
                         ['patch', 'rect', self.tagMatStFlange, self.NfibeY, 1, *crdsI2, *crdsJ2], #2 Top Flange
                         ['patch', 'rect', self.tagMatCtConf,   N_CZone,     1, *crdsI6, *crdsJ6], #6 Bot BZ Concrete (Confined)
                         ['patch', 'rect', self.tagMatCtConf,   N_CZone,     1, *crdsI7, *crdsJ7], #7 Top BZ Concrete (Confined)
                         ['patch', 'rect', self.tagMatCtUnconf, N_UZone,     1, *crdsI8, *crdsJ8], #8 Mid Concrete (Unconfined)
                         ['patch', 'rect', self.tagMatStWeb,    N_CZone,     1, *crdsI3, *crdsJ3], #3 Bot BZ Web
                         ['patch', 'rect', self.tagMatStWeb,    N_CZone,     1, *crdsI4, *crdsJ4], #4 Top BZ Web
                         ['patch', 'rect', self.tagMatStWeb,    N_UZone,     1, *crdsI5, *crdsJ5], #5 Mid Web
                         ]
            else:
                fib_sec = [['section', 'Fiber', self.tagSec, '-GJ', GJ],
        
                         ['patch', 'rect', self.tagMatStFlange, self.NfibeY, 1, *crdsI1, *crdsJ1], #1 Bot Flange
                         ['patch', 'rect', self.tagMatStFlange, self.NfibeY, 1, *crdsI2, *crdsJ2], #2 Top Flange
                         ['patch', 'rect', self.tagMatStWeb,    N_UZone,     1, *crdsI5, *crdsJ5], #5 Mid Web
                         ['patch', 'rect', self.tagMatCtConf,   N_UZone,     1, *crdsI8, *crdsJ8], #8 Mid Concrete (Confined)
                         ]
        
            matcolor = ['y', 'b', 'r', 'g', 'y', 'b', 'r', 'g', 'm', 'k', 'y', 'b', 'r', 'g', 'm', 'k']
            opv.plot_fiber_section(fib_sec, matcolor=matcolor)
            plt.axis('equal')
            plt.show()
        
        
    def printVar(self):
        print(f"ALR\t\t= {self.ALR*100:.3f}%")
        print(f"Es\t\t= {self.Es/GPa:.0f} GPa")
        print(f"Esh\t\t= {self.Esh/GPa:.0f} GPa")
        print(f"Fy\t\t= {self.Fy/MPa:.2f} MPa")
        print(f"Fu\t\t= {self.Fu/MPa:.2f} MPa")
        print(f"eps_sh\t= {self.eps_sh:.5f}")
        print(f"eps_ult\t= {self.eps_ult:.5f}")
        print(f"fts\t\t= {self.Ct_unconf.fts/MPa:.2f} MPa")
        print(f"Ets\t\t= {self.Ct_unconf.Ets/MPa:.2f} MPa")
        print(f"rhos\t= {self.rhos:.5f}")
        print(f"R\t\t= {self.R:.5f}")
        print(f"frp\t\t= {self.frp/MPa:.2f} MPa")
        print(f"fpc\t\t= {self.Ct_unconf.fpc/MPa:.2f} MPa")
        print(f"epsc0\t= {self.Ct_unconf.epsc0:.5f}")
        print(f"Ec\t\t= {self.Ct_unconf.Ec/GPa:.0f} GPa")
        print(f"fpcu\t= {self.Ct_unconf.fpcu/MPa:.2f} MPa")
        if (self.Ct_unconf.fpcu/self.Ct_unconf.fpc < 0.45):
            print(f"Warning!!! fpcu/fpc = {self.Ct_unconf.fpcu/self.Ct_unconf.fpc:.2f} < 0.45")
        if (self.Ct_unconf.fpcu/self.Ct_unconf.fpc > 0.80):
            print(f"Warning!!! fpcu/fpc = {self.Ct_unconf.fpcu/self.Ct_unconf.fpc:.2f} > 0.80")
        print(f"epscU\t= {self.Ct_unconf.epscU:.5f}")
        if (self.Ct_unconf.epscU/self.Ct_unconf.epsc0 < 1.6):
            print(f"Warning!!! epscU/epsc0 = {self.Ct_unconf.epscU/self.Ct_unconf.epsc0:.2f} < 1.6")
        if (self.Ct_unconf.epscU/self.Ct_unconf.epsc0 > 5.5):
            print(f"Warning!!! epscU/epsc0 = {self.Ct_unconf.epscU/self.Ct_unconf.epsc0:.2f} > 5.5")
        print(f"fpcc\t= {self.Ct_conf.fpc/MPa:.2f} MPa")
        print(f"epscc0\t= {self.Ct_conf.epsc0:.5f}")
        print(f"Ecc\t\t= {self.Ct_conf.Ec/GPa:.0f} GPa")
        print(f"fpccu\t= {self.Ct_conf.fpcu/MPa:.2f} MPa")
        if (self.Ct_conf.fpcu/self.Ct_conf.fpc < 0.45):
            print(f"Warning!!! fpccu/fpcc = {self.Ct_conf.fpcu/self.Ct_conf.fpc:.2f} < 0.45")
        if (self.Ct_conf.fpcu/self.Ct_conf.fpc > 0.90):
            print(f"Warning!!! fpccu/fpcc = {self.Ct_conf.fpcu/self.Ct_conf.fpc:.2f} > 0.90")
        print(f"epsccU\t= {self.Ct_conf.epscU:.5f}")
        if (self.Ct_conf.epscU/self.Ct_conf.epsc0 < 4.0):
            print(f"Warning!!! epsccU/epscc0 = {self.Ct_conf.epscU/self.Ct_conf.epsc0:.2f} < 4.0")
        if (self.Ct_conf.epscU/self.Ct_conf.epsc0 > 8.7):
            print(f"Warning!!! epsccU/epscc0 = {self.Ct_conf.epscU/self.Ct_conf.epsc0:.2f} > 8.7")
        print(f"r\t\t= {self.r:.5f}")
        if (self.r < 0.08):
            print(f"Warning!!! r = {self.r:.2f} < 0.08")
        if (self.r > 0.7):
            print(f"Warning!!! r = {self.r:.2f} > 0.7")
        print(f"Cd\t\t= {self.Cd:.5f}")
        if (self.Cd < 0.2):
            print(f"Warning!!! Cd = {self.Cd:.2f} < 0.2")
        if (self.Cd > 0.75):
            print(f"Warning!!! Cd = {self.Cd:.2f} > 0.75")
        
        
        
# #tags       = [tagSec, tagMatStFlange, tagMatStWeb, tagMatCtUnconf, tagMatCtConf]
# tags        = [1,      1,              2,           3,              4           ]
# #propStPart = [B,         H,         Es,      Fy,      Fu,      eps_sh, eps_ult, nu,   alpha, beta, gamma, Cf,  a1,  limit] 
# propWeb     = [3/16*inch, 0.24 *m,   200*GPa, 422*MPa, 473*MPa, 0.007,  0.12,    0.28, 0.65,  1.0,  1.0,   0.5, 4.3, 0.01]
# propFlange  = [11*inch,   3/16*inch, 200*GPa, 422*MPa, 473*MPa, 0.007,  0.12,    0.28, 0.65,  1.0,  1.0,   0.5, 4.3, 0.01]
# #propCore   = [tc,     fpc,    wc,     lamConf, lamUnconf]
# propCore    = [9*inch, 50*MPa, 0.2*mm, 0.05,    0.25     ]
# lsr         = 48.
# b           = 114*mm
# NfibeY      = 1
# #beam= compo(*tags, P,      lsr, b, NfibeY, *propWeb, *propFlange, *propCore)
# beam = compo(*tags, 200*kN, lsr, b, NfibeY, *propWeb, *propFlange, *propCore)

# compo.printVar(beam)

# print(f"Es of flange material is {beam.St_flange.Es}")
# print(f"Es of web    material is {beam.St_web.Es}")
# print(f"Average Es the beam   is {beam.Es}")
# print(f"Depth of the beam is {beam.d}")
# print(f"EIeff = {beam.EIeff}")
# print(f"EAeff = {beam.EAeff}")

# ops.wipe()


class I_Shaped:
    def __init__(self, tagSec, tagMat, NfibeY, Fy, E, b1, R0, cR1, cR2, a1, a2, a3, a4, 
                 tw, Hw, bf, tf, linearity=False):
        class steel:
            def __init__(self, E, Fy, Asw, nho_St):
                self.E      = E
                self.Fy     = Fy
                self.G      = G = E /(2 *(1 +nho_St))
                self.GAveff = G *Asw
        
        self.tagMatStFlange     = tagMat
        self.tagMatStWeb        = tagMat
        
        self.tagSec     = tagSec
        self.tagMat     = tagMat
        self.NfibeY     = NfibeY
        self.Fy         = Fy
        self.E          = E
        self.b1         = b1
        self.R0         = R0
        self.cR1        = cR1
        self.cR2        = cR2
        self.a1         = a1
        self.a2         = a2
        self.a3         = a3
        self.a4         = a4
        
        self.tw         = tw
        self.Hw         = Hw
        self.bf         = bf
        self.tf         = tf
        self.d          = d     = Hw +2 *tf
        self.A          = A     = tw *Hw + 2 *(bf *tf)
        self.I          = I     = 1/12 *(bf *d **3) -1/12 *((bf -tw) *Hw **3)
        self.Z          = Z     = (bf * tf) *(Hw +tf) + (tw *Hw /2) *(Hw /2)
        self.St_Asw     = Asw = tw *Hw
        
        self.St_web     = steel(E, Fy, Asw, nho_St)
        self.GAveff     = self.St_web.GAveff
        self.EIeff      = E *I
        self.EAeff      = E *A
        self.Mp         = Mp    = Z *Fy
        self.Vp         = Vp    = Asw *(0.6*Fy)
        self.eMax       = 2.6*Mp/Vp
        self.eMin       = 1.6*Mp/Vp
        self.eAve       = 2.0*Mp/Vp
        ops.uniaxialMaterial('Steel02', tagMat, Fy, E, b1, *[R0,cR1,cR2], *[a1, a2, a3, a4])
    
    def defFibSection(self, plot_section=False):
        tw, Hw, bf, tf          = self.tw, self.Hw, self.bf, self.tf
        tagSec, tagMat, NfibeY  = self.tagSec, self.tagMat, self.NfibeY
        
        # Define Fiber Section
        GJ = 1e8
        # Section Geometry
        ##  Bottom Flange
        crdsI1 = [-(Hw /2 +tf), -bf /2]
        crdsJ1 = [- Hw /2     ,  bf /2]
        ##  Top Flange      
        crdsI3 = [  Hw /2     , -bf /2]
        crdsJ3 = [ (Hw /2 +tf),  bf /2]
        ##  Web
        crdsI2 = [-Hw /2      , -tw   ]
        crdsJ2 = [ Hw /2      ,  tw   ]
        
        NfibeZ = 1
        ops.section('Fiber', tagSec, '-GJ',   GJ)
        ops.patch  ('rect',  tagMat,  NfibeY, NfibeZ, *crdsI1, *crdsJ1) #Bot Flange
        ops.patch  ('rect',  tagMat,  NfibeY, NfibeZ, *crdsI3, *crdsJ3) #Top Flange
        ops.patch  ('rect',  tagMat,  NfibeY, NfibeZ, *crdsI2, *crdsJ2) #Web
        
        if plot_section == True: # True , False
            # this part of the code is just to sent out a varibale for plotting the fiber section
            fib_sec =[['section', 'Fiber', tagSec, '-GJ', GJ],
    
                      ['patch', 'rect', tagMat, NfibeY, NfibeZ, *crdsI1, *crdsJ1], #Bot Flange
                      ['patch', 'rect', tagMat, NfibeY, NfibeZ, *crdsI3, *crdsJ3], #Top Flange
                      ['patch', 'rect', tagMat, NfibeY, NfibeZ, *crdsI2, *crdsJ2], #Web
                      ]
        
            matcolor = ['y', 'b', 'r', 'g', 'y', 'b', 'r', 'g', 'm', 'k', 'y', 'b', 'r', 'g', 'm', 'k']
            opv.plot_fiber_section(fib_sec, matcolor=matcolor)
            plt.axis('equal')
            plt.show()




