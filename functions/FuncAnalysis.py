exec(open("MAIN.py").readlines()[18]) # It SHOULD read and execute exec(open("Input/units    .py").read())
exec(open("MAIN.py").readlines()[19]) # It SHOULD read and execute exec(open("Input/inputData.py").read())
import time, os, eqsig, winsound, sys
import openseespy.opensees     as ops
import numpy                   as np
import functions.FuncRecorders as fr
import functions.FuncPlot      as fp
import matplotlib.pyplot       as plt
# from colorama import Fore, Style # print(f"{Fore.YELLOW} your text {Style.RESET_ALL}")



def analyzeEigen(nEigen, printIt=True):
    omega2List  = sorted(ops.eigen(nEigen))
    Periods     = []
    for omega2 in omega2List:
        T       = 2 * np.pi/omega2**0.5
        Periods.append(T)
    if printIt == True:
        for index, Period in enumerate(Periods):
            print(f"Period{index:02} = {Period}")
    return Periods

def rayleighDamping(zeta):
    eigenList   = ops.eigen(2)
    omegaI2     = eigenList[0]
    omegaJ2     = eigenList[1]
    omegaI      = omegaI2 **0.5
    omegaJ      = omegaJ2 **0.5
    alphaM      = 2.0*zeta*(omegaI*omegaJ)/(omegaI+omegaJ)
    betaKinit   = 2.0*zeta/(omegaI+omegaJ)
    # rayleigh(alphaM, betaK, betaKinit, betaKcomm)
    ops.rayleigh(alphaM, 0.0, betaKinit, 0.0)

def Sa(T,S_MS=1.5, S_M1=0.9):
    Method  = 1
    S_DS    = 2/3 *S_MS
    S_D1    = 2/3 *S_M1
    Ts      = S_D1/S_DS
    T0      = 0.2 *Ts
    TL      = 8
    if Method == 1:
        if 0 <= T < T0:
            Sa  = (0.4 +0.6 *T /T0) *S_DS
        elif T0 <= T < Ts:
            Sa  = S_DS
        elif Ts <= T < TL:
            Sa  = S_D1 /T
        elif T >= TL:
            Sa  = S_D1 *TL /T **2
        else:
            print(f"T = {T} which shows there is something wrong!!!\nProgram exits here."); sys.exit()
    elif Method == 2:
        Sa      = min(S_DS, S_D1 /T)
    return Sa

def verDistFact(We, T, h_1, h_typ, n_story):
    if 0 <= T < 0.5:
        k = 1
    elif T >= 2.5:
        k = 2
    else:
        k = 0.5 *(T -0.5) +1
    def h(n):
        if n == 1:
            return h_1
        else:
            return (h_1 + (n-1) *h_typ)
    wx = We /n_story
    
    sumWH = 0
    for i in range(1, n_story+1):
        sumWH += wx *h(i) **k
    Cvx = [0]
    for x in range(1, n_story+1):
        factor = wx *h(x) **k /sumWH
        Cvx.append(factor)
    return Cvx

def gravity(load, tagNodeLoad):
    
    tagTSGravity    = 100000
    ops.timeSeries('Linear', tagTSGravity)
    # ops.timeSeries('Constant', tagTSGravity)
    
    tagPtnGravity   = tagTSGravity
    ops.pattern('Plain', tagPtnGravity, tagTSGravity)
    print(f"Type of tagNodeLoad is {type(tagNodeLoad)}")
    if type(tagNodeLoad) == int: 
        print("Loading is based on Cantilever Column Structure.")
        ops.load(tagNodeLoad, 0.0, load, 0.0)
    else:
        print("Loading is based on Shear Wall Structure.")
        for element, tagNodes in tagNodeLoad.items():
            if element == "wall":
                for tagNode in tagNodes:
                    ops.load(tagNode, 0.0, -abs(load["wallG"]), 0.0)
            elif element == "leaningColumn":
                for tagNode in tagNodes:
                    ops.load(tagNode, 0.0, -abs(load["leaningColumnG"]), 0.0)
            else:
                print("In GravityLoading element type was unknown!"); sys.exit()
    
    
    # Gravity Analysis:
    tol = 1.0e-5
    iteration = 100
    
    ops.wipeAnalysis()
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGen')
    ops.test('NormDispIncr', tol, iteration)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 0.1)
    # ops.integrator('LoadControl', 1)
    ops.analysis('Static')
    ops.analyze(10)
    # ops.analyze(1)
    ops.loadConst('-time', 0.0)

def convergeIt(typeAnalysis, tagNodeLoad, tagNodeBase, dofNodeControl, incrFrac, numFrac, disp, dispIndex, dispList, dispTarget, t_beg, Beams, numIncrInit=5):
    if type(tagNodeLoad) == list:
        tagNodeControl  = tagNodeLoad[-1]
    else:
        tagNodeControl  = tagNodeLoad
    
    def dispShear(Beams):
        rotation={}
        deltaY={}
        for tagEle, tagNodes in Beams.items():
            coordI          = ops.nodeCoord(tagNodes[0])
            coordJ          = ops.nodeCoord(tagNodes[1])
            L               = ((coordJ[0] -coordI[0]) **2 +(coordJ[1] -coordI[1]) **2) **0.5
            dispYI          = ops.nodeDisp(tagNodes[0], 2)
            dispYJ          = ops.nodeDisp(tagNodes[1], 2)
            deltaY[tagEle]  = abs(dispYJ - dispYI)
            rotation[tagEle]= abs(dispYJ - dispYI) /L
        return deltaY, rotation
    
    def curD():
        if typeAnalysis=="NTHA":
            
            t = ops.getTime()
            # print(f" = {}")
            drift       = []
            tagNodePrv  = tagNodeBase
            # print(f"tagNodePrv = {tagNodePrv}")
            if type(tagNodeLoad) == list:
                for i, tagNode in enumerate(tagNodeLoad):
                    if i >0:
                        height  = ops.nodeCoord(tagNode)[1] - ops.nodeCoord(tagNodePrv)[1]
                        # print(f"height = {height}")
                        dispTop = ops.nodeDisp(tagNode,    dofNodeControl)
                        # print(f"dispTop = {dispTop}")
                        dispBot = ops.nodeDisp(tagNodePrv, dofNodeControl)
                        # print(f"dispBot = {dispBot}")
                        drift.append(abs(dispTop - dispBot)/height)
                        # print(f"drift[{i}] = {drift[i-1] *100} %")
                    tagNodePrv = tagNode
                    # print(f"tagNodePrv = {tagNodePrv}")
                driftMax = max(drift)
                # print(f"driftMax = {driftMax *100} %")
            else:
                height  = ops.nodeCoord(tagNodeLoad)[1] - ops.nodeCoord(tagNodeBase)[1]
                dispTop = ops.nodeDisp(tagNodeLoad, dofNodeControl)
                dispBot = ops.nodeDisp(tagNodeBase, dofNodeControl)
                driftMax=driftAve= abs(dispTop - dispBot)/height
            return t, driftMax
        else:
            d = ops.nodeDisp(tagNodeControl, dofNodeControl)
            return d, -1
    
    def msgCurrentState(algorithm, tester):
        if typeAnalysis=="NTHA":
            remD    = dispTarget - curD()[0]
            print(f"\n{'~'*100}")
            print(f"{'-'*60}\nAlgorithm:\t{algorithm}")
            print(f"{'-'*60}\ntester:\t\t{tester}\n{'-'*60}")
            print(f"======>>> durGM\t\t\t\t= {dispTarget}")
            print(f"======>>> dT({iii:05}/{numFrac:05})\t= {dispTar}")
            print(f"======>>> Current   Time\t= {curD()[0]}")
            print(f"======>>> Remaining Time\t= {remD}")
            print(f"numIncr\t\t\t= {numIncr}")
            print(f"dt\t\t\t\t= {incr}")
        else:
            remD    = dispTarget - curD()[0]
            print(f"\n{'~'*100}")
            print(f"{'-'*60}\nAlgorithm:\t{algorithm}")
            print(f"{'-'*60}\ntester:\t\t{tester}\n{'-'*60}")
            print(f"======>>> disp({dispIndex+1:02}/{len(dispList):02})\t\t\t\t= {disp}")
            print(f"======>>> dispTarget\t\t\t\t= {dispTarget}")
            print(f"======>>> dispTar({iii:02}/{numFrac:02})\t\t\t= {dispTar}")
            print(f"======>>> Current   Displacement\t= {curD()[0]}")
            print(f"======>>> Remaining Displacement\t= {remD}")
            print(f"numIncr\t\t\t= {numIncr}")
            print(f"Incr\t\t\t= {incr}")
    
    def msgReducingIncrSize():
        if typeAnalysis=="NTHA":
            print(f"\n{'#'*65}\nAnalysis Failed!!\nReducing the dt size:\n{'#'*65}")
            print(f"\n>>>>>>>>> tolerance\t\t\t= {tol}")
            print(f"======>>> Current   Time\t= {curD()[0]}")
            print(f"======>>> Remaining Time\t= {remD}")
            print(f"numIncr\t\t\t= {numIncr}")
            print(f"dt\t\t\t\t= {incr}")
            
        else:
            print(f"\n{'#'*65}\nAnalysis Failed!!\nReducing the Incr size:\n{'#'*65}")
            print(f"\n>>>>>>>>> tolerance\t\t\t= {tol}")
            print(f"======>>> Current   Displacement\t= {curD()[0]}")
            print(f"======>>> Remaining Displacement\t= {remD}")
            print(f"numIncr\t\t\t= {numIncr}")
            print(f"Incr\t\t\t= {incr}")
    
    ##########################################################
    def core():
        algorithm   = algorithmList [algorithmIndex]
        tester      = testerList    [testerIndex]
        if tester == 'NormUnbalance':
            tol = tolForce
            ops.test(tester, tol, numIter, pFlag, nType)
        elif tester == 'NormDispIncr' or tester == 'EnergyIncr':
            tol = tolDisp
            ops.test(tester, tol, numIter, pFlag, nType)
        elif tester == 'NormDispAndUnbalance' or tester == 'NormDispOrUnbalance':
            tol = tolDisp
            ops.test(tester, tolDisp, tolForce, numIter, pFlag, nType)
        elif tester == 'RelativeNormUnbalance' or tester == 'RelativeNormDispIncr' or tester == 'RelativeEnergyIncr':
            tol = tolDisp
            ops.test(tester, 1.0000001, numIter, pFlag, nType)
        elif tester == 'RelativeTotalNormDispIncr':
            tol = tolDisp
            ops.test(tester, 0.0001, numIter, pFlag, nType)
        else:
            tol = tolDisp
            ops.test(tester, tol, numIter, pFlag, nType)
        # ops.test(tester, tol, numIter)
        
        secant      = False # False True
        if algorithm == 'Newton' or algorithm == 'ModifiedNewton' or algorithm == 'BFGS' or algorithm == 'Broyden':
            ops.algorithm(algorithm, secant)
        elif algorithm == 'NewtonLineSearch':
            # secant              = False # False True
            Bisection           = False
            RegulaFalsi         = False
            InitialInterpolated = False
            tol=0.1; maxIter=100; minEta=0.1; maxEta=10.0
            ops.algorithm('NewtonLineSearch', Bisection, secant, RegulaFalsi, InitialInterpolated, tol, maxIter, minEta, maxEta)  
        elif algorithm == 'KrylovNewton' or algorithm == 'SecantNewton' or algorithm == 'RaphsonNewton' or algorithm == 'PeriodicNewton' :
            # iterate     ='current' # current initial noTangent
            # increment   ='current' # current initial noTangent
            # maxDim      = 3
            ops.algorithm(algorithm) #, iterate, increment, maxDim)  
            
        if typeAnalysis == "NTHA":
            ops.integrator('Newmark', gamma, beta)
            ops.analysis('Transient')
        else:
            ops.integrator('DisplacementControl', tagNodeControl, dofNodeControl, incr)
            # numIter     = 1500
            # dUmin       = 1e-13
            # dUmax       = 1e-8
            # ops.integrator('ParallelDisplacementControl', tagNodeControl, dofNodeControl, incr, numIter, dUmin, dUmax)
            ops.analysis('Static')
        msgCurrentState(algorithm, tester)
        
        # Run Analysis
        if typeAnalysis == "NTHA": 
            OK = ops.analyze(numIncr, incr)
        else:
            OK = ops.analyze(numIncr)
        print(f"AnalyzeOutput\t= {OK}")
        t_now=time.time(); elapsed_time=t_now-t_beg; mins=int(elapsed_time/60); secs=int(elapsed_time%60)
        print(f"\nElapsed time: {mins} min + {secs} sec")
        if Beams != 0:
            deltaYBeam, rotationBeam = dispShear(Beams)
            removedElements = []
            for tagEle, tagNodes in Beams.items():
                print(f"deltaYBeam {tagEle} = {deltaYBeam[tagEle]*1000:.1f} mm")
                if rotationBeam[tagEle] > 0.08:
                    if tagEle not in removedElements:
                        ops.remove('ele', tagEle)
                        ops.element('elasticBeamColumn', tagEle, *[tagNodes[0], tagNodes[1]], 1, 1, 1, 1)
                        removedElements.append(tagEle)
                    print(f"==--> Rotation of Beam {tagEle} is {rotationBeam[tagEle]:.3f} which exceeded 0.08 rad, thus removed!")
        return OK, algorithm, tester, tol
    ##########################################################
    
    
    for iii in range(1, numFrac+1):
        dispTar         = iii * incrFrac
        testerList      = [
            # 'NormDispAndUnbalance',
            # 'RelativeNormUnbalance',
            # 'RelativeEnergyIncr',
            # 'RelativeNormDispIncr',
            # 'NormDispOrUnbalance',
            'EnergyIncr', 
            'NormUnbalance', 
            'NormDispIncr', 
            # 'RelativeTotalNormDispIncr',
            ]
        algorithmList   = [*(1*[
            'KrylovNewton', 
            'Newton', 
            'NewtonLineSearch', 
            # 'ModifiedNewton', 
            # 'RaphsonNewton', 
            # 'SecantNewton', 
            # 'PeriodicNewton', 
            # 'BFGS', 
            # 'Broyden', 
            ])]
        gamma = 0.5; beta = 0.25
        numIter     = 500
        numIncrMax  = 30000
        incrMin     = 1e-5
        
        tolForce    = 1e-8 *N
        tolDisp     = 1e-8 *m
        
        numIncr     = numIncrInit
        incr        = incrFrac/numIncrInit
        j = 1; jj = 1
        for i in range(100000000):
            
            ###########################################################################
            # testerCounter   = 0
            # algorithmCounter= 0
            pFlag           = 0 # 0, 1, 2, 4, 5
            nType           = 0 # Type of norm, (0 = max-norm, 1 = 1-norm, 2 = 2-norm)
            for tester_ in testerList:
                if "testerIndex" not in locals():
                    testerIndex = 0
            ###########################################################################
                
                ###########################################################################
                for algorithm_ in algorithmList:
                    if "algorithmIndex" not in locals():
                        algorithmIndex = 0
                ###########################################################################
                    
                    ###########################################################################
                    OK, algorithm, tester, tol = core()
                    ###########################################################################     
                
                ###########################################################################     
                    if OK == 0: break
                    elif OK != 0:
                        print(f"\n=============== THE ALGORITHM {algorithm} FAILED TO CONVERGE!!! ===============")
                        algorithmIndex  += 1
                        if algorithmIndex/len(algorithmList) >=1:
                            algorithmIndex   = algorithmIndex %len(algorithmList)
                            break
                            # testerIndex += 1
                            # if testerIndex/len(testerList) ==1:
                            #     testerIndex = 0
                            
                        # algorithmIndex   = algorithmIndex %len(algorithmList)
                        
                        # if algorithmIndex == len(algorithmList):
                        #     # algorithmCounter    =   0
                        #     testerIndex         +=  1
                        #     if testerIndex == len(testerList): 
                        #         testerIndex = 0
                        #         break
                            
                        # if testerCounter == len(testerList):
                        #     algorithmIndex += 1
                        #     testerCounter   = 0
                        #     if algorithmIndex == len(algorithmList)-1: 
                        #         algorithmIndex = 0
                        #         break
                ###########################################################################     

            ###########################################################################     
                if OK == 0: break
                elif OK != 0:
                    t_now=time.time(); elapsed_time=t_now-t_beg; mins=int(elapsed_time/60); secs=int(elapsed_time%60)
                    print(f"\nElapsed time: {mins} min + {secs} sec")
                    print(f"\n=============== The tester {tester} failed to converge!!! ===============")
                    
                    testerIndex += 1
                    if testerIndex/len(testerList) >=1:
                        testerIndex   = testerIndex %len(testerList)
                        break
                    
                        
                    # algorithmCounter  += 1
                    # testerIndex       += 1
                    # if algorithmIndex == len(algorithmList): 
                    #     algorithmIndex = 0
                    #     break
                    
                    # testerIndex    += 1
                    # testerCounter  += 1
                    # if testerIndex == len(testerList): 
                    #     testerIndex = 0
                    #     break
            ###########################################################################     

            if OK == 0: break
            else:
                # if tester == 'NormUnbalance':
                #     # tol = min(1.5*tol, 5 *N)
                #     tolForce    = min(2 *tolForce, 10 *kN)
                # else:
                #     # tol = min(1.5*tol, 1 *mm)
                #     tolDisp     = min(2 *tolDisp,  5 *mm)
                tolForce    = min(10 *tolForce, 10 *kN) # To decrease precision replace 1 with a higher number
                tolDisp     = min(10 *tolDisp,  5 *mm)
                remD    = dispTar - curD()[0]; print(f"{remD = }")
                if abs(remD) >= 0.00001:
                    numIncr = int(numIncr*1.001**i + j)
                    j       += 1; print(f"{j = }")
                    incr    = remD/numIncr
                else:
                    OK      = 1
                    break
                    numIncr = 2
                    incr    = remD/numIncr
                    jj      += 1; print(f"{jj = }")
                msgReducingIncrSize()
                if numIncr >= numIncrMax or abs(incr) <= incrMin or jj>200:
                    print("\nIncrement size is too small!!!")
                    t_now=time.time(); elapsed_time=t_now-t_beg; mins=int(elapsed_time/60); secs=int(elapsed_time%60)
                    print(f"\nElapsed time: {mins} min + {secs} sec")
                    # winsound.Beep(440, 1000)  # generate a 440Hz sound that lasts 500 milliseconds
                    text = " pushover" if typeAnalysis!="NTHA" else ""
                    print("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                    print(f"*!*!*!*!*!* The {typeAnalysis}{text} analysis failed to converge!!! *!*!*!*!*!*")
                    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                    break
        if OK < 0: break
        if typeAnalysis == "NTHA":
            driftMaxAllowed = 0.2
            if curD()[1] >= driftMaxAllowed: 
                print(f"driftMax = {curD()[1] *100} % >= {driftMaxAllowed *100} % ==> the next scaleFactor will be applied now!")
                break
    return OK

def pushoverDCF(dispTarget, incrInit, numIncrInit, tagNodeLoad, tagNodeLoad2, Beams, distributeOnWalls=True): 
    
    t_beg           = time.time()
    T1              = analyzeEigen(1)[0]
    dofNodeControl  = 1
    tagTSLinear     = 1
    ops.timeSeries('Linear',   tagTSLinear)
    tagPatternPlain = 1
    ops.pattern('Plain', tagPatternPlain, tagTSLinear)
    #   load(nodeTag,     *loadValues)
    if type(tagNodeLoad) == list:
        if distributeOnWalls == False:
            tagNodeControl  = tagNodeLoad[-1]
            n_story         = len(tagNodeLoad)-1
            Cvx             = verDistFact(We, T1, h_1, h_typ, n_story)
            for i, tagNode in enumerate(tagNodeLoad):
                # ops.load(tagNode, *[i/n_story, 0, 0])
                ops.load(tagNode, *[Cvx[i], 0, 0])
        else:
            def tagNodeLoadStory(tagNodeLoad, n_story):
                tagNodeStory        = {}
                for i in range(1, n_story+1):
                    tagNodeStory[i] = []
                    for tagNode in tagNodeLoad:
                        tagCoordYI  = f"{tagNode}"[1:-3]
                        if tagCoordYI == f"{i:02}":
                            tagNodeStory[i].append(tagNode)
                return tagNodeStory
            
            tagNodeControl      = tagNodeLoad[-1]
            nWalls              = 2
            n_story             = int(len(tagNodeLoad2)/nWalls)
            tagNodeLoadStories  = tagNodeLoadStory(tagNodeLoad2, n_story)
            Cvx                 = verDistFact(We, T1, h_1, h_typ, n_story)
            for story, tagNodeList in tagNodeLoadStories.items():
                for tagNode in tagNodeList:
                    ops.load(tagNode, *[Cvx[story]/nWalls, 0, 0])
    else:
        tagNodeControl  = tagNodeLoad
        ops.load(tagNodeControl, *[1, 0, 0])
    
    #  Define Analysis Options
    ops.wipeAnalysis()
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('FullGeneral')   # 'FullGeneral', 'UmfPack', 'SparseSYM', 
    disp        = dispTarget; dispIndex   = 0
    delta       = dispTarget
    numFrac     = int(delta/incrInit)
    incrFrac    = delta/numFrac
    asTagNodeBase = 1 #it is not going to be used at all in this analysis. it is just to fill a positional argument
    OK          = convergeIt("Monotonic", tagNodeLoad, asTagNodeBase, dofNodeControl, incrFrac, numFrac, disp, dispIndex, ['This is a list'], dispTarget, t_beg, Beams, numIncrInit)
    return OK

def calcDrift(tagNodeLoad, tagNodeBase, dofNodeControl):
    driftList   = []
    if type(tagNodeLoad) == list:
        tagNodePrv  = tagNodeBase[0]
    if type(tagNodeLoad) == list:
        for i, tagNode in enumerate(tagNodeLoad):
            if i>0:
                height  = ops.nodeCoord(tagNode)[1] - ops.nodeCoord(tagNodePrv)[1]
                dispTop = ops.nodeDisp(tagNode,    dofNodeControl)
                dispBot = ops.nodeDisp(tagNodePrv, dofNodeControl)
                driftList.append(abs(dispTop - dispBot)/height)
            tagNodePrv = tagNode
        driftMax = max(driftList)
        driftAve = sum(driftList)/len(driftList)
    else:
        height  = ops.nodeCoord(tagNodeLoad)[1] - ops.nodeCoord(tagNodeBase)[1]
        dispTop = ops.nodeDisp(tagNodeLoad, dofNodeControl)
        dispBot = ops.nodeDisp(tagNodeBase, dofNodeControl)
        driftMax=driftAve= abs(dispTop - dispBot)/height
    return driftMax, driftAve

def pushoverLCF(tagNodeLoad, tagNodeBase, tagEleList, S_MS=1.5, S_M1=0.9):
    t_beg           = time.time()
    T1              = analyzeEigen(3, True)[0]
    Cvx             = verDistFact(We, T1, h_1, h_typ, n_story)
    C_V_base        = Sa(T1, S_MS, S_M1) /(R /Ie)
    V_base          = C_V_base *We
    dofNodeControl  = 1
    tagTSLinear     = 2
    ops.timeSeries('Linear',   tagTSLinear)
    tagPatternPlain = 2
    ops.pattern('Plain', tagPatternPlain, tagTSLinear)
    if type(tagNodeLoad) == list:
        tagNodeControl  = tagNodeLoad[-1]
        for i, tagNode in enumerate(tagNodeLoad):
            ops.load(tagNode, *[Cvx[i] *V_base, 0, 0])
    else:
        tagNodeControl  = tagNodeLoad
        ops.load(tagNodeControl, *[force, 0, 0])
        
    #  Define Analysis Options
    ops.wipeAnalysis()
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('FullGeneral')   # 'FullGeneral', 'UmfPack', 'SparseSYM',
    ops.test('NormUnbalance', 1e-6, 100)
    ops.algorithm('Linear')
    ops.integrator('LoadControl', 1)
    ops.analysis('Static')
    ops.analyze(1)
    # modalProp           = ops.modalProperties('-print', '-return')
    driftMax, driftAve  = calcDrift(tagNodeLoad, tagNodeBase, dofNodeControl)
    response            = {}
    shear               = []
    for tagEle in tagEleList:
        response[tagEle]= ops.eleResponse(tagEle, 'force')
        shearforce      = max(abs(response[tagEle][1]), abs(response[tagEle][4]))
        shear.append(shearforce)
    shearAverage        = sum(shear)/len(shear)
    t_fin               = time.time()
    duration            = t_fin - t_beg
    mins = int(duration/60)
    secs = int(duration%60)
    print(f"\n\n\n{'|LCF|'*12}")
    print(f"Monotonic LC Pushover Analysis Finished in {mins}min+{secs}sec.")
    print(f"{'|LCF|'*12}\n\n\n")
    return T1, driftMax, V_base, shearAverage
    

def cyclicAnalysis(dispList, incrInit, tagNodeLoad, Beams, numIncrInit=2):
    asTagNodeBase   = 1 #it is not going to be used at all in this analysis. it is just to fill a positional argument
    t_beg           = time.time()
    dofNodeControl  = 1
    tagTSLinear     = 1
    ops.timeSeries('Linear',   tagTSLinear)
    tagPatternPlain = 1
    ops.pattern('Plain', tagPatternPlain, tagTSLinear)
    if type(tagNodeLoad) == list:
        tagNodeControl  = tagNodeLoad[-1]
        n_story = len(tagNodeLoad)-1
        for i, tagNode in enumerate(tagNodeLoad):
            ops.load(tagNode, *[i/n_story, 0, 0])
    else:
        tagNodeControl  = tagNodeLoad
        ops.load(tagNodeControl, *[1, 0, 0])
    
    #  Define Analysis Options
    ops.wipeAnalysis()
    ops.constraints('Transformation')
    ops.numberer('RCM') # Plain, RCM, AMD, ParallelPlain, ParallelRCM
    ops.system('FullGeneral') # BandGen, BandSPD, ProfileSPD, SuperLU, UmfPack, FullGeneral, SparseSYM, ('Mumps', '-ICNTL14', icntl14=20.0, '-ICNTL7', icntl7=7)
    
    # Run Analysis
    for dispIndex, disp in enumerate(dispList):
        print(f"\n\ndisp({dispIndex+1}/{len(dispList)})\t= {disp}")
        dispTargetList  = [disp, -disp]
        # dispTargetList  = [disp, 0, -disp, 0]
        for index, dispTarget in enumerate(dispTargetList):
            curD        = ops.nodeDisp(tagNodeControl, dofNodeControl)
            delta       = dispTarget - curD
            numFrac     = int(abs(delta)/incrInit)
            if numFrac == 0: numFrac=1
            incrFrac    = delta/numFrac
            OK          = convergeIt('Cyclic', tagNodeLoad, asTagNodeBase, dofNodeControl, 
                                     incrFrac, numFrac, disp, dispIndex, dispList, 
                                     dispTarget, t_beg, Beams, numIncrInit)
            if OK < 0: break
        if OK < 0: break
    return OK


def NTHA1(tagNodeLoad, tagNodeBase, filePath, scaleFactor, dtGM, NPTS, Tmax, tag, Beams, numIncrInit=2):
    # ops.wipeAnalysis()
    t_beg           = time.time()
    rayleighDamping(0.05)
    dofNodeControl  = 1
    tagTSPath       = tag
    ops.timeSeries('Path', tagTSPath, '-dt', dtGM, '-filePath', filePath, '-factor', scaleFactor *g)
    tagPatternNTHA  = tag
    direction       = 1
    ops.pattern('UniformExcitation', tagPatternNTHA, direction,'-accel', tagTSPath)
    
    #  Define Analysis Options
    ops.wipeAnalysis()
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('FullGeneral')
    
    # Run Analysis
    OK = convergeIt('NTHA', tagNodeLoad, tagNodeBase, dofNodeControl, dtGM, NPTS, Tmax, 0, 
                    ["No list required!"], Tmax, t_beg, Beams, numIncrInit)
    ops.wipeAnalysis()
    return OK


def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num - 1] = text + '\n'  # array index starts at 0, subtract 1
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()

def get_file_names(directory): # This functions returns a list containing the file names in the given directory
    return [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]

def read_ground_motion_record(filename):
    with open(filename, 'r') as file:
        data = file.readlines()
    # Read all numbers from each line
    ground_motion_record = [float(num) for line in data for num in line.split()]
    return np.array(ground_motion_record)

def get_spectral_acceleration(filePath, dt, T1, SDC, RAS_average, outputDirIDA, approach=1):
    Ti          = 0.2 *T1
    Tf          = 1.5 *T1
    def scale_Interval():
        arrSF       = np.empty((0, 2))
        for i, T in enumerate(periods):
            if Ti <= T <= Tf:
                Sa_MCE90    = 0.9 *Sa(T, S_MS, S_M1)
                indexT      = np.where(RAS_average[:,0]==T)[0][0]
                Sa_Avera    = RAS_average[indexT, 1]
                SF          = Sa_MCE90 /Sa_Avera
                arrSF = np.append(arrSF, np.array([[T, SF]]), axis=0)
        # print(f"{arrSF =}")
        SFmax       = max(arrSF[:,1])
        return SFmax
    
    def Sa(T,S_MS=1.5, S_M1=0.9):
        Ts = S_M1/S_MS
        T0 = 0.2*Ts
        TL = 8
        if 0 <= T < T0:
            Sa  = (0.4 +0.6 *T /T0) *S_MS
        elif T0 <= T < Ts:
            Sa  = S_MS
        elif Ts <= T < TL:
            Sa  = S_M1 /T
        elif T >= TL:
            Sa  = S_M1 *TL /T **2
        return Sa
    S_MS, S_M1 = (1.50, 0.9) if SDC == "Dmax" else (0.75, 0.30)

    periods = RAS_average[:,0]
    a       = read_ground_motion_record(filePath)
    record  = eqsig.AccSignal(a * 1, dt)
    record.generate_response_spectrum(response_times=periods)
    
    Sa_MCE  = [Sa(T, S_MS, S_M1) for T in periods]
    Sa_MCE90= [0.9*Sa(T, S_MS, S_M1) for T in periods]
    Sa_GMR  = record.s_a # This gives the array of points corresponding to periods
    Sa_AVE  = RAS_average[:,1]
    
    # Find the corresponding value on the vertical axis
    S_MT    = np.interp(T1, periods, Sa_MCE)
    S_GT    = np.interp(T1, periods, Sa_GMR)
    S_AT    = np.interp(T1, periods, Sa_AVE)
    SaMax   = max(S_MT, S_GT, S_AT)
    Sai     = Sa(Ti, S_MS, S_M1)
    Saf     = Sa(Tf, S_MS, S_M1)
    rec     = filePath[9:-4]
    if 1:
        fig, ax = plt.subplots(figsize=(7, 5), dpi=150)
        ax.set_ylabel('Sa [g]')
        ax.set_xlabel('Period [sec]')
        plt.plot(periods, Sa_MCE,           'r-',  label=f"MCE_AS    [SDC: {SDC}]")
        plt.plot(periods, Sa_MCE90,         'r-.', label=f"90%MCE_AS [SDC: {SDC}]")
        plt.plot(periods, Sa_AVE,           'k--', label="Average R_AS")
        plt.plot(periods, Sa_GMR,           'b--', label=f"GMR R_AS: {rec}")
        plt.plot([T1, T1], [0,    SaMax],   'g-',  label=f"T1      ={T1:.4f} sec")
        plt.plot([Ti, Ti], [0,    Sai  ],   'y--', label=f"0.2T1 ={Ti:.4f} sec")
        plt.plot([Tf, Tf], [0,    Saf  ],   'y--', label=f"1.5T1 ={Tf:.4f} sec")
        plt.plot([0,  T1], [S_MT, S_MT],    'r--', label=f"S_MT ={S_MT:.4f} g")
        plt.plot([0,  T1], [S_AT, S_AT],    'k--', label=f"S_AT  ={S_AT:.4f} g")
        plt.plot([0,  T1], [S_GT, S_GT],    'b--', label=f"S_GT  ={S_GT:.4f} g")
        if approach == 1:
            SF_MCE = S_MT /S_AT
            Sa_AVEs= SF_MCE *Sa_AVE
            plt.plot(periods, Sa_AVEs,      'k-',  label=f"Average R_AS [SF_MCE={SF_MCE:.3f}]")
        elif approach == 2:
            SF_MCE = S_MT /S_GT
            Sa_GMRs= SF_MCE *Sa_GMR
            plt.plot(periods, Sa_GMRs,      'b-',  label=f"GMR R_AS [SF_MCE={SF_MCE:.3f}]")
        elif approach == 3:
            SF_MCE = scale_Interval()
            Sa_AVEs= SF_MCE *Sa_AVE
            plt.plot(periods, Sa_AVEs,      'k-',  label=f"Average R_AS [SF_MCE={SF_MCE:.3f}]")
            
        fig.suptitle("Acceleration Response Spectrum", fontsize=16)
        plt.legend()
        plt.tight_layout()
        os.makedirs(f"{outputDirIDA}/{rec}", exist_ok=True)
        plt.savefig(f"{outputDirIDA}/{rec}/Spectrum-{rec}.png")
        # plt.show()
    
    return S_MT, SF_MCE, S_GT


def plotMomCurv(outputDir, tagEle, section, typeBuild):
    if typeBuild == "CantileverColumn":
        momenIndex = 1
    else:
        momenIndex = 2
    zero        = np.array([0.])
    # Store the section's Moment in an array
    momentTXT   = np.loadtxt(f"{outputDir}/moment{tagEle}.txt", delimiter= ' ')
    moment      = np.append(zero, momentTXT[:, momenIndex])
    Mpeak       = max(moment)
    Mpeak60perc = 0.6*Mpeak
    # Store the section's top and bottom strains in an array
    SSListTop   = np.loadtxt(f"{outputDir}/SS_top{tagEle}.txt", delimiter= ' ')
    SSListBot   = np.loadtxt(f"{outputDir}/SS_bot{tagEle}.txt", delimiter= ' ')
    # Calculate the curvature 
    StrainTop   = np.append(zero, SSListTop[:,1])
    StrainBot   = np.append(zero, SSListBot[:,1])
    h           = section.Hw + section.tw
    curvature   = ((StrainTop -StrainBot)
                   /h)
    curAtM60per = np.interp(Mpeak60perc, moment, curvature)
    EI          = Mpeak60perc /curAtM60per
    curAtMpeakE = 1/EI *Mpeak
    MomCurv     = np.column_stack((curvature, moment))
    np.savetxt(f"{outputDir}/MomentCurvature_{tagEle}.txt", MomCurv)
    
    fig, ax = plt.subplots()
    fig.suptitle(f"Momemnt-Curvature: {tagEle}")
    ax.set_xlabel(f'curvature (m^-1)')
    ax.set_ylabel('Moment (kN.m)')
    plt.plot(curvature, moment*N/kN, linewidth=0.8, label=f"M-c: Mpeak={Mpeak*N/kN:.1f} kN.m")
    plt.plot([curAtM60per, curAtM60per], [0, Mpeak60perc*N/kN], 'r--')        # Vertical Line
    plt.plot([0, curAtM60per], [Mpeak60perc*N/kN, Mpeak60perc*N/kN], 'r--')   # Horizontal Line
    plt.plot([0, curAtMpeakE], [0, Mpeak*N/kN], 'g--', label = f" EI = {EI/(kN*m**2):.1f} kN.m^2")                        # Horizontal Line
    plt.legend()
    # plt.show()
    
    return(EI)

















