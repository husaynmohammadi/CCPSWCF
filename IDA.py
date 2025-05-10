import os, sys, time, shutil
import openseespy.opensees     as ops
import functions.FuncPlot      as fp
import functions.FuncAnalysis  as fa
import numpy                   as np
# import matplotlib.pyplot       as plt


#=============================================================================
#    Options
#=============================================================================
exec(open("Input/inputData.py").readlines()[16])  # Assigns SDC from input
# SDC     = "Dmax"  # "Dmax", "Dmin"
recordToLogIDA  = True

AlgorithmIDA    = "Hunt_Fill2"       # Hunt_Fill, list__SF_CLP, manual__SF_CLP
# SF_CLP          = 20                # If "manual__SF_CLP" is selected above

approach        = 3                 # "1" for scaling the Average RAS to MCE_AS at T1, 
                                    # "2" for scaling GMR RAS to MCE_AS at T1, 
                                    # "3" for scaling the Average RAS to MCE_AS between 0.2T1 to 1.5T1

extraTime       = 0

#=============================================================================
#    Pre-Settings for IDA
#=============================================================================
g = 9.80665
start_timeIDA = time.time()
fa.replace_line('MAIN.py', 27, "recordToLog     = False                      # True, False")
fa.replace_line('MAIN.py', 31, "linearity       = False")
fa.replace_line('MAIN.py', 32, "typeBuild       = 'coupledWalls'            # 'CantileverColumn', 'coupledWalls', 'buildBeam', 'ShearCritBeam'")
fa.replace_line('MAIN.py', 34, "typeAnalysis    = ['NTHA']             # 'monotonic', 'cyclic', 'NTHA'")
fa.replace_line('MAIN.py', 78, "plot_MomCurv    = False")
fa.replace_line('MAIN.py', 83, "Pu_1wall        = -load['wallG']")

recList         = fa.get_file_names("Input/GM")
RAS_average     = np.loadtxt("Input/GM/RAS_average.txt")
outputDirIDA    = "Output/IDA"; os.makedirs(outputDirIDA, exist_ok=True)

def storeData(filePath, array):
    """
    Writes a 2D array to a text file. If the file exists, it reads the existing data,
    appends new data while avoiding duplicates, and then writes the combined data back to the file.
    
    Parameters:
    - filePath: Path to the text file
    - array: 2D numpy array with two columns
    """
    if os.path.exists(filePath):
        existing_data = np.loadtxt(filePath, delimiter=',')
        if existing_data.size == 0:
            existing_data = array
        else:
            # Combine existing data and new data, avoiding duplicates
            combined_data = np.vstack((existing_data, array))
            combined_data = np.unique(combined_data, axis=0)
    else:
        combined_data = array
    sorted_data = combined_data[combined_data[:, 1].argsort()]
    # with open(filePath, 'w') as file:
    #     for row in combined_data:
    #         file.write(f"{row[0]}\t{row[1]}\n")
    np.savetxt(filePath, sorted_data, delimiter=',', fmt='%s')
#=============================================================================
#    Start IDA
#=============================================================================
# Find T1
with open("MAIN.py") as file:
    lines = file.readlines()[:126]
code_to_exec = ''.join(lines)
exec(code_to_exec)
T1 = Periods[0]

# Record to log option check
if recordToLogIDA == True: sys.stdout = open('logIDA.txt', 'w') 

# Begin
rec_i = 1
rec_f = rec_i+1
recList     = recList[rec_i-1:rec_f-1]
numRecords  = len(recList)
list_S_CT   = []
for i_rec, rec in enumerate(recList):
    t_begIDA1   = time.time()
    filePath    = f"Input/GM/{rec}" 
    dtGM        = float(rec[3:10])
    NPTS        = int(rec[11:16])
    duration    = NPTS *dtGM +extraTime
    S_MT,SF_MCE, S_GT = fa.get_spectral_acceleration(filePath, dtGM, T1, SDC, RAS_average, outputDirIDA, approach) #Here it should scale the RAS per SDC (Dmax or Dmin) and T1
    
    if AlgorithmIDA == "list__SF_CLP":
        SF_CLPList      = [ 
                            0.1, #1
                            1.00, #2
                            2.00, #3
                            3.00, #4
                            4.00, #5
                            5.00, #6
                            6.00, #7
                            7.00, #8
                            8.00, #9
                            9.00, #10
                            10., #11
                            ]
        list_driftMax   = [0]
        list_SCTtest    = [0]
        for tag, SF_CLP in enumerate(SF_CLPList):
            S_CTtest    = SF_CLP *S_MT
            list_SCTtest.append(S_CTtest)
            print(f"\n\n\n\n\n{'#'*65}")
            print(f"Running record {i_rec+1:02}/{numRecords:02}: {rec}")
            print(f"Sa = {S_MT:.3f}*g SF = {SF_CLP:.3f} ==> S_CTtest = {S_CTtest:.3f}")
            print(f"{'#'*65}\n\n")
            
            
            ops.wipe(); exec(open("MAIN.py").read())
            list_driftMax.append(100 *abs(driftMax))
            
            durIDA1         = time.time() - t_begIDA1; mins = int(durIDA1 /60)
            SF_CLPi, IDRi   = fp.plotIDA(list_driftMax, list_SCTtest, outputDirIDA, rec, mins, SF_CLP)
            S_CTi           = SF_CLPi *S_MT
            arrIDA          = np.array([[IDRi, S_CTi]])
            storeData(f"{outputDirIDA}/{rec[:-4]}/IDA-{rec[:2]}.csv", arrIDA)
            if abs(driftMax) >= 0.2:
                break
    
    elif AlgorithmIDA == "Hunt_Fill":
        SF_CLP20        = 21 #as min
        Sa1             = 5
        SF_CLP          = Sa1 #Since Main.py puts SF_CLP in as a variable
        tag             = 1
        
        # Find Sa20min
        list_driftMax   = [0]
        list_SCTtest    = [0]
        flagLessThan20  = 0
        # SF_CLP20        = 1000
        while True:
            # First try
            ops.wipe();  exec(open("MAIN.py").read())
            driftMAX    = 100 *abs(driftMax)
            
            # Case #) If IDR >=20
            if driftMAX >= 20:
                print(f"driftMAX > 20")
                if flagLessThan20 == 0:
                    print("flagLessThan20 == 0")
                    if Sa1 <= SF_CLP20:
                        SF_CLP20 = Sa1
                        print("SF_CLP20 = {SF_CLP20:.4f}")    
                    Sa1     = Sa1 *0.7
                    SF_CLP  = Sa1
                    print(f"Sa1 = {Sa1:.4f}")
                else:
                    print("flagLessThan20 == 1")
                    if SF_CLP <= SF_CLP20:
                        SF_CLP20 = SF_CLP
                        print("SF_CLP20 = {SF_CLP20:.4f}") 
                    Sa1     = Sa2 + (Sa1 - Sa2)/3
                    SF_CLP  = Sa1
                    print(f"Sa1 = {Sa1:.4f}")
                    print(f"Sa2 = {Sa2:.4f}")
                arrayNumber = np.array([[20,SF_CLP20*S_MT]])
                storeData(f"{outputDirIDA}/{rec[:-4]}/SF_CLP20min.csv", arrayNumber)
            
            # Case #) If IDR <= 12
            elif driftMAX <= 12:
                print(f"driftMAX < 12")
                print(f"driftMAX = {driftMAX:.2f}")
                flagLessThan20 = 1
                list_SCTtest.append(SF_CLP *S_MT)
                list_driftMax.append(driftMAX)
                
                Sa2     = Sa1
                Sa1     = Sa1 *1.5
                if Sa1 <= SF_CLP20:
                    SF_CLP  = Sa1
                else:
                    Sa1 = Sa2 + (SF_CLP20 - Sa2)/2
                    SF_CLP  = Sa1
                print(f"Sa1 = {Sa1:.4f}")
                print(f"Sa2 = {Sa2:.4f}")
                durIDA1         = time.time() - t_begIDA1; mins = int(durIDA1 /60)
                SF_CLPi, IDRi   = fp.plotIDA(list_driftMax, list_SCTtest, outputDirIDA, rec, mins, SF_CLP*S_MT)
                S_CTi           = SF_CLPi *S_MT
                arrIDA          = np.array([[IDRi, S_CTi]])
                storeData(f"{outputDirIDA}/{rec[:-4]}/IDA.csv", arrIDA)
            
            
            # Case 1) If 12 < IDR < 20:
            elif 12 < driftMAX < 20:
                SF_CLPtarget    = SF_CLP
                list_SCTtest.append(SF_CLP *S_MT)
                list_driftMax.append(driftMAX)
                fractionizerList= [ 
                                    0.95,
                                    0.9,
                                    0.85,
                                    0.7,
                                    0.6,
                                    0.5,
                                    0.1,
                                    ]
                for tag, fractionizer in enumerate(fractionizerList):
                    SF_CLP      = fractionizer *SF_CLPtarget
                    S_CTtest    = SF_CLP *S_MT
                    list_SCTtest.append(S_CTtest)
                    print(f"\n\n\n\n\n{'#'*65}")
                    print(f"Running record {i_rec+1:02}/{numRecords:02}: {rec}")
                    print(f"Sa = {S_MT:.3f}*g SF = {SF_CLP:.3f} ==> S_CTtest = {S_CTtest:.3f}")
                    print(f"{'#'*65}\n\n")
                    
                    
                    ops.wipe(); exec(open("MAIN.py").read())
                    list_driftMax.append(100 *abs(driftMax))
                                        
                    durIDA1         = time.time() - t_begIDA1; mins = int(durIDA1 /60)
                    # Sort and update lists before plotting
                    arrayPLOT = np.array([list_driftMax, list_SCTtest]).T
                    sorted_arrayPLOT = arrayPLOT[arrayPLOT[:, 1].argsort()].T
                    list_driftMax = sorted_arrayPLOT[0].tolist()
                    list_SCTtest = sorted_arrayPLOT[1].tolist()
                    SF_CLPi, IDRi   = fp.plotIDA(list_driftMax, list_SCTtest, outputDirIDA, rec, mins, SF_CLP*S_MT)
                    S_CTi           = SF_CLP *S_MT
                    arrIDA          = np.array([[IDRi, S_CTi]])
                    storeData(f"{outputDirIDA}/{rec[:-4]}/IDA-{rec[:2]}.csv", sorted_arrayPLOT.T)
                break
                    
    elif AlgorithmIDA == "Hunt_Fill2":
        print('AlgorithmIDA == "Hunt_Fill2"')
        SF_CLPList      = [ 
                            0.1,
                            1.00,
                            2.00,
                            3.00,
                            4.00,
                            5.00,
                            6.00,
                            7.00,
                            8.00,
                            9.00,
                            10.,
                            ]
        list_driftMax   = [0]
        list_SCTtest    = [0]
        for tag, SF_CLP in enumerate(SF_CLPList):
            S_CTtest    = SF_CLP *S_MT
            print(f"\n\n\n\n\n{'#'*65}")
            print(f"Running record {i_rec+1:02}/{numRecords:02}: {rec}")
            print(f"Sa = {S_MT:.3f}*g SF = {SF_CLP:.3f} ==> S_CTtest = {S_CTtest:.3f}")
            print(f"{'#'*65}\n\n")
            ops.wipe(); exec(open("MAIN.py").read())
            
            if abs(driftMax) < 0.2:
                SF_CLPmaxBefore10 = SF_CLP
                print(f"SF_CLPmaxBefore10 = {SF_CLPmaxBefore10}")
                
                durIDA1         = time.time() - t_begIDA1; mins = int(durIDA1 /60)
                list_SCTtest.append(S_CTtest)
                list_driftMax.append(100 *abs(driftMax))
                SF_CLPi, IDRi   = fp.plotIDA(list_driftMax, list_SCTtest, outputDirIDA, rec, mins, SF_CLP)
                S_CTi           = SF_CLPi *S_MT
                arrIDA          = np.array([[IDRi, S_CTi]])
                storeData(f"{outputDirIDA}/{rec[:-4]}/IDA-{rec[:2]}.csv", arrIDA)
            elif abs(driftMax) >= 0.2:
                SF_CLPmin20 = SF_CLP
                print(f"SF_CLPmin20 = {SF_CLPmin20}")
                div = 5
                SF_CLPList = [(SF_CLPmaxBefore10 + i/div * (SF_CLPmin20 - SF_CLPmaxBefore10)) for i in range(div+1)]
                for tag, SF_CLP in enumerate(SF_CLPList):
                    S_CTtest    = SF_CLP *S_MT
                    print(f"\n\n\n\n\n{'#'*65}")
                    print(f"Running record {i_rec+1:02}/{numRecords:02}: {rec}")
                    print(f"Sa = {S_MT:.3f}*g SF = {SF_CLP:.3f} ==> S_CTtest = {S_CTtest:.3f}")
                    print(f"{'#'*65}\n\n")
                    ops.wipe(); exec(open("MAIN.py").read())
                    
                    if abs(driftMax) < 0.2:
                        SF_CLPmaxBefore10 = SF_CLP
                        print(f"SF_CLPmaxBefore10 = {SF_CLPmaxBefore10}")
                        
                        durIDA1         = time.time() - t_begIDA1; mins = int(durIDA1 /60)
                        list_SCTtest.append(S_CTtest)
                        list_driftMax.append(100 *abs(driftMax))
                        SF_CLPi, IDRi   = fp.plotIDA(list_driftMax, list_SCTtest, outputDirIDA, rec, mins, SF_CLP)
                        S_CTi           = SF_CLPi *S_MT
                        arrIDA          = np.array([[IDRi, S_CTi]])
                        storeData(f"{outputDirIDA}/{rec[:-4]}/IDA-{rec[:2]}.csv", arrIDA)
                    elif abs(driftMax) >= 0.2:
                        SF_CLPmin20 = SF_CLP
                        print(f"SF_CLPmin20 = {SF_CLPmin20}")
                        SF_CLPList = [(SF_CLPmaxBefore10 + i/div * (SF_CLPmin20 - SF_CLPmaxBefore10)) for i in range(div+1)]
                        for tag, SF_CLP in enumerate(SF_CLPList):
                            S_CTtest    = SF_CLP *S_MT
                            print(f"\n\n\n\n\n{'#'*65}")
                            print(f"Running record {i_rec+1:02}/{numRecords:02}: {rec}")
                            print(f"Sa = {S_MT:.3f}*g SF = {SF_CLP:.3f} ==> S_CTtest = {S_CTtest:.3f}")
                            print(f"{'#'*65}\n\n")
                            ops.wipe(); exec(open("MAIN.py").read())
                            list_SCTtest.append(S_CTtest)
                            list_driftMax.append(100 *abs(driftMax))
                            durIDA1         = time.time() - t_begIDA1; mins = int(durIDA1 /60)
                            SF_CLPi, IDRi   = fp.plotIDA(list_driftMax, list_SCTtest, outputDirIDA, rec, mins, SF_CLP)
                            S_CTi           = SF_CLPi *S_MT
                            arrIDA          = np.array([[IDRi, S_CTi]])
                            storeData(f"{outputDirIDA}/{rec[:-4]}/IDA-{rec[:2]}.csv", arrIDA)
                            if abs(driftMax) >= 0.2:
                                break
                        break
                break
        
    # Save the log to record directory
    shutil.copy('logIDA.txt', f'{outputDirIDA}/{rec[:-4]}/logIDA.txt')
        
print(f"list_S_CT = {list_S_CT}") 

finish_timeIDA  = time.time()
elapsedTime     = finish_timeIDA - start_timeIDA
mins = int(elapsedTime /60); secs = int(elapsedTime %60)
print("\n\n\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
print(f"IDA Finished in {mins}min+{secs}sec.")
print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n\n")

if recordToLogIDA == True:
    sys.stdout.close()
    sys.stdout = sys.__stdout__




fa.replace_line('MAIN.py', 27, "recordToLog     = True                      # True, False")




























