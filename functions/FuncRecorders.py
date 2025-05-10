import openseespy.opensees as ops




def recordPushover(tagNodeControl, tagNodeBase, outputDir):
    
    #   recorder('Node', '-file', filename,                        '-time', '-node', *nodeTags=[], '-dof', *dofs=[], respType)
    if type(tagNodeControl) == list:
        ops.recorder('Node', '-file', f"{outputDir}/top_disp.txt", '-time', '-node', tagNodeControl[-1], '-dof', 1,        'disp')
    else:
        ops.recorder('Node', '-file', f"{outputDir}/top_disp.txt", '-time', '-node', tagNodeControl,     '-dof', 1,        'disp')
    if type(tagNodeBase) == list:
        ops.recorder('Node', '-file', f"{outputDir}/reaction.txt", '-time', '-node', *tagNodeBase,       '-dof', 1,        'reaction')
    else: 
        ops.recorder('Node', '-file', f"{outputDir}/reaction.txt", '-time', '-node',  tagNodeBase,       '-dof', 1,        'reaction')

def recordStressStrain(outputDir, tagEleList, section):
    
    # coords_Fiber
    
    Hw = section.Hw
    tw = section.tw
    H  = Hw + 2*tw
    
    coordTop = [ H/2, 0]
    coordBot = [-H/2, 0]
    coordMid = [0, 0.01]
    coords   = {
        'top': coordTop,
        'mid': coordMid,
        'bot': coordBot
        }

    #   recorder('Element', '-file', filename,                          '-ele', *eleTags=[], *args)
    #                                                                                        *args for dispBeamColumn elements are  when creating an ElementRecorder object are 'force,' and 'section $secNum secArg1 secArg2...'
    #                                                                                        'section $secNum secArg1 secArg2...'  Where $secNum refers to the integration point whose data is to be output valid entries being 1 through $numIntgrPts.
    #                                                                                                         'fiber', y, z,                tagMat, 'stressStrain'        #(It can also be used as 'fiber', y, z, tagMat, 'stressStrain' so as to record the closest fiber to [y,z] with material tag tagMat)
    ##  Top Flange
    ops.recorder('Element', '-file', f"{outputDir}/fiberSt_top.txt", '-ele', *tagEleList,           
                 'section', section.tagSec,    'fiber', *coordTop, section.tagMatStFlange, 'stressStrain'); print(f"tagMatStFlange \t= {section.tagMatStFlange}")
    ##  Bottom Flange
    ops.recorder('Element', '-file', f"{outputDir}/fiberSt_bot.txt", '-ele', *tagEleList,           
                 'section', section.tagSec,    'fiber', *coordBot, section.tagMatStFlange, 'stressStrain')
    ##  Mid Web
    ops.recorder('Element', '-file', f"{outputDir}/fiberSt_mid.txt", '-ele', *tagEleList,           
                 'section', section.tagSec,    'fiber', *coordMid, section.tagMatStWeb, 'stressStrain'); print(f"tagMatStWeb \t= {section.tagMatStWeb}")
    ##  Top Confined Concrete
    ops.recorder('Element', '-file', f"{outputDir}/fiberCt2_top.txt", '-ele', *tagEleList,           
                  'section', section.tagSec,    'fiber', *coordTop, section.tagMatCtConf, 'stressStrain'); print(f"tagMatCtConf \t= {section.tagMatCtConf}")
    ##  Bottom Confined Concrete
    ops.recorder('Element', '-file', f"{outputDir}/fiberCt2_bot.txt", '-ele', *tagEleList,           
                  'section', section.tagSec,    'fiber', *coordBot, section.tagMatCtConf, 'stressStrain')
    ##  Mid Confined Concrete
    ops.recorder('Element', '-file', f"{outputDir}/fiberCt2_mid.txt", '-ele', *tagEleList,           
                  'section', section.tagSec,    'fiber', *coordMid, section.tagMatCtConf, 'stressStrain')
    ##  Top Unconfined Concrete
    ops.recorder('Element', '-file', f"{outputDir}/fiberCt1_top.txt", '-ele', *tagEleList,           
                 'section', section.tagSec,    'fiber', *coordTop, section.tagMatCtUnconf, 'stressStrain'); print(f"tagMatCtUnconf \t= {section.tagMatCtUnconf}")
    ##  Bottom Unconfined Concrete
    ops.recorder('Element', '-file', f"{outputDir}/fiberCt1_bot.txt", '-ele', *tagEleList,           
                 'section', section.tagSec,    'fiber', *coordBot, section.tagMatCtUnconf, 'stressStrain')
    ##  Mid Unconfined Concrete
    ops.recorder('Element', '-file', f"{outputDir}/fiberCt1_mid.txt", '-ele', *tagEleList,           
                 'section', section.tagSec,    'fiber', *coordMid, section.tagMatCtUnconf, 'stressStrain')
    
    return(coords)
    

def recordDataNTHA(tagNodeBase, tagNodeLoad, outputDir, tag):
    #01 Record Displacements, Velocities, and Accelerations
    if type(tagNodeLoad) != list:
        tagNodeLoading = [tagNodeLoad]
    else:
        tagNodeLoading = tagNodeLoad #[1:]
    ops.recorder('Node',         '-file', f"{outputDir}/disp{tag}.txt",   '-time', '-node', *tagNodeLoading, '-dof', 1, 'disp')
    ops.recorder('Node',         '-file', f"{outputDir}/velo{tag}.txt",   '-time', '-node', *tagNodeLoading, '-dof', 1, 'vel')
    ops.recorder('Node',         '-file', f"{outputDir}/acce{tag}.txt",   '-time', '-node', *tagNodeLoading, '-dof', 1, 'accel')
    ops.recorder('EnvelopeNode', '-file', f"{outputDir}/acceEnv{tag}.txt",'-time', '-node', *tagNodeLoading, '-dof', 1, 'accel')
    
    #02 Record Reactions
    if type(tagNodeBase) != list:
        tagNodeBase = [tagNodeBase]
    ops.recorder('Node',         '-file', f"{outputDir}/R{tag}.txt",      '-time', '-node', *tagNodeBase, '-dof', *[1, 2, 3], 'reaction')
    # ops.recorder('Node',         '-file', f"{outputDir}/Vy.txt",      '-time', '-node', *tagNodeBase, '-dof', 2, 'reaction')
    # ops.recorder('Node',         '-file', f"{outputDir}/Mz.txt",      '-time', '-node', *tagNodeBase, '-dof', 3, 'reaction')
    
    
    

def recordMomCurv(tagNode, tagEle, section, outputDir):
    
    # Record at each step: the moment of the element at a section at the given node
    ops.recorder('Node', '-file', f"{outputDir}/moment{tagEle}.txt", '-time', '-node',  tagNode, '-dof', 3, 'reaction')
    
    # Record at each step: the strains in the fibers farthest from the neutral axis of the section; one at top and another at bottom
    Hw = section.Hw
    tp = section.tw
    H  = Hw +2 *tp
    
    coordTop = [ H /2, 0]
    coordBot = [-H /2, 0]
    
    ##  Top Flange
    ops.recorder('Element', '-file', f"{outputDir}/SS_top{tagEle}.txt", '-ele', tagEle,           
                 'section', section.tagSec,    'fiber', *coordTop, section.tagMatStFlange, 'stressStrain')
    ##  Bottom Flange
    ops.recorder('Element', '-file', f"{outputDir}/SS_bot{tagEle}.txt", '-ele', tagEle,           
                 'section', section.tagSec,    'fiber', *coordBot, section.tagMatStFlange, 'stressStrain')
































