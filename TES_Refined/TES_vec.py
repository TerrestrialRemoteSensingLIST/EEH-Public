# Temperature and Emissivity Separation algorithm
import numpy as np

def nem(e_max, ls2, ls4, ls5, ld2, ld4, ld5):
    print('NEM starting...')
    
    N = 15   
    threshold1 = 0.007 # NET
    threshold2 = 0.007 # NET
    
    c1 = 1.19104E+8 # unit: W.m-2
    c2 = 14388 # unit: um.K
    lamda = np.array([8.8161,10.4919,12.0876])
    
    # Output variables 
    qa = np.zeros(ls2.shape)
    t_nem = np.zeros(ls2.shape)
    eb2_nem = np.zeros(ls2.shape)
    eb4_nem = np.zeros(ls2.shape)
    eb5_nem = np.zeros(ls2.shape)
    
    # Flag variables
    flag = np.zeros(ls2.shape) #0: continue, 1: break
    count = np.zeros(ls2.shape)
    
    eb2_max = np.ones(ls2.shape)*e_max
    eb4_max = np.ones(ls2.shape)*e_max
    eb5_max = np.ones(ls2.shape)*e_max
    
    # Surface-leaving radiance
    radiance2 = []
    radiance4 = []
    radiance5 = []
    
    radiance2.append(ls2 - (1-eb2_max)*ld2) 
    radiance4.append(ls4 - (1-eb4_max)*ld4) 
    radiance5.append(ls5 - (1-eb5_max)*ld5) 
    
    for n in range(N):
        # Determine whether to continue
        if np.all(flag):
            print('NEM module completed')
            break          
        
        print('NEM Iteration: ', n+1)
        
        print('Undecided pixels: ',qa[flag == 0].shape[0])
                
        i1 = np.logical_and(np.logical_or.reduce((radiance2[n] <= 0, radiance4[n] <= 0, radiance5[n] <= 0)), flag == 0)
        qa[i1] = -5
        flag[i1] = 1
        count[i1] = n
        
        i2 = np.logical_and.reduce((radiance2[n] > 0, radiance4[n] > 0, radiance5[n] > 0, flag == 0))
        t2 = c2/lamda[0]/(np.log(1+c1*eb2_max[i2]/(lamda[0]**5*radiance2[n][i2])))
        t4 = c2/lamda[1]/(np.log(1+c1*eb4_max[i2]/(lamda[1]**5*radiance4[n][i2])))
        t5 = c2/lamda[2]/(np.log(1+c1*eb5_max[i2]/(lamda[2]**5*radiance5[n][i2])))
        t_nem[i2] = np.amax(np.array([t2,t4,t5]), axis = 0)
        
        eb2_nem[i2] = radiance2[n][i2]/(c1/(lamda[0]**5*(np.exp(c2/t_nem[i2]/lamda[0])-1)))
        eb4_nem[i2] = radiance4[n][i2]/(c1/(lamda[1]**5*(np.exp(c2/t_nem[i2]/lamda[1])-1)))
        eb5_nem[i2] = radiance5[n][i2]/(c1/(lamda[2]**5*(np.exp(c2/t_nem[i2]/lamda[2])-1)))
        
        j1 = np.logical_and.reduce((radiance2[n] > 0, radiance4[n] > 0, radiance5[n] > 0, flag == 0,
                                   np.logical_or.reduce(( 
                                   eb2_nem <= 0.5, eb2_nem >= 1, eb4_nem <= 0.5, 
                                   eb4_nem >= 1, eb5_nem <= 0.5, eb5_nem >= 1))
                          ))
        qa[j1] = -3
        flag[j1] = 1
        count[j1] = n
        
#         j2 = np.logical_and.reduce((radiance2[n] > 0, radiance4[n] > 0, radiance5[n] > 0, flag == 0,
#                                     eb2_nem > 0.5, eb2_nem < 1, eb4_nem > 0.5, 
#                                     eb4_nem < 1, eb5_nem > 0.5, eb5_nem < 1, flag == 0))
        
        radiance2.append(ls2 - (1-eb2_nem)*ld2)
        radiance4.append(ls4 - (1-eb4_nem)*ld4)
        radiance5.append(ls5 - (1-eb5_nem)*ld5)
        
        eb2_max = eb2_nem.copy()
        eb4_max = eb4_nem.copy()
        eb5_max = eb5_nem.copy()
        
        if n > 1:
            k1 = np.logical_and(np.amax([np.abs(radiance2[n+1]-radiance2[n]),np.abs(radiance4[n+1]-radiance4[n]),
                                   np.abs(radiance5[n+1]-radiance5[n])],axis=0) < threshold2, flag == 0)
            qa[k1] = 1
            flag[k1] = 1
            count[k1] = n
            
            k2 = np.logical_and(np.amax([np.abs(radiance2[n+1]-radiance2[n]),np.abs(radiance4[n+1]-radiance4[n]),
                                   np.abs(radiance5[n+1]-radiance5[n])],axis=0) >= threshold2, flag == 0)
            qa[k2] = 0
            flag[k2] = 0
            count[k2] = n
            
            k3 = np.logical_and.reduce((np.amin([np.abs(radiance2[n+1]+radiance2[n-1]-2*radiance2[n]),
                                   np.abs(radiance4[n+1]+radiance4[n-1]-2*radiance4[n]),
                                   np.abs(radiance5[n+1]+radiance5[n-1]-2*radiance5[n])],axis=0) >= threshold1, 
                                   flag == 0, qa == 0))
            qa[k3] = -1
            flag[k3] = 1
            count[k3] = n
            
    l = np.logical_and.reduce((qa == 0, count == N-1, flag == 0))
    qa[l] = -2
    flag[l] = 1
    
    if np.all(qa):
        print('All pixels completed')
    else:
        print('Further check needed!!')
     
    
    return (qa, t_nem, eb2_nem, eb4_nem, eb5_nem)


#Inputs into the TES module include the at-surface radiance and atmospheric downwelling radiance
def ecostress_tes(ls,ld,alpha1,alpha2,alpha3):
    print('TES starting...')
    
    c1 = 1.19104E+8 #unit: W.m-2
    c2 = 14388 #unit: um.K
    lamda = np.array([8.8161,10.4919,12.0876])
         
    v1 = 1.7E-4
    
    #Coefficients for EEH TES
    #alpha1 = 0.9895
    #alpha2 = 0.7994
    #alpha3 = 0.8572
    
    #Coefficients for input samples for SAIL271
    # alpha1 = 0.9692
    # alpha2 = 0.8117
    # alpha3 = 0.9957
    
    #Coefficients for SAIL271
    # alpha1 = 0.9824
    # alpha2 = 0.8931
    # alpha3 = 0.9757
    
    out_emib2 = np.zeros((ls.shape[0], ls.shape[1]))
    out_emib4 = np.zeros((ls.shape[0], ls.shape[1]))
    out_emib5 = np.zeros((ls.shape[0], ls.shape[1]))
    out_t = np.zeros((ls.shape[0], ls.shape[1]))
    out_qa = np.zeros((ls.shape[0], ls.shape[1]))
    out_mmd = np.zeros((ls.shape[0], ls.shape[1]))
    
    eb2_nem = np.zeros((ls.shape[0], ls.shape[1]))
    eb4_nem = np.zeros((ls.shape[0], ls.shape[1]))
    eb5_nem = np.zeros((ls.shape[0], ls.shape[1]))
    t_nem = np.zeros((ls.shape[0], ls.shape[1]))
    qa_nem = np.zeros((ls.shape[0], ls.shape[1]))
    
    lemi_b2 = np.zeros((ls.shape[0], ls.shape[1]))
    lemi_b4 = np.zeros((ls.shape[0], ls.shape[1]))
    lemi_b5 = np.zeros((ls.shape[0], ls.shape[1]))
    
    eb2_tes = np.zeros((ls.shape[0], ls.shape[1]))
    eb4_tes = np.zeros((ls.shape[0], ls.shape[1]))
    eb5_tes = np.zeros((ls.shape[0], ls.shape[1]))
    t_tes = np.zeros((ls.shape[0], ls.shape[1]))
    mmd = np.zeros((ls.shape[0], ls.shape[1]))
    
    ## Cloud-contaminated pixels
    i1 = np.logical_or.reduce((ls[:,:,0] <= 0, ls[:,:,1] <= 0, ls[:,:,2] <= 0, 
                               np.isnan(ls[:,:,0]), np.isnan(ls[:,:,1]), np.isnan(ls[:,:,2])
                      ))
    out_qa[i1] = -4
    
    ## Normal pixels
    i2 = np.logical_and.reduce((ls[:,:,0] > 0, ls[:,:,1] > 0, ls[:,:,2] > 0))
    (qa_nem1,t_nem1,eb2_nem1,eb4_nem1,eb5_nem1) = nem(0.99,ls[:,:,0][i2],ls[:,:,1][i2],ls[:,:,2][i2],ld[:,:,0][i2],ld[:,:,1][i2],ld[:,:,2][i2])
    qa_nem[i2] = qa_nem1.copy()
    t_nem[i2] = t_nem1.copy()
    eb2_nem[i2] = eb2_nem1.copy()
    eb4_nem[i2] = eb4_nem1.copy()
    eb5_nem[i2] = eb5_nem1.copy()   
    
    # Unreasonable values occurred in NEM module
    j1 = np.logical_and.reduce((ls[:,:,0] > 0, ls[:,:,1] > 0, ls[:,:,2] > 0,
                       np.logical_or.reduce((qa_nem == -5, qa_nem == -3, qa_nem == -2, qa_nem == -1))
                      ))   
#     out_emib2[j1] = 0.
#     out_emib4[j1] = 0.
#     out_emib5[j1] = 0.
#     out_t[j1] = 0.
    out_qa[j1] = qa_nem[j1].copy()
    
#     print('out_qa[j1].shape',out_qa[j1].shape)
    
    # Pixel containing rock or soil
    j2 = np.logical_and.reduce((ls[:,:,0] > 0, ls[:,:,1] > 0, ls[:,:,2] > 0,
                       qa_nem != -5, qa_nem != -3, qa_nem != -2, qa_nem != -1, 
                       np.var(np.array([eb2_nem,eb4_nem,eb5_nem]), axis = 0) > v1
                      ))
    
    (qa_nem2,t_nem2,eb2_nem2,eb4_nem2,eb5_nem2) = nem(0.96,ls[:,:,0][j2],ls[:,:,1][j2],ls[:,:,2][j2],
                                                 ld[:,:,0][j2],ld[:,:,1][j2],ld[:,:,2][j2])
    qa_nem[j2] = qa_nem2.copy()
    t_nem[j2] = t_nem2.copy()
    eb2_nem[j2] = eb2_nem2.copy()
    eb4_nem[j2] = eb4_nem2.copy()
    eb5_nem[j2] = eb5_nem2.copy()
    
#     print('qa_nem[j2].shape',qa_nem[j2].shape)
    
    j11 = np.logical_and.reduce((ls[:,:,0] > 0, ls[:,:,1] > 0, ls[:,:,2] > 0,
                        np.logical_or.reduce((qa_nem == -5, qa_nem == -3, qa_nem == -2, qa_nem == -1))
                       ))
#     out_emib2[i2][j11] = 0
#     out_emib4[i2][j11] = 0
#     out_emib5[i2][j11] = 0
#     out_t[i2][j11] = 0
    out_qa[j11] = qa_nem[j11].copy()
    
#     print('out_qa[j11].shape',out_qa[j11].shape)
    
    # RATIO module
    j3 = np.logical_and.reduce((ls[:,:,0] > 0, ls[:,:,1] > 0, ls[:,:,2] > 0,
                       qa_nem != -5, qa_nem != -3, qa_nem != -2, qa_nem != -1))
    ratiob2 = eb2_nem[j3]/np.mean(np.array([eb2_nem[j3],eb4_nem[j3],eb5_nem[j3]]), axis = 0)
    ratiob4 = eb4_nem[j3]/np.mean(np.array([eb2_nem[j3],eb4_nem[j3],eb5_nem[j3]]), axis = 0)
    ratiob5 = eb5_nem[j3]/np.mean(np.array([eb2_nem[j3],eb4_nem[j3],eb5_nem[j3]]), axis = 0)
    print('RATIO module completed')
    
#     print('ratiob2.shape', ratiob2.shape)
    
    # MMD module
    mmd[j3] = np.amax(np.array([ratiob2,ratiob4,ratiob5]), axis = 0) - np.amin(np.array([ratiob2,ratiob4,ratiob5]), axis = 0)  
    e_min = alpha1 - alpha2*mmd[j3]**alpha3
    eb2_tes[j3] = ratiob2*e_min/np.amin(np.array([ratiob2,ratiob4,ratiob5]), axis = 0)
    eb4_tes[j3] = ratiob4*e_min/np.amin(np.array([ratiob2,ratiob4,ratiob5]), axis = 0)
    eb5_tes[j3] = ratiob5*e_min/np.amin(np.array([ratiob2,ratiob4,ratiob5]), axis = 0)    
    
    lemi_b2[j3] = ls[:,:,0][j3] - (1-eb2_tes[j3])*ld[:,:,0][j3]
    lemi_b4[j3] = ls[:,:,1][j3] - (1-eb4_tes[j3])*ld[:,:,1][j3]
    lemi_b5[j3] = ls[:,:,2][j3] - (1-eb5_tes[j3])*ld[:,:,2][j3]   
    
    k1 = np.logical_and.reduce((ls[:,:,0] > 0, ls[:,:,1] > 0, ls[:,:,2] > 0,
                       qa_nem != -5, qa_nem != -3, qa_nem != -2, qa_nem != -1, 
                       np.logical_or.reduce((lemi_b2 <= 0, lemi_b4 <= 0, lemi_b5 <= 0))
                      ))
    out_qa[k1] = -5
    
#     print('out_qa[k1].shape',out_qa[k1].shape)
    
    k2 = np.logical_and.reduce((ls[:,:,0] > 0, ls[:,:,1] > 0, ls[:,:,2] > 0, 
                       qa_nem != -5, qa_nem != -3, qa_nem != -2, qa_nem != -1,     
                       lemi_b2 > 0, lemi_b4 > 0, lemi_b5 > 0))   
    
    t_b2 = c2/lamda[0]/(np.log(1+c1*eb2_tes[k2]/(lamda[0]**5*lemi_b2[k2])))
    t_b4 = c2/lamda[1]/(np.log(1+c1*eb4_tes[k2]/(lamda[1]**5*lemi_b4[k2])))
    t_b5 = c2/lamda[2]/(np.log(1+c1*eb5_tes[k2]/(lamda[2]**5*lemi_b5[k2])))
    emi_array = np.array([eb2_tes[k2], eb4_tes[k2], eb5_tes[k2]])
    t_array = np.array([t_b2, t_b4, t_b5])
    index = np.argmax(emi_array,axis = 0)   
    t_tes[k2] = t_array[index,np.arange(t_array.shape[1])]
    
    print('MMD module completed')
    
    out_emib2[k2] = eb2_tes[k2].copy()
    out_emib4[k2] = eb4_tes[k2].copy()
    out_emib5[k2] = eb5_tes[k2].copy()    
    out_t[k2] = t_tes[k2].copy()    
    out_mmd[k2] = mmd[k2].copy() 
    out_qa[k2] = qa_nem[k2].copy()
    
#     print(out_t[k2])
                        
    return (out_qa,out_emib2,out_emib4,out_emib5,out_t,out_mmd)


def LST_Estimate(r2, r4, r5, upclear_f, dnclear_f, trans_f, alpha1, alpha2, alpha3):
    #Calculate at-surface radiance using the calculated transmittance and atmosphere upwelling radiance from at-sensor radiance

    ls2 = (r2 - upclear_f[:,:,0])/trans_f[:,:,0]
    ls4 = (r4 - upclear_f[:,:,1])/trans_f[:,:,1]
    ls5 = (r5 - upclear_f[:,:,2])/trans_f[:,:,2]

    #at-surface radiance
    ls = np.empty((r2.shape[0],r2.shape[1],3),dtype=np.float32)
    ls[:,:,0] = ls2
    ls[:,:,1] = ls4
    ls[:,:,2] = ls5

    (qa, emib2, emib4, emib5, lst, mmd) = ecostress_tes(ls,dnclear_f, alpha1, alpha2, alpha3) 
    
    print('TES module completed')
    
    return (lst, emib2, emib4, emib5, mmd, qa)
