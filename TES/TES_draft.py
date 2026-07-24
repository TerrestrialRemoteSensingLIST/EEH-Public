# Temperature and Emissivity Separation algorithm
import numpy as np

def nem(e_max,ls,ld):
    N = 13   
    t1 = 0.007 #NET
    t2 = 0.007 #NET
    
    c1 = 1.19104E+8 #unit: W.m-2
    c2 = 14388 #unit: um.K
    lamda = np.array([[8.8161,10.4919,12.0876]])

    radiance = []
    radiance.append(ls - (1-e_max)*ld) #radiance[k]
    for k in range(N): 
        if np.any(radiance[k] <= 0):
            qa = -5
            t_nem = 0.
            e_nem = np.asarray([0.,0.,0.])
            break
        #if np.any(radiance[k] <= 0)

        t = c2/lamda/(np.log(1+c1*e_max/(lamda**5*radiance[k])))
        t_nem = np.amax(t)
        e_nem = radiance[k]/(c1/(lamda**5*(np.exp(c2/t_nem/lamda)-1)))

        if not(np.amin(e_nem)>0.5 and np.amax(e_nem)<1):
            qa = -3 #emissivity out of reasonale limits
            t_nem = 0.
            e_nem = np.asarray([0.,0.,0.])
            break
        #if not(np.amin(e)>0.5 and np.amax(e)<1)

        radiance.append(ls - (1-e_nem)*ld) #radiance[k+1]
        e_max = e_nem

        if (k > 1) and np.amax(np.abs(radiance[k+1] - radiance[k])) < t2:
            qa = 1 #convergence met
            break
        else:
            qa = 0 #convergence not met, continue to iterate
        #if np.amax(np.abs(radiance[k+1] - radiance[k])) < t2

        if (qa == 0) and (k > 1) and (np.amin(np.abs(radiance[k+1] + radiance[k-1] - 2*radiance[k])) > t1):
            qa = -1 #divergence met
            t_nem = 0.
            e_nem = np.asarray([0.,0.,0.])
            break
        #if (qa == 0) and (k>=1) and (np.amin(np.abs(radiance[k+1] + radiance[k-1] - 2*radiance[k])) > t1)
    #for k in range(N)

    if (k == N-1) and (qa == 0):
        qa = -2
        t_nem = 0.
        e_nem = np.asarray([0.,0.,0.])
    #if (k == N-1) and (qa == 0)

    return qa,t_nem,e_nem
#def nem(e_max,ls,ld)

#Inputs into the TES module include the at-surface radiance and atmospheric downwelling radiance
def ecostress_tes(ls,ld):
    c1 = 1.19104E+8 #unit: W.m-2
    c2 = 14388 #unit: um.K
    lamda = np.array([[8.8161,10.4919,12.0876]])
         
    v1 = 1.7E-4
    alpha1 = 0.9895
    alpha2 = 0.7994
    alpha3 = 0.8572
    
    if ls[0] > 0 and ls[1] > 0 and ls[2] > 0:
        qa,t_nem,e_nem = nem(0.99,ls,ld)
        
        if (qa == -5) or (qa == -3) or (qa == -2) or (qa == -1):
            out_emi = e_nem
            out_t = t_nem
            out_qa = qa
            mmd = 0.
        else:
            #Adjust the maximum emissivity
            if np.var(e_nem) > v1: #pixel consisting of rock or soil 
                dorefine = 1
            else: #pixel is near graybody
                dorefine = 0
            #if np.var(e) > v1

            if dorefine:
                qa,t_nem,e_nem = nem(0.96,ls,ld)
            #if dorefine
            
            if (qa == -5) or (qa == -3) or (qa == -2) or (qa == -1):
                out_emi = e_nem
                out_t = t_nem
                out_qa = qa
                mmd = 0.
            else:               
                #RATIO module
                ratio = e_nem/np.mean(e_nem)

                #MMD module
                mmd = np.amax(ratio) - np.amin(ratio)
                e_min = alpha1 - alpha2*mmd**alpha3
                e_tes = ratio*e_min/np.amin(ratio)

                lemi = ls - (1-e_tes)*ld

                if np.any(lemi <= 0):
                    out_emi = np.asarray([0.,0.,0.])
                    out_t = 0.
                    out_qa = -5
                    mmd = 0.
                else:
                    t = c2/lamda/(np.log(1+c1*e_tes/(lamda**5*lemi)))

                    index = np.where(e_tes == np.amax(e_tes))
                    t_tes = t[index]

                    out_emi = e_tes
                    out_t = t_tes
                    out_qa = qa
    else:
        out_emi = np.asarray([0.,0.,0.])
        out_t = 0.
        out_qa = -4
        mmd = 0.
    #if ls[i,j,0] > 0 and ls[i,j,1] > 0 and ls[i,j,1] > 0           
                        
    return out_qa,out_emi,out_t,mmd
#def ecostress_tes(ls,ld)

def LST_Estimate(r2, r4, r5, upclear_f, dnclear_f, trans_f):
    #Calculate at-surface radiance using the calculated transmittance and atmosphere upwelling radiance from at-sensor radiance

    ls2 = (r2_sub - upclear_f[:,:,0])/trans_f[:,:,0]
    ls4 = (r4_sub - upclear_f[:,:,1])/trans_f[:,:,1]
    ls5 = (r5_sub - upclear_f[:,:,2])/trans_f[:,:,2]

    #at-surface radiance
    ls = np.empty((row,col,3),dtype=np.float64)
    ls[:,:,0] = ls2
    ls[:,:,1] = ls4
    ls[:,:,2] = ls5

    lst = np.empty((row,col),dtype=np.float64)
    emi = np.empty((row,col,3),dtype=np.float64)
    mmd = np.empty((row,col),dtype=np.float64)
    qa = np.empty((row,col),dtype=np.int16)
    for i in range(row):
        for j in range(col):
            qa0,emi0,lst0,mmd0 = ecostress_tes(ls[i,j,:],dnclear_f[i,j,:]) 
            lst[i,j] = lst0
            emi[i,j,:] = emi0
            mmd[i,j] = mmd0
            qa[i,j] = qa0
            
    return (lst, emi, mmd, qa)
