from PyAstronomy import pyasl
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

def yes_no(question_to_be_answered):
    while True:
        choice = input(question_to_be_answered).lower()
        if choice[:1] == 'y': 
            return True
        elif choice[:1] == 'n':
            return False
        else:
            print("Please respond with 'Yes' or 'No'\n")

def rangee(lim_i,lim_s,default,question_to_be_ans):
    while True:
        text_q = question_to_be_ans

        user_number=default
        user_number = input (text_q)
        if not user_number:
                user_number=default
        try:
            val = float(user_number)
            if((val >=lim_i) & (val<=lim_s)):
                return False,val
            else:
                print("Out of range")
     
        except ValueError:
            print("It's a string")

def gauss_mod(c):

    alfa=np.empty(np.size(lline))*0
    beta=np.empty(np.size(lline))*0
    gammap=np.empty(np.size(lline))*0
    delta=np.empty(np.size(lline))*0
    me1m=np.empty(np.size(lline))*0

    me1=100
    scme1=100
    A=c-ygaussmin
    xo = xgaussmin
    si = sigline
    contr=1
    
    while ((scme1>0.001) & (contr<20)):
        for i in range(0, (np.size(lline)-1)):
            alfa[i]=fline[i]-(c-A*(np.exp(-(lline[i]-xo)**2/(2*si**2))))
            beta[i]=-(np.exp(-(lline[i]-xo)**2/(2*si**2)))
            gammap[i]=-(A/si**2)*(lline[i]-xo)*np.exp(-(lline[i]-xo)**2/(2*si**2))
            delta[i]=-((A*(lline[i]-xo)**2)/(si**3))*np.exp(-(lline[i]-xo)**2/(2*si**2))
        
        M=np.empty([3,3])*0
        V=np.empty([3])*0
        M[0,0]=np.sum(beta**2)
        M[0,1]=np.sum(beta*gammap)
        M[0,2]=np.sum(beta*delta)
        M[1,0]=np.sum(beta*gammap)
        M[1,1]=np.sum(gammap**2)
        M[1,2]=np.sum(gammap*delta)
        M[2,0]=np.sum(beta*delta)
        M[2,1]=np.sum(gammap*delta)
        M[2,2]=np.sum(delta**2)
        V[0]=np.sum(beta*alfa)
        V[1]=np.sum(gammap*alfa)
        V[2]=np.sum(delta*alfa)
        
        Minv=np.linalg.pinv(M) 
        C=Minv.dot(V)
    
        A=A+C[0]
        xo=xo+C[1]
        si=si+C[2]
        for j in range(0, (np.size(lline)-1)):
            me1m[j]=fline[j]-(c-A*np.exp(-(lline[j]-xo)**2/(2*si**2)))
        somma=0
        for k in range(0, (np.size(lline)-1)):
            somma=(somma+(me1m[k])**2)
        me1park=me1
        me1=(np.sqrt(1/((np.size(lline)**2)-5)*somma))
        scme1=(np.abs(me1park-me1))
        contr=contr+1
    return  A,xo,si

#inputs
file_name_spectra_fit='comb_24.fits'     #spectra fits format in angstroems
filename_line_list_Fe='linelist_fe.dat'	  #line list Fe 
filename_line_list_other='linelist_other.dat'  #line list other elements

name_out=file_name_spectra_fit.split('.fits')[0]    #name for   outputs

#Insturments parameters 
Rs=0.25 #resolution
#A= 1    #FWHM
lint= 3
lyi=0.0
lys =1.55  
l, f = pyasl.read1dFitsSpec(file_name_spectra_fit)  #reading spectrum
minw=l[0]     #Min wavelenght 
maxw=l[-1]   #Max  wavelenght

####################Reading linelists##########################################
#### format linelists -- lambda col 0, numb element col 1,ep col 2,loggf col 3,damp col4
cond=True
while (cond==True):
    choice_list = input('What line-list want to measure: Fe (1) or others(2)? ').lower()

    if choice_list==str(1):
        linelistFE=pd.read_table(filename_line_list_Fe,delim_whitespace=True, names = ['lambda','el','ep', 'loggf', 'damp'],header=1)
        #lambda=0, el=1,ep=2,loggf=3,damp=4
        linelist_array=np.array(linelistFE)
        el_name='Fe'
        cond= False
    elif choice_list==str(2):
        linelist_other=pd.read_table(filename_line_list_other,delim_whitespace=True, names = ['lambda','el','ep', 'loggf', 'damp'],header=1)
        #lambda=0, el=1,ep=2,loggf=3,damp=4
        linelist_array=np.array(linelist_other)
        el_name='Other'
        cond= False    
    else:
        cond=True
#######################################################################################     
indice_range=np.where ((linelist_array[:,0]> minw) &  ((linelist_array[:,0])<maxw))
lambda_new=linelist_array[indice_range,0][0]
el_new=linelist_array[indice_range,1][0]
ep_new=linelist_array[indice_range,2][0]
loggf_new=linelist_array[indice_range,3][0]
damp_new=linelist_array[indice_range,4][0]
 
c=1  #continuum inital value
infoeqw=np.empty(np.size(lambda_new))*0
si_back=np.empty(np.size(lambda_new))*0  
cont_default=1       #Counter inital value
rng=0.5              #range capture default
k=0                  #Counter number of number on line from linelist
sav=0                #counter number of lines measured and ready to save  

######################loop to measure eqw lines from linelist and spectra#####################
while (k!=(np.size(lambda_new)-1)):
     
    boole, rng=rangee(0.3,4,rng,'Enter the range to capture (0.3 - 4), prev.:'+str(rng)+', press enter to keep') 
    boole2, c=rangee (-2,2,c,'Enter the continuum, prev.:'+str(c)+' press enter to keep')

    #######################
    lline_ind= np.where((l> (lambda_new[k]-rng/2)) &  (l<(lambda_new[k]+rng/2)))
    lline=l[lline_ind]
    fline=f[lline_ind]

    xgaussmin=lline[np.argmin(fline)]
    ygaussmin=fline[np.argmin(fline)]

    lxi=lambda_new[k]-lint  #lower limit  y axis (flux)
    lxs=lambda_new[k]+lint  #Upper limit  y axis (flux)
    ########################  
    sigline=(Rs/(2*1.1774))

    lline_ind = np.where((l> (lambda_new[k]-rng/2)) &  (l<(lambda_new[k]+rng/2)))
    lline=l[lline_ind]
    fline=f[lline_ind]

    lxi=lambda_new[k]-lint  #Lower limit  x axis (wavelength)
    lxs=lambda_new[k]+lint  #Upper limit  y axis (wavelength)

    ###capturing data using the range givenn######
    indice_range_2=np.where ((l> lxi) &  (l<lxs))
    ll=l[indice_range_2]
    fl=f[indice_range_2]
    lineavert_y=np.arange(lyi,lys+2,1)
    lineavert_x=np.empty(np.size(lineavert_y))*0+lambda_new[k]
    
    ##################################
    A,xo,si=gauss_mod(1)
    xgauss=np.arange(lxi,lxs,0.001)
    ygauss=c-A*np.exp(-(xgauss-xo)**2/(2*si**2))
    ##################################

    infoeqw[k]= np.sum(((c-ygauss)/c)*0.001)
    si_back[k]=si

    plt.ion()
    plt.figure(figsize=(18,6),dpi=70)
    plt.subplot2grid((1,2), (0,0), colspan=1, rowspan=1)
    a=plt.plot (l,f,color='blue')
    b=plt.plot (lineavert_x,lineavert_y,c='red')
    plt.plot(xgauss, ygauss,c='orange')
    plt.xlim(lxi,lxs)
    plt.ylim(lyi,lys)
    plt.xlabel ('$\lambda [\AA$]',fontsize = 13)
    plt.ylabel('Flux',fontsize = 13)
    text='Element: '+str(el_new[k]) +'   W: '+ str(lambda_new[k])+'[$\AA$]'
    text2='Contin.:' + str(c)
    text4=str(k+1)+'/'+ str(np.size(lambda_new)+1)
    textsav='Lines saved: ' + str(sav)
    plt.text(lxi,lys+0.1,text, fontsize=14)
    plt.text(lxs,lys+0.1,text2, fontsize=14,horizontalalignment='right')
    plt.text(lxi,lyi-0.15,text4, fontsize=14)
    plt.text(lxs,lyi-0.15,textsav, fontsize=14,horizontalalignment='right')
    plt.yticks(np.arange(lyi,lys, 0.1)) 

    #  subplot 1
    plt.subplot2grid((1,2), (0,1))
    text3='Range.:' + str(rng)
    plt.text(lxi+lint-rng,lys+0.1,text3, fontsize=14)
    plt.locator_params(axis='x', nbins=5)
    plt.locator_params(axis='y', nbins=5)
    plt.title('Range Capture ')
    plt.scatter (l,f,color='blue')
    plt.scatter(lline,fline,c='red')
    plt.plot (lineavert_x,lineavert_y,c='red',)
    plt.plot(xgauss, ygauss,c='orange')
    plt.xlim(lxi+lint-rng,lxs-lint+rng)
    plt.ylim(lyi,lys)
    plt.xlabel ('$\lambda [\AA$]',fontsize = 13)
    plt.yticks(np.arange(lyi,lys, 0.1)) 
    plt.show()

    loop = yes_no(
        'Do you want  to repeat  the process for this line? (yes/no)')
        
    if loop == False:
        safe=  yes_no('Do you want  save this line? (yes/no)')
        if safe== True:
            sav=sav+1
        if safe== False:
            infoeqw[k]= 9999  
        k=k+1
    plt.close() 

######################## Save measurements for Fe or other #############################################################
index_val=np.where ((infoeqw<2) & (infoeqw>0.00001) & (si_back<100))
formato="%4.3f %2.1f %2.3f %2.3f  %2.3f %1.1f %3.1f" 
infoeqwNew=infoeqw*1000
np.savetxt(name_out+'.'+el_name+'.eqw.dat',np.c_[lambda_new[index_val],el_new[index_val],ep_new[index_val],loggf_new[index_val],damp_new[index_val],infoeqwNew[index_val],si_back[index_val]],fmt=formato)  #output file
#####################################################################################
lambda_new_s=lambda_new[index_val]
ep_new_s=ep_new[index_val]
el_new_s=el_new[index_val]
loggf_new_s=loggf_new[index_val]
damp_new_s=damp_new[index_val]
infoeqwNew_s=infoeqwNew[index_val]
si_back_s=si_back[index_val]
#####################################################################################
exit_q=yes_no('Do you want follow the next step (Sigma Rejection)?  (y/n)')
if exit_q==False:
    sys.exit()
########################################Sigma Rejection#############################################
loop_sigma=True
while (loop_sigma==True):
    n = 3  # degree of polynomial
    p, C_p = np.polyfit(lambda_new[index_val], si_back[index_val], n, cov=True)  # C_z is estimated covariance matrix
   
    t=lambda_new[index_val]
    t=np.array(t)
    ss=si_back[index_val]

    # Matrix with rows 1, t, t**2, ...:
    TT = np.vstack([t**(n-i) for i in range(n+1)]).T
    yi = np.dot(TT, p)  # matrix multiplication calculates the polynomial values
    C_yi = np.dot(TT, np.dot(C_p, TT.T)) # C_y = TT*C_z*TT.T
    sig_yi = np.sqrt(np.diag(C_yi))  # Standard deviations are sqrt of diagonal
    boole3, nsig=rangee(0,1000,1,'sigma rejection (positive number)')
    # Do the plotting:

    plt.figure(figsize=(9, 6), dpi=100)
    plt.title("Fit for Polynomial (degree {}) with $\pm {} \sigma$-interval".format(n,nsig))
    plt.fill_between(t, yi+nsig*sig_yi, yi-nsig*sig_yi, alpha=.25)
    plt.plot(t, yi,'-')
    plt.axis('tight')
    plt.xlim(5000,5300)
        
    poly_m=0
    for i in range (0,n+1,1):
        #poly_m=poly_m+p[i]*lambda_new[index_val]**(n-i)
        poly_m=poly_m+p[i]*t**(n-i)

    diffsip=ss-poly_m

    indicesig=np.where((np.abs(diffsip)<nsig*sig_yi))

    lambnew2=t[indicesig]
    sipnew2=ss[indicesig]

    plt.scatter(lambda_new[index_val], si_back[index_val],marker='x',c='orange')
    plt.scatter(lambnew2,sipnew2,c='green',marker='x')
    plt.xlabel ('$\lambda [\AA$]',fontsize = 13)
    plt.ylabel('$\sigma$',fontsize = 13)
    plt.xlim(np.min(lambnew2)-5,np.max(lambnew2)+5)
    plt.ylim(np.min(sipnew2)-0.03,np.max(sipnew2)+0.03)
    plt.show()

    loop_sigma = yes_no('Do you want  change the sigma? (yes/no) ')
    plt.close() 

ep_new2=ep_new_s[indicesig]
loggf_new2=loggf_new_s[indicesig]
damp_new2=damp_new_s[indicesig]
infoeqwNew2=infoeqwNew_s[indicesig]
el_new2=el_new_s[indicesig]
lambnew2=lambda_new[indicesig]
si_back_s=si_back[indicesig]

if choice_list==str(1):
    index_sel=np.where ((el_new2==26.0) | (el_new2==26.1))  
elif choice_list==str(2):
    index_sel=np.where ((el_new2!=26.0) & (el_new2!=26.1))

lambnew3=lambnew2[index_sel]
sipnew3=sipnew2[index_sel]
ep_new3=ep_new2[index_sel]
loggf_new3=loggf_new2[index_sel]
damp_new3=damp_new2[index_sel]
infoeqwNew3=infoeqwNew2[index_sel]
el_new3=el_new2[index_sel]

ind_sort=np.argsort(el_new3)  #sort for element

lambnew3=lambnew3[ind_sort]
sipnew3=sipnew3[ind_sort]
ep_new3=ep_new3[ind_sort]
loggf_new3=loggf_new3[ind_sort]
damp_new3=damp_new3[ind_sort]
infoeqwNew3=infoeqwNew3[ind_sort]
el_new3=el_new3[ind_sort]
disp3=el_new3*0

########################   Save measurements with MOOG format ##############################################################
formato="%4.3f %2.1f %2.3f %2.3f  %2.3f %1.1f %3.1f" 
np.savetxt(name_out+'.'+el_name+'.MOOG.dat',np.c_[lambnew3,el_new3,ep_new3,loggf_new3,damp_new3,disp3,infoeqwNew3],fmt=formato)  #output file
######################################################################################
######################################################################################
#####################################END##############################################





