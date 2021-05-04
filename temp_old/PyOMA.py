# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 15:49:48 2021

@author: dagpa
"""

import numpy as np
import pandas as pd
from scipy import signal
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
import matplotlib.patches as patches
import seaborn as sns
import mplcursors

# =============================================================================
# =============================================================================


def MaC(Fi1,Fi2):
    '''This function return the Modal Assurance Criterion (MAC) for two modal 
    shape vectors.'''
    return np.abs(Fi1.conj().T@Fi2)**2 / ((Fi1.conj().T@Fi1)*(Fi2.conj().T@Fi2))


# =============================================================================
# =============================================================================
    

def SSIdatStaDiag(data, fs, br, ordmax=None, lim=(0.01,0.05,0.02,0.1), method='1'):
    '''This function return the Stabilization Diagram (Plot) for the given
    data calculated according the data driven SSI algorithm.'''
    
    ndat=data.shape[0] # Number of data points
    nch=data.shape[1] # Number of channel
    
    # If the maximum order is not given (default) it is set as the maximum
    # allowable model order which is: number of block rows * number of channels
    if ordmax == None:
        ordmax = br*nch
        
    freq_max = fs/2 # Nyquist Frequency
    
    # unpack the limits used for the construction of the Stab Diag
    _lim_f, _lim_s, _lim_ms, _lim_s1 = lim[0], lim[1], lim[2], lim[3]

    Yy=np.array(data.transpose()) # Matrice dati trasposta
    
    j=ndat-2*br+1; # DIMENSIONE MATRICE DI HANKEL

    H=np.zeros((nch*2*br,j)) # PREALLOCA LA MATRICE DI HANKEL
    for _k in range(0,2*br):
     	H[_k*nch:((_k+1)*nch),:]=(1/j**0.5)*Yy[:,_k:_k+j] # CALCOLO MATRICE DI HANKEL
    
    # FATTORIZZAZIONE LQ DELLA MATRICE DI HANKEL
    Q , L = np.linalg.qr(H.T)
    L = L.T
    Q = Q.T
    
    _a = nch*br
    _b = nch
    
    L21 = L[_a:_a+_b,:_a]
    L22 = L[_a:_a+_b,_a:_a+_b]
    L31 = L[_a+_b:,:_a]
    L32 = L[_a+_b:,_a:_a+_b]
    
    Q1 = Q[:_a,:]
    Q2 = Q[_a:_a+_b,:]
    
    P_i = np.vstack((L21,L31)) @ Q1 # Projection Matrix P_i
    P_im1 = np.hstack((L31,L32)) @ np.vstack((Q1, Q2)) # Projection P_(i-1)
    Y_i = np.hstack((L21,L22)) @ np.vstack((Q1, Q2)) # Output sequence
    
    # SINGULAR VALUE DECOMPOSITION
    U1, S1, V1_t = np.linalg.svd(P_i,full_matrices=False)
    S1 = np.diag(S1)
    S1rad=np.sqrt(S1)
    
    # Ciclo per ordine del sistema crescente
    Fr=np.full((ordmax, ordmax+1), np.nan) # inizializzo la matrice che conterrà le frequenze
    Fr_lab=np.full((ordmax, ordmax+1), np.nan)  # inizializzo la matrice che conterrà le labels per i poli(frequenze) da stampare nel diagramma di stabilizzazione
    Sm=np.full((ordmax, ordmax+1), np.nan) # inizializzo la matrice che conterrà gli smorzamenti
    Ms = []  # inizializzo la matrice (3D) che conterrà le forme modali
    for z in range(0, int((ordmax)/2+1)):
        Ms.append(np.zeros((nch, z*(2))))

    for _ind in range(0, ordmax+1, 2):
        # Inizializzo le matrici che mi servono
        S11 = np.zeros((_ind, _ind)) # Inizializzo
        U11 = np.zeros((br*nch, _ind)) # Inizializzo
        V11 = np.zeros((_ind, br*nch)) # Inizializzo
        O_1 = np.zeros((br*nch - nch, _ind)) # Inizializzo
        O_2 = np.zeros((br*nch - nch, _ind)) # Inizializzo
        
        # Estrazione di sottomatrici per ordine crescente del sistema
        S11[:_ind, 0:_ind] = S1rad[0:_ind,0:_ind] # ESTRAZIONE DI UNA SOTTOMATRICE DEI SINGULAR VALUES DA ORDINE MIN A ORD MAX
        U11[:br*nch, 0:_ind] = U1[0:br*nch,0:_ind] # ESTRAZIONE DI UNA SOTTOMATRICE DEI LEFT SINGULAR VECTORS DA ORDINE MIN A ORD MAX
        V11[0:_ind, :br*nch] = V1_t[0:_ind,0:br*nch] # 

        O = U11 @ S11 # Observability matrix
        S = np.linalg.pinv(O) @ P_i # Kalman filter state sequence

        O_1[:,:] = O[:O.shape[0] - nch,:]
        O_2[:,:] = O[nch:,:]

        # STIMA DELLA MATRICE DINAMICA A TEMPO DISCRETO
        if method == '2': # Method 2 
            A = np.linalg.pinv(O_1) @ O_2 
            C = O[:nch,:]     
            # Ci sarebbero da calcolare le matrici G e R0 

        else:  # METODO 1
            Sp1 = np.linalg.pinv(O_1) @ P_im1 # kalman state sequence S_(i+1)
        
            AC = np.vstack((Sp1,Y_i)) @ np.linalg.pinv(S) 
            A = AC[:Sp1.shape[0]]
            C = AC[Sp1.shape[0]:]
            # Ci sarebbero da calcolare le matrici G e R0 

      
        [_AuVal, _AuVett] = np.linalg.eig(A) # CALCOLO AUTOVALORI ED AUTOVETTORI
        Lambda =(np.log(_AuVal))*fs # CALCOLO Lambda
        fr = abs(Lambda)/(2*np.pi) # CALCOLO FRQUENZE
        smorz = -((np.real(Lambda))/(abs(Lambda))) # CALCOLO SMORZAMENTO
        
        # calcolo la matrice C per definizione dei modi
        # le deformate modali M in forma complessa sono calcolate attraverso
        # l'espressione M=CL dove L è la matrice le cui colonne sono gli autovettori
        # di A (_AuVett)

        Mcomp = C@_AuVett
        Mreal =np.real(C@_AuVett)
        
        Fr[:len(fr),_ind] = fr # SALVA LE FREQUENZE    
        Sm[:len(fr),_ind] = smorz # SALVA gli smorzamenti
        Ms[int(_ind/2)] = Mcomp # SALVA le forme modali
# =============================================================================
# Controllo stabilità dei poli 
        for idx, (_freq, _smor) in enumerate(zip(fr,smorz)):
            if _ind == 0 or _ind == 2: # alla prima iterazione/primo ordine sono tutti nuovi poli
                Fr_lab[:len(fr),_ind] = 0 # 0 è l'etichetta per i nuovi poli e poli instabili

            elif np.isnan(_freq) == True:
                _freq=0
                Fr[idx,_ind] = 0
    
            else:
                # Trovo l'indice del polo che minimizza la differenza alla iterazione/ordine n-1
                ind2 = np.nanargmin(abs(_freq - Fr[:,(_ind) - 2]) 
                                    - min(abs(_freq - Fr[:,(_ind) - 2])))
                        
                Fi_n = Mcomp[:, idx] # forma modale all'iterazione = _iteraz (attuale)
                Fi_nmeno1 = Ms[int(_ind/2-1)][:,ind2] # forma modale all'iterazione = _iteraz-1 (precedente)
                
                # aMAC = np.abs(Fi_n@Fi_nmeno1)**2 / ((Fi_n@Fi_n)*(Fi_nmeno1@Fi_nmeno1)) # autoMAC
                aMAC =  MaC(Fi_n.real,Fi_nmeno1.real)
                
                if min(abs(_freq - Fr[:,(_ind) - 2])/_freq) < _lim_f: # CONDIZIONE 1 SULLA FREQUENZA
                    
        
                    if (_smor - Sm[ind2, (_ind) - 2])/_smor < _lim_s: # CONDIZIONE 2 SULLO SMORZAMENTO
                        if 1 - aMAC < _lim_ms: # CONDIZIONE 3 SULLE FORME MODALI
                            Fr_lab[idx,_ind] = 4 # Se il polo è stabile sia per smorzamento che per forma modale (viene classificato come stabile/verde)
                        else:
                            Fr_lab[idx,_ind] = 2 # Stabile per smorzamento ma non per forma modale
        
                    elif 1 - aMAC < _lim_ms: # CONDIZIONE 3 SULLE FORME MODALI
                        Fr_lab[idx,_ind] = 3 # Stabile per forma modale ma non per smorzamento
        
                    else:
                        Fr_lab[idx,_ind] = 1 # Stabile solo per frequenza
                else:
                    Fr_lab[idx,_ind] = 0  # Nuovo polo o polo instabile

# ============================================================================= 
# PLOT DIAGRAMMA DI STABILIZZAZZIONE 
# =============================================================================
    _x = Fr.flatten(order='f')
    _y = np.array([_i//len(Fr) for _i in range(len(_x))])
    _l = Fr_lab.flatten(order='f')
    _d = Sm.flatten(order='f')
    
    _df = pd.DataFrame(dict(Frequency=_x, Order=_y, Label=_l, Damp=_d))
    _df1 = _df.copy() # DataFrame che mi serve per plot diagramma stabilizzazione
    
# =============================================================================
# Qui creo un dataframe "ridotto" (senza nan) dove sono salvate: frequenze,
# smorzamenti e indici della forma modale associata al polo (order, emme)
    df2 = _df.copy()
    df2 = df2.dropna()
    emme = []
    for effe,order in zip(df2.Frequency,df2.Order):
        emme.append(np.nanargmin(abs(effe - Fr[:, order]))) # trovo l'indice 
    emme = np.array(emme)
    df2['Emme'] = emme
    df2 = df2.drop_duplicates()
    df2 = df2.drop(columns='Label')
# =============================================================================

    _df1.Frequency = _df1.Frequency.where(_df.Damp < _lim_s1) # qui rimuovo tutti i poli che hanno smorzamento maggiore di lim_s1
    _df1.Frequency = _df1.Frequency.where(_df.Damp > 0) # qui rimuovo tutti i poli che hanno smorzamento negativo (< 0)

    _colors = {0:'Red', 1:'darkorange', 2:'gold', 3:'yellow', 4:'Green'} # assegno i colori alle etichette dei poli
    
    fig1, ax1 = plt.subplots()
    ax1 = sns.scatterplot(_df1['Frequency'], _df1['Order'], hue=_df1['Label'], palette=_colors)
    
    ax1.set_xlim(left=0, right=freq_max)
    ax1.set_ylim(bottom=0, top=ordmax)
    ax1.xaxis.set_major_locator(MultipleLocator(freq_max/10))
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax1.xaxis.set_minor_locator(MultipleLocator(freq_max/100))
    ax1.set_title('''{0} - shift: {1}'''.format('Stabilization Diagram', br))
    ax1.set_xlabel('Frequency [Hz]')
    mplcursors.cursor()
    plt.show()

    Results={}
    Results['Data'] = {'Data': data}
    Results['Data']['Samp. Freq.'] = fs
    Results['All Poles'] = _df1
    Results['Reduced Poles'] = df2
    Results['Modes'] = Ms
    
    return fig1, Results


# =============================================================================
# =============================================================================


def SSIcovStaDiag(data, fs, br, ordmax=None, lim=(0.01,0.05,0.02,0.1), method='1'):
    ndat=data.shape[0] #NUMERO DI DATI CAMPIONATI
    nch=data.shape[1] #NUMERO DI CANALI ACQUISITI
    if ordmax == None:
        ordmax = br*nch
    freq_max = fs/2 # Frequenza di Nyquist

    _lim_f, _lim_s, _lim_ms, _lim_s1 = lim[0], lim[1], lim[2], lim[3]

    Yy=np.array(data.transpose())
# =============================================================================
# METODO LIBRO (Rainieri & Fabbrocino)
# =============================================================================
# Mi calcolo tutti i R[i] (con i da 0 a 2*br)
    R_is = np.array([1/(ndat - _s)*(Yy[:, : ndat - _s]@Yy[:, _s:].T) for _s in range(br*2+1)]) 
# Mi costruisco la matrice di Toepliz
    Tb = np.vstack([np.hstack([R_is[_o,:,:] for _o in range(br+_l, _l,-1)]) for _l in range(br)])
    
# One-lag shifted matrice di Toeplitz per calcolo di A secondo "NExT-ERA"
    Tb2 = np.vstack([np.hstack([R_is[_o,:,:] for _o in range(br+_l, _l,-1)]) for _l in range(1,br+1)])

# =============================================================================
# SINGULAR VALUE DECOMPOSITION
# =============================================================================
    U1, S1, V1_t = np.linalg.svd(Tb)
    S1 = np.diag(S1)
    S1rad=np.sqrt(S1)
# =============================================================================
# Ciclo per ordine del sistema crescente
# =============================================================================
    Fr=np.full((ordmax, ordmax+1), np.nan) # inizializzo la matrice che conterrà le frequenze
    Fr_lab=np.full((ordmax, ordmax+1), np.nan)  # inizializzo la matrice che conterrà le labels per i poli(frequenze) da stampare nel diagramma di stabilizzazione
    Sm=np.full((ordmax, ordmax+1), np.nan) # inizializzo la matrice che conterrà gli smorzamenti
    Ms = []  # inizializzo la matrice (3D) che conterrà le forme modali
    for z in range(0, int((ordmax)/2+1)):
        Ms.append(np.zeros((nch, z*(2))))
# =============================================================================
    for _ind in range(0, ordmax+1, 2):
# Inizializzo le matrici che mi servono
        S11 = np.zeros((_ind, _ind)) # Inizializzo
        U11 = np.zeros((br*nch, _ind)) # Inizializzo
        V11 = np.zeros((_ind, br*nch)) # Inizializzo
        O_1 = np.zeros((br*nch - nch, _ind)) # Inizializzo
        O_2 = np.zeros((br*nch - nch, _ind)) # Inizializzo
# =============================================================================
# Estrazione di sottomatrici per ordine crescente del sistema
        S11[:_ind, 0:_ind] = S1rad[0:_ind,0:_ind] # ESTRAZIONE DI UNA SOTTOMATRICE DEI SINGULAR VALUES DA ORDINE MIN A ORD MAX
        U11[:br*nch, 0:_ind] = U1[0:br*nch,0:_ind] # ESTRAZIONE DI UNA SOTTOMATRICE DEI LEFT SINGULAR VECTORS DA ORDINE MIN A ORD MAX
        V11[0:_ind, :br*nch] = V1_t[0:_ind,0:br*nch] # 
# =============================================================================
        O = U11 @ S11 # CALCOLO MATRICE DI OSSERVABILITA
        _GAM = S11 @ V11 # CALCOLO MATRICE DI CONTROLLABILITA
        
        O_1[:,:] = O[:O.shape[0] - nch,:]
        O_2[:,:] = O[nch:,:]
# =============================================================================
# STIMA DELLA MATRICE DINAMICA A TEMPO DISCRETO
        if method == '2':
            A = np.linalg.inv(S11)@U11.T@Tb2@V11.T@np.linalg.inv(S11) # METODO "NExT-ERA"
        else:
            A = np.linalg.pinv(O_1)@O_2 # METODO 2 (BALANCED_REALIZATION)
        
        [_AuVal, _AuVett] = np.linalg.eig(A) # CALCOLO AUTOVALORI ED AUTOVETTORI
        Lambda =(np.log(_AuVal))*fs # CALCOLO Lambda
        fr = abs(Lambda)/(2*np.pi) # CALCOLO FRQUENZE
        smorz = -((np.real(Lambda))/(abs(Lambda))) # CALCOLO SMORZAMENTO
        
        # calcolo la matrice C per definizione dei modi
        # le deformate modali M in forma complessa sono calcolate attraverso
        # l'espressione M=CL dove L è la matrice le cui colonne sono gli autovettori
        # di A (_AuVett)
        C = O[:nch,:]
        
        Mcomp = C@_AuVett
        Mreal =np.real(C@_AuVett)
        
        Fr[:len(fr),_ind] = fr # SALVA LE FREQUENZE    
        Sm[:len(fr),_ind] = smorz # SALVA gli smorzamenti
        Ms[int(_ind/2)] = Mcomp # SALVA le forme modali
# =============================================================================
# Controllo stabilità dei poli 
        for idx, (_freq, _smor) in enumerate(zip(fr,smorz)):
            if _ind == 0 or _ind == 2: # alla prima iterazione/primo ordine sono tutti nuovi poli
                Fr_lab[:len(fr),_ind] = 0 # 0 è l'etichetta per i nuovi poli e poli instabili

            elif np.isnan(_freq) == True:
                _freq=0
                Fr[idx,_ind] = 0
    
            else:
                # Trovo l'indice del polo che minimizza la differenza alla iterazione/ordine n-1
                ind2 = np.nanargmin(abs(_freq - Fr[:,(_ind) - 2]) 
                                    - min(abs(_freq - Fr[:,(_ind) - 2])))
                        
                Fi_n = Mcomp[:, idx] # forma modale all'iterazione = _iteraz (attuale)
                Fi_nmeno1 = Ms[int(_ind/2-1)][:,ind2] # forma modale all'iterazione = _iteraz-1 (precedente)
                
                # aMAC = np.abs(Fi_n@Fi_nmeno1)**2 / ((Fi_n@Fi_n)*(Fi_nmeno1@Fi_nmeno1)) # autoMAC
                aMAC =  MaC(Fi_n.real,Fi_nmeno1.real)
                    
                if min(abs(_freq - Fr[:,(_ind) - 2])/_freq) < _lim_f: # CONDIZIONE 1 SULLA FREQUENZA
                    
        
                    if (_smor - Sm[ind2, (_ind) - 2])/_smor < _lim_s: # CONDIZIONE 2 SULLO SMORZAMENTO
                        if 1 - aMAC < _lim_ms: # CONDIZIONE 3 SULLE FORME MODALI
                            Fr_lab[idx,_ind] = 4 # Se il polo è stabile sia per smorzamento che per forma modale (viene classificato come stabile/verde)
                        else:
                            Fr_lab[idx,_ind] = 2 # Stabile per smorzamento ma non per forma modale
        
                    elif 1 - aMAC < _lim_ms: # CONDIZIONE 3 SULLE FORME MODALI
                        Fr_lab[idx,_ind] = 3 # Stabile per forma modale ma non per smorzamento
        
                    else:
                        Fr_lab[idx,_ind] = 1 # Stabile solo per frequenza
                else:
                    Fr_lab[idx,_ind] = 0  # Nuovo polo o polo instabile

# ============================================================================= 
# PLOT DIAGRAMMA DI STABILIZZAZZIONE 
# =============================================================================
    _x = Fr.flatten(order='f')
    _y = np.array([_i//len(Fr) for _i in range(len(_x))])
    _l = Fr_lab.flatten(order='f')
    _d = Sm.flatten(order='f')
    
    _df = pd.DataFrame(dict(Frequency=_x, Order=_y, Label=_l, Damp=_d))
    _df1 = _df.copy() # DataFrame che mi serve per plot diagramma stabilizzazione
    
# =============================================================================
# Qui creo un dataframe "ridotto" (senza nan) dove sono salvate: frequenze,
# smorzamenti e indici della forma modale associata al polo (order, emme)
    df2 = _df.copy()
    df2 = df2.dropna()
    emme = []
    for effe,order in zip(df2.Frequency,df2.Order):
        emme.append(np.nanargmin(abs(effe - Fr[:, order]))) # trovo l'indice 
    emme = np.array(emme)
    df2['Emme'] = emme
    df2 = df2.drop_duplicates()
    df2 = df2.drop(columns='Label')
# =============================================================================

    _df1.Frequency = _df1.Frequency.where(_df.Damp < _lim_s1) # qui rimuovo tutti i poli che hanno smorzamento maggiore di lim_s1
    _df1.Frequency = _df1.Frequency.where(_df.Damp > 0) # qui rimuovo tutti i poli che hanno smorzamento negativo (< 0)

    _colors = {0:'Red', 1:'darkorange', 2:'gold', 3:'yellow', 4:'Green'} # assegno i colori alle etichette dei poli
    
    fig1, ax1 = plt.subplots()
    ax1 = sns.scatterplot(_df1['Frequency'], _df1['Order'], hue=_df1['Label'], palette=_colors)
    
    ax1.set_xlim(left=0, right=freq_max)
    ax1.set_ylim(bottom=0, top=ordmax)
    ax1.xaxis.set_major_locator(MultipleLocator(freq_max/10))
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax1.xaxis.set_minor_locator(MultipleLocator(freq_max/100))
    ax1.set_title('''{0} - shift: {1}'''.format('Stabilization Diagram', br))
    ax1.set_xlabel('Frequency [Hz]')
    mplcursors.cursor()
    plt.show()
    
    Results={}
    Results['Data'] = {'Data': data}
    Results['Data']['Samp. Freq.'] = fs
    Results['All Poles'] = _df1
    Results['Reduced Poles'] = df2
    Results['Modes'] = Ms
    
    
    return fig1, Results


# =============================================================================
# =============================================================================


def SSIModEX(FreQ, Results, deltaf=0.05, aMaClim=0.9):
    
    df2 = Results['Reduced Poles']
    Ms = Results['Modes'] 
    
    Freq = []
    Damp = []
    Fi = []
    for _x in FreQ:
        xmeno1, xpiu1 = _x-deltaf, _x+deltaf
        
        df3 = df2.where((df2.Frequency < xpiu1) & (df2.Frequency > xmeno1))
        df3 = df3.dropna()
        
        npoli = len(df3['Frequency'].values)
        
        AutoMacche = np.zeros((npoli, npoli))
        for _b in range(npoli):
            _zuno = int(df3['Order'].values[_b]/2)
            fiuno = Ms[_zuno][:,int(df3['Emme'].values[_b])]
            for _k in range(len(df3['Frequency'].values)):
                _zdue = int(df3['Order'].values[_k]/2)
                fidue = Ms[_zdue][:,int(df3['Emme'].values[_k])]
                    
                AutoMacche[_b, _k] = MaC(fiuno.real,fidue.real)
        
        SAmaC = np.sum(AutoMacche, axis=1)
        idxmax = np.argmax(SAmaC)
        MSrefidx1 = int(df3['Order'].values[idxmax]/2)
        MSrefidx2 = int(df3['Emme'].values[idxmax])
        firef = Ms[MSrefidx1][:,MSrefidx2]
        
        AMaC = AutoMacche[idxmax]
        df3['AMaC'] = AMaC
        df3 = df3.where(df3.AMaC > aMaClim)
        df3 = df3.dropna()
        FrMean = df3.Frequency.mean()
        FrStd = df3.Frequency.std()
        DampMean = df3.Damp.mean()
        DampStd = df3.Damp.std()
        
        FI1 = np.array([Ms[int(df3['Order'].values[_i]/2)][:,int(df3['Emme'].values[_i])] for _i in range(len(df3))])
        FI1real = FI1.real
        sgn = []
        for _l in FI1real[0]:
            if _l>0:
               sgn.append(1)
            else:
                sgn.append(-1)
        Fimean = np.mean(np.abs(FI1), axis=0)*sgn
        _idx = np.argmax(abs(Fimean))
        Fi_mean_norm = Fimean/Fimean[_idx]   
        
        Freq.append(FrMean)
        Damp.append(DampMean)
        Fi.append(Fi_mean_norm)

    Freq = np.array(Freq)
    Damp = np.array(Damp)
    Fi = np.array(Fi)

    Results={}
    Results['Frequencies'] = Freq
    Results['Damping'] = Damp
    Results['Mode Shapes'] = Fi.T
        
    return Results


# =============================================================================
# =============================================================================
    

def FDDsvp(data, fs, df=0.01, pov=0.5, window='hann'):
    
    ndat=data.shape[0] #NUMERO DI DATI CAMPIONATI
    nch=data.shape[1] #NUMERO DI CANALI ACQUISITI
    freq_max = fs/2 # Frequenza di Nyquist
    nxseg = fs/df # numero di punti per segmenti (su cui mediare)
#    nseg = ndat // nxseg # numero di segmenti (su cui mediare)
    noverlap = nxseg // (1/pov) # Numero di punti che si sovrappongono tra i segmenti (Default 50%)
    
    PSD_matr = np.zeros((nch, nch, int((nxseg)/2+1)), dtype=complex) # Inizializzo la SD matrix
    S_val = np.zeros((nch, nch, int((nxseg)/2+1))) # Inizializzo la matrice dove salverò i Singular Values
    S_vec = np.zeros((nch, nch, int((nxseg)/2+1)), dtype=complex) # Inizializzo la matrice dove salverò i Singular Vectors
    
    # loop dove mi calcolo le Auto e Cross-Spectral Density
    # (si passa al dominio della frequenza)
    for _i in range(0, nch):
        for _j in range(0, nch):
            _f, _Pxy = signal.csd(data[:, _i],data[:, _j], fs=fs, nperseg=nxseg, noverlap=noverlap, window=window)
            PSD_matr[_i, _j, :] = _Pxy
            
    # loop dove mi calcolo i singular value      
    for _i in range(np.shape(PSD_matr)[2]):
        U1, S1, _V1_t = np.linalg.svd(PSD_matr[:,:,_i])
        U1_1=np.transpose(U1) 
        S1 = np.diag(S1)
        S1rad=np.sqrt(S1)
        S_val[:,:,_i] = S1rad
        S_vec[:,:,_i] = U1_1
    
    # Plot dei singular values (in scala logaritmica)
    fig, ax = plt.subplots()
    for _i in range(nch):
    #    ax.semilogy(_f, S_val[_i, _i]) # scala log
        ax.plot(_f[:], 10*np.log10(S_val[_i, _i])) # decibel
    ax.grid()
    ax.set_xlim(left=0, right=freq_max)
    ax.xaxis.set_major_locator(MultipleLocator(freq_max/10))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.xaxis.set_minor_locator(MultipleLocator(freq_max/100))
    ax.set_title("Singular values plot - (Freq. res. ={0})".format(df))
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel(r'dB $[g^2/Hz]$')    
    # ax.set_ylabel(r'dB $\left[\frac{\left(\frac{m}{s^2}\right)^2}{Hz}\right]$')    
    mplcursors.cursor()
    
    Results={}
    Results['Data'] = {'Data': data}
    Results['Data']['Samp. Freq.'] = fs
    Results['Data']['Freq. Resol.'] = df
    Results['Singular Values'] = S_val
    Results['Singular Vectors'] = S_vec
    Results['PSD Matrix'] = PSD_matr
    
    return fig, Results


# =============================================================================
# =============================================================================
    

def FDDmodEX(FreQ, Results, ndf=2):
#    data = Results['Data']['Data']
    fs = Results['Data']['Samp. Freq.']
    df = Results['Data']['Freq. Resol.']
    S_val = Results['Singular Values']
    S_vec = Results['Singular Vectors']
    deltaf=ndf*df
#    ndat=data.shape[0] #NUMERO DI DATI CAMPIONATI
#    nch=data.shape[1] #NUMERO DI CANALI ACQUISITI
    freq_max = fs/2 # Frequenza di Nyquist
#    nxseg = fs/df # numero di punti per segmenti (su cui mediare)
#    nseg = ndat // nxseg # numero di segmenti (su cui mediare)
    f = np.linspace(0, int(freq_max), int(freq_max*(1/df)+1)) # vettore delle linee spettrali
 
    Freq = []
    index = []
    Fi = []

    for _x in FreQ:
#        idx = np.argmin(abs(f-_x)) 
        lim = (_x - deltaf, _x + deltaf) # banda di frequenza in cui cercare il picco
        idxlim = (np.argmin(abs(f-lim[0])), np.argmin(abs(f-lim[1])))
        # vettore dei rapporti tra primo e secondo valore singolare all interno della banda di ricerca
        diffS1S2 = S_val[0,0,idxlim[0]:idxlim[1]]/S_val[1,1,idxlim[0]:idxlim[1]]
        maxDiffS1S2 = np.max(diffS1S2)
        idx1 = np.argmin(abs(diffS1S2 - maxDiffS1S2))
        idxfin = idxlim[0] + idx1
        # =============================================================================
        # Estrazione Forma Modale
        fr_FDD = f[idxfin]
        fi_FDD = S_vec[0,:,idxfin]
        idx3 = np.argmax(abs(fi_FDD))
        fi_FDDn = fi_FDD/fi_FDD[idx3] # Forma Modale (Normalizzata spostamento unitario)
        fiFDDn = np.array(fi_FDDn)
        
        Freq.append(fr_FDD)
        Fi.append(fiFDDn)
        index.append(idxfin)
        
    Freq = np.array(Freq)
    Fi = np.array(Fi)
    index = np.array(index)   
        
    Results={}
    Results['Frequencies'] = Freq
    Results['Mode Shapes'] = Fi.T
    Results['Freq. index'] = index

    return Results


# =============================================================================
# =============================================================================


def EFDDmodEX(FreQ, Results, ndf=2, MAClim=0.8, sppk=1, npmax=15, method='FSDD'):
    data = Results['Data']['Data']
    fs = Results['Data']['Samp. Freq.']
    df = Results['Data']['Freq. Resol.']
    S_val = Results['Singular Values']
    S_vec = Results['Singular Vectors']
    PSD_matr = Results['PSD Matrix']*2
    
    Res = FDDmodEX(FreQ, Results, ndf=ndf)
    Freq, Fi, index = Res['Frequencies'], Res['Mode Shapes'], Res['Freq. index']
    
#    ndat=data.shape[0] #NUMERO DI DATI CAMPIONATI
    nch=data.shape[1] #NUMERO DI CANALI ACQUISITI
    freq_max = fs/2 # Frequenza di Nyquist
#    nxseg = fs/df # numero di punti per segmenti (su cui mediare)
#    nseg = ndat // nxseg # numero di segmenti (su cui mediare)
#    dur=ndat/fs # [sec] Durata acquisizione
    tlag = 1/df # Durata 
    Nf = freq_max/df+1 # numero di linee spettrali
    f = np.linspace(0, int(freq_max), int(Nf)) # vettore con tutte le linee spettrali
    
    nIFFT = (int(Nf))*20 # num punti per trasformata inversa (se > di len(_f), viene zero-paddata)
    
    # Normalizzazione (spostamento max = 1) di tutte le forme modali
    S_vec_n = S_vec.copy()
    for _n in range(nch):
        Fis = S_vec[_n,:,:]
        idxFimax = np.array([np.argmax(abs(Fis[:,_i])) for _i in range(int(Nf))])
        S_vec_n[_n,:,:] = np.array([Fis[:,_j]/Fis[idxFimax[_j],_j] for _j in range(int(Nf))]).T
    
    Freq_E = []
    Fi_E = []
    Damp_E = []
    
    for _l,_f in enumerate(Freq):
        # Pre-allocazione
        _fi=Fi[:,_l]
        SDOFsval = np.zeros(int(Nf))
        SDOFsval1 = np.zeros(int(Nf),dtype=complex)
        SDOFsvec = np.zeros((nch,int(Nf)),dtype=complex)
                
        _p = 0
        while MaC(_fi, S_vec_n[0,:, index[_l] + _p]) > MAClim:
                SDOFsval[(index[_l] + _p)] = S_val[0,0, index[_l] + _p] # si aggiungono i relativi Singular values alla SFOF bell
                SDOFsval1[(index[_l] + _p)] = _fi.conj().T@PSD_matr[:,:, index[_l] + _p]@_fi # si aggiungono i relativi Singular values alla SFOF bell
                SDOFsvec[:,(index[_l] + _p)] = S_vec_n[0,:, index[_l] + _p] # e anche i relativi Singular vectors alla SFOF bell
                _p +=1
        while MaC(_fi, S_vec_n[1,:, index[_l] + _p]) > MAClim:
                SDOFsval[(index[_l] + _p)] = S_val[1,1, index[_l] + _p] # si aggiungono i relativi Singular values alla SFOF bell
                SDOFsval1[(index[_l] + _p)] = _fi.conj().T@PSD_matr[:,:, index[_l] + _p]@_fi # si aggiungono i relativi Singular values alla SFOF bell
                SDOFsvec[:,(index[_l] + _p)] = S_vec_n[1,:, index[_l] + _p] # e anche i relativi Singular vectors alla SFOF bell
                _p +=1
                
        _p = 0
        while MaC(_fi, S_vec_n[0,:, index[_l] - _p]) > MAClim:
                SDOFsval[(index[_l] - _p)] = S_val[0,0, index[_l] - _p] # si aggiungono i relativi Singular values alla SFOF bell
                SDOFsval1[(index[_l] - _p)] = _fi.conj().T@PSD_matr[:,:, index[_l] - _p]@_fi # si aggiungono i relativi Singular values alla SFOF bell
                SDOFsvec[:,(index[_l] - _p)] = S_vec_n[0,:, index[_l] - _p] # e anche i relativi Singular vectors alla SFOF bell
                _p +=1     
        while MaC(_fi, S_vec_n[1,:, index[_l] - _p]) > MAClim:
                SDOFsval[(index[_l] - _p)] = S_val[1,1, index[_l] - _p] # si aggiungono i relativi Singular values alla SFOF bell
                SDOFsval1[(index[_l] - _p)] = _fi.conj().T@PSD_matr[:,:, index[_l] - _p]@_fi # si aggiungono i relativi Singular values alla SFOF bell
                SDOFsvec[:,(index[_l] - _p)] = S_vec_n[1,:, index[_l] - _p] # e anche i relativi Singular vectors alla SFOF bell
                _p +=1  
        
        if method == 'EFDD':
            SDOFsval = SDOFsval
        elif method == 'FSDD':
            SDOFsval = SDOFsval1
        
        # indice dei SV in SDOFsval        
        idSV = np.array(np.where(SDOFsval)).T
        fsval = f[idSV]
        
        # Forme modali associate a tutti i Singular Values e mediate rispetto ai Singular Values stessi
        FIs = [SDOFsval[idSV[_u]]*SDOFsvec[:,idSV[_u]] for _u in range(len(idSV))]
        FIs = np.squeeze(np.array(FIs))
        
        # media e devstd
        stdFIs = np.std(FIs.real,axis=0)
        meanFi = np.mean(FIs,axis=0)
 
        _idx6 = np.argmax(abs(meanFi))
        # Forma Modale (Normalizzata spostamento unitario)
#        meanFireal = meanFi.real/meanFi[_idx6].real 
        meanFi = meanFi/meanFi[_idx6] 
       
        # PLOT 1 - Plotto la SDOF Bell Function estratta 
        _fig, ((_ax1,_ax2),(_ax3,_ax4)) = plt.subplots(nrows=2,ncols=2)
        _ax1.plot(f, 10*np.log10(S_val[0,0]), c='b')
        _ax1.plot(fsval, 10*np.log10(SDOFsval[idSV].real), c='r',label='SDOF bell')
        _ax1.set_title("SDOF Bell function")
        _ax1.set_xlabel('Frequency [Hz]')
        _ax1.set_ylabel(r'dB $[V^2/Hz]$')
        _ax1.legend()
        
        # Funzione di Autocorrelazione (Free Decay)
        SDOFcorr = np.fft.ifft(SDOFsval,n=nIFFT,axis=0,norm='ortho').real # y
        timeLag = np.linspace(0,tlag,len(SDOFcorr)) # t
        
        # PLOT 2 (AUTOCORRELATION)
        idxmax = np.argmax(SDOFcorr)
        normSDOFcorr = SDOFcorr[:len(SDOFcorr)//2]/SDOFcorr[idxmax]
        _ax2.plot(timeLag[:len(SDOFcorr)//2], normSDOFcorr)
        _ax2.set_title("Auto-correlation Function")
        _ax2.set_xlabel('Time lag[s]')
        _ax2.set_ylabel('Normalized correlation')    
            
        # Trovo dove incrocio asse x = 0
        sgn = np.sign(normSDOFcorr).real # trovo il segno
        sgn1 = np.diff(sgn,axis=0) # trovo dove cambia il segno (incrocio asse x=0)
        zc1 = np.where(sgn1)[0] # indici di Zero Crossing
    
        # trovo i massimi e i minimi (picchi) della autocorrelazione
        maxSDOFcorr = [np.max(normSDOFcorr[zc1[_i]:zc1[_i+2]]) for _i in range(0,len(zc1)-2,2)]
        minSDOFcorr = [np.min(normSDOFcorr[zc1[_i]:zc1[_i+2]]) for _i in range(0,len(zc1)-2,2)]
        if len(maxSDOFcorr) > len(minSDOFcorr):
            maxSDOFcorr = maxSDOFcorr[:-1]
        elif len(maxSDOFcorr) < len(minSDOFcorr):
            minSDOFcorr = minSDOFcorr[:-1]
        minmax = np.array((minSDOFcorr, maxSDOFcorr))
        minmax = np.ravel(minmax, order='F')
        
        # trovo gli indici dei massimi e i minimi (picchi) della autocorrelazione
        maxSDOFcorr_idx = [np.argmin(abs(normSDOFcorr-maxx)) for maxx in maxSDOFcorr]
        minSDOFcorr_idx = [np.argmin(abs(normSDOFcorr-minn)) for minn in minSDOFcorr]
        minmax_idx = np.array((minSDOFcorr_idx, maxSDOFcorr_idx))
        minmax_idx = np.ravel(minmax_idx, order='F')
    
        minmax_fit = np.array([minmax[_a] for _a in range(sppk,sppk+npmax)])
        minmax_fit_idx = np.array([minmax_idx[_a] for _a in range(sppk,sppk+npmax)])
        
        # PLOT 3 (PORTION for FIT)
        _ax3.plot(timeLag[:minmax_fit_idx[-1]], normSDOFcorr[:minmax_fit_idx[-1]])
        _ax3.scatter(timeLag[minmax_fit_idx], normSDOFcorr[minmax_fit_idx])
        _ax3.set_title("Portion for fit")
        _ax3.set_xlabel('Time lag[s]')
        _ax3.set_ylabel('Normalized correlation')    
        # plt.show()
        
        Td = np.diff(timeLag[minmax_fit_idx])*2
        Td_EFDD = np.mean(Td)
        
        fd_EFDD = 1/Td_EFDD
    # =============================================================================
        # build a rectangle in axes coords
        left, width = .25, .5
        bottom, height = .25, .5
        right = left + width
        top = bottom + height
        # axes coordinates are 0,0 is bottom left and 1,1 is upper right

        if method == 'EFDD':
            delta = np.array([np.log(np.abs(minmax[0])/np.abs(minmax[_i])) for _i in range(len(minmax_fit))])
        elif method == 'FSDD':
            delta = np.array([2*np.log(np.abs(minmax[0])/np.abs(minmax[_i])) for _i in range(len(minmax_fit))])
            
        _fit = lambda x,m:m*x
        m, _ = curve_fit(_fit, np.arange(len(minmax_fit)), delta)
    
        xi_EFDD = m/np.sqrt(4*np.pi**2 + m**2)
    
        fn_EFDD = fd_EFDD/np.sqrt(1-xi_EFDD**2)
    
    # =============================================================================
         # PLOT 4 (FIT)
        _ax4.scatter(np.arange(len(minmax_fit)), delta)
        _ax4.plot(np.arange(len(minmax_fit)), m*np.arange(len(minmax_fit)))
     
        _ax4.text(left, top, r'''$f_n$ = %.2f
$\xi$ = %.2f%s'''% (fn_EFDD, float(xi_EFDD)*100,"%"),transform=_ax4.transAxes)
         
     
         # _ax4.text(0.6, 0.6, r'''$f_n$ = {0:.2f} $Hz$
         #           $\xi$ = {1:.2f} %'''.format(fn_EFDD, float(xi_EFDD)*100))
        _ax4.set_title("Fit - Frequency and Damping")
        _ax4.set_xlabel(r'counter $k^{th}$ extreme')
        _ax4.set_ylabel(r'$2ln\left(r_0/|r_k|\right)$')    
    # =============================================================================
    # =============================================================================
        plt.tight_layout()
    
        Freq_E.append(fn_EFDD)
        Damp_E.append(xi_EFDD)
        Fi_E.append(meanFi)
    
    Freq = np.array(Freq_E)
    Damp = np.array(Damp_E)
    Fi = np.array(Fi_E)

    Results={}
    Results['Frequencies'] = Freq
    Results['Damping'] = Damp
    Results['Mode Shapes'] = Fi.T

    return Results

