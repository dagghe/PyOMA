# -*- coding: utf-8 -*-
"""
Created on Sun Oct 20 21:59:26 2019

@author: dagpa
"""
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
from mpl_toolkits.mplot3d import Axes3D
import glob
import numpy as np
from scipy import signal
from scipy.optimize import curve_fit
import pandas as pd
import seaborn as sns
import mplcursors
from sklearn.preprocessing import StandardScaler
import threading

# =============================================================================
# IMPORTA DATI
# =============================================================================
# Salva il percorso per la cartella contenente i file d'interesse

_file = str(sys.argv[1]) # Selezionare l'indice del file d'interesse
antype = str(sys.argv[2]) #seleziono il tipo: fdd efdd SSI
stage = str(sys.argv[3]) #stage: 1 o 2

# =============================================================================
_title= _file.rsplit('\\')[-1]  # crea una variable per il titolo (nome del file) 
_title= _title.rsplit('.')[0]
_title= _title.rsplit('_')[0]
df = pd.read_csv(_file, header=0, sep="\t", index_col=False) # crea un dataframe dal file selezionato (attenzione a header e separatore)
dati = np.array(df)
#dati = np.delete(dati, (0), axis=1) # colonne da rimuovere
dati = np.delete(dati, (0, 5, 12), axis=1) # colonne da rimuovere

Data= {'Title':_title}

# Pre-Processing (Standardizzazzione, De-trend, Decimazione)
# =============================================================================
_StandScaler = StandardScaler(with_mean=True, with_std=False)
_StandScaler.fit(dati)
_dati_std = _StandScaler.transform(dati)

#Standardizzazione dati (RIMOZIONE MEDIA)
dati = _dati_std

_fs=int(sys.argv[4]) # [Hz] - Frequenza di campionamento ACQUISITORE

_q = int(sys.argv[5]) # Fattore di decimazione
_tipo = 'fir' # tipo filtro 
dati = signal.detrend(dati, axis=0) # Rimozione trend

if _q != 1:
    dati = signal.decimate(dati, _q, ftype=_tipo, axis=0) # Decimazione segnale
    _fs=_fs/_q # [Hz] - Frequenza di campionamento DECIMATA

_freq_max = _fs/2 # Frequenza di Nyquist
Data['Decimated Sampling Frequency'] = _fs
Data['Nyquist Frequency'] = _freq_max


_n_dat=dati.shape[0] #NUMERO DI DATI CAMPIONATI
_n_can=dati.shape[1] #NUMERO DI CANALI ACQUISITI
_dur=_n_dat/_fs # [sec] Durata acquisizione


Data['number of samples N'] = _n_dat
Data['number of channels'] = _n_can
Data['Acquisition lenght [sec]'] = _dur





def fDD():


    # Singular Value Plot - Frequency Domain Decomposition   
    # =============================================================================
    # Per ora non vengono salvate le forme modali ma solo le singular values
    _nseg = int(sys.argv[6]) # numero di segmenti in cui dividere le storie temporali per mediare le PSD (Hanning, 50% overlap di default non modificati)
    _nxseg = (_n_dat)//_nseg 
    _PSD_matr = np.zeros((_n_can, _n_can, int((_nxseg)/2+1)), dtype=complex) # Inizializzo la SD matrix
    _S_val = np.zeros((_n_can, _n_can, int((_nxseg)/2+1))) # Inizializzo la matrice dove salverò i Singular Values
    _S_vec = np.zeros((_n_can, _n_can, int((_nxseg)/2+1)), dtype=complex) # Inizializzo la matrice dove salverò i Singular Vectors
    
    # loop dove mi calcolo le Auto e Cross-Spectral Density
    # (si passa al dominio della frequenza)
    for _i in range(0, _n_can):
        for _j in range(0, _n_can):
            _f, _Pxy = signal.csd(dati[:, _i],dati[:, _j], fs=_fs, nperseg=_nxseg)
            _PSD_matr[_i, _j, :] = _Pxy
            
    # loop dove mi calcolo i singular value      
    for _i in range(np.shape(_PSD_matr)[2]):
        _U1, _S1, _V1_t = np.linalg.svd(_PSD_matr[:,:,_i])
        _U1_1=np.transpose(_U1) 
        _S1 = np.diag(_S1)
        _S1rad=np.sqrt(_S1)
        _S_val[:,:,_i] = _S1rad
        _S_vec[:,:,_i] = _U1_1
    
    # Plot dei singular values (in scala logaritmica)
    _fig, _ax = plt.subplots()
    for _i in range(_n_can):
    #    ax.semilogy(_f, _S_val[_i, _i]) # scala log
        _ax.plot(_f, 10*np.log10(_S_val[_i, _i])) # decibel
    _ax.grid()
    _ax.set_xlim(left=0, right=_freq_max)
    _ax.xaxis.set_major_locator(MultipleLocator(_freq_max/10))
    _ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    _ax.xaxis.set_minor_locator(MultipleLocator(_freq_max/100))
    _ax.set_title('''Singular Values FDD 
                 {0}'''.format(_title))


    
    _tlag = _dur/_nseg # Durata 
    
    Data['Frequency Domain Decomposition'] = {'INPUT': {'number of segments':_nseg}}
    Data['Frequency Domain Decomposition']['INPUT']['n of point per segment'] = _nxseg
    Data['Frequency Domain Decomposition']['OUTPUT'] = {'PSD Matrix': _PSD_matr}
    Data['Frequency Domain Decomposition']['OUTPUT']['Singular Values matrix'] = _S_val
    Data['Frequency Domain Decomposition']['OUTPUT']['Singular Vector matrix'] = _S_vec
    Data['Frequency Domain Decomposition']['OUTPUT']['Frequency'] = _f
        

    mplcursors.cursor() 
    plt.plot()

    if antype == "fdd" and stage == "2": 

        _x = float(sys.argv[7])

        
        # Frequency Domain Decomposition 
        # =============================================================================
        # Estreazione Frequenza d'interesse
        _idx = np.argmin(abs(_f-_x))
        _lim = (_x-_x*0.01,_x+_x*0.01)
        _idxlim = (np.argmin(abs(_f-_lim[0])), np.argmin(abs(_f-_lim[1])))
        # da controlloare (S1-S2 o S1/S2 ????)!!!!
        _diffS1S2 = _S_val[0,0,_idxlim[0]:_idxlim[1]] - _S_val[1,1,_idxlim[0]:_idxlim[1]]
        _maxDiffS1S2 = np.max(_diffS1S2)
        _idx1 = np.argmin(abs(_diffS1S2 - _maxDiffS1S2))
        _idxfin = _idxlim[0] + _idx1
        # =============================================================================
        # Estrazione Forma Modale
        _fr_FDD = _f[_idxfin]
        _fi_FDD = _S_vec[0,:,_idxfin]
        _idx3 = np.argmax(abs(_fi_FDD))
        _fi_FDDn = _fi_FDD/_fi_FDD[_idx3] # Forma Modale (Normalizzata spostamento unitario)
        
        _fig, _ax = plt.subplots()
        _ax.plot(_fi_FDDn.real)
        plt.show()
        
        Data['Frequency Domain Decomposition']['RISULTATI'] = {'Initial Frequency': _x}
        Data['Frequency Domain Decomposition']['RISULTATI']['FDD Frequency'] = _fr_FDD
        Data['Frequency Domain Decomposition']['RISULTATI']['FDD Mode Shape'] = _fi_FDDn


    elif antype == "efdd" and stage == "2": 
        # Enhanced Frequency Domain Decomposition 
        # =============================================================================
        
        _x = float(sys.argv[7])
        _npmax = int(sys.argv[8]) # numero di massimi da considerare per decremento log
        _skipPeak = int(sys.argv[9]) # numero di picchi (massimi) da saltare
        _limMAC = float(sys.argv[10]) # limite MAC per selezionamento SDOF Bell Function
        
        # Frequency Domain Decomposition 
        # =============================================================================
        # Estreazione Frequenza d'interesse
        _idx = np.argmin(abs(_f-_x))
        _lim = (_x-_x*0.01,_x+_x*0.01)
        _idxlim = (np.argmin(abs(_f-_lim[0])), np.argmin(abs(_f-_lim[1])))
        # da controlloare (S1-S2 o S1/S2 ????)!!!!
        _diffS1S2 = _S_val[0,0,_idxlim[0]:_idxlim[1]] - _S_val[1,1,_idxlim[0]:_idxlim[1]]
        _maxDiffS1S2 = np.max(_diffS1S2)
        _idx1 = np.argmin(abs(_diffS1S2 - _maxDiffS1S2))
        _idxfin = _idxlim[0] + _idx1
        # =============================================================================
        # Estrazione Forma Modale
        _fr_FDD = _f[_idxfin]
        _fi_FDD = _S_vec[0,:,_idxfin]
        _idx3 = np.argmax(abs(_fi_FDD))
        _fi_FDDn = _fi_FDD/_fi_FDD[_idx3] # Forma Modale (Normalizzata spostamento unitario)
        
        _nIFFT = (len(_f)-1)*20 # num punti per trasformata inversa (se > di len(_f), viene zero-paddata)
        _idxdasaltare = _skipPeak*2+1 # numero di massimi da considerare per decremento log
        
        _Fi = _S_vec[0,:,:]
        _idx4 = np.array([np.argmax(abs(_Fi[:,_i])) for _i in range(_Fi.shape[1])])
        _Fi_FDDsn = np.array([_Fi[:,_j]/_Fi[_idx4[_j],_j] for _j in range(len(_idx4))]).T
        
        
        _SDOFsval = np.zeros(_Fi.shape[1])
        for _i in range(len(_f)):
            _aMAC= np.abs(_fi_FDDn.real@_Fi_FDDsn[:,_i].real)**2 / \
                ((_fi_FDDn.real@_fi_FDDn.real)*(_Fi_FDDsn[:,_i].real@_Fi_FDDsn[:,_i].real)) # autoMAC
            if _aMAC > _limMAC:
                _SDOFsval[_i] = _S_val[0,0,_i]
            else:
                _SDOFsval[_i] = 0
        
        aaaa = np.array(np.where(_SDOFsval)).T
        BBB = [_S_val[0,0,aaaa[_u]]*_S_vec[0,:,aaaa[_u]] for _u in range(len(aaaa))]
        
        
        #_fig, _ax = plt.subplots()
        #_ax.plot(BBB[109].ravel())
        #plt.show()
        
        # Plotto la SDOF Bell Function estratta (da sola)
        _fig, _ax = plt.subplots()
        _ax.plot(_f, _SDOFsval)
        mplcursors.cursor() 
        plt.show()
        
        # Funzione di Autocorrelazione (Free Decay)
        _SDOFcorr = np.fft.ifft(_SDOFsval,n=_nIFFT,axis=0,norm='ortho').real # y
        _timeLag = np.linspace(0,_tlag,len(_SDOFcorr)) # t
        
        # Plotto la funzione di autocorrelazione 
        _fig, _ax = plt.subplots()
        _ax.plot(_timeLag[:len(_SDOFcorr)//2], _SDOFcorr[:len(_SDOFcorr)//2])
        mplcursors.cursor() 
        plt.show()
        
        # Trovo dove incrocio asse x = 0
        _asd = np.sign(_SDOFcorr).real # trovo il segno
        _asd1 = np.diff(_asd,axis=0) # trovo dove cambia il segno (incrocio asse x=0)
        _zc1 = np.where(_asd1)[0] # indici di Zero Crossing
        
        # Estraggo una porzione (data dal numero di picchi) della autocorrelazione
        _SDOFcorr_perFit = _SDOFcorr[_zc1[_idxdasaltare]:_zc1[_idxdasaltare+_npmax*2]]
        _tlag_fit = np.arange(_timeLag[_zc1[_idxdasaltare]],_timeLag[_zc1[_idxdasaltare+_npmax*2]],_timeLag[1])
        
        # trovo i massimi e i minimi (picchi) della autocorrelazione
        _maxSDOFcorr_perFit = np.array([np.max(_SDOFcorr[_zc1[_i]:_zc1[_i+1]]) for _i in range(_idxdasaltare,_idxdasaltare+_npmax*2,2)])
        _minSDOFcorr_perFit = np.array([np.min(_SDOFcorr[_zc1[_i]:_zc1[_i+1]]) for _i in range(_idxdasaltare+1,_idxdasaltare+_npmax*2,2)])
        _minmax = np.array((_maxSDOFcorr_perFit,_minSDOFcorr_perFit))
        _minmax = np.ravel(_minmax, order='F')
        
        # trovo gli indici dei massimi e i minimi (picchi) della autocorrelazione
        _maxSDOFcorr_perFit_idx1 = [np.argmin(abs(_SDOFcorr_perFit-_maxSDOFcorr)) for _maxSDOFcorr in _maxSDOFcorr_perFit]
        _minSDOFcorr_perFit_idx1 = [np.argmin(abs(_SDOFcorr_perFit-_minSDOFcorr)) for _minSDOFcorr in _minSDOFcorr_perFit]
        _minmax_idx = np.array((_maxSDOFcorr_perFit_idx1, _minSDOFcorr_perFit_idx1))
        _minmax_idx = np.ravel(_minmax_idx, order='F')
        
        _fig, _ax = plt.subplots()
        _ax.plot(_tlag_fit[:], _SDOFcorr_perFit)
        #_ax.scatter(_timeLag[_maxSDOFcorr_perFit_idx1], _maxSDOFcorr_perFit)
        #_ax.scatter(_timeLag[_minSDOFcorr_perFit_idx1], _minSDOFcorr_perFit)
        _ax.scatter(_tlag_fit[_minmax_idx], _minmax)
        plt.show()
        
        _Td = np.diff(_timeLag[_maxSDOFcorr_perFit_idx1])
        _Td_EFDD = np.mean(_Td)
        _omegaD = 2*np.pi/_Td_EFDD
        _fd_EFDD = 1/_Td_EFDD
        
        _delta = np.array([2*np.log(_minmax[0]/np.abs(_minmax[_i])) for _i in range(len(_minmax))])
        
        _fit = lambda x,m:m*x
        _m, _ = curve_fit(_fit, np.arange(len(_minmax)), _delta)
        
        
        _fig, _ax = plt.subplots()
        _ax.scatter(np.arange(len(_minmax)), _delta)
        _ax.plot(np.arange(len(_minmax)), _fit(np.arange(len(_minmax)), _m))
        # ax = sns.regplot(np.arange(len(minmax)),delta)
        _ax.legend()
        plt.show()
        
        _xi_EFDD = _m/np.sqrt(4*np.pi**2 + _m**2)
        _fn_EFDD = _fd_EFDD/np.sqrt(1-_xi_EFDD**2)
        
        
        Data['Enhanced FDD'] = {'INPUT':{'number of peaks(max)':_npmax}}
        Data['Enhanced FDD']['INPUT']['Number of peak to skip'] = _skipPeak
        Data['Enhanced FDD']['INPUT']['AutoMAC limit Bell Function'] = _limMAC
        Data['Enhanced FDD']['INPUT']['Number of points IFT'] = _nIFFT
        Data['Enhanced FDD']['OUTPUT'] = {'Sing. Val. SDOF Bell func.':_SDOFsval}
        Data['Enhanced FDD']['OUTPUT']['AutoCorr function (Free Decay)'] = \
            _SDOFcorr[:len(_SDOFcorr)//2]
        Data['Enhanced FDD']['OUTPUT']['Time Lag AutoCorr f.'] = \
            _timeLag[:len(_SDOFcorr)//2]
        Data['Enhanced FDD']['OUTPUT']['Extracted Free Decay to fit'] = _SDOFcorr_perFit
        Data['Enhanced FDD']['OUTPUT']['Time Lag AutoCorr f.'] = _timeLag
        Data['Enhanced FDD']['OUTPUT']['All Peaks (min,max)'] = _minmax
        Data['Enhanced FDD']['OUTPUT']['All Peaks index (min,max)'] = _minmax_idx
        Data['Enhanced FDD']['RISULTATI'] = {'Damped Frequency': _fd_EFDD}
        Data['Enhanced FDD']['RISULTATI']['xi'] = float(_xi_EFDD)
        Data['Enhanced FDD']['RISULTATI']['Natural Frequency'] = float(_fn_EFDD)

def SSI():
    
    # INPUT
    # =============================================================================
    _br=int(sys.argv[6])#SHIFT TEMPORALE (numero di blocchi riga)
    _ordine_max = int(sys.argv[7]) # ORDINE MASSIMO PER DIAGRAMMA DI STABILIZZAZIONE
    _ordine_min = int(sys.argv[8]) # L'ordine minimo (da cui parte il ciclo) è fissato a priori uguale a questo valore
    
    _lim_f = 0.01 # limite CONDIZIONE 1 sulla frequenza: ( f(n,m) - f(n-1,m) < lim_f; dove f è la frequenza, n è l'ordine considerato e m è il modo considerato)
    _lim_s = 0.05 # limite CONDIZIONE 2 sullo smorzamento: ( d(n,m) - d(n-1,m) < lim_s; dove f è lo smorzamento, n è l'ordine considerato e m è il modo considerato)
    _lim_ms = 0.02 # limite CONDIZIONE 1 sulla forma modale
    
    _lim_s1 = 0.1 # limite di smorzamento oltre il quale rimuovere i poli dalla matrice Fr quando si plotta il diagramma di stabilizzazione
    
    _n_dat=dati.shape[0] #NUMERO DI DATI CAMPIONATI
    _n_can=dati.shape[1] #NUMERO DI CANALI ACQUISITI
    _dur=_n_dat/_fs # [sec] Durata acquisizione
    
    
    Data['number of samples N'] = _n_dat
    Data['number of channels'] = _n_can
    Data['Acquisition lenght [sec]'] = _dur
    
    Data['Cov. Stochastic Subspace Identification'] = \
        {'INPUT': {'number of block rows ii': _br}}
    Data['Cov. Stochastic Subspace Identification']['INPUT']['Stable freq. limit'] = _lim_f
    Data['Cov. Stochastic Subspace Identification']['INPUT']['Stable damp. limit'] = _lim_s
    Data['Cov. Stochastic Subspace Identification']['INPUT']['Stable mode shape limit'] = _lim_ms    
        
    if stage == "1" or stage == "2":
        
        
        # SSI-COV INDENTIFICAZIONE FORME E FREQUENZE DI VIBRARE
        # =============================================================================
        _Yy=np.array(dati.transpose())
        _j=_n_dat-2*_br+1; # DIMENSIONE MATRICE DI HANKEL
        
        _H=np.zeros((_n_can*2*_br,_j)) # PREALLOCA LA MATRICE DI HANKEL
        for _k in range(0,2*_br):
        	_H[_k*_n_can:((_k+1)*_n_can),:]=_Yy[:,_k:_k+_j] # CALCOLO MATRICE DI HANKEL
        
        _Hp = _H[:_br*_n_can,:] # MATRICE HANKEL PASSATO
        _Hf = _H[_br*_n_can:,:] # MATRICE HANKEL FUTURO
        
        _fatt=1/(_j-_br); # FATTORE MATRICE COVARIANZA DEGLI OUTPUT
        Tb=_fatt*(_Hf@_Hp.T) # MATRICE DI COVARIANZA DEGLI OUTPUT 
        # =============================================================================
        # SINGULAR VALUE DECOMPOSITION
        # =============================================================================
        _U1, _S1, _V1_t = np.linalg.svd(Tb)
        _U1_1=np.transpose(_U1) 
        _S1 = np.diag(_S1)
        _S1rad=np.sqrt(_S1)
        # =============================================================================
        # Ciclo per ordine del sistema crescente
        # =============================================================================
        Fr=np.full((_ordine_max, _ordine_max), np.nan) # inizializzo la matrice che conterrà le frequenze
        Fr_lab=np.full((_ordine_max, _ordine_max), np.nan)  # inizializzo la matrice che conterrà le labels per i poli(frequenze) da stampare nel diagramma di stabilizzazione
        Sm=np.full((_ordine_max, _ordine_max), np.nan) # inizializzo la matrice che conterrà gli smorzamenti
        Ms = []  # inizializzo la matrice (3D) che conterrà le forme modali
        for z in range(0, int((_ordine_max - _ordine_min)/2+1)):
            Ms.append(np.zeros((_n_can, _ordine_min + z*(2))))
        
        for _iteraz, _indice in enumerate(range(_ordine_min, _ordine_max+1, 2)): # INIZIA IL CICLO PER IL CALCOLO DELLE FREQUENZE AL CRESCERE DELL'ORDINE DEL SISTEMA
                
            _S11 = np.zeros((_indice, _indice)) # Inizializzo
            _U11 = np.zeros((_ordine_max, _indice)) # Inizializzo
            _O_1 = np.zeros((_ordine_max - _n_can, _indice)) # Inizializzo
            _O_2 = np.zeros((_ordine_max - _n_can, _indice)) # Inizializzo
        
            _S11[:_indice, 0:_indice] = _S1rad[0:_indice,0:_indice] # ESTRAZIONE DI UNA SOTTOMATRICE DEI SINGULAR VALUES DA ORDINE MIN A ORD MAX
        
            _U11[:_ordine_max+1, 0:_indice] = _U1[0:_ordine_max,0:_indice] # ESTRAZIONE DI UNA SOTTOMATRICE DEI LEFT SINGULAR VECTORS DA ORDINE MIN A ORD MAX
        
            _O = _U11 @ _S11 # CALCOLO MATRICE DI OSSERVABILITA
        
            _O_1[:,:] = _O[:_O.shape[0] - _n_can,:]
            _O_2[:,:] = _O[_n_can:,:]
        
            _A = np.dot(np.linalg.pinv(_O_1), _O_2) # STIMA DELLA MATRICE DINAMICA A TEMPO DISCRETO
            [_AuVal, _AuVett] = np.linalg.eig(_A) # CALCOLO AUTOVALORI ED AUTOVETTORI
            _Lambda =(np.log(_AuVal))*_fs # CALCOLO _Lambda
            _fr = abs(_Lambda)/(2*np.pi) # CALCOLO FRQUENZE
            _smorz = -((np.real(_Lambda))/(abs(_Lambda))) # CALCOLO SMORZAMENTO
            
            # calcolo la matrice C per definizione dei modi
            # le deformate modali M in forma complessa sono calcolate attraverso
            # l'espressione M=CL dove L è la matrice le cui colonne sono gli autovettori
            # di A (_AuVett)
            _C = _O[:_n_can,:]
            _M =np.real(_C@_AuVett)
        
            Fr[:len(_fr),_indice-1] = _fr # SALVA LE FREQUENZE    
            Sm[:len(_fr),_indice-1] = _smorz # SALVA gli smorzamenti
            Ms[_iteraz] = _M # SALVA le forme modali
        
            # Quì vengono fatti i controlli per vedere se i poli sono stabili.
            # Il funzionamento è il seguente: si fa un ciclo su tutti i poli (frequenze)
            # trovati ad ogni nuovo ordine/iterazione (vettore _fr). Si fa la differenza 
            # tra questo polo e il vettore di tutti i poli alla iterazione/ordine precedente (_fr alla iterazione n-1).
            # Si ottiene quindi un vettore con le suddette differenze. Di questo si 
            # prende il minimo e si controlla che rispetti il limite prefissato;
            # se viene superato questo controllo, si procede a controllare la differenza
            # tra gli smorzamenti e tra le forme modali (autoMAC). Si può quindi 
            # classificare un polo come:
            # 0 = nuovo o instabile ; 1 = Frequenza stabile; 2 = stabile per Frequenza + Smorzamento
            # 3 = stabile per Frequenza + Forma modale; 4 = Stabile per Frequenza + Smorzamento + Forma modale
            for idx, (_freq, _smor) in enumerate(zip(_fr,_smorz)):
                if _iteraz == 0: # alla prima iterazione/primo ordine sono tutti nuovi poli
                    Fr_lab[:len(_fr),_indice-1] = 0 # 0 è l'etichetta per i nuovi poli e poli instabili
                elif _iteraz != 0:
                    if min(abs(_freq - Fr[:,(_indice-1) - 2])/_freq) < _lim_f: # CONDIZIONE 1 SULLA FREQUENZA
                        _indice2 = np.nanargmin(  # Trovo l'indice del polo che minimizza la differenza alla iterazione/ordine n-1
                                                abs(_freq - Fr[:,(_indice-1) - 2])
                                                - min(abs(_freq - Fr[:,(_indice-1) - 2])))
                            
                        _Fi_n = _M[:, idx] # forma modale all'iterazione = _iteraz (attuale)
                        _Fi_nmeno1 = Ms[_iteraz-1][:,_indice2] # forma modale all'iterazione = _iteraz-1 (precedente)
                        
                        _aMAC = np.abs(_Fi_n@_Fi_nmeno1)**2 / ((_Fi_n@_Fi_n)*(_Fi_nmeno1@_Fi_nmeno1)) # autoMAC
                        
        
                        if (_smor - Sm[_indice2, (_indice-1) - 2])/_smor < _lim_s: # CONDIZIONE 2 SULLO SMORZAMENTO
                            if 1 - _aMAC < _lim_ms: # CONDIZIONE 3 SULLE FORME MODALI
                                Fr_lab[idx,_indice-1] = 4 # Se il polo è stabile sia per smorzamento che per forma modale (viene classificato come stabile/verde)
                            else:
                                Fr_lab[idx,_indice-1] = 2 # Stabile per smorzamento ma non per forma modale
        
        
                        elif 1 - _aMAC < _lim_ms: # CONDIZIONE 3 SULLE FORME MODALI
                            Fr_lab[idx,_indice-1] = 3 # Stabile per forma modale ma non per smorzamento
        
                        else:
                            Fr_lab[idx,_indice-1] = 1 # Stabile solo per frequenza
                    else:
                        Fr_lab[idx,_indice-1] = 0  # Nuovo polo o polo instabile
        
        
        
        
        # PLOT DIAGRAMMA DI STABILIZZAZZIONE + PLOT DEI SINGULAR VALUES SOVRAPPOSTO
        # =============================================================================
        _x = Fr.flatten(order='f')
        _y = np.array([_i//len(Fr) + 1 for _i in range(len(_x))])
        _l = Fr_lab.flatten(order='f')
        _d = Sm.flatten(order='f')
        
        _df = pd.DataFrame(dict(Frequency=_x, Order=_y, Label=_l, Damp=_d))
        _df1 = _df.copy()
        _df1.Frequency = _df1.Frequency.where(_df.Damp < _lim_s1) # qui rimuovo tutti i poli che hanno smorzamento maggiore di lim_s1
        _colors = {0:'Red', 1:'darkorange', 2:'gold', 3:'yellow', 4:'Green'} # assegno i colori alle etichette dei poli
        
        fig1, ax1 = plt.subplots(figsize=(12,6))
        ax1 = sns.scatterplot(_df1['Frequency'], _df1['Order'], hue=_df1['Label'], palette=_colors)
        
        ax1.set_xlim(left=0, right=_freq_max)
        ax1.set_ylim(bottom=0, top=_ordine_max)
        ax1.xaxis.set_major_locator(MultipleLocator(_freq_max/10))
        ax1.xaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax1.xaxis.set_minor_locator(MultipleLocator(_freq_max/100))
        ax1.set_title('''{0}
        shift: {1}'''.format(_title, _br))
        # =============================================================================
        # Creo un secondo ax che condivide lo stesso asse-x con quello definito
        # precedentemente; poi plotto le singolar value trovate prima con la FDD su 
        # scala logaritmica (semilogy)
        
        #ax2 = ax1.twinx()  
        #
        #for _i in range(_n_can):
        #    ax2.semilogy(f, _S_val[_i, _i], c='LightGrey')
        mplcursors.cursor()
        plt.show()
        
        Risultati = _df1

    if stage == "2":
        # ESTRAZIONE FORME MODALI D'INTERESSE
        # =============================================================================
        # Istruzioni:
        # 1. Rimuovere (# ctrl-1) plot dei singular value da plot sovrapposto (secondo asse_y è in scala logaritmica)
        # 2. Plottare diagramma di stabilizzazzione e individuare i poli stabili 
        #    relativi ai modi che si desidera plottare (segnarsi le coordinate x,y).
        # 3. Andare nella matrice Fr delle frequenze, trovare la colonna = a y(ordine),
        #    trovare la frequenza x nel vettore colonna e segnarsi l'indice di riga = m;
        # 4. Andare nella lista dei modi Ms l'indice da scegliere è: 
        #    z = (y(ordine) - ordine_min)/2 (l'ordine si può anche trovare al secondo indice nella lista Ms)
        #    nella matrice scegliere la colonna m;
        #    questa corrisponderà alla forma modale cercata: fi = Ms[z][:,m]
         
        # scrivere qui la x e la y del polo individuato sul diagramma di stabilizzazzione
        # (se non ci sono 2 poli troppo vicini la frequenza può essere approssimata 
        # ai termini interi) 
        _x = float(sys.argv[9])
        _y = int(sys.argv[10])
        
        _z = int((_y - _ordine_min)/2) # indice/ della matrice FI di interesse (relativa all'ordine del sistema)
        # Trovo l'indice del polo 
        _m = np.nanargmin(abs(_x - Fr[:, _y-1])) # trovo l'indice 
        # N.B. y-1 a causa dell'indicizzazione a 0 nelle matrici Fr e Sm
        freq = Fr[_m, _y-1] # estraggo frequenza d'interesse
        xi = Sm[_m, _y-1] # estraggo smorzamento associato
        fi = Ms[_z][:,_m] # estraggo forma modale associata
        
        print('''
        Frequenza = {0:.3f} Hz
        Smorzamento = {1:.2f} %'''.format(freq, xi*100))
        
        _idx = np.argmax(abs(fi))
        fi_norm = fi/fi[_idx]
        
        
        fig, ax = plt.subplots()
        ax.plot(fi_norm)
        plt.show()
        
        
        

if antype == "fdd" or antype == "efdd":
    fDD()
else:
    SSI()        