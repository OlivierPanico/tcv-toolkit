#%%
import MDSplus as mds
from MDSplus.mdsExceptions import TreeNNF
import numpy as np
import sys
from matplotlib import pyplot as plt
import scipy.io
from scipy import interpolate
import warnings
from matlabtools import Struct

# CLASSE PER CARICARE VARI SEGNALI E FARE QUALCHE MEDIA


# GET_THOMSON: temperatura, densita, tempo, r, z, del Thomson
# MEDIE_THOMSON_Z: profili medi T(r), n(r) tra 2 tempi scelti
# GET_THOMSON_FIT: carica i fit del thomson. NON usare, fanno schifo
# GET_CRXS: carica i fit del scambio carica. NON usare, fanno schifo
# GET_CRXS_NOFIT: profili dello scambio carica (velocita, Ti, ni) grezzi
# GET_PRESSURE: Carica pressione neutri baratron equatore esterno (NON SO SE FUNZIONA PER SPARI DEL 2022)
# GET_FLUX: Flusso del fuelling (NON SO SE FUNZIONA PER SPARI DEL 2022)
# GET_PARAM: parametri principali dello sparo (corrente, tensione...)
# GET_IP: carica solo la corrente di plasma
# GET_NBI: quasi tutti i parametri del HBI 1
# GET_ECRH: poche potenza totale del ECRH (2022: non funziona piu)
# GET_ECRH2: carica i 3 segnali di L1, L4 e L5
# GET_HA: tutti i filtri con anche il nome della specie
# GET_HA_ONE: carica un filtro Ha, basta mettere il numero di quello che vuoi caricare
# MEDIE_CX_NOFIT: medie tra 2 istanti tempoarali delle varie quantita del CRXS
# GET_THB: segnali grezzi del THB
# GET_THB_SINGLE: carica 1 segnale (667, canale 3). Veloce se serve solo 1 segnale
# GET_THB_PUFFING: segnali del puffing della GPI, che e anche il puffing usato dal THB
# GET_RCP: dati della reciprocating Langmuir probe (equatore esterno)
# GET_BOLO: carica i dati raw del bolometro. Se si da come input la los, carica solo quella line di vista
# GET_WALL_PROBES: carica i dati delle sonde di Langmuir alla prima parete. Jsat, ne, Te
# GET_FUELLING: carica i segnali in volts delle 3 valvole di fuelling, sia le misure che le richieste
# GET_FIR: carica il segnali del FIR (tutte le corde e la media)




class TCVShot(Struct):

    def __init__(self, shot):
        self.shot = shot
        self.r_vessel = np.array([1.1360, 1.1360, 1.1120, 1.0880, 1.0640, 1.0399, 1.0159,
                                 0.9919, 0.9679, 0.6707, 0.6240, 0.6240, 0.6240, 0.6724, 0.9679, 1.1360, 1.1360])
        self.z_vessel = np.array([0, 0.5494, 0.5781, 0.6067, 0.6354, 0.6640, 0.6927,
                                 0.7213, 0.7500, 0.7500, 0.7033, 0, -0.7033, -0.7500, -0.7500, -0.5494, 0.0000])

        print('')
        print('Opening shot ', shot)
        try:
            self.tree = mds.Tree('tcv_shot', self.shot)
            self.tag = True
        except:
            self.tag = False
            pass

    # CARICA Thomson

    def get_thomson(self):
        print(' Loading Thomson  \n')
        self.th_r = self.tree.getNode('\diagz::thomson_set_up:radial_pos').data()
        self.th_z = self.tree.getNode('\diagz::thomson_set_up:vertical_pos').data()
        try:
            self.time_th = self.tree.getNode('\RESULTS::thomson:times').data()
        except:
            self.time_th = self.tree.getNode('\RESULTS::thomson:te:foo').dim_of().data()
        try:
            self.th_t = self.tree.getNode('\RESULTS::thomson:te:foo').data()
        except:
            self.th_t = self.tree.getNode('\RESULTS::thomson:te').data()
        try:
            self.th_n = self.tree.getNode('\RESULTS::thomson:ne').data()
        except:
            self.th_n = self.tree.getNode('\RESULTS::thomson:ne:foo').data()
        self.th_error_t = self.tree.getNode('\RESULTS::thomson:te:error_bar').data()
        self.th_error_n = self.tree.getNode('\RESULTS::thomson:ne:error_bar').data()

    # calcola le medie con errori dei profili grezzi in Z del THOMSON da t1 a t2
    def medie_thomson_z(self, t1, t2):
        nzeta = len(self.th_z)

        media_te1 = np.zeros([nzeta, 2])
        media_ne1 = np.zeros([nzeta, 2])
        for i in range(nzeta):
            ind1 = (self.time_th < t2)*(self.time_th > t1) * \
                    (self.th_t[i, :] > 1.)*(self.th_n[i, :] > 1e15)
            media_te1[i, 0] = np.nanmean(self.th_t[i, ind1])
            media_te1[i, 1] = np.nanstd(self.th_t[i, ind1])
            media_ne1[i, 0] = np.nanmean(self.th_n[i, ind1])
            media_ne1[i, 1] = np.nanstd(self.th_n[i, ind1])

        media_te = media_te1
        media_ne = media_ne1
        return media_te, media_ne

    # CARICA Thomson FIT (FA SCHIFO)

    def get_thomson_fit(self):
        print('Loading Thomson FIT  \n')
        self.time_th = self.tree.getNode(
            '\RESULTS::thomson.profiles.auto:time').data()
        self.th_t = self.tree.getNode('\RESULTS::thomson.profiles.auto:te').data()
        self.th_n = self.tree.getNode('\RESULTS::thomson.profiles.auto:ne').data()
        self.th_rho = self.tree.getNode('\RESULTS::thomson.profiles.auto:rho').data()
        self.ne_axis = self.tree.getNode(
            '\RESULTS::thomson.profiles.auto:ne_axis').data()

    # CARICA CRXS (FA SCHIFO)
    def get_crxs(self):
        print('Loading CRXS FIT  \n')
        self.VtorNode = self.tree.getNode(r'\results::cxrs.proffit:vi_tor')
        self.rho_vtor = self.VtorNode.getDimensionAt(0).data()
        self.time_vtor = self.VtorNode.getDimensionAt(1).data()
        self.vtor = self.VtorNode.data()
        self.vtor_err = self.tree.getNode('\RESULTS::CXRS:PROFFIT:VI_TOR:ERR').data()

        self.VpolNode = self.tree.getNode(r'\results::cxrs.proffit:vi_pol')
        self.vpol = self.VpolNode.data()
        self.rho_vpol = self.VpolNode.getDimensionAt(0).data()
        self.time_vpol = self.VpolNode.getDimensionAt(1).data()
        self.vpol_err = self.tree.getNode('\RESULTS::CXRS:PROFFIT:VI_POL:ERR').data()
  

        try:
            self.niNode = self.tree.getNode('\RESULTS::CXRS.PROFFIT:NI')
            self.rho_ni = self.niNode.getDimensionAt(0).data()
            self.time_ni = self.niNode.getDimensionAt(1).data()
            self.ni = self.niNode.data()
            self.ni_err = self.tree.getNode('\RESULTS::CXRS.PROFFIT:NI:ERR').data()
        except:
            print(f'No ni for shot {self.shot}')
            self.ni = np.nan * np.ones_like(self.vtor)
            self.ni_err = np.nan * np.ones_like(self.vtor)
            self.rho_ni = np.nan * np.ones_like(self.vtor)
            self.time_ni = np.nan * np.ones_like(self.vtor)
        
        self.tiNode = self.tree.getNode('\RESULTS::CXRS.PROFFIT:TI')
        self.rho_ti = self.tiNode.getDimensionAt(0).data()
        self.time_ti = self.tiNode.getDimensionAt(1).data()
        self.ti = self.tiNode.data()
        self.ti_err = self.tree.getNode('\RESULTS::CXRS.PROFFIT:TI:ERR').data()

    def get_crxs_nofit(self):
        print('Loading CRXS WITHOUT FIT  \n')
        try:
            self.VtorNode_raw = self.tree.getNode(r'\results::cxrs:vi_tor')
            self.VpolNode_raw = self.tree.getNode(r'\results::cxrs:vi_pol')
            self.tag_cx_raw = True
        except:
            self.tag_cx_raw = False
            return

        self.rho_vtor_raw = self.VtorNode_raw.getDimensionAt(0).data()
        self.time_vtor_raw = self.VtorNode_raw.getDimensionAt(1).data()
        self.vtor_raw = self.VtorNode_raw.data()
        self.vtor_err_raw = self.tree.getNode('\RESULTS::CXRS:VI_TOR:ERR').data()
        self.ni_raw = self.tree.getNode('\RESULTS::CXRS.NI').data()
        self.ti_raw = self.tree.getNode('\RESULTS::CXRS.TI').data()
        self.ni_err_raw = self.tree.getNode('\RESULTS::CXRS:NI:ERR').data()
        self.ti_err_raw = self.tree.getNode('\RESULTS::CXRS:TI:ERR').data()
        self.time_vpol_raw = self.VpolNode_raw.getDimensionAt(1).data()
        self.vpol_err_raw = self.tree.getNode('\RESULTS::CXRS:VI_POL:ERR').data()
        self.vpol_raw = self.VpolNode_raw.data()
        self.rho_vpol_raw = self.VpolNode_raw.getDimensionAt(0).data()

    # CARICA Pressione neutri baratron equatore esterno

    def get_pressure(self):
        print('Loading pressure from baratron  \n')

        if self.shot > 55358:
            data_baratron = - \
                self.tree.getNode('\BASE::trch_baratron:channel_001').data()*0.267
            self.t_press = self.tree.getNode(
                '\BASE::trch_baratron:channel_001').dim_of().data()

        if self.shot <= 55358:
            data_baratron = - \
                self.tree.getNode('\BASE::trch_ece_pols:channel_002').data()*0.267
            self.t_press = self.tree.getNode(
                '\BASE::trch_ece_pols:channel_002').dim_of().data()

        ioff = (self.t_press > -6) & (self.t_press < -4)
        offset = np.mean(data_baratron[ioff])
        self.press = data_baratron-offset  # pressione in Pa

    # Flusso di particelle del fueling

    def get_flux(self):
        print('Loading particle fluxes (fueling)  \n')
        flux_d1 = self.tree.getNode('\diagz::flux_gaz:piezo_1:flux ').data()
        self.flux_t = self.tree.getNode(
            '\diagz::flux_gaz:piezo_1:flux ').dim_of().data()
        # da mbar*l/sec a particelle/sec (2 perche deuterio)
        self.flux_d = 2*flux_d1*0.1/(1.38e-23*293)

    # CARICA IP e altri parametri interessanti della scarica (IP lentissimo!)

    def get_param(self, noc=0):
        if noc == 0:
            print('Loading Ip and others...  \n')
            print('Lentissimo per caricare la corrente... \n')
            try:
                self.ip = self.tree.getNode('\magnetics::iplasma:trapeze').data()
                self.tip = self.tree.getNode('\magnetics::iplasma:trapeze').dim_of().data()
            except:
                self.ip = 0
                self.tip = 0
            self.ip2 = self.tree.getNode(r'\results::i_p').data()
            self.tip2 = self.tree.getNode(r'\results::i_p').dim_of().data()

        self.q95 = self.tree.getNode('\RESULTS::Q_95').data()
        self.betap = self.tree.getNode('\RESULTS::BETA_POL').data()
        self.betat = self.tree.getNode('\RESULTS::BETA_TOR').data()
        try:
            self.filling = self.tree.getNode('\diagz::flux_gaz:piezo_1:flux').data()
            self.t_filling = self.tree.getNode(
                '\diagz::flux_gaz:piezo_1:flux').dim_of().data()
        except:
            pass
            try:
                self.tvloop = self.tree.getNode('\magnetics::vloop').dim_of().data()
                self.vloop = self.tree.getNode('\magnetics::vloop').data()
            except:
                pass

        try:
            self.taue = self.tree.getNode(r'\results::conf:taue').data()
            self.t_taue = self.tree.getNode(r'\results::conf:taue').dim_of().data()
            self.tau = self.tree.getNode(r'\results::conf:tau').data()
            self.t_tau = self.tree.getNode(r'\results::conf:tau').dim_of().data()
        except:
            pass

        try:
            self.we = self.tree.getNode(r'\results::conf:we').data()
            self.t_we = self.tree.getNode(r'\results::conf:we').dim_of().data()
            self.t_wtot = self.tree.getNode(r'\results::total_energy').dim_of().data()
            self.wtot = self.tree.getNode(r'\results::total_energy').data()
            self.wi = self.tree.getNode(r'\results::conf:wi').data()
            self.t_wi = self.tree.getNode(r'\results::conf:wi').dim_of().data()
            self.h98 = self.tree.getNode(r'\results::conf:h_scal').data()
            self.t_h98 = self.tree.getNode(r'\results::conf:h_scal').dim_of().data()

        except:
            pass
        self.t_zaxis = self.tree.getNode(r'\results::z_axis').dim_of().data()
        self.zaxis = self.tree.getNode(r'\results::z_axis').data()
        self.t_kappa = self.tree.getNode(r'\results::kappa_edge').dim_of().data()
        self.kappa = self.tree.getNode(r'\results::kappa_edge').data()

    def get_ip(self):
        self.ip2 = self.tree.getNode(r'\results::i_p').data()
        self.tip2 = self.tree.getNode(r'\results::i_p').dim_of().data()

    # CARICA Beam

    def get_nbi(self):
        print('Loading NBI \n')
        try:
            self.time_beam = self.tree.getNode(
                '\ATLAS::NBH.DATA.MAIN_ADC:DATA').dim_of().data()
            self.tag_nbi = True
        except:
            self.tag_nbi = False
            return

        self.neu_eff = self.tree.getNode('\ATLAS::NBH.DATA.MODEL:NEU_EFF').data()

        self.t_neu_eff = self.tree.getNode(
            '\ATLAS::NBH.DATA.MODEL:NEU_EFF').dim_of().data()
        self.ref_energy = self.tree.getNode('\ATLAS::NBH.DATA.MODEL:ENERGY').data()
        self.ref_neutral_p = self.tree.getNode(
            '\ATLAS::NBH.DATA.MODEL:POWR_NEUTRAL').data()
        self.time_ref = self.tree.getNode(
            '\ATLAS::NBH.DATA.MODEL:POWR_NEUTRAL').dim_of().data()
        sig = self.tree.getNode('\ATLAS::NBH.DATA.MAIN_ADC:DATA ').data()
        self.neutral_p = sig[36, :]
        self.ion_p = sig[35, :]
        self.energy = sig[32, :]
        self.current = sig[34, :]
        self.ion_current = sig[33, :]
        self.magnet_i = sig[0, :]
        self.second_grid = sig[11, :]
        self.rf = sig[7, :]

    # CARICA ECRH

    def get_ecrh(self):
        print('Loading ECRH \n')
        try:
            self.time_ecrh = self.tree.getNode(
                '\RESULTS::TORAY.INPUT:P_GYRO').dim_of().data()
            self.ecrh = self.tree.getNode('\RESULTS::TORAY.INPUT:P_GYRO').data()
            self.tag_ecrh = True
        except:
            self.tag_ecrh = False

    # CARICA ECRH

    def get_ecrh2(self):
        print('Loading ECRH \n')
        try:
            conn = mds.Connection('tcvdata')
            t = conn.openTree('tcv_shot', self.shot)
            a1 = conn.get('\ECRH::ECRH_POWER:FORWARD:X2_1')
            a4 = conn.get('\ECRH::ECRH_POWER:FORWARD:X2_4')
            a5 = conn.get('\ECRH::ECRH_POWER:FORWARD:X2_5')
            self.l1_pow = a1.data()
            self.l4_pow = a4.data()
            self.l5_pow = a5.data()

            self.l1_t = conn.get('data(dim_of(\ATLAS::SYSTEM.ECRH.POWER:FORWARD:X2_1))')
            self.l4_t = conn.get('data(dim_of(\ATLAS::SYSTEM.ECRH.POWER:FORWARD:X2_4))')
            self.l5_t = conn.get('data(dim_of(\ATLAS::SYSTEM.ECRH.POWER:FORWARD:X2_5))')
            self.tag_ecrh2 = True

            fit4 = interpolate.interp1d(
                self.l4_t, self.l4_pow, fill_value='extrapolate')
            fit5 = interpolate.interp1d(
                self.l5_t, self.l5_pow, fill_value='extrapolate')

            self.ecrh_tot = fit4(self.l1_t)+fit5(self.l1_t)+self.l1_pow

        except:
            self.tag_ecrh2 = False

    # CARICA FOTODIODI PER HALFA E ALTRI INFLUSSI

    def get_ha(self):
        print('Loading Ha: \n')
        print('And correcting the offset at t<0')
        for i in range(18):
            call = 'pd_calibrated('+str(i+1)+')'
            print(call)
            sig = mds.Data.execute(call).data()

            if i == 0:
                self.time_ha = mds.Data.execute(call).dim_of().data()
                ntha = len(self.time_ha)
                sig_ha1 = np.zeros([ntha, 18])
                sig_ha2 = np.zeros([ntha, 18])

            sig_ha2[:, i] = mds.Data.execute(call).data()
            ind_zero = self.time_ha < -0.01

            sig_ha1[:, i] = sig_ha2[:, i]-np.mean(sig_ha2[ind_zero, i])

            self.sig_ha = sig_ha1
            self.ha_element = list(range(18))
            for jj in range(10):
                self.ha_element[jj] = 'H'
            self.ha_element[6] = 'H Hor'
            self.ha_element[10] = 'OIII'
            self.ha_element[11] = 'B+Bre'
            self.ha_element[12] = 'CIII'
            self.ha_element[13] = 'CIII+CIV'
            self.ha_element[14] = 'B+Bre'
            self.ha_element[15] = 'CII'
            self.ha_element[16] = 'OII'
            self.ha_element[17] = 'HeII'
            # 0 e 10-17 sono alto basso al centro della camera
            # 6 e all'equatore

        print('')
        print('')

    # CARICA FOTODIODO SINGOLO

    def get_ha_one(self, numero):
        print('Loading Ha: \n')
        call = 'pd_calibrated('+str(numero)+')'
        print(call)
        sig = mds.Data.execute(call).data()

        self.time_ha = mds.Data.execute(call).dim_of().data()
        self.sig_ha = mds.Data.execute(call).data()

    # CARICA FOTODIODO SINGOLO

    def get_ha_one2(self, numero):
        print('Loading Ha '+str(numero))
        conn = mds.Connection('tcvdata')
        tree_f = conn.openTree('tcv_shot', self.shot)
        if numero < 10:
            self.ha = conn.get('\BASE::PD:PD_00'+str(numero)).data()
            p = self.tree.getNode('\BASE::PD:PD_00'+str(numero))
            self.tha = p.getDimensionAt(0).data()
        if numero >= 10:
            self.ha = conn.get('\BASE::PD:PD_0'+str(numero)).data()
            p = self.tree.getNode('\BASE::PD:PD_0'+str(numero))
            self.tha = p.getDimensionAt(0).data()

    # calcola le medie con errori dei profili FIT del THOMSON da t1 a t2
    def medie_thomson(self, t1, t2):
        ind1 = (self.time_th < t2)*(self.time_th > t1)
        nrho = len(self.th_rho)

        media_te1 = np.zeros([nrho, 2])
        media_ne1 = np.zeros([nrho, 2])
        for i in range(nrho):
            media_te1[i, 0] = np.nanmean(self.th_t[i, ind1])
            media_te1[i, 1] = np.nanstd(self.th_t[i, ind1])
            media_ne1[i, 0] = np.nanmean(self.th_n[i, ind1])
            media_ne1[i, 1] = np.nanstd(self.th_n[i, ind1])

        media_te = media_te1
        media_ne = media_ne1
        return media_te, media_ne

    # calcola le medie con errori dei profili CXRS (TI e Vtor)
    def medie_cx(self, t1, t2):
        ind1 = (self.time_v < t2)*(self.time_v > t1)
        nrho = len(self.rho_v[0, :])

        media_ti1 = np.zeros([nrho, 2])
        media_vtor1 = np.zeros([nrho, 2])

        for i in range(nrho):
            media_ti1[i, 0] = np.nanmean(self.ti[ind1, i])
            media_ti1[i, 1] = np.nanstd(self.ti[ind1, i])
            media_vtor1[i, 0] = np.nanmean(self.vtor[ind1, i])
            media_vtor1[i, 1] = np.nanstd(self.vtor[ind1, i])

        media_ti = media_ti1
        media_vtor = media_vtor1
        return media_ti, media_vtor

    # calcola le medie con errori dei profili CXRS senza fit (TI e Vtor)
    def medie_cx_nofit(self, t1, t2, rho1, rho2):
        ind1 = (self.time_v_raw < t2)*(self.time_v_raw > t1)
        tion_r1 = self.ti_raw[ind1, :]
        vtor_r1 = self.vtor_raw[ind1, :]
        nion_r1 = np.array(self.ni_raw[ind1, :])
        rho_r1 = np.array(self.rho_v_raw[ind1, :])  # [t,r]
        time1 = self.time_v_raw[ind1]

        cc = np.argwhere(ind1 == True)
        nt = len(cc)
        erre = [0]
        v = [0]
        n = [0]
        t = [0]
        for i in range(nt):
            erre.extend(rho_r1[i, :])
            v.extend(vtor_r1[i, :])
            t.extend(tion_r1[i, :])
            n.extend(nion_r1[i, :])

        ntot = len(erre)
        erre = np.array(erre[1:ntot-1])
        v = np.array(v[1:ntot-1])
        n = np.array(n[1:ntot-1])
        t = np.array(t[1:ntot-1])

        ind_rho = (erre < rho2)*(erre > rho1)

        media_ti = np.zeros(2)
        media_vtor = np.zeros(2)
        media_rho = np.zeros(2)
        media_ni = np.zeros(2)

        media_ti[0] = np.nanmean(t[ind_rho])
        media_ti[1] = np.nanstd(t[ind_rho])
        media_rho[0] = np.nanmean(erre[ind_rho])
        media_rho[1] = np.nanstd(erre[ind_rho])
        media_ni[0] = np.nanmean(n[ind_rho]*1e-19)*1e19
        media_ni[1] = np.nanstd(n[ind_rho]*1e-19)*1e19
        media_vtor[0] = np.nanmean(v[ind_rho])
        media_vtor[1] = np.nanstd(v[ind_rho])

        return media_rho, media_vtor, media_ni, media_ti

    def get_fir(self):
        print('Loading FIR  \n')
        self.time_fir = self.tree.getNode(
            '\Results::fir_lin_int_dens_array').dim_of().data()
        self.fir = self.tree.getNode('\Results::fir_lin_int_dens_array').data()
        self.fir_r = self.tree.getNode('\diagz::fir_array:radii').data()
        self.fir_avg = self.tree.getNode(r'\results::fir:n_average').data()  # m^-3
        self.t_fir_avg = self.tree.getNode(
            r'\results::fir:n_average').dim_of().data()

    def get_cece(self):
        print('Loading CECE 6 channels  \n')
        self.time_cece = self.tree.getNode(
            r'\atlas::dt100_northwest_003:channel_001').dim_of().data()
        nsam_cece = len(self.time_cece)
        self.cece_sig = np.zeros((nsam_cece, 6))
        self.cece_sig[:, 0] = self.tree.getNode(
            r'\atlas::dt100_northwest_003:channel_001').data()
        self.cece_sig[:, 1] = self.tree.getNode(
            r'\atlas::dt100_northwest_003:channel_002').data()
        self.cece_sig[:, 2] = self.tree.getNode(
            r'\atlas::dt100_northwest_003:channel_003').data()
        self.cece_sig[:, 3] = self.tree.getNode(
            r'\atlas::dt100_northwest_003:channel_004').data()
        self.cece_sig[:, 4] = self.tree.getNode(
            r'\atlas::dt100_northwest_003:channel_005').data()
        self.cece_sig[:, 5] = self.tree.getNode(
            r'\atlas::dt100_northwest_003:channel_006').data()

    def get_thb(self):
        print('Loading THB RAW signals \n')
        conn = mds.Connection('tcvdata')
        tree_thb = conn.openTree('tcv_shot', self.shot)
        self.time_thb = np.array(
            conn.get(r'dim_of(\ATLAS::SYSTEM.GPI_THB.CHANNEL_005:LINE_501)'))

        nsam = len(self.time_thb)
        ll = ['501', '667', '706', '728']
        ch = np.array(range(8))
        self.sig_thb = np.zeros((nsam, 8, 4))
        for i in range(8):
            for j, line in enumerate(ll):
                self.sig_thb[:, i, j] = conn.get(
                    r'\ATLAS::SYSTEM.GPI_THB.CHANNEL_00'+str(i+1)+':LINE_'+line)

    def get_thb_single(self):
        print('Loading THB RAW signal 667 - 3 \n')
        conn = mds.Connection('tcvdata')
        tree_thb = conn.openTree('tcv_shot', self.shot)
        self.time_thb_s = np.array(
            conn.get(r'dim_of(\ATLAS::SYSTEM.GPI_THB.CHANNEL_003:LINE_667)'))
        self.sig_thb_s = conn.get(r'(\ATLAS::SYSTEM.GPI_THB.CHANNEL_007:LINE_667)')

    def get_thb_puffing(self):
        print('Loading THB puffing from GPI \n')
        c = mds.Connection('gpi_07')
        tree_puff = c.openTree('gpic', self.shot)
        self.f_ref = c.get('flowref').data()
        self.f_time = c.get('flowref').dim_of().data()*.5*1e-3
        self.f_meas = c.get('flowmeas').data()

    def get_toray(self):
        print('Loading TORAY \n')

        self.toray_x_xray = self.tree.getNode(
            r'\results::toray.output_x:x_ray').data()
        self.toray_x_yray = self.tree.getNode(
            r'\results::toray.output_x:y_ray').data()
        self.toray_x_zray = self.tree.getNode(
            r'\results::toray.output_x:z_ray').data()
        # W/m^3/Wtot Density of absorbed power normalized on total injected power
        self.toray_x_pdens = self.tree.getNode(
            r'\results::toray.output_x:pdens').data()
        # W/Wtot	Integral of power deposition normalized on total injected power
        self.toray_x_pint = self.tree.getNode(
            r'\results::toray.output_x:pint').data()
        # time,numero gyro da 0 a 8. Il 9 e il totale
        self.toray_x_ptot = self.tree.getNode(
            r'\results::toray.output_x:ptot').data()
        # A/Wtot	Integral of current drive normalized on total injected power as a function of rho
        self.toray_x_icdint = self.tree.getNode(
            r'\results::toray.output_x:icdint').data()
        self.toray_x_icddend = self.tree.getNode(
            r'\results::toray.output_x:icddens').data()
        self.toray_nbrays = self.tree.getNode(r'\results::toray.input:nbrays').data()

        self.torayNode = self.tree.getNode(r'\results::toray.output_x:pdens')
        self.toray_rho = self.torayNode.getDimensionAt(0).data()
        self.toray_time = self.torayNode.getDimensionAt(2).data()

        self.torayNode2 = self.tree.getNode(r'\results::toray.output_x:x_ray')
        self.toray_a = self.torayNode2.getDimensionAt(0).data()
        self.toray_b = self.torayNode2.getDimensionAt(1).data()
        self.toray_c = self.torayNode2.getDimensionAt(2).data()

    def get_crxs_sys(self):
        # print('Loading CRXS WITHOUT ANYTHING  \n')
        
        
        # for nsys in range(4):
        #     for quantity in ['vi', 'ni', 'ti']:
                
        #         node_path = r'\results::cxrs:acq_00' + f'{nsys+1}:{quantity.upper()}'
        #         node_path_err = r'\results::cxrs:acq_00' f'{nsys+1}:{quantity.upper()}:ERR'
        #         key = f'{quantity}_sys{nsys+1}'
        #         try:
        #             val = self.tree.getNode(node_path)
        #             val_err = self.tree.getNode(node_path_err)
                    
        #             setattr(self, f'{key}_node', val)
        #             setattr(self, f'{key}', val.data())
        #             # setattr(self, f'{key}_err', np.nan * np.ones_like(val.data()))
        #             setattr(self, f'{key}_err', val_err.data())
        #             setattr(self, f'time_{key}', val.getDimensionAt(1).data())
        #             setattr(self, f'rho_{key}', val.getDimensionAt(0).data())
                    
        #         except Exception as e:
        #             print(e)
        #             print(f'No "{quantity}_sys{nsys}" (on node path "f{node_path}") node for shot {self.shot}')
                    
  
        try:

            self.VtorNode_raw1 = self.tree.getNode(r'\results::cxrs:acq_001:VI')
            self.VtorNode_raw2 = self.tree.getNode(r'\results::cxrs:acq_002:VI')
            self.VtorNode_raw3 = self.tree.getNode(r'\results::cxrs:acq_003:VI')
            self.VtorNode_raw4 = self.tree.getNode(r'\results::cxrs:acq_004:VI')

            self.time_c1 = self.VtorNode_raw1.getDimensionAt(1).data()
            self.time_c2 = self.VtorNode_raw2.getDimensionAt(1).data()
            self.time_c3 = self.VtorNode_raw3.getDimensionAt(1).data()
            self.time_c4 = self.VtorNode_raw4.getDimensionAt(1).data()

            self.rho_c1 = self.VtorNode_raw1.getDimensionAt(0).data()
            self.rho_c2 = self.VtorNode_raw2.getDimensionAt(0).data()
            self.rho_c3 = self.VtorNode_raw3.getDimensionAt(0).data()
            self.rho_c4 = self.VtorNode_raw4.getDimensionAt(0).data()

            self.vtor_c1 = self.VtorNode_raw1.data()  # LFS
            self.vtor_c2 = self.VtorNode_raw2.data()  # HFS
            self.vtor_c3 = self.VtorNode_raw3.data()  # POLODAL VELOCITY
            self.vtor_c4 = self.VtorNode_raw4.data()  # POLODAL EDGE VELOCITY
            
            self.vi_c1_err = self.tree.getNode(r'\CXRS:ACQ_001:VI:ERR').data()
            self.vi_c2_err = self.tree.getNode(r'\CXRS:ACQ_002:VI:ERR').data()
            self.vi_c3_err = self.tree.getNode(r'\CXRS:ACQ_003:VI:ERR').data()
            self.vi_c4_err = self.tree.getNode(r'\CXRS:ACQ_004:VI:ERR').data()

            self.tion_c1 = self.tree.getNode(r'\results::cxrs:acq_001:TI').data()
            self.tion_c2 = self.tree.getNode(r'\results::cxrs:acq_002:TI').data()
            self.tion_c3 = self.tree.getNode(r'\results::cxrs:acq_003:TI').data()
            self.R1 = self.tree.getNode(r'\results::cxrs:acq_001:R').data()
            self.R2 = self.tree.getNode(r'\results::cxrs:acq_002:R').data()
            self.R3 = self.tree.getNode(r'\results::cxrs:acq_003:R').data()
            self.Z1 = self.tree.getNode(r'\results::cxrs:acq_001:Z').data()
            self.Z2 = self.tree.getNode(r'\results::cxrs:acq_002:Z').data()
            self.Z3 = self.tree.getNode(r'\results::cxrs:acq_003:Z').data()
   
        except TreeNNF as e:
            warnings.warn(f'{e}')

    # calcola le medie con errori dei profili CXRS senza fit (TI e Vtor)
    def medie_cx_sys(self, t1, t2, rho1, rho2):
        ind1 = (self.time_c1 < t2)*(self.time_c1 > t1)
        ind2 = (self.time_c2 < t2)*(self.time_c2 > t1)

        tion_r1 = self.tion_c1[ind1, :]
        vtor_r1 = self.vtor_c1[ind1, :]
        rho_r1 = np.array(self.rho_c1[ind1, :])  # [t,r]
        time1 = self.time_c1[ind1]

        tion_r2 = self.tion_c2[ind2, :]
        vtor_r2 = self.vtor_c2[ind2, :]
        rho_r2 = np.array(self.rho_c2[ind2, :])  # [t,r]
        time2 = self.time_c2[ind2]

        cc = np.argwhere(ind1 == True)
        nt = len(cc)
        erre1 = [0]
        v1 = [0]
        n1 = [0]
        t1 = [0]
        for i in range(nt):
            erre1.extend(rho_r1[i, :])
            v1.extend(vtor_r1[i, :])
            t1.extend(tion_r1[i, :])

        cc = np.argwhere(ind2 == True)
        nt = len(cc)
        erre2 = [0]
        v2 = [0]
        n2 = [0]
        t2 = [0]
        for i in range(nt):
            erre2.extend(rho_r2[i, :])
            v2.extend(vtor_r2[i, :])
            t2.extend(tion_r2[i, :])

        ntot = len(erre1)
        erre1 = np.array(erre1[1:ntot-1])
        v1 = np.array(v1[1:ntot-1])
        n1 = np.array(n1[1:ntot-1])
        t1 = np.array(t1[1:ntot-1])
        ind_rho = (erre1 < rho2)*(erre1 > rho1)
        media_ti1 = np.zeros(2)
        media_vtor1 = np.zeros(2)
        media_rho1 = np.zeros(2)
        media_ti1[0] = np.nanmean(t1[ind_rho])
        media_ti1[1] = np.nanstd(t1[ind_rho])
        media_rho1[0] = np.nanmean(erre1[ind_rho])
        media_rho1[1] = np.nanstd(erre1[ind_rho])
        media_vtor1[0] = np.nanmean(v1[ind_rho])
        media_vtor1[1] = np.nanstd(v1[ind_rho])

        ntot = len(erre2)
        erre2 = np.array(erre2[1:ntot-1])
        v2 = np.array(v2[1:ntot-1])
        n2 = np.array(n2[1:ntot-1])
        t2 = np.array(t2[1:ntot-1])
        ind_rho2 = (erre2 < rho2)*(erre2 > rho1)
        media_ti2 = np.zeros(2)
        media_vtor2 = np.zeros(2)
        media_rho2 = np.zeros(2)
        media_ti2[0] = np.nanmean(t2[ind_rho2])
        media_ti2[1] = np.nanstd(t2[ind_rho2])
        media_rho2[0] = np.nanmean(erre2[ind_rho2])
        media_rho2[1] = np.nanstd(erre2[ind_rho2])
        media_vtor2[0] = np.nanmean(v2[ind_rho2])
        media_vtor2[1] = np.nanstd(v2[ind_rho2])
        # plt.plot(erre2[ind_rho])
        return media_rho1, media_vtor1, media_ti1, media_rho2, media_vtor2, media_ti2

    def get_xtomo(self, camera):

        print('Loading XTOMO CAMERA '+str(camera))

        Channels = self.tree.getNode(
            r'\base::xtomo:array_{:03}:source'.format(camera))
        mat = scipy.io.loadmat('/home/agostini/tcv/cat_defaults2008.mat')

        index = np.arange(20)
        self.x_chord = mat['xchord'][:, (camera - 1) * 20 + index] / 100.
        self.y_chord = mat['ychord'][:, (camera - 1) * 20 + index] / 100.

        for ind, channel in enumerate(Channels):
            if ind == 0:
                self.time_xtomo = self.tree.getNode(channel).dim_of().data()
                nsam = len(self.time_xtomo)
                nsig = len(Channels)
                matrice = np.zeros((nsam, nsig))

            sig = self.tree.getNode(channel).data()
            matrice[:, ind] = sig

        self.sig_xtomo = matrice

    # Reciprocating probe equatore esterno

    def get_rcp(self):
        print('Loading Langmuir RCP from pulsfile')

        self.rcp_ne1 = self.tree.getNode('\RESULTS::FP:NE_1').data()
        self.rcp_rho1 = self.tree.getNode('\RESULTS::FP:RHO_1').data()
        self.rcp_ne2 = self.tree.getNode('\RESULTS::FP:NE_2').data()
        self.rcp_rho2 = self.tree.getNode('\RESULTS::FP:RHO_2').data()
        self.rcp_te1 = self.tree.getNode('\RESULTS::FP:TE_1').data()
        self.rcp_te2 = self.tree.getNode('\RESULTS::FP:TE_2').data()
        self.rcp_time1 = self.tree.getNode('\RESULTS::FP:T_1').data()
        self.rcp_time2 = self.tree.getNode('\RESULTS::FP:T_2').data()

    def get_bolo(self, los=0):

        self.xbolo = np.loadtxt('/home/agostini/tcv/lib/bolo_x.txt')
        self.ybolo = np.loadtxt('/home/agostini/tcv/lib/bolo_y.txt')
        if los == 0:
            print(' Loading BOLO all data')
            self.bolo_sig = self.tree.getNode('\BASE::BOLO_U:AMPLITUDE:RAW').data()
            self.bolo_time = self.tree.getNode('\BASE::BOLO_U:DIM').dim_of().data()

        if los != 0:
            print(' Loading BOLO LOS '+str(los))
            self.bolo_time = self.tree.getNode('\BASE::BOLO_U:DIM').dim_of().data()
            if los <= 99:
                self.bolo_sig = self.tree.getNode(
                    '\ATLAS::SYSTEM.BOLO.AMPLITUDE:CHANNEL_0'+str(los)).data()
            if los > 99:
                self.bolo_sig = self.tree.getNode(
                    '\ATLAS::SYSTEM.BOLO.AMPLITUDE:CHANNEL_'+str(los)).data()

    # Sonde langmuir della prima parete
    # Jsat, te, ne, posizione(R,Z)
    # Jsat e ne,Te hanno campionamenti e asse dei tempi diversi

    def get_wall_probes(self):

        print('')
        print(' Loading wall probes')

        try:
            self.prob = self.tree.getNode('\RESULTS::LANGMUIR:PROBES').data()
            self.pos = self.tree.getNode('\RESULTS::LANGMUIR:POS').data()
            p = self.tree.getNode('\RESULTS::LANGMUIR:RHO_PSI')
            self.trhopsi = p.getDimensionAt(1).data()
            self.rhopsi = self.tree.getNode('\RESULTS::LANGMUIR:RHO_PSI').data()

            self.jsat = self.tree.getNode('\RESULTS::LANGMUIR:JSAT2').data()
            p = self.tree.getNode('\RESULTS::LANGMUIR:JSAT2')
            self.tjsat = p.getDimensionAt(1).data()
            self.rhopsi2 = self.tree.getNode('\RESULTS::LANGMUIR:RHO_PSI2').data()

            self.dens = self.tree.getNode('\RESULTS::LANGMUIR:FOUR_PAR_FIT:dens').data()
            p = self.tree.getNode('\RESULTS::LANGMUIR:FOUR_PAR_FIT:dens')
            self.tdens = p.getDimensionAt(1).data()
            self.tel = self.tree.getNode('\RESULTS::LANGMUIR:FOUR_PAR_FIT:te').data()
        except:
            print('')
            print(' Wall probes not found ')
            self.prob = np.array([0])

    # richiesa e misura di flusso (in Volt) del fuelling per le 3 valvole

    def get_fuelling(self):
        print(' Loading fuelling...')
        conn = mds.Connection('tcvdata')
        tree_f = conn.openTree('tcv_shot', self.shot)
        self.piezo1_r = conn.get(
            r'\ATLAS::DT4G_MIX_002:CHANNEL_050').data() 	# request Piezo1
        self.piezo1_m = conn.get(
            r'\ATLAS::DT4G_MIX_002:CHANNEL_052').data()  # measure
        self.piezo2_r = conn.get(
            r'\ATLAS::DT4G_MIX_002:CHANNEL_053') .data() 	# request P2
        self.piezo2_m = conn.get(
            r'\ATLAS::DT4G_MIX_002:CHANNEL_055').data()  # measure
        self.piezo3_r = conn.get(
            r'\ATLAS::DT4G_MIX_002:CHANNEL_054').data()  	# request P2
        self.piezo3_m = conn.get(
            r'\ATLAS::DT4G_MIX_002:CHANNEL_056').data()  # measure

        p = self.tree.getNode(r'\ATLAS::DT4G_MIX_002:CHANNEL_050')
        self.time_piezo = p.getDimensionAt(0).data()



# %%
