"""
zhao_core.py — Python port of the φC31 integrase/RDF recombination model
from Zhao, Pokhilko, Ebenhöh, Rosser & Colloms, NAR 2019 (47:4896–4909)
"A single-input binary counting module based on serine integrase
 site-specific recombination" (Supplementary Data, eqs. 1–38).

Faithful translation of Model_tetR_110219.m (github.com/alex297/
model-of-binary-counter-based-on-recombination-with-serine-integrase).

State vector y (38 ODEs), units: µM, hours.
  y[0]  LR                         y[19] int2-rdf
  y[1]  int                        y[20] PB-int2-rdf
  y[2]  int2-rdf2                  y[21] LR-int2-rdf
  y[3]  int2                       y[22] rdf
  y[4]  PB-int2                    y[23] PB
  y[5]  LR-int2                    y[24] PB-int4-rdf
  y[6]  LR-int4                    y[25] PB-int4-rdf2
  y[7]  PB-int4                    y[26] PB-int4-rdf3
  y[8]  LR-int4 synapse second     y[27] LR-int4-rdf
  y[9]  PB-int4 synapse            y[28] LR-int4-rdf2
  y[10] LR-int4 synapse first      y[29] LR-int4-rdf3
  y[11] PB-int2-rdf2               y[30] PB-int6i
  y[12] LR-int2-rdf2               y[31] PB-int6-rdf4i
  y[13] PB-int4-rdf4               y[32] REP mRNA   (tetR or bm3r1 mRNA; was y(33) RDF mRNA — see below)
  y[14] LR-int4-rdf4               y[33] PB-int6-rdfi
  y[15] PB-int4-rdf4 synapse 2nd   y[34] PB-int6-rdf2i
  y[16] PB-int4-rdf4 synapse 1st   y[35] PB-int6-rdf3i
  y[17] LR-int4-rdf4 synapse       y[36] REP mRNA
  y[18] int-rdf                    y[37] REP protein

NOTE on indices: MATLAB y(33) is RDF mRNA and y(37)/y(38) are tetR mRNA/protein.
Here: y[32] = RDF mRNA, y[36] = repressor mRNA, y[37] = repressor protein.
"""

import numpy as np

# ---------------------------------------------------------------- parameters
def default_params():
    """Default parameter set from tetR_110219.m / Supplementary Table S1."""
    P = dict(
        # dissociation/equilibrium constants (µM), Table S1 & Pokhilko 2016 (ref 8)
        Kr1=1.0, Kr2=1.0, Ki=0.02, Kii=0.3, Kir=0.05,
        Kmod=3.4, Kmodr=1.9, Ks01=0.001, Ks02=0.007,
        Ks1=0.1, Ks2=0.12, Ks3=0.1, Ks4=0.013,
        Kb1=0.02, Kb2=0.01, Kb3=0.025, Kb4=0.05,
        Dtot=0.017,           # total switch-plasmid DNA, µM (10 copies after division)
        # v49c A-module growth curve: nominal Td=50 min.
        # k_dil = ln(2) / (50/60 h) = 0.8317766 h^-1.
        # Keep every B-module growth-dilution term on this same time basis.
        k_dil=0.8317766,      # dilution by growth, h^-1 (50 min doubling)
        k_tscr=120.0,         # transcription rate constant, h^-1
        k_rna=4.0,            # mRNA degradation, h^-1 (~10 min half-life)
        # --- delay-circuit repressor (TetR in the paper; BM3R1 in TEMPO) ---
        K_rep=0.01,           # repression threshold Ktet = 0.01 µM (fitted)
        n_rep=2.0,            # Hill coefficient (n=2 in the paper's eq. 36)
        krep_tsl=0.3,         # repressor translation rate, h^-1 (ktet_tsl, fitted)
        leak_rep=0.0,         # explicit leak floor of repressed promoter (paper: 0)
        krdf_tsl=4.0,         # RDF translation rate, h^-1 (fitted)
        k_tag_int=0.0,        # extra degradation of Int-containing species, h^-1
                              # (ssrA degradation tag on the integrase; tagged
                              #  Int in DNA complexes releases free DNA, same
                              #  bookkeeping as dilution)
        # --- integrase input (square arabinose pulse, for verification) ---
        k_int=3.0,            # Int production rate during pulse, µM/h (fitted)
        ara_on=0.5, ara_off=0.7, period=24.0, kt=0.3,
    )
    return P


def _rate_constants(P):
    """Derived rate constants, exactly as in Model_tetR_110219.m lines 48-76."""
    kp = 60 * 60
    C = dict(
        kps01=40 * 60, kps02=40 * 60, kp=kp,
        kmii=P['Kii'] * kp, kmir=P['Kir'] * kp,
        kmb1=P['Kb1'] * kp, kmb2=10 * 60,
        kpb2=P['Kb2'] * (10 * 60), kmb3=P['Kb3'] * kp, kpb4=P['Kb4'] * kp,
        kpi=3 * 60, kmi=P['Ki'] * (3 * 60),
        kps1=0.8 * 60, kps3=0.8 * 60,
        kms01=P['Ks01'] * (40 * 60), kms02=P['Ks02'] * (40 * 60),
        kms1=P['Ks1'] * (0.8 * 60), kms3=P['Ks3'] * (0.8 * 60),
        kms2=0.0001 * 60, kps2=P['Ks2'] * (0.0001 * 60),
        kms4=0.005 * 60, kps4=P['Ks4'] * (0.005 * 60),
        kpr=1 * 60, kmr1=(1 * 60) / P['Kr1'], kmr2=(1 * 60) / P['Kr2'],
        kpmod=1 * 60, kpmodr=1 * 60,
        kmmod=(1 * 60) / P['Kmod'], kmmodr=(1 * 60) / P['Kmodr'],
    )
    return C


# ---------------------------------------------------------------- input pulse
def square_pulse(t, P):
    """Dimensionless arabinose pulse, eq. for ara(t) in Supplementary Data."""
    tm = t - P['period'] * np.floor(t / P['period'])
    return 0.5 * (np.tanh((tm - P['ara_on']) / P['kt'])
                  - np.tanh((tm - P['ara_off']) / P['kt']))


# ---------------------------------------------------------------- RHS
def rhs(t, y, P, C, int_production):
    """
    int_production: callable(t) -> µM/h, production rate of free integrase.
    """
    F = np.zeros(38)
    (kps01, kps02, kp, kmii, kmir, kmb1, kmb2, kpb2, kmb3, kpb4, kpi, kmi,
     kps1, kps3, kms01, kms02, kms1, kms3, kms2, kps2, kms4, kps4,
     kpr, kmr1, kmr2, kpmod, kpmodr, kmmod, kmmodr) = (
        C['kps01'], C['kps02'], C['kp'], C['kmii'], C['kmir'], C['kmb1'],
        C['kmb2'], C['kpb2'], C['kmb3'], C['kpb4'], C['kpi'], C['kmi'],
        C['kps1'], C['kps3'], C['kms01'], C['kms02'], C['kms1'], C['kms3'],
        C['kms2'], C['kps2'], C['kms4'], C['kps4'], C['kpr'], C['kmr1'],
        C['kmr2'], C['kpmod'], C['kpmodr'], C['kmmod'], C['kmmodr'])
    k_dil = P['k_dil']
    kd_int = k_dil + P.get('k_tag_int', 0.0)  # Int-containing species

    Lr_t = (y[0] + y[5] + y[7] + y[8] + y[10] + y[12] + y[14] + y[17]
            + y[21] + y[27] + y[28] + y[29])
    Bp_t = (y[23] + y[4] + y[6] + y[9] + y[11] + y[13] + y[15] + y[16]
            + y[20] + y[24] + y[25] + y[26] + y[30] + y[31] + y[33]
            + y[34] + y[35])

    F[0] = (kpb2 * y[5] - kmb2 * y[0] * y[3] + kmb3 * (y[12] + y[21])
            - kp * y[0] * (y[2] + y[19]) + kd_int * (Lr_t - y[0]))
    F[1] = (int_production(t) + kmii * (2 * y[3] + y[19]) + kmir * y[18]
            - kp * (2 * y[1] * y[1] + y[1] * y[18] + y[1] * y[22])
            - kd_int * y[1])
    F[2] = (kp * (y[19] * y[22] + y[18] * y[18]) - (kmir + kmii) * y[2]
            - kp * y[2] * (y[23] + y[11] + y[0] + y[12])
            - kps01 * y[2] * (y[5] + y[21]) - kps02 * y[2] * (y[4] + y[20])
            + kpb4 * (y[11] + y[13]) + kms01 * (y[28] + y[29])
            + kms02 * (y[25] + y[26]) + kmb3 * (y[12] + y[14]) - kd_int * y[2])
    F[3] = (kp * y[1] * y[1] - kmii * y[3] - kp * y[3] * (y[22] + y[23] + y[4])
            - kmb2 * y[3] * (y[0] + y[5]) + kpb2 * (y[5] + y[7])
            - kps01 * y[3] * (y[12] + y[21]) - kps02 * y[3] * (y[11] + y[20])
            + kmir * y[19] + kmb1 * (y[4] + y[6]) + kms01 * (y[27] + y[28])
            + kms02 * (y[24] + y[25])
            - (kpi * y[3] * (y[6] + y[13] + y[24] + y[25] + y[26])
               - kmi * (y[30] + y[31] + y[33] + y[34] + y[35])) - kd_int * y[3])
    F[4] = (kp * y[23] * y[3] - kmb1 * y[4] - kp * y[4] * y[3]
            - kps02 * y[4] * (y[19] + y[2]) + kmb1 * y[6]
            + kms02 * (y[24] + y[25]) - kp * y[22] * y[4] + kmir * y[20]
            - kd_int * y[4])
    F[5] = (kmb2 * y[3] * (y[0] - y[5]) - kpb2 * (y[5] - y[7])
            - kps01 * y[5] * (y[19] + y[2]) + kms01 * (y[27] + y[28])
            - kp * y[22] * y[5] + kmir * y[21] - kd_int * y[5])
    F[6] = (kp * y[4] * y[3] - kmb1 * y[6] - kps1 * y[6] + kms1 * y[9]
            - (kpi * y[3] * y[6] - y[30] * kmi) - kd_int * y[6])
    F[7] = (kmb2 * y[5] * y[3] - kpb2 * y[7] - kms2 * y[7] + kps2 * y[8]
            - kd_int * y[7])
    F[8] = (kpmod * y[10] - kmmod * y[8] - kps2 * y[8] + kms2 * y[7]
            - kd_int * y[8])
    F[9] = (kps1 * y[6] - kms1 * y[9] - kpr * y[9] + kmr1 * y[10]
            - kd_int * y[9])
    F[10] = (kpr * y[9] - kmr1 * y[10] - kpmod * y[10] + kmmod * y[8]
             - kd_int * y[10])
    F[11] = (kp * y[23] * y[2] - kpb4 * y[11] + kp * y[20] * y[22]
             - kmir * y[11] - kp * y[11] * y[2] - kps02 * y[11] * (y[3] + y[19])
             + kms02 * (y[25] + y[26]) + kpb4 * y[13] - kd_int * y[11])
    F[12] = (kp * y[0] * y[2] - kmb3 * y[12] + kp * y[21] * y[22]
             - kmir * y[12] - kp * y[12] * y[2] - kps01 * y[12] * (y[3] + y[19])
             + kms01 * (y[28] + y[29]) + kmb3 * y[14] - kd_int * y[12])
    F[13] = (kp * y[2] * y[11] - kpb4 * y[13] + kps4 * y[15] - kms4 * y[13]
             - (kpi * y[3] * y[13] - y[31] * kmi) - kd_int * y[13])
    F[14] = (kp * y[2] * y[12] - kmb3 * y[14] - kps3 * y[14] + kms3 * y[17]
             - kd_int * y[14])
    F[15] = (kpmodr * y[16] - kmmodr * y[15] - kps4 * y[15] + kms4 * y[13]
             - kd_int * y[15])
    F[16] = (kpr * y[17] - kmr2 * y[16] + kmmodr * y[15] - kpmodr * y[16]
             - kd_int * y[16])
    F[17] = (kps3 * y[14] - kms3 * y[17] - kpr * y[17] + kmr2 * y[16]
             - kd_int * y[17])
    F[18] = (kp * (y[1] * y[22] - y[1] * y[18] - 2 * y[18] * y[18])
             - kmir * y[18] + kmii * (y[19] + 2 * y[2]) - kd_int * y[18])
    F[19] = (kp * (y[3] * y[22] + y[1] * y[18]) - (kmir + kmii) * y[19]
             - kp * y[22] * y[19] + kmir * y[2] - kp * y[19] * (y[23] + y[0])
             - kps01 * y[19] * (y[5] + y[21] + y[12])
             - kps02 * y[19] * (y[4] + y[20] + y[11]) + kpb4 * y[20]
             + kms01 * (y[27] + y[28] + y[29]) + kms02 * (y[24] + y[25] + y[26])
             + kmb3 * y[21] - kd_int * y[19])
    F[20] = (kp * y[23] * y[19] - kpb4 * y[20]
             - kps02 * y[20] * (y[3] + y[19] + y[2])
             + kms02 * (y[24] + y[25] + y[26]) + kp * y[22] * y[4]
             - kmir * y[20] - kp * y[20] * y[22] + kmir * y[11] - kd_int * y[20])
    F[21] = (kp * y[0] * y[19] - kmb3 * y[21]
             - kps01 * y[21] * (y[19] + y[2] + y[3])
             + kms01 * (y[27] + y[28] + y[29]) + kp * y[22] * y[5]
             - kmir * y[21] - kp * y[21] * y[22] + kmir * y[12] - kd_int * y[21])
    F[22] = (P['krdf_tsl'] * y[32]
             + kmir * (y[18] + y[19] + y[2] + y[20] + y[21] + y[11] + y[12])
             - kp * y[22] * (y[1] + y[3] + y[19] + y[4] + y[5] + y[20] + y[21])
             - k_dil * y[22])
    F[23] = (kmb1 * y[4] + kpb4 * (y[11] + y[20])
             - kp * y[23] * (y[3] + y[19] + y[2]) + kd_int * (Bp_t - y[23]))
    F[24] = (kps02 * (y[3] * y[20] + y[19] * y[4]) - y[24] * 2 * kms02
             - (kpi * y[3] * y[24] - y[33] * kmi) - kd_int * y[24])
    F[25] = (kps02 * (y[3] * y[11] + y[19] * y[20] + y[2] * y[4])
             - y[25] * 3 * kms02 - (kpi * y[3] * y[25] - y[34] * kmi)
             - kd_int * y[25])
    F[26] = (kps02 * (y[19] * y[11] + y[2] * y[20]) - y[26] * 2 * kms02
             - (kpi * y[3] * y[26] - y[35] * kmi) - kd_int * y[26])
    F[27] = (kps01 * (y[3] * y[21] + y[19] * y[5]) - y[27] * 2 * kms01
             - kd_int * y[27])
    F[28] = (kps01 * (y[3] * y[12] + y[19] * y[21] + y[2] * y[5])
             - y[28] * 3 * kms01 - kd_int * y[28])
    F[29] = (kps01 * (y[19] * y[12] + y[2] * y[21]) - y[29] * 2 * kms01
             - kd_int * y[29])
    F[30] = kpi * y[3] * y[6] - y[30] * kmi - kd_int * y[30]
    F[31] = kpi * y[3] * y[13] - y[31] * kmi - kd_int * y[31]
    # y[32]: RDF mRNA — transcribed from LR DNA, repressed by REP (Hill)
    rep = y[37]
    hill = P['leak_rep'] + (1.0 - P['leak_rep']) / (
        1.0 + (rep / P['K_rep']) ** P['n_rep'])
    F[32] = P['k_tscr'] * Lr_t * hill - P['k_rna'] * y[32]
    F[33] = kpi * y[3] * y[24] - y[33] * kmi - kd_int * y[33]
    F[34] = kpi * y[3] * y[25] - y[34] * kmi - kd_int * y[34]
    F[35] = kpi * y[3] * y[26] - y[35] * kmi - kd_int * y[35]
    # y[36]: repressor (TetR/BM3R1) mRNA — transcribed from PB DNA
    F[36] = P['k_tscr'] * Bp_t - P['k_rna'] * y[36]
    # y[37]: repressor protein
    F[37] = P['krep_tsl'] * y[36] - k_dil * y[37]
    return F


# ---------------------------------------------------------------- observables
_LR_IDX = [0, 5, 7, 8, 10, 12, 14, 17, 21, 27, 28, 29]
_PB_IDX = [4, 6, 9, 11, 13, 15, 16, 20, 23, 24, 25, 26, 30, 31, 33, 34, 35]

def LR_total(Y):
    return Y[:, _LR_IDX].sum(axis=1)

def PB_total(Y):
    return Y[:, _PB_IDX].sum(axis=1)

def int_total(Y):
    return (Y[:, 1] + Y[:, 18]
            + 2 * (Y[:, 2] + Y[:, 19] + Y[:, 20] + Y[:, 21] + Y[:, 3]
                   + Y[:, 4] + Y[:, 5] + Y[:, 11] + Y[:, 12])
            + 4 * (Y[:, 6] + Y[:, 7] + Y[:, 8] + Y[:, 9] + Y[:, 10]
                   + Y[:, 13] + Y[:, 14] + Y[:, 15] + Y[:, 16] + Y[:, 17]
                   + Y[:, 27] + Y[:, 28] + Y[:, 29] + Y[:, 24] + Y[:, 25]
                   + Y[:, 26])
            + 6 * (Y[:, 30] + Y[:, 31] + Y[:, 33] + Y[:, 34] + Y[:, 35]))

def rdf_total(Y):
    return (Y[:, 22] + Y[:, 18] + Y[:, 19] + Y[:, 20] + Y[:, 21] + Y[:, 24]
            + Y[:, 27] + Y[:, 33]
            + 2 * (Y[:, 2] + Y[:, 11] + Y[:, 12] + Y[:, 25] + Y[:, 28]
                   + Y[:, 34])
            + 3 * (Y[:, 26] + Y[:, 29] + Y[:, 35])
            + 4 * (Y[:, 13] + Y[:, 14] + Y[:, 15] + Y[:, 16] + Y[:, 17]
                   + Y[:, 31]))


# ---------------------------------------------------------------- initial states
def y0_PB(P, rep_mrna=0.4, rep=0.06):
    """100 % PB, repressor near steady state (README, Fig. 5E init)."""
    y = np.zeros(38)
    y[23] = P['Dtot']
    y[36] = rep_mrna
    y[37] = rep
    return y

def y0_PB_ss(P):
    """100 % PB, closer to steady state (README, Fig. S11A init)."""
    y = np.zeros(38)
    y[22] = 0.01          # rdf
    y[23] = P['Dtot']
    y[32] = 0.005         # rdf mRNA
    y[36] = 0.4           # rep mRNA
    y[37] = 0.05          # rep
    return y

def y0_LR_ss(P):
    """100 % LR, closer to steady state (README, Fig. S11B init)."""
    y = np.zeros(38)
    y[0] = P['Dtot']
    y[22] = 0.5           # rdf
    y[32] = 0.25          # rdf mRNA
    y[36] = 0.06
    y[37] = 0.01
    return y
