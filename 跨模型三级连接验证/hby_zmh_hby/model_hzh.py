"""HB Y bit0 -> ZMH carry0/bit1 -> HBY carry1/bit2. Time in hours."""
from pathlib import Path
import sys
import numpy as np

PARENT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(PARENT))
from hybrid_model import HbyReceiver,donor_namespace,hill,BIT_NAMES,TAIL_NAMES

NAMES=tuple('b0_'+n for n in BIT_NAMES)+('A0_zmh','F0_zmh','pb1_zmh','I1_zmh','T1_zmh','RDF1_zmh')+TAIL_NAMES
IDX={n:i for i,n in enumerate(NAMES)}


class HZHModel:
    names=NAMES
    def __init__(self,mode,han=None):
        if mode not in ('square','han'):raise ValueError(mode)
        if mode=='han' and han is None:raise ValueError('Han upstream required')
        self.mode,self.han=mode,han
        self.tail=HbyReceiver()
        self.p=self.tail.p
        self.z=donor_namespace()['p'].copy()
        self.square=dict(start_h=5.,period_h=10.,width_h=100/60,amplitude_au_per_h=6.)
        if han is not None:
            hp=han.p
            # Prescribed upstream promoter; M_I0 is the shared C31 mRNA itself.
            self.h31=np.array([han.module.promoter_activities(row[1],row[3],row[5],hp)[3]
                               for row in han.sol.y.T])
            self.copies_per_au=602.214076*self.tail.c.receiver_uM_per_au

    def input_source(self,t):
        """Target immature-protein source, used for reporting/input comparison."""
        if self.mode=='square':
            q=self.square
            on=t>=q['start_h'] and ((t-q['start_h'])%q['period_h'])<q['width_h']
            return q['amplitude_au_per_h'] if on else 0.
        return 60*self.han.interp(60*t,self.han.flux_copies_min)/self.copies_per_au

    def bit0_rhs(self,t,b):
        p,c=self.p,self.tail.c
        m,iu,i,mt,tu,tr,mr,ru,r,comp,s=b
        lm,km=self.tail.lmb,self.tail.kmb
        d=np.zeros(11)
        if self.mode=='square':
            d[0]=self.input_source(t)*lm/c.translation_h-lm*m
            translation=c.translation_h*m
        else:
            hp=self.han.p
            tx=60*hp.oscillator_plasmid_copies*hp.c31_tx_per_plasmid_per_min*self.han.interp(60*t,self.h31)/self.copies_per_au
            loss=60*(hp.c31_mrna_intrinsic_loss_per_min+hp.mu)
            d[0]=tx-loss*m
            translation=60*hp.c31_translation_per_mrna_per_min*m
        d[1]=translation-km*iu
        d[2]=km*iu-p['gamma_int'][0]*i
        d[3:6]=self.tail.expression(mt,tu,tr,p['alpha_rep']*(1-s),p['gamma_rep'],lm,km)
        d[6:9]=self.tail.expression(mr,ru,r,p['alpha_rdf']*s*(1-hill(tr,p['K_rep'],p['n_rep'])),p['gamma_rdf'],lm,km)
        bind,un=c.kon*i*r,c.koff*comp
        d[2]+=-bind+un;d[8]+=-bind+un
        d[9]=bind-un-c.complex_decay*comp
        vf=p['k_fwd']*hill(i,p['K_D_int'][0],2)*p['K_inh']/(p['K_inh']+max(r,0))
        vr=p['k_rev']*hill(comp,self.tail.k_complex,2)
        d[10]=vf*(1-s)-vr*s
        return d

    def middle_rhs(self,b0s,y):
        p=self.z
        a,f,pb,i,tr,r=y
        act=a**p['n_A0']/(p['K_A0']**p['n_A0']+a**p['n_A0'])
        gate=p['K_R0']**p['n_R0']/(p['K_R0']**p['n_R0']+f**p['n_R0'])
        ia=i**p['n_int1']/(p['K_D_int1']**p['n_int1']+i**p['n_int1'])
        vf=p['k_fwd']*pb*ia*(p['Kinh']/(p['Kinh']+r))
        vr=p['k_rev']*(1-pb)*((i*r)**2/(p['K_D_comp']**2+(i*r)**2))
        return 60*np.array([p['alpha_A0']*(1-b0s)-p['gamma_A0']*a,
             p['alpha_R0']*act-p['gamma_R0']*f,-vf+vr,
             p['alpha_Int1']*act*gate-p['gamma_int1']*i,
             p['alpha_rep1']*pb-p['gamma_rep1']*tr,
             p['alpha_rdf1']*(1-pb)*(1.0/(1.0+(tr/p['K_rep'])**p['n']))-p['gamma_rdf1']*r])

    def initial_state(self):
        b0=self.tail.initial_state()[6:].copy() # same Rep initial condition, zero Int/RDF/S
        if self.mode=='han':b0[0]=self.han.sol.y[6,0]/self.copies_per_au
        p=self.z;a=p['alpha_A0']/p['gamma_A0']
        act=a**p['n_A0']/(p['K_A0']**p['n_A0']+a**p['n_A0'])
        middle=[a,p['alpha_R0']*act/p['gamma_R0'],1.,0.,p['alpha_rep1']/p['gamma_rep1'],0.]
        return np.r_[b0,middle,self.tail.initial_state()]

    def rhs(self,t,y):
        if np.shape(y)!=(34,):raise ValueError('New HZH layout requires 34 states')
        return np.r_[self.bit0_rhs(t,y[:11]),self.middle_rhs(y[10],y[11:17]),
                     self.tail.rhs(t,y[17:],pb1=y[13],int0=y[2])]

    def diagnostic_projection(self,y):
        """Named conversion to the preliminary readout's old 27-state contract."""
        from hybrid_model import INDEX
        out=np.zeros((27,y.shape[1]))
        out[INDEX['pb0']]=1-y[IDX['b0_S']]
        out[INDEX['pb1']]=y[IDX['pb1_zmh']]
        out[INDEX['int0']]=y[IDX['b0_I']]
        out[INDEX['b2_S']]=y[IDX['b2_S']]
        return out
