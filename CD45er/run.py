import sys
import pysam
import pandas as pd
import numpy as np
from scipy.stats import norm, t

import matplotlib.pyplot as plt
import seaborn as sns

import pickle

isoforms = np.array(['RX','RAX','RBX','RABX','RBCX','RABCX'])
segments = np.array(['R','A','B','C','X'])

# Iterates over a 10x bam file and compiles isoform compatibility and position for every read
def annotate_reads(bam,E):

    # Maps genomic position to exon number
    pos_min = E['st'].min()
    pos_max = E['en'].max()
    chrom = str(E['chr'].unique()[0])
    exons = np.zeros(pos_max-pos_min).astype(int)
    for ind,g in E.iterrows():
        exons[(g['st']-pos_min):(g['en']-pos_min)]=ind+1

    # Begin read annotation
    samfile = pysam.AlignmentFile(bam, "rb")
    res = list()

    for read in samfile.fetch(chrom,pos_min,pos_max):
        tags = dict(read.get_tags())
        # Ignore reads without an assigned cell barcode or UMI
        if (not 'CB' in tags.keys()) or (not 'UB' in tags.keys()):
            continue
        else:
            cb = tags['CB']
            ub = tags['UB']
    
        exon_seen = np.zeros(33)

        # Record last position in an exon supported by a read
        last_pos = -1

        for pos in read.get_reference_positions():
            # If not currently in an exon, continue
            if (pos<(pos_min-1)) or (pos>=(pos_max-1)):
                continue

            # Look up which exon we're in and record
            exon=exons[pos-(pos_min-1)]        
            if exon>0:
                exon_seen[exon-1]=1
                last_pos=pos

        # If no exons covered, ignore
        if sum(exon_seen)==0:
            continue
        
        # Create a string of blocks R(exon3)-A-B-C-X(exon7) of exons supported by read
        segments_seen = "".join(segments[exon_seen[2:7]==1])
    

        # Record in which exon the read ended, and at what position
        last_exon = np.where(exon_seen)[0][-1]
        exon_pos=last_pos-E.loc[last_exon,'st']
    
        res.append(pd.Series({'readname':read.qname,'CB':cb,'UB':ub,'seg':segments_seen,'end':last_pos,
                          'last_exon':last_exon,'exon_pos':exon_pos}))
    
    samfile.close()
    
    df = pd.concat(res,axis=1).T
    
    return(df)

# Giving a dataframe of annotated reads, performs EM to estimate 
#  - the isoform probability of every transcript (UMI)
#  - the overall frequency of different isoforms (pi)
#  - The fragment length distribution parameters mu/sig2 (only if not provided)
def isoform_quantification(df,E,mu=None,sig2=None,max_fraglen=500):

    exon_map = {'A':4,'B':5,'C':6}

    # Binary matrix of exon inclusion
    # exon x isoform matrix
    exon_inc = np.ones((E.shape[0],len(isoforms)))
    exon_inc[3:6] = 0
    for i in range(0,len(isoforms)):
        iso = isoforms[i].replace("R","").replace("X","")
        for exon in iso:
            exon_inc[exon_map[exon]-1,i] = 1

    # Pre-compute exon x isoform matrix of transcript lengths to each exon
    shifts = np.zeros((E.shape[0],len(isoforms)))

    for i in range(0,len(isoforms)): 
        elen = E['len'].values.copy()
        elen[exon_inc[:,i]==0]=0
        shifts[:,i] = [0] + list(elen.cumsum()[0:-1])
    
    exon_last_use = np.where(shifts.min(axis=1)>=max_fraglen)[0][0]

    idx = (df['last_exon']>=exon_last_use)

    # Fraction reads too far 3' 
    print(f"Filtering {sum(idx)}/{len(idx)} reads with minimum possible length >={max_fraglen}")

    # Reads ending before the A segment will have equal likelihood to all isoforms
    splicing_informative=~idx&(df['last_exon']>=4)
    print(f"{sum(splicing_informative)}/{len(splicing_informative)} reads informative for splicing")

    df = df[~idx]
    df = df.sort_values(['CB','UB']).reset_index(drop=True)

    # Calculate fragment length implied by each isoform for every read
    frag_len = shifts[df['last_exon'].values.astype(int),:] + df['exon_pos'].values.reshape(-1,1)    

    # Plot fragment length distribution if assuming all exons are included
    plt.hist(frag_len[:,-1],50);
    plt.title('Fragment length assuming RABC',fontsize=16)
    plt.ylabel('#Reads',fontsize=16)
    plt.xlabel('Length',fontsize=16)
    plt.savefig('RABC_fraglen_histogram.pdf')

    # Determine compatibility matrix
    compat = np.array([[seg in iso for iso in isoforms] for seg in df['seg']])

    # Seed initial values for mu/sig2 (to be later optimized)
    if (mu is None) and (sig2 is None):
        mu = np.median(frag_len.min(axis=1))
        sig2 = frag_len[frag_len<=600].std()**2
        
        estimate_fraglen_dist = True
    else:
        estimate_fraglen_dist = False

    # Initialize isoform class weights uniformly
    pi = np.ones(len(isoforms))/len(isoforms)

    # Run EM
    mu,sig2,pi,dfup = EM(df,frag_len,compat,mu,sig2,pi,estimate_fraglen_dist=estimate_fraglen_dist)

    print('Estimated parameters:')

    if estimate_fraglen_dist:
        print(f'mu:{mu},var:{sig2}')

    pis = pd.Series(pi,index=isoforms)
    print(pis)

    # Deduplicate using UMIs
    p_iso = dfup[['CB','UB'] + list(isoforms)].drop_duplicates().reset_index(drop=True)

    return(mu,sig2,pis,p_iso)

# Run EM engine
def EM(df,frag_len,compat,mu,sig2,pi,nit=100,estimate_fraglen_dist=True):

    for i in range(0,nit):

        ## E-step

        # Calculate likelihood of each read
        lik = norm.logpdf(frag_len.astype(float),mu,np.sqrt(sig2))
        lik[~compat] = -1*np.Inf

        # Sum across reads to get unnormalized posterior probability for each UMI
        lik = pd.DataFrame(lik,columns=isoforms)
        lik[['CB','UB']]=df[['CB','UB']]
        lik_sum = lik.groupby(['CB','UB']).sum()+np.log(pi)
        lik_sum = lik_sum.subtract(lik_sum.max(axis=1),axis=0)
    
        # Translate into probabilities
        p = np.exp(lik_sum)
        p = p.divide(p.sum(axis=1),axis=0)
    
        # Cases where reads of same UMI are incompatible with all isoforms (rare)
        p[p.isna()]=0

        # Map probabilities back to read level
        dfup = df.join(p,on=['CB','UB'])

        ## M-step

        # Update params
        if estimate_fraglen_dist:
            est_fraglen = (dfup[isoforms]*frag_len).sum(axis=1)
            mu=est_fraglen.mean()
            var=((dfup[isoforms]*(frag_len - mu)**2)).sum().sum()/dfup[isoforms].sum().sum()

        pi = p.sum()/p.sum().sum()

    return(mu,sig2,pi,dfup)

## Takes two inputs:
### (1) A 10x bam with CB/UB tags consistent with cellranger processing
### (2) A gtf-format file with the exon locations of CD45

def main():
    bam_path = sys.argv[1]
    exon_path = sys.argv[2] 
    
    E=pd.read_csv(exon_path,sep='\t',header=None,
              names=['chr','havana','exon','st','en','dot','strand','dot2','info'])
    E = E.sort_values(['st','en']).reset_index(drop=True)
    E['len'] = E['en']-E['st']+1
    
    # Annotate reads with consistency to each isoform position
    df = annotate_reads(bam_path,E)
    df.to_csv('CD45_annotated_reads.csv',index=False)
    
    # Run EM algorithm to estimate isoforms
    mu,sig2,pi,p_iso = isoform_quantification(df,E,max_fraglen=500)
    p_iso.to_csv('isoform_probabilities.csv',index=False)

    # Save output
    with open("parameters.pkl","wb") as f:
        pickle.dump({'mu':mu,'sig2':sig2,'pi':pi,'p_iso':p_iso},f)


if __name__ == "__main__":
    main()


