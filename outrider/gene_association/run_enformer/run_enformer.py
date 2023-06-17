#!/usr/bin/env python
# coding: utf-8

# In[14]:


import pandas as pd
import numpy as np
import pysam
import tensorflow.compat.v2 as tf

import os
import sys
from collections import defaultdict
import pickle


# In[15]:


#import tensorflow_hub as hub
#enformer_model = hub.load("https://tfhub.dev/deepmind/enformer/1").model


# In[16]:

promoter_length=5000
promoter_dir='symm'
max_gnomAD_AF='5e-4'

SEQ_LENGTH = 393216 #Enformer input sequences length
N_bins = 896 #Number of Enformer output bins


# In[36]:


enformer_model_dir = '/s/project/mll/sergey/effect_prediction/tools/enformer/model/'

fasta_fa = '/s/project/mll/sergey/effect_prediction/outrider/aml_association/ass_genes.fa'

tss_tsv = '/s/project/mll/sergey/effect_prediction/outrider/gene_coords/transcripts_canonocal_GRCh37.bed'


# In[18]:


variants_tsv=f"/s/project/mll/sergey/effect_prediction/outrider/aml_association/promoter_variants/max_gnomAD_{max_gnomAD_AF}/promoters_{promoter_length}_{promoter_dir}_short.tsv"

output_dir = f'/s/project/mll/sergey/effect_prediction/outrider/aml_association/enformer/max_gnomAD_{max_gnomAD_AF}/promoters_{promoter_length}_{promoter_dir}/'

# In[19]:


#targets_idx = np.array(np.arange(0,674)) #DNASE


# In[20]:


variants_df = pd.read_csv(variants_tsv, sep='\t', names=['chrom','pos','ref','alt','info_','vcf_name','geneID','rstart','rstop'])
variants_df['pos'] = variants_df['pos']-1 #to 0-based


# In[21]:


tss_df = pd.read_csv(tss_tsv, sep="\t", usecols = [1,3], names=['tss','geneID']) #0-based TSS coordinates for each gene


# In[22]:


regions = variants_df[['geneID','rstart','rstop']].drop_duplicates()
regions = regions.merge(tss_df).set_index('geneID')


# In[23]:


def insert_variant(ref, alt, pos, seq, seq_pos):
    '''
    insert a variant into an existing sequence
    seq - array of 'A', 'T', 'C', 'G' or 'N'
    seq_pos - absolute positions of sequence bp in the genome
    '''
    varpos = seq_pos.index(pos) #index inside the sequence of variant position (relative position)
    if len(alt)==len(ref):
        assert seq[varpos]==ref, 'Wrong reference allele'
        seq[varpos] = alt
    elif len(alt)>len(ref): #insertion
        assert seq[varpos]==ref, 'Wrong reference allele'
        seq = seq[:varpos] + list(alt) + seq[varpos+1:]
        seq_pos = seq_pos[:varpos] + [seq_pos[varpos]]*len(alt) + seq_pos[varpos+1:] #assign all inserted bases the same position
    else: #deletion
        assert seq[varpos:varpos+len(ref)]==list(ref), 'Wrong reference allele'
        seq = seq[:varpos+1] + seq[varpos+len(ref):]
        seq_pos = seq_pos[:varpos+1] + seq_pos[varpos+len(ref):]
    return seq, seq_pos

def center_around_tss(seq, seq_pos, tss_pos):
    '''
    center the sequence around the TSS
    seq - array of 'A', 'T', 'C', 'G' or 'N'
    seq_pos - absolute positions of sequence bp in the genome
    tss_pos - absolute position of TSS in the genome
    '''

    centered_seq = ['N']*SEQ_LENGTH #initialize centered sequence

    tss_idx = seq_pos.index(tss_pos) #TSS index in the input sequence

    left_seq = seq[max(0,tss_idx-SEQ_LENGTH//2):tss_idx] #part of the input sequence to the left of TSS
    right_seq = seq[tss_idx:tss_idx+SEQ_LENGTH//2] #part of the input sequence to the right of TSS

    #insert left and right parts of the input sequence to the centered sequence
    centered_seq[SEQ_LENGTH//2:SEQ_LENGTH//2+len(right_seq)] =  right_seq
    centered_seq[SEQ_LENGTH//2-len(left_seq):SEQ_LENGTH//2] = left_seq

    return centered_seq

def reverse_complement(seq):
    '''
    reverse complement of a given sequence
    '''
    s = list(map(lambda x:{'A':'T','C':'G','T':'A','G':'C'}.get(x,'N'),seq))
    return s[::-1]

def roll_seq(seq, shift):
    '''
    shift a sequence to right (positive shift) or to left (negative shift)
    pad with 'N'
    '''
    if shift>0:
        return ['N']*shift + seq[:-shift]
    else:
        return seq[-shift:] + ['N']*(-shift)

def one_hot(seq):
    '''
    One-hot encoding in order 'ACGT'
    '''
    seq = np.array(seq)
    s = np.vstack((seq=='A',seq=='C',seq=='G',seq=='T')).astype(int).T
    return np.expand_dims(s,0)

def enformer_predict(refseq_c, altseq_c):
    '''
    get enformer predictions for centered reference and alternative sequences
    '''
    all_pred = []

    for seq in refseq_c, reverse_complement(refseq_c), altseq_c, reverse_complement(altseq_c):
        for subseq in one_hot(seq), one_hot(roll_seq(seq,47)),one_hot(roll_seq(seq,-47)),one_hot(roll_seq(seq,25)),one_hot(roll_seq(seq,-25)):
            pred = enformer_model.predict_on_batch(subseq)['human'].numpy()
            all_pred.append(pred[:,N_bins//2,:]) #only the central bin

    all_pred = np.vstack(all_pred)

    ref_pred = all_pred[:10,:]#.mean(axis=0) #average for seq, shifted seq (right), shifted seq (left) and reverse complement
    alt_pred = all_pred[10:,:]#.mean(axis=0)

    #log2fc = np.log2(alt_pred[targets_idx]/ref_pred[targets_idx]).mean()

    return ref_pred, alt_pred


# In[37]:


os.makedirs(output_dir, exist_ok=True)

fasta = pysam.FastaFile(fasta_fa)

enformer_model = tf.keras.models.load_model(enformer_model_dir).model


# In[51]:


N_genes = len(variants_df.geneID.unique())

for gene_idx, geneID in enumerate(variants_df.geneID.unique()):

    print(f'gene {geneID} ({gene_idx+1}/{N_genes})')

    refseq = fasta.fetch(reference = geneID).upper() # reference sequence around TSS

    refseq = list(refseq) # string to list

    region_coords = regions.loc[geneID] # start/end coordinates of reference sequence and TSS position

    refseq_pos = list(range(region_coords.rstart,region_coords.rstop)) # all positions in the reference sequence

    refseq_c = center_around_tss(refseq, refseq_pos, region_coords.tss) # place the reference sequence in the center

    gene_df = variants_df[variants_df.geneID==geneID] # only variants for given gene

    gene_preds = defaultdict(dict) # predictions for given gene, each entry is a vcf name

    N_samples = len(gene_df.vcf_name.unique())

    for sample_idx, sample in enumerate(gene_df.vcf_name.unique()):

        print(f'sample {sample} ({sample_idx+1}/{N_samples})')

        altseq, altseq_pos = list(refseq), list(refseq_pos) # initialize altseq with refseq

        for pos,ref,alt in gene_df.loc[gene_df.vcf_name==sample,['pos','ref','alt']].values:

            print(pos,ref,alt)

            try:
                altseq, altseq_pos = insert_variant(ref,alt,pos,altseq, altseq_pos)
            except Exception as ex:
                print(f'ERROR:{ex}')

        altseq_c = center_around_tss(altseq, altseq_pos, region_coords.tss) # place the alternative sequence in the center

        #mutation_at_tss = (gene_df.pos-region_coords.tss==0).sum()

        #if not mutation_at_tss:
        #    assert refseq_c[SEQ_LENGTH//2]==altseq_c[SEQ_LENGTH//2], 'ref/alt mismatch at TSS'


        ref_pred, alt_pred =  enformer_predict(refseq_c, altseq_c)

        gene_preds[sample]['ref'] = ref_pred
        gene_preds[sample]['alt'] = alt_pred

    with open(output_dir + geneID + '.pickle', 'wb') as f:
        pickle.dump(gene_preds, f)

#enformer_log2fc = get_log2fc(refseq_c, altseq_c)
#predictions.append((geneID, sample, enformer_log2fc))

#pd.DataFrame(predictions, columns=['geneID','vcf_name','enformer_log2fc'])


# In[78]:


#with open(output_dir + 'geneID.pickle', 'rb') as f:
#        a=pickle.load(f)
