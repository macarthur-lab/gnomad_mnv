# -*- coding: utf-8 -*-
__author__ = 'QingboWang'
import numpy as np
import copy as cp
import sys, math, cmath
import time as tm
import pandas as pd
import os
import datetime
import random as rd
from matplotlib import pyplot as plt
from scipy import stats
import re
import mmap
import glob
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import seaborn as sns
from matplotlib.colors import ListedColormap
from collections import OrderedDict
from scipy.stats import *

mut_table = pd.read_csv("/Users/qingbowang/Desktop/MNV/fordist_1KG_mutation_rate_table.txt", sep="\s+")


def draw_heatmap(crstb, title, outdir, num_style="d"):
    mask = crstb.applymap(lambda x: x == 0)
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 10)
    ax.set_aspect('equal')
    ax2 = sns.heatmap(crstb, linewidths=.5, annot=True, fmt=num_style, mask=mask, linecolor="black")
    ax.set_title(title)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    plt.savefig(outdir, dpi=300)

def prob_dNV_null(refspm1b, altspm1b, verbose=False):#probability under the null (independent) model
    #takes ref2base +- 1base and that of alt
    if ((refspm1b[0] != altspm1b[0]) | (refspm1b[3] != altspm1b[3]) | (refspm1b[1] == altspm1b[1]) | (refspm1b[2] == altspm1b[2])):
        return 0 #if not MNV, return 0 (might be better to return NA? but will go on. can mask later.)
    intm1_4b = refspm1b[0] + altspm1b[1] + refspm1b[2:] #intermediate 1, where ref1-> alt1 already happened
    intm2_4b = refspm1b[:2] + altspm1b[2] + refspm1b[3] #intermediate 2, where ref2-> alt2 already happened
    path1 = mut_table.loc[(mut_table["from"]==refspm1b[:3]) & (mut_table["to"]==intm1_4b[:3]), "mu_snp"].values[0] * \
            mut_table.loc[(mut_table["from"] == intm1_4b[1:]) & (mut_table["to"] == altspm1b[1:]), "mu_snp"].values[0]
    path2 = mut_table.loc[(mut_table["from"]==refspm1b[1:]) & (mut_table["to"]==intm2_4b[1:]), "mu_snp"].values[0] * \
            mut_table.loc[(mut_table["from"] == intm2_4b[:3]) & (mut_table["to"] == altspm1b[:3]), "mu_snp"].values[0]
    if verbose:
        print("prob. of {0} to {1}\n".format(refspm1b, altspm1b))
        print("1st path: {0} -> {1} -> {2}\n".format(refspm1b, intm1_4b, altspm1b))
        print("p1 = p({0}->{1}) * p({2}->{3})\n".format(refspm1b[:3], intm1_4b[:3], intm1_4b[1:], altspm1b[1:]))
        print("   ={0} * {1} = {2}\n".format((mut_table.loc[(mut_table["from"]==refspm1b[:3]) & (mut_table["to"]==intm1_4b[:3]), "mu_snp"].values[0]), \
                                             (mut_table.loc[(mut_table["from"] == intm1_4b[1:]) & (mut_table["to"] == altspm1b[1:]), "mu_snp"].values[0]), (path1)))
        print("2nd path: {0} -> {1} -> {2}\n".format(refspm1b, intm2_4b, altspm1b))
        print("p2 = p({0}->{1}) * p({2}->{3})\n".format(refspm1b[1:], intm2_4b[1:], intm2_4b[:3], altspm1b[:3]))
        print("   ={0} * {1} = {2}\n".format((mut_table.loc[(mut_table["from"]==refspm1b[1:]) & (mut_table["to"]==intm2_4b[1:]), "mu_snp"].values[0]), \
                                             (mut_table.loc[(mut_table["from"] == intm2_4b[:3]) & (mut_table["to"] == altspm1b[:3]), "mu_snp"].values[0]), (path2)))
        print ("sum of probability = p1 + p2 = {0}".format((path1+path2)))
    return (path1 + path2)


def calc_ratio(afs):
    if afs[0]=="." or afs[1]==".": return (0)
    elif float(afs[0])*float(afs[1])==0: return (0)
    else: return (min(float(afs[1])/float(afs[0]), float(afs[0])/float(afs[1])))
def log2_adjusted(x):
    if x==0: return (0) #non significant
    else: return (np.log2(x))
def fisher_OR_and_pval(x1, x2, y1, y2): #1 case, 1 alt, 2 case, 2 alt
    oddsratio, pvalue = stats.fisher_exact([[x1, x2], [y1, y2]])
    return (oddsratio, pvalue)
def log2OR_adjusted(OR, P):
    if P>0.05 / (8**2): return (0) #non significant
    else: return (log2_adjusted(OR))
def plot_heatmap_fisher(table1, table2, title, dir, only_signif=True): #16x16 tables.
    #odds ratio for each entry of the table
    #if mask=True: mask all the non significant ones
    out = table1.copy() #for output
    pval_table = table1.copy()  # for output
    X2 = sum(table1.sum(axis=0)) #=x2
    Y2 = sum(table2.sum(axis=0)) #=y2
    for i in range(table1.shape[0]):
        for j in range(table1.shape[1]):
            (OR, P) = fisher_OR_and_pval(table1.iloc[i,j], X2, table2.iloc[i,j], Y2)
            log2OR = log2_adjusted(OR)
            out.iloc[i,j] = log2OR
            pval_table.iloc[i, j] = P
    if (only_signif): mask = pval_table.applymap(lambda x: x > 0.05/(16*9))
    else: mask = out.applymap(lambda x: x == 0)
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 10)
    ax.set_aspect('equal')
    ax2 = sns.heatmap(out, linewidths=.5, cbar_kws={"label": "log2(OR)"}, linecolor="black", fmt='.2f', annot=True, mask=mask)
    ax.set_title("{0} \n n= {1}, {2}".format(title, X2, Y2))
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    plt.savefig(dir, dpi=300)

def calc_ratio_zeroadjusted(v1,v2):
    if v2==0: return 0
    else: return (float(v1)/v2)
    
def plot_heatmap_ratio(table1, table2, title, dir): #16x16 tables.
    #the ratio of table2 compare to table1, for each cell entry
    out = table1.copy() #for output
    X2 = sum(table1.sum(axis=0)) #=x2
    Y2 = sum(table2.sum(axis=0)) #=y2
    for i in range(table1.shape[0]):
        for j in range(table1.shape[1]):
            out.iloc[i,j] = calc_ratio_zeroadjusted(table1.iloc[i,j],table2.iloc[i,j])
    mask = out.applymap(lambda x: x == 0)
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 10)
    ax.set_aspect('equal')
    ax2 = sns.heatmap(out, linewidths=.5, cbar_kws={"label": "fraction"}, linecolor="black", fmt='.2f', annot=True, mask=mask)
    ax.set_title("{0} \n n= {1}, {2}".format(title, X2, Y2))
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    plt.savefig(dir, dpi=300)


def draw_null_matrix_dnv(obs_refs, cols, cov): 
    #returns a matrix of probability of each entry, given the number of reference 2bp as obs_refs
    #and the column names as cols, 
    #and the coverage as cov
    null_table = pd.DataFrame(np.zeros([16, 16]))  # just to put the input
    null_table.columns = cols
    null_table.index = cols
    for fourbs in obs_refs.index:
        tb = pd.DataFrame(np.zeros([16, 16]))  # just to put the input
        tb.columns = null_table.columns
        tb.index = null_table.index
        refs = fourbs[1] + "," + fourbs[2]  # refs, as M,N
        for alts in tb.columns:
            fourbs_alt = fourbs[0] + alts[0] + alts[2] + fourbs[3]
            tb.loc[refs, alts] = prob_dNV_null(fourbs, fourbs_alt)
        tb = tb * obs_refs[fourbs] * cov.loc[fourbs, "FofC"] #might not need to be squared
        null_table = null_table + tb
    return (null_table)

def max_repeat(context, mer):
    #mer needs to be smaller than 4
    r = ["A","T","G","C"]
    if mer==2:
        r2 = []
        for i in r:
            for j in r:
                r2.append(i+j)
        r = r2
    if mer==3:
        r3 = []
        for i in r:
            for j in r:
                for k in r:
                    r3.append(i+j+k)
        r = r3
    cnt_max = 0
    for unit in r:
        cnt = 0
        unit_now = unit
        while unit_now in context:
            cnt = cnt + 1
            unit_now = unit_now + unit #add a repeat unit count if it is continuing
        if cnt_max<cnt: cnt_max=cnt
    return (cnt_max)

def revcomp(seq):
    comp = {}
    comp["A"] = "T"
    comp["T"] = "A"
    comp["G"] = "C"
    comp["C"] = "G"
    comp["N"] = "N"
    out = ""
    for i in seq[::-1]:
        out = out + comp[i]
    return (out)

def collapse_crstb_to_revcomp(crstb):
    # collaspse to a table, instead of another matrix
    flt = crstb.stack().reset_index()
    flt.columns = ['refs', 'alts', 'cnt']
    for i in range(flt.shape[0]):
        if i in flt.index:  # if it hasn't been deleted yet
            refs_revcomp = revcomp(flt.refs[i])
            alts_revcomp = revcomp(flt.alts[i])
            ix_revcomp = flt[(flt.refs == refs_revcomp) & (flt.alts == alts_revcomp)].index[0]
            if i != ix_revcomp:  # if revcomp is not yourself
                flt.loc[i, "cnt"] = flt.loc[i, "cnt"] + flt.loc[ix_revcomp, "cnt"]
                flt.drop(ix_revcomp, inplace=True)
    flt = flt[flt.cnt > 0]  # deleting the no dNVs
    flt.reset_index(inplace=True)
    del flt["index"]
    return (flt)


def calc_symmetry(crosstab, out_dir):
    #calculate how symmetric they are, return as the ratio(that are closer to one)
    sym = cp.copy(crosstab)
    for i in range(crosstab.shape[0]):
        for j in range(crosstab.shape[1]):
            refs = crosstab.index[i]
            alts = crosstab.index[j]
            if refs==0:
                sym.iloc[i, j] = 0
            else:
                revcomp_refs = str(Seq(refs).reverse_complement())
                revcomp_alts = str(Seq(alts).reverse_complement())
                if ((refs==revcomp_refs) & (alts==revcomp_alts)):
                    sym.iloc[i, j] = 0 #if revcomp is yourself, return 0
                else:
                    sym.iloc[i, j] = crosstab.loc[refs, alts] / crosstab.loc[revcomp_refs, revcomp_alts]
    # plot
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 10)
    ax.set_aspect('equal')
    ax2 = sns.heatmap(sym, linewidths=.5, annot=True, mask = sym.applymap(lambda x: x == 0))
    ax.set_title("number of MNV divided by \n that of its complementary pattern")
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    plt.savefig(out_dir, dpi=300)
    #and finally return the matrix itself
    return (sym)



