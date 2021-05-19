#!/usr/bin/env python

import re, os, sys, pandas, pathlib, CytoSig

from statsmodels.stats.multitest import multipletests
from scipy import stats

base_path = pathlib.Path(__file__).parent.absolute()
base_path = os.path.dirname(base_path)

data_path = os.path.join(base_path, 'data')


def load_gene_map():
    """ load ligand and receptor gene members of each cytokine """
    
    spliter = re.compile('[,+]')
    
    included = set()
    gene_map = []
    
    for title in ['Cytokine', 'Chemokine', 'Growth_Factor', 'Inhibitory']:
        xls = pandas.ExcelFile(os.path.join(data_path, 'annotation', title + '.xlsx'), engine='openpyxl')
        print(title)
        
        for family in xls.sheet_names:
            if family == 'Extra': continue
            
            data = xls.parse(family, index_col=0)[['Gene', 'Receptor']]
            data = data.loc[~data.index.isnull()]
            
            for v in data:
                flag = ~data.loc[:,v].isnull()
                data.loc[flag,v] = data.loc[flag,v].apply(lambda v: set(spliter.split(v.replace(' ','').replace('*',''))) )
            
            data = data.loc[data.index.difference(included)]
            
            gene_map.append(data)
            included.update(data.index)
        xls.close()
    
    gene_map = pandas.concat(gene_map)
    
    return gene_map



def load_run_list():
    """ load file list from cohort folder """
    run_lst = []
    
    cohorts = ['TCGA', 'GTEx']
    
    for cohort in cohorts:
        fpath = os.path.join(data_path, cohort)
        
        files = pandas.read_csv(os.path.join(fpath, 'file.list'), header=None, sep='\t').iloc[:,0]
        for f in files: run_lst.append([cohort, os.path.join(fpath, f)])
    
    return run_lst



def self_correlation_test_all(run_lst, gene_map):
    """ compute the correlation between cytokine ligand or receptor genes and cytokine target activities"""
    # read parallel running indices 
    inx, Nnode = int(sys.argv[1]), int(sys.argv[2])
    
    if len(run_lst) != Nnode:
        sys.stderr.write('Please use %d as N node\n' % len(run_lst))
        sys.exit(1)
    
    cohort, expression = run_lst[inx]
    output = os.path.join(data_path, 'output', cohort + '.' + os.path.basename(expression))
    
    expression = pandas.read_csv(expression + '.gz', sep='\t', index_col=0)
    signature = pandas.read_csv(os.path.join(data_path, 'diff.merge.gz'), sep='\t', index_col=0)
    
    # alpha = 0: simple linear regression, without any ridge penalty
    # nrand = 0: use student t-test, no permutation test
    signature = signature.apply(
        lambda v: CytoSig.ridge_significance_test(
            v.dropna(), expression, alpha=0, alternative="greater", nrand=0, flag_normalize=False, flag_const=True, verbose = False
            )[2].iloc[0] # first variable of zscore matrix
        )
    
    
    # compute correlations
    merge = []
    
    for title in signature:
        gid = title.split('@')[0].split('&')[0]
        
        if gid not in gene_map.index:
            print('Cannot find %s\n' % gid)
            continue
        
        pivots = gene_map.loc[gid].dropna()
        
        if pivots is None:
            print('Cannot find %s\n' % gid)
            continue
        
        rmap = {}
        
        for cat, s in pivots.iteritems():
            if len(expression.index.intersection(s)) == 0:
                print('Complete missing pivots %s\t%s\t%s\n' % (title, cat, ','.join(s)))
                continue
            
            pivot = expression.loc[expression.index.intersection(s)].mean()
            pivot.name = title
            
            if pivot.std() == 0:
                print('jump constant pivot vector %s\n' % title)
                continue
            
            try:
                r_pearson = signature[title].corr(pivot)
                r_spearman = signature[title].corr(pivot, method='spearman')
            
            except:
                print('Correlation failure on %s\n' % title)
                continue
            
            rmap['pearson.' + cat] = r_pearson
            rmap['spearman.' + cat] = r_spearman
        
        merge.append(pandas.Series(rmap, name = title))
    
    data = pandas.concat(merge, axis=1, sort=False).transpose()
    
    data = data.reindex(columns=['pearson.Gene', 'pearson.Receptor', 'spearman.Gene', 'spearman.Receptor'])
    data.to_csv(output, sep='\t', index_label=False)




def merge_files(run_lst):
    """ merge correlation results from individual cancer or tissue types to generate cohort level scores """
    result_path = os.path.join(data_path, 'output')
    
    for pivot in ['Gene', 'Receptor']:
        cohort_map = {}
        
        for cohort, f in run_lst:
            # jump blood cancer, only focus on solid tumors or tissues
            if f.find('ALL') >= 0 or f.find('AML') >= 0: continue
                
            f = os.path.basename(f)
                
            data = pandas.read_csv(os.path.join(result_path, cohort + '.' + f), sep='\t', index_col=0)['pearson.' + pivot].dropna()
            data.name = f
                
            lst = cohort_map.get(cohort)
            if lst is None: lst = cohort_map[cohort] = []
            lst.append(data)    
        
        for cohort, lst in cohort_map.items():
            out = os.path.join(result_path, 'merge.' + pivot + '.' + cohort)
            
            data = pandas.concat(lst, axis=1, join='outer', sort=False)
            data = data.loc[data.isnull().mean(axis=1) < 0.1]
                
            data = data.apply(
                lambda v: pandas.Series([v.dropna().median(), stats.wilcoxon(v.dropna(), alternative="greater")[1]], index=['med', 'p']), axis=1
                )
                            
            data['FDR'] = multipletests(data['p'], method='fdr_bh')[1]
            data.to_csv(out + '.stat', sep='\t', index_label=False)
            
    


def self_correlation_filter(cnt_thres = 5, qthres = 0.05):
    signal_rename_map = {
        'IFNA': 'IFN1',
        'IFNB': 'IFN1',
        
        'IL36A': 'IL36',
        'IL36B': 'IL36',
        'IL36G': 'IL36',
    }

    output = os.path.join(data_path, 'output', 'diff.centroid')

    included = None
    
    for cohort in ['TCGA', 'GTEx']:
        data_Gene = pandas.read_csv(os.path.join(data_path, 'output', 'merge.Gene.' + cohort + '.stat') , sep='\t', index_col=0)
        data_Receptor = pandas.read_csv(os.path.join(data_path, 'Output', 'merge.Receptor.' + cohort + '.stat') , sep='\t', index_col=0)
        
        s = data_Gene.index[data_Gene['FDR'] < qthres].union(data_Receptor.index[data_Receptor['FDR'] < qthres])
        
        if included is None:
            included = s
        else:
            included = included.intersection(s)
    
    data = pandas.read_csv(os.path.join(data_path, 'diff.merge.gz'), sep='\t', index_col=0)
    data = data[included]
    
    # family merge, merge highly similar signals
    flag_family = [v.split('@')[0].split('&')[0] for v in data.columns]
    flag_family = [signal_rename_map[v] if signal_rename_map.get(v) is not None else v for v in flag_family]
    
    data_group = data.groupby(flag_family, axis=1)
    
    #merge = []
    merge_centroid = []
    
    for gid, data in data_group:    
        if data.shape[1] < cnt_thres:
            print('jump', gid, 'by low dataset count', data.shape[1])
            continue
        
        #merge.append(data)
        
        # every element must have effect value
        data = data.loc[(~data.isnull()).sum(axis=1) >= cnt_thres]
        arr = data.median(axis=1)
        arr.name = gid
        
        merge_centroid.append(arr)
    
    #data = pandas.concat(merge, axis=1, join='inner')
    #data.to_csv(output + '.filter.gz', sep='\t', index_label=False, compression='gzip')
    
    data_centroid = pandas.concat(merge_centroid, axis=1, join='outer', sort=False)
    data_centroid_compact = data_centroid.dropna()
    
    print(data_centroid.shape, data_centroid_compact.shape)
    
    # here is the signature matrix in CytoSig prediction model
    data_centroid_compact.to_csv(output, sep='\t', index_label=False)



def main():
    # map between cytokine name to ligand and receptor genes
    gene_map = load_gene_map()
    
    # load file list from TCGA and GTEx cohorts
    run_lst = load_run_list()
    
    if len(sys.argv) == 1:
        # merge result files
        merge_files(run_lst)
        self_correlation_filter()
    else:
        # compute correlations
        self_correlation_test_all(run_lst, gene_map)

if __name__ == '__main__': main()
