import numpy as np
import pandas as pd
import scanpy as sc
sc.settings.verbosity = 3
import matplotlib.pyplot as plt
import seaborn as sns

import subprocess
import re
import os 
import sys
import gc
from multiprocessing import Pool, cpu_count

import sklearn.metrics, sklearn.cluster
from sklearn.decomposition import PCA
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import kneighbors_graph, radius_neighbors_graph

from time import ctime
import scipy.sparse as sps
import igraph as ig
import leidenalg as la
import sklearn

import sklearn.metrics, sklearn.cluster
from sklearn.decomposition import PCA
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import kneighbors_graph, radius_neighbors_graph
import sklearn.neighbors
from tqdm import tqdm
# from scanpy import _utils as utils



def rank_weight(csr_mtx, rank_cutoff, transpose=False, ret='value', include_self=True, alpha=1, decay='distance2'):
    from scipy.stats import rankdata
    if transpose:
        tmp = csr_mtx.tocoo().T.tocsr()
    else:
        tmp = csr_mtx.tocoo().tocsr()

    data = tmp.data[tmp.data != 0]
    ptr = tmp.indptr
    row, col = np.nonzero(tmp)
    
    value = np.array([])
    data2 = np.array([])
    row2 = np.array([])
    col2 = np.array([])
    for i in range(tmp.shape[0]):
        row_data = data[ptr[i]:ptr[i+1]].copy()
        row_ind = row[ptr[i]:ptr[i+1]]
        col_ind = col[ptr[i]:ptr[i+1]]
        
        if not include_self:
            indices = row_ind != col_ind
            row_data = row_data[indices]
            row_ind = row_ind[indices]
            col_ind = col_ind[indices]

        rank_data = rankdata(- row_data, method='min')

        indices = rank_data <= rank_cutoff

        
        value = np.append(value, row_data[indices])
        data2 = np.append(data2, rank_data[indices])
        col2 = np.append(col2, col_ind[indices])
        row2 = np.append(row2, row_ind[indices])
        
    if ret == 'rank':
        return sps.csr_matrix((data2+1, (row2, col2)),shape=tmp.shape)
    if ret == 'value':
        return sps.csr_matrix((value, (row2, col2)),shape=tmp.shape)
    if ret == 'weight':
        w = sps.csr_matrix((value, (row2, col2)),shape=tmp.shape)
        w = dist_decay(w, kernel=decay, alpha=alpha)
        # w = weight_norm(w, alpha=alpha)
        return w

def weight_norm(csr_mtx, alpha): ### normalize the total weights per spot to 1 
    con = csr_mtx.copy()
    con.data = con.data ** alpha
    scalers = np.array(con.sum(axis=1))[:, 0]
    scalers = 1 / scalers
    sklearn.utils.sparsefuncs.inplace_row_scale(con, scalers)
    return con
        
def dist_decay(con, kernel='Gaussion', alpha=2):
    
    r = con.data.mean()/alpha
    
    if kernel == 'Gaussian':
        con.data = np.exp(-con.data**2/(2*r**2))
    
    if kernel == 'Exponential':
        con.data = np.exp(-con.data/(2*r**2))
    
    if kernel == 'Laplacian':
        con.data = np.exp(-con.data/(2*r))
    
    if kernel == 'Cauchy':
        con.data = np.exp(1/(con.data**2/r + 1))
        
    if kernel == 'distance':
        con.data = 1/(1+con.data)
    
    if kernel == 'distance2':
        con.data = 1/(1+con.data**2)
    return con

def build_connect(adata, radius, include_self=False,scaling=True, norm=True, decay='Gaussian', alpha=2, rank_cutoff=100, rank2=True, n_neighbors=15, method='radius', use_raw=False):
    
    if method == 'radius':
        con  = radius_neighbors_graph(adata.obsm['spatial'], radius=radius, include_self=include_self, mode='distance')
               
    
    if method == 'knn':
        con = kneighbors_graph(adata.obsm['spatial'], n_neighbors=n_neighbors, include_self=include_self, mode='distance')
    
    con = dist_decay(con, kernel=decay, alpha=alpha)
    
    if norm:   ### normalize the total weights per spot to 1 
        con = weight_norm(con, alpha=1)           
    print('Step1: cell-cell connectivity calculation finished,', ctime()) 
    
    if use_raw:
        adata = adata.raw.to_adata()

    reads = adata.X.T.copy()
    
    v = mul_rank_weight(reads, rank_cutoff=rank_cutoff, ret='rank', alpha=1, transpose=False, nthreads=10)
   
    v = reads.multiply(v > 0)
    v_min = mul_rank_weight(-v, rank_cutoff=1, ret='value', alpha=1, transpose=False, nthreads=10)
    v_min = (-v_min).max(axis=1).toarray()[:, 0]
    sklearn.utils.sparsefuncs.inplace_row_scale(v, 1/v_min)

    print('Step2: ranking SNV expression values finished,', ctime())

    s_con = mul_sparse_dot(con, v.T, nthreads=20)
    
    counts = (con != 0).astype('float64')
    value = (v.T != 0).astype('float64')

    mode = np.array(counts.sum(axis=1))[:, 0]
    p = np.array(value.sum(axis=0))[0, :]/value.shape[0]

    scaler_mode = 1/mode
    scaler_p = 1/p
    
    counts = counts.dot(value)
    # return counts, value, scaler_mode, scaler_p
    sklearn.utils.sparsefuncs.inplace_row_scale(counts.T, scaler_p)
    sklearn.utils.sparsefuncs.inplace_row_scale(counts, scaler_mode)
    counts.data = np.log10(1 + counts.data)
    
    s_con = s_con.multiply(counts)

    print('Step3: cell-SNV connectivity calculation finished,', ctime())

    if rank2:
        s_con = mul_rank_weight(s_con, rank_cutoff=100, ret='value', nthreads=20)
        print('Step4: ranking cell-SNV connectivity finished,', ctime())

    ss_con = mul_sparse_dot(v, s_con, nthreads=10)
    print('Step5: SNV-SNV connectivity calculation finished,', ctime())
    
    return ss_con





def create_graph(con, directed=True):
    import igraph as ig
    sources, targets = con.nonzero()
    weights = con[sources, targets]
    if isinstance(weights, np.matrix):
        weights = weights.A1
        
    g = ig.Graph(directed=directed)
    g.add_vertices(con.shape[0]) 
    g.add_edges(list(zip(sources, targets)))
    g.es['weight'] = weights
    
    return g


def partition_leiden(time, g, ptype, max_comm_size,snv_var):
    part = la.find_partition(g, partition_type=ptype, seed=123, resolution_parameter=time, max_comm_size=max_comm_size, weights=np.array(g.es['weight']).astype(np.float64))
    y_predict = np.array(part.membership)
    peak_labels_df = pd.DataFrame(y_predict, index = snv_var, columns=['group']) 
    accum = np.sum(peak_labels_df.value_counts() > 1) / len(peak_labels_df.value_counts())
    return accum



def build_sg(snv, con, save=None,resolution=None, rank_cutoff=50, ptype=la.RBConfigurationVertexPartition, max_comm_size=0, syn=False):
    con = mul_rank_weight(con, rank_cutoff=rank_cutoff, ret='value',transpose=False, include_self=False, nthreads=20)
    if syn:
        con = (con + con.T) * 0.5
        con = weight_norm(con, alpha=1)

        g = create_graph(con, directed=False)
        
    else:
        con = weight_norm(con, alpha=1)
        g = create_graph(con, directed=True)
        
    if resolution == None:
        snv_var = snv.var_names.copy()
        args = [(time, g, ptype, max_comm_size,snv_var) for time in range(1,10)]

        
        with Pool(processes = 2) as pool:
            accum_list = list(tqdm(pool.starmap(partition_leiden, args), total=len(range(1,10))))
    
        plt.figure(figsize=(5, 5))
        plt.plot(accum_list, marker='o', linestyle='--', color='black', markerfacecolor='black', markeredgecolor='black')
        plt.title("Proportion of Communities with Multiple Members using Leiden Algorithm")
        plt.xlabel("Resolution Parameter (Time)")
        plt.ylabel("Proportion of Multi-member Communities")
        plt.grid(False)
        plt.show()
    else:
        part = la.find_partition(g, partition_type=ptype,  seed=123, resolution_parameter=resolution, max_comm_size=max_comm_size, weights=np.array(g.es['weight']).astype(np.float64))
        y_predict = np.array(part.membership)
        
        peak_labels_df = pd.DataFrame(y_predict, index = snv.var_names, columns = ['group'])
        
        snv.var['snv_group'] = y_predict
        
        counts = snv.var['snv_group'].value_counts()
        
        groups = list(set(y_predict))
        sg_matrix = np.array([snv.X.T[np.where(y_predict==x)[0],:].sum(axis=0) for x in groups])
        sg_matrix = sg_matrix[:,0,:].T
        
        sg_df = pd.DataFrame(sg_matrix, index = snv.obs_names, columns = groups)
        snv.obsm['snv_group'] = sg_df.values
        peak_values = peak_labels_df.value_counts().reset_index()
        peak_values['group'] = 'sg' + peak_values['group'].map(str)
        peak_values = peak_values.set_index('group')
        peak_values.columns = ['snv_count']
        ac = sc.AnnData(snv.obsm['snv_group'].copy(), obs=snv.obs)
        ac.var_names = ['sg'+str(x) for x in ac.var_names]
        ac.var['snv_count'] = peak_values['snv_count']
        ac.obsm['spatial'] = snv.obsm['spatial']
        return snv,ac


def snp_norm(snp, adata, method='umi'):
    if method == 'umi':
        scalers = 10000/adata.obs['total_counts'].values
        sklearn.utils.sparsefuncs.inplace_row_scale(snp.X, scalers)
    if method == 'gene':
        scalers = 10000/adata.obs['n_genes_by_counts'].values
        sklearn.utils.sparsefuncs.inplace_row_scale(snp.X, scalers)
    if method == 'total':
        sc.pp.calculate_qc_metrics(snp, percent_top=None, log1p=False, inplace=True)
        scalers = np.array(10000/snp.X.sum(axis=0))[0]
        sklearn.utils.sparsefuncs.inplace_column_scale(snp.X, scalers)
    snp.X.data = np.log2(snp.X.data+1)
    

import threading

class sparse_dot(threading.Thread):
    def __init__(self, csr_a, b, count, ret_dict):
        threading.Thread.__init__(self)
        self.csr_a = csr_a
        self.b = b
        self.count = count
        self.ret_dict = ret_dict
    
    def run(self):
        ret = self.csr_a.dot(self.b)
        self.ret_dict[self.count] = ret
    
    # def __del__(self):
    #     print('thread %d is ending' % self.count)

class sparse_rank(threading.Thread):
    def __init__(self, csr_mtx, rank_cutoff, ret, include_self, count, ret_dict, alpha):
        threading.Thread.__init__(self)
        self.csr_mtx = csr_mtx
        self.rank_cutoff = rank_cutoff
        self.ret = ret
        self.include_self = include_self
        self.count = count
        self.ret_dict = ret_dict
        self.alpha = alpha
    
    def run(self):
        tmp = rank_weight(self.csr_mtx, 
                          rank_cutoff=self.rank_cutoff, 
                          transpose=False, 
                          ret=self.ret, 
                          include_self=self.include_self,
                          alpha=self.alpha)
        self.ret_dict[self.count] = tmp


def mul_rank_weight(csr_mtx, rank_cutoff, nthreads, transpose=False, ret='value', include_self=True, alpha=1):

    if transpose:
        csr_mtx = csr_mtx.tocoo().T.tocsr()
        
    nrow = csr_mtx.shape[0]

    step = nrow // nthreads
    if (step * nthreads) != nrow:
        step = nrow // nthreads + 1
    slices = list(range(0, nrow, step))
    slices.append(nrow)
    
    ret_dict = dict()
    threads = list()
    for i in range(nthreads):
        a = csr_mtx[slices[i]:slices[i+1], :]
        count = i
        thread = sparse_rank(a, 
                             rank_cutoff=rank_cutoff, 
                             ret=ret, 
                             include_self=include_self, 
                             count=count, 
                             ret_dict=ret_dict,
                             alpha=alpha)
        threads.append(thread)
        thread.start()
    
    for thread in threads:
        thread.join()
    
    return sps.vstack([ret_dict[i] for i in range(nthreads)])
        

def mul_sparse_dot(csr_a, b, nthreads):
    
    if not sps.isspmatrix_csr(csr_a):
        csr_a = csr_a.tocsr()
        
    nrow = csr_a.shape[0]

    step = nrow // nthreads
    if (step * nthreads) != nrow:
        step = nrow // nthreads + 1
    slices = list(range(0, nrow, step))
    slices.append(nrow)
    
    ret_dict = dict()
    threads = list()
    for i in range(nthreads):
        a = csr_a[slices[i]:slices[i+1], :]
        count = i
        thread = sparse_dot(a, b, count, ret_dict)
        threads.append(thread)
        thread.start()
    
    for thread in threads:
        thread.join()
    
    return sps.vstack([ret_dict[i] for i in range(nthreads)])
        

def sg_cluster(adata, res=0.5, npcs=15, neighbors=15, norm=False, subset=None):
    sc.settings.verbosity = 1
    
    ac = sc.AnnData(adata.obsm['snv_group'].copy(), obs=adata.obs)
    if subset != None:
        t = adata.var['snv_group'].value_counts()
        t = t[t>subset].index.tolist()
        t = list(map(str,t))
        ac = ac[:,ac.var_names.isin(t)]
    ac.var_names = ['sg'+str(x) for x in ac.var_names]
    ac.layers['raw_counts'] = ac.X.copy()
    if norm:
        sc.pp.normalize_total(ac, target_sum=1e4)
    sc.pp.log1p(ac)
    sc.pp.scale(ac, max_value=10)
    sc.pp.pca(ac, svd_solver='full')
    sc.pl.pca_variance_ratio(ac, n_pcs=50)
    sc.pp.neighbors(ac, n_neighbors=neighbors, n_pcs=npcs)
    sc.tl.umap(ac)
    sc.tl.leiden(ac, resolution=res)
    # sc.pl.umap(ac, color=['leiden'])
    adata.obs['leiden_sg'] = ac.obs['leiden']
    adata.obsm['X_umap'] = ac.obsm['X_umap']
    return ac

# process

def snv2avinput(
    snv,
    sample_name,
    annovar,
    annovar_ref,
    gtf,
    spe,
    ref_name,
    outdir,
    overwrite = False
):

    if overwrite:
        avi = pd.DataFrame(snv.var_names)
        avi.columns=['snv']
        avi['chrom'] = avi['snv'].apply(lambda x : x.split('_')[0])
        avi['site1'] = avi['snv'].apply(lambda x : x.split('_')[1].split(':')[0])
        avi['site2'] = avi['site1']
        avi['ref_base'] = avi['snv'].apply(lambda x : x.split('_')[1].split(':')[1].split('>')[0])
        avi['alt_base'] = avi['snv'].apply(lambda x : x.split('_')[1].split(':')[1].split('>')[1])
        avi = avi[['chrom','site1','site2','ref_base','alt_base']]
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        avi.to_csv(f'{outdir}/{sample_name}.avinput',index=None,sep='\t',header=None)
        with open(f'{outdir}/{sample_name}.anno.sh','w') as f:
            f.write(f'perl {annovar} {outdir}/{sample_name}.avinput {annovar_ref} --out {outdir}/{sample_name} --buildver {spe} --remove --protocol {ref_name} --operation g --nastring .\n')
        subprocess.call(f"bash {outdir}/{sample_name}.anno.sh",shell=True)
    
    annotation = pd.read_csv(f'{outdir}/{sample_name}.{spe}_multianno.txt',sep='\t')
    annotation['snp'] = annotation['Chr']+'_'+annotation['Start'].map(str) + ':' + annotation['Ref'] + '>' +annotation['Alt']
    annotation = annotation.set_index('snp')
    gene_info = pd.read_csv(gtf, sep='\t', comment='#', usecols=[0,2,3,4, 8], names=['chr','type','start','end', 'info'])
    gene_info = gene_info[gene_info['type'] == 'gene']
    gene_info['gene_id'] = gene_info['info'].str.extract(r'gene_id "([^"]+)";')
    gene_info = gene_info.drop_duplicates(subset='gene_id').reset_index(drop=True)
    gene_info['gene_name'] = gene_info['info'].str.extract(r'gene_name "([^"]+)";')
    gene_info['gene_lenth'] = gene_info['end'] - gene_info['start']
    gene_annotation_dict = dict(zip(gene_info['gene_id'],gene_info['gene_name']))
    gene_len_dict = dict(zip(gene_info['gene_id'],gene_info['gene_lenth']))
    def map_gene(genecode):
        return gene_annotation_dict.get(genecode,'None')
    def map_len(genecode):
        return str(gene_len_dict.get(genecode,'None'))
    annotation['gene_name'] = annotation[f'Gene.{ref_name}'].apply(lambda x : ';'.join(list(map(map_gene,x.split(';')))))
    annotation['gene_lenth'] = annotation[f'Gene.{ref_name}'].apply(lambda x : ';'.join(list(map(map_len,x.split(';')))))
    annotation = annotation[[f'Func.{ref_name}',f'ExonicFunc.{ref_name}','gene_name','gene_lenth']]
    snv.var[f'Func'] = annotation[f'Func.{ref_name}']
    snv.var[f'ExonicFunc'] = annotation[f'ExonicFunc.{ref_name}']
    snv.var['gene_name'] = annotation['gene_name']
    snv.var['gene_lenth'] = annotation['gene_lenth']
    return snv

def summary(region):
    if re.match(r'^(ncRNA_)?exonic(;splicing)?$', region):
        part = "Exonic"
    elif re.match(r'^(ncRNA_)?intronic$', region) or re.match(r'^(ncRNA_)?splicing$', region):
        part = "Intronic"
    elif re.match(r'^UTR', region):
        part = "UTR"
    else:
        part = "Intergenic"
    return part

def processsnv(
    sample_name,
    snv,
    snv_depth,
    gtf,
    annovar_ref,
    annovar_spe,
    annovar_ref_name,
    annovar,
    outdir,
    thrshold = 20,
    min_cells = 5,
):
    snv.var['SNVDepth'] = np.sum(snv_depth.X,axis = 0).reshape(-1,1)
    snv.obs['TotalDepth'] = np.sum(snv_depth.X,axis = 1)
    snv.var['SNVCount'] = np.sum(snv.X,axis = 0).reshape(-1,1)
    sc.pp.filter_genes(snv,min_cells = 1)
    sc.pp.filter_genes(snv_depth,min_cells = 1)
    thrshold = 20
    sub_snv = snv[:,snv.var['SNVDepth'] >= thrshold].copy()
    sub_snv_depth = snv_depth[:,snv.var['SNVDepth'] >= thrshold].copy()
    sc.pp.filter_genes(sub_snv,min_cells = min_cells)
    sub_snv = snv2avinput(
        sub_snv,
        sample_name,
        annovar = annovar,
        annovar_ref = annovar_ref,
        gtf = gtf,
        spe = annovar_spe,
        ref_name = annovar_ref_name,
        outdir = outdir,
        overwrite=True
    )
    sub_snv.var['Func_L0'] = sub_snv.var['Func'].apply(lambda x : summary(x))
    return sub_snv

def _normalize_data(X, counts, after= None):
    from sklearn.utils.sparsefuncs import inplace_row_scale
    X = X.copy()
    after = np.median(counts_greater_than_zero, axis=0) if after is None else after
    counts += counts == 0
    counts = counts / after
    inplace_row_scale(X, 1 / counts)
    return X

def normalize_with_rna(snv,rna,target_sum = 1e4):
    snv.obs['UMI_counts'] = rna.obs['total_counts']
    umi_per_cell = np.ravel(snv.obs['UMI_counts'])
    snv.layers['norm'] = np.log1p(_normalize_data(snv.X,umi_per_cell,after = target_sum))
    snv.obs['snvperumi'] = np.sum(snv.layers['norm'],axis = 1)
    matrix = snv.layers['norm'].todense()
    matrix = np.where(matrix>0,1,0)
    snv.obs['SNVtypes'] = np.sum(matrix,axis = 1 ).reshape(-1,1)

from anndata import AnnData
from scipy.sparse import csc_matrix
def bulid_windows(snv,window_size = 100000,basis = 'basis'):
    snv_df = snv.to_df().T
    snv_df['chrom'] = snv_df.index.map(lambda x: x.split('_')[0])
    snv_df['num'] = snv_df.index.map(lambda x: str(int(x.split('_')[1].split(':')[0])//window_size))
    snv_df['windows'] = snv_df['chrom'] +'@'+ snv_df['num']
    del snv_df['chrom']
    del snv_df['num']
    snv_df = snv_df.groupby('windows').sum().T
    snv_gene_adata = AnnData(snv_df)
    snv_gene_adata.X = csc_matrix(snv_gene_adata.X)
    snv_gene_adata.obsm[basis] = snv.obsm[basis]
    return snv_gene_adata


def get_min_distance(adata,basis = 'spatial'):
    locations = adata.obsm[basis].copy()
    distances = np.linalg.norm(locations[:, np.newaxis] - locations, axis=2)
    np.fill_diagonal(distances, np.inf)
    min_distances = np.min(distances, axis=1)
    average_min_distance = np.mean(min_distances)
    return average_min_distance

import igraph as ig
def netplot(snv_gene_adata,snv_group,ss_con,topn = 5,rank_cutoff=50,ret='value',transpose=False, include_self=False, nthreads=20,save = None):
    procon = mul_rank_weight(ss_con, rank_cutoff = rank_cutoff, ret=ret,transpose=transpose, include_self=include_self, nthreads=nthreads)
    procon = weight_norm(procon, alpha=1)
    inlist = snv_gene_adata[:,snv_gene_adata.var['snv_group'] == snv_group].var_names

    indices = [np.where(snv_gene_adata.var_names == gene)[0][0] for gene in inlist if gene in snv_gene_adata.var_names]
    sub_con = procon[indices, :][:, indices]
    g = create_graph(sub_con, directed=True)
    g.vs['name'] = inlist
    visual_style = {}
    visual_style["bbox"] = (300, 300) 
    visual_style["margin"] = 20
    visual_style["edge_color"] = "gray"
    degrees = g.degree()
    n = len(g.vs)
    top_indices = sorted(range(n), key=lambda i: degrees[i], reverse=True)[:topn]
    top_degrees = [degrees[i] for i in top_indices]
    max_degree_indices = sorted(range(len(degrees)), key=lambda i: degrees[i], reverse=True)[:topn]
    labels = ["" for _ in range(len(g.vs))]
    for idx in max_degree_indices:
        labels[idx] = g.vs[idx]["name"]
        
    min_top_degree, max_top_degree = min(top_degrees), max(top_degrees)
    node_sizes = [topn] * n 
    for idx in top_indices:
        if max_top_degree > min_top_degree: 
            size = ((degrees[idx] - min_top_degree) / (max_top_degree - min_top_degree) * 15) + 10  
        else:
            size = 25
        node_sizes[idx] = max(size, 5)
        
    layout = g.layout("mds")   #"grid"  "graphopt" "davidson_harel"  "drl" "mds" "kk"
    visual_style = {}
    visual_style["bbox"] = (350, 350)
    visual_style["edge_color"] = "gray"
    visual_style['edge_arrow_size'] = 0.1  
    if save != None:
        target = save
    
    ig.plot(g, layout=layout,
            vertex_size=node_sizes,edge_width = 0.2,
            vertex_color=['#E84B50' if idx in top_indices else '#00A1DF' for idx in range(n)],
            vertex_label=labels,
            **visual_style,target=save)
