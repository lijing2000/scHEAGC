import torch
import random
import numpy as np
import scanpy as sc
from sklearn import metrics
from munkres import Munkres
import torch.nn.functional as F
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score as ari_score
from sklearn.metrics.cluster import normalized_mutual_info_score as nmi_score
from sklearn.metrics import adjusted_mutual_info_score
import opt


def setup_seed(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.manual_seed(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True


def numpy_to_torch(a, sparse=False):
    if sparse:
        a = torch.sparse.Tensor(a)
        a = a.to_sparse()
    else:
        a = torch.FloatTensor(a)
    return a


# the reconstruction function
def reconstruction_loss(X, A_norm, X_hat, Z_hat, A_hat):
    loss_ae = F.mse_loss(X_hat, X)
    loss_w = F.mse_loss(Z_hat, torch.spmm(A_norm, X))
    loss_a = F.mse_loss(A_hat, A_norm.to_dense())
    loss_igae = loss_w + opt.args.alpha_value * loss_a
    loss_rec = loss_ae + loss_igae
    return loss_rec


def target_distribution(Q):
    weight = Q ** 2 / Q.sum(0)
    P = (weight.t() / weight.sum(1)).t()
    return P


# clustering guidance
def distribution_loss(Q, P):
    loss = F.kl_div((Q[0].log() + Q[1].log() + Q[2].log()) / 3, P, reduction='batchmean')
    # loss = F.kl_div(Q[0].log(), P, reduction='batchmean')
    return loss


def off_diagonal(x):
    n, m = x.shape
    assert n == m

    return x.flatten()[:-1].view(n - 1, n + 1)[:, 1:].flatten()


def cross_correlation(Z_v1, Z_v2):

    return torch.mm(F.normalize(Z_v1, dim=1), F.normalize(Z_v2, dim=1).t())


def correlation_reduction_loss(S):

    return torch.diagonal(S).add(-1).pow(2).mean() + off_diagonal(S).pow(2).mean()


def drr_loss(cons):

    S_N = cross_correlation(cons[0], cons[1])
    L_N = correlation_reduction_loss(S_N)

    S_F = cross_correlation(cons[2], cons[3])
    L_F = correlation_reduction_loss(S_F)

    loss_drr = opt.args.lambda1 * L_N + opt.args.lambda2 * L_F

    return loss_drr


def clustering(Z, y):
    model = KMeans(n_clusters=opt.args.n_clusters, n_init=10)
    cluster_id = model.fit_predict(Z.data.cpu().numpy())

    ari, nmi, ami, acc = eva(y, cluster_id, show_details=True)

    return ari, nmi, ami, acc, model.cluster_centers_


def assignment(Q, y):
    y_pred = torch.argmax(Q, dim=1).data.cpu().numpy()
    ari, nmi, ami, acc = eva(y, y_pred, show_details=False)
    return ari, nmi, ami, acc, y_pred


def cluster_acc(y_true, y_pred):
    y_true = y_true - np.min(y_true)
    l1 = list(set(y_true))
    num_class1 = len(l1)
    l2 = list(set(y_pred))
    num_class2 = len(l2)
    ind = 0

    if num_class1 != num_class2:
        for i in l1:
            if i in l2:
                pass
            else:
                y_pred[ind] = i
                ind += 1
    l2 = list(set(y_pred))
    numclass2 = len(l2)
    if num_class1 != numclass2:
        print('error')
        return
    cost = np.zeros((num_class1, numclass2), dtype=int)
    
    for i, c1 in enumerate(l1):
        mps = [i1 for i1, e1 in enumerate(y_true) if e1 == c1]
        for j, c2 in enumerate(l2):
            mps_d = [i1 for i1 in mps if y_pred[i1] == c2]
            cost[i][j] = len(mps_d)
    m = Munkres()
    cost = cost.__neg__().tolist()
    indexes = m.compute(cost)
    new_predict = np.zeros(len(y_pred))
    for i, c in enumerate(l1):
        c2 = l2[indexes[i][1]]
        ai = [ind for ind, elm in enumerate(y_pred) if elm == c2]
        new_predict[ai] = c
    acc = metrics.accuracy_score(y_true, new_predict)

    return acc


def eva(y_true, y_pred, show_details=True):

    acc = cluster_acc(y_true, y_pred)
    nmi = nmi_score(y_true, y_pred, average_method='arithmetic')
    ari = ari_score(y_true, y_pred)
    ami = adjusted_mutual_info_score(y_true, y_pred)
    if show_details:
        print("\n","ARI: {:.4f},".format(ari), "NMI: {:.4f},".format(nmi), "AMI: {:.4f}".format(ami), "ACC: {:.4f},".format(acc))
        
    return ari, nmi, ami, acc


def normalize2(adata, copy=True, highly_genes = None, filter_min_counts=True, size_factors=True, normalize_input=True, logtrans_input=True):
    if isinstance(adata, sc.AnnData):
        if copy:
            adata = adata.copy()
    elif isinstance(adata, str):
        adata = sc.read(adata)
    else:
        raise NotImplementedError
    # norm_error = 'Make sure that the dataset (adata.X) contains unnormalized count data.'
    # assert 'n_count' not in adata.obs, norm_error
    # if adata.X.size < 50e6: # check if adata.X is integer only if array is small
    #     if sp.sparse.issparse(adata.X):
    #         assert (adata.X.astype(int) != adata.X).nnz == 0, norm_error
    #     else:
    #         assert np.all(adata.X.astype(int) == adata.X), norm_error

    if filter_min_counts:
        sc.pp.filter_genes(adata, min_counts=1)
        sc.pp.filter_cells(adata, min_counts=1)
    if size_factors or normalize_input or logtrans_input:
        adata.raw = adata.copy()
    else:
        adata.raw = adata
    if size_factors:
        adata.X = adata.X.astype(float)
        sc.pp.normalize_per_cell(adata)
        adata.obs['size_factors'] = adata.obs.n_counts / np.median(adata.obs.n_counts)
    else:
        adata.obs['size_factors'] = 1.0
    if logtrans_input:
        sc.pp.log1p(adata)
    if highly_genes != None:
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes = highly_genes, subset=True)
    if normalize_input:
        sc.pp.scale(adata)
    return adata



def estimate_clusters(data, method='silhouette', max_clusters=50):
    """
    估计最佳聚类数量
    
    Parameters:
    -----------
    data : array-like
        特征矩阵
    method : str
        估计方法，可选 'silhouette', 'elbow', 'gap', 'leiden'
    max_clusters : int
        最大聚类数量
    
    Returns:
    --------
    n_clusters : int
        估计的最佳聚类数量
    """
    from sklearn.metrics import silhouette_score
    from sklearn.cluster import KMeans
    import numpy as np
    
    if method == 'silhouette':
        # 轮廓系数法
        silhouette_scores = []
        cluster_range = range(2, min(max_clusters, data.shape[0]//2))
        
        for n_clusters in cluster_range:
            kmeans = KMeans(n_clusters=n_clusters, random_state=42)
            cluster_labels = kmeans.fit_predict(data)
            if len(np.unique(cluster_labels)) < 2:
                continue
            silhouette_avg = silhouette_score(data, cluster_labels)
            silhouette_scores.append(silhouette_avg)
        
        if silhouette_scores:
            best_n = cluster_range[np.argmax(silhouette_scores)]
        else:
            best_n = 10  # 默认值
        
        return best_n
    
    elif method == 'elbow':
        # 肘部法则
        inertias = []
        cluster_range = range(1, min(max_clusters, data.shape[0]//2))
        
        for n_clusters in cluster_range:
            kmeans = KMeans(n_clusters=n_clusters, random_state=42)
            kmeans.fit(data)
            inertias.append(kmeans.inertia_)
        
        # 计算二阶差分寻找拐点
        diff = np.diff(inertias)
        diff2 = np.diff(diff)
        if len(diff2) > 0:
            best_n = np.argmax(np.abs(diff2)) + 3  # 加3是因为起始点
        else:
            best_n = 10
        
        return min(best_n, max_clusters)
    
    else:
        return 10  # 默认值
