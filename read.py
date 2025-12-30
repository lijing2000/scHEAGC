#import anndata as ad
#import numpy as np
#
## 读取数据
#adatarna = ad.read_h5ad("./input/PBMC-10k/Pbmc10k-RNA.h5ad")

## 按出现顺序生成标签
#cell_types = adatarna.obs['cell_type'].values
#type_to_label = {}
#labels = []
#next_label = 0
#
#for cell_type in cell_types:
#    if cell_type not in type_to_label:
#        type_to_label[cell_type] = next_label
#        next_label += 1
#    labels.append(type_to_label[cell_type])
#
#true_labels = np.array(labels)
#
## 保存到label.npy文件
#np.save("./input/Ma/label.npy", true_labels)
#
#print(f"标签已保存到: ./input/Ma/label.npy")
#print(f"标签形状: {true_labels.shape}")
#print(f"标签范围: {np.min(true_labels)} 到 {np.max(true_labels)}")
#print(f"唯一标签数: {len(np.unique(true_labels))}")
#print(f"前10个标签示例: {true_labels[:10]}")



#import scanpy as sc
#import pandas as pd
#
## ---------------------- 1. 加载AnnData数据 ----------------------
## 替换为你的h5ad文件路径
#adata = sc.read_h5ad("./input/PBMC-10k/Pbmc10k-RNA.h5ad")
#
## 查看obs列名，确认cell_type和数字标签列（关键！）
#print("Obs列名列表：")
#print(adata.obs.columns.tolist())
#
## ---------------------- 2. 确认数字标签列 ----------------------
## 你需要先确认：数字标签存在哪个列？（比如你可能命名为cell_type_label/seurat_clusters/sub.cluster）
## 👇 替换为你的数字标签列名（比如seurat_clusters/sub.cluster，或自定义的cell_type_num）
#num_label_col = "seurat_clusters"  # 示例：如果数字标签存在seurat_clusters列
#
## ---------------------- 3. 提取cell_type与数字标签的映射 ----------------------
## 方法1：按cell_type首次出现顺序生成映射（符合你“出现先后”的规则）
## 提取所有唯一的cell_type（按首次出现顺序）
#unique_cell_types = adata.obs["cell_type"].drop_duplicates().tolist()
## 生成“细胞类型→数字标签”字典（按首次出现顺序从0开始编号）
#cell_type_to_num = {ct: idx for idx, ct in enumerate(unique_cell_types)}
#
## 方法2：如果数字标签已存在（比如seurat_clusters），提取现有映射（去重）
## 生成cell_type和数字标签的对应表（去重）
#mapping_df = adata.obs[["cell_type", num_label_col]].drop_duplicates()
## 转换为字典（数字标签→细胞类型）
#num_to_cell_type = dict(zip(mapping_df[num_label_col], mapping_df["cell_type"]))
## 转换为字典（细胞类型→数字标签）
#cell_type_to_num_exist = dict(zip(mapping_df["cell_type"], mapping_df[num_label_col]))
#
## ---------------------- 4. 输出映射结果 ----------------------
#print("\n=== 细胞类型 → 数字标签 映射表 ===")
## 方法1的结果（按首次出现顺序）
#print("（按cell_type首次出现顺序生成的映射）")
#for ct, num in cell_type_to_num.items():
#    print(f"{num}\t→ {ct}")
#
## 方法2的结果（基于现有数字标签列）
#print("\n（基于现有数字标签列的映射）")
#for num, ct in num_to_cell_type.items():
#    print(f"{num}\t→ {ct}")
#
## 可选：保存映射表到CSV文件
#mapping_df.to_csv("cell_type_num_mapping.csv", index=False)
#print("\n映射表已保存到 cell_type_num_mapping.csv")






#import pandas as pd
#import numpy as np
#import os
#
#def run_cellmarker_annotation_xlsx(deg_dict, database_path, species='Human'):
#    """
#    针对 PBMC-10k 优化的 CellMarker 2.0 自动匹配脚本
#    """
#    if not os.path.exists(database_path):
#        print(f"错误：找不到文件 {database_path}")
#        return None
#
#    print(f"正在加载数据库: {database_path}...")
#    
#    # 1. 加载文件
#    if database_path.endswith('.xlsx'):
#        df_db = pd.read_excel(database_path, engine='openpyxl')
#    else:
#        df_db = pd.read_csv(database_path)
#
#    # 2. 优化筛选条件：针对 PBMC 扩大组织搜索范围，确保找回 HSPC 和 DC 亚型
#    # 搜索范围包括：血液、外周血、血清、以及骨髓（很多免疫细胞起源于此）
#    tissue_keywords = 'Blood|Peripheral blood|Serum|Bone marrow|PBMC'
#    
#    mask = (df_db['species'].str.contains(species, case=False, na=False)) & \
#           (df_db['tissue_type'].str.contains(tissue_keywords, case=False, na=False))
#    df_db_filtered = df_db[mask].copy()
#
#    if df_db_filtered.empty:
#        print("警告：未找到匹配的物种或组织，请检查数据库列名是否为 'species' 和 'tissue_type'。")
#        return None
#
#    # 3. 执行基因匹配
#    annotation_results = []
#
#    for cluster_id, my_genes in deg_dict.items():
#        # 清洗查询基因
#        my_genes_set = set([str(g).strip().upper() for g in my_genes if pd.notna(g)])
#        
#        for _, row in df_db_filtered.iterrows():
#            cell_type = row['cell_name']
#            raw_markers = str(row['marker'])
#            
#            if raw_markers == 'nan' or not raw_markers:
#                continue
#                
#            # 清洗数据库基因（处理多种分隔符）
#            db_markers = set([m.strip().upper() for m in raw_markers.replace(',', ' ').replace(';', ' ').split()])
#            
#            intersection = my_genes_set.intersection(db_markers)
#            
#            if len(intersection) > 0:
#                # 记录每一个命中的记录
#                for gene in intersection:
#                    annotation_results.append({
#                        'Cluster': cluster_id,
#                        'CellMarker_CellType': cell_type,
#                        'Gene': gene
#                    })
#
#    # 4. 数据聚合统计
#    raw_res = pd.DataFrame(annotation_results)
#    if raw_res.empty:
#        print("未发现匹配项。")
#        return None
#
#    # 核心步骤：按 Cluster 和 细胞类型分组，统计不重复基因的个数
#    summary = raw_res.groupby(['Cluster', 'CellMarker_CellType']).agg(
#        Match_Score=('Gene', 'nunique'),
#        Matched_Genes=('Gene', lambda x: ', '.join(sorted(set(x))))
#    ).reset_index()
#
#    # 计算该簇总共输入了多少个基因供评审参考 (这里假设每个簇输入数量一致)
#    summary['Total_Query_Genes'] = summary['Cluster'].map(lambda x: len(deg_dict[x]))
#
#    # 5. 排序：每个 Cluster 匹配得分最高的细胞类型排在前面
#    summary = summary.sort_values(['Cluster', 'Match_Score'], ascending=[True, False])
#    
#    # 每个簇仅保留前 5 个最可能的匹配结果（按分数排）
#    final_report = summary.groupby('Cluster').head(5)
#    
#    return final_report
#
## --- 数据输入 ---
#my_degs = {
#    'Cluster 0': ['TCF7L2', 'MTSS1', 'FCGR3A', 'PSAP', 'IFITM3'],
#    'Cluster 1': ['INPP4B', 'ANK3', 'IL32', 'ITGB1', 'CDC14A'],
#    'Cluster 2': ['LEF1', 'NELL2', 'BACH2', 'PDE3B', 'THEMIS'],
#    'Cluster 4': ['BANK1', 'EBF1', 'MS4A1', 'RALGPS2', 'OSBPL10'],
#    'Cluster 5': ['LTB', 'RORA', 'IL7R', 'RPS2', 'EEF1A1'],
#    'Cluster 6': ['CCL5', 'NKG7', 'GZMA', 'A2M', 'IL32'],
#    'Cluster 7': ['DPYD', 'PLXDC2', 'ARHGAP26', 'NEAT1', 'LRMDA'],
#    'Cluster 8': ['GNLY', 'PRF1', 'NKG7', 'CD247', 'KLRD1'],
#    'Cluster 9': ['PSAP', 'CST3', 'IFI30', 'CALHM6', 'HLA-DRB1'],
#    'Cluster 10': ['LEF1', 'FHIT', 'CAMK4', 'BCL11B', 'ANK3'],
#    'Cluster 11': ['RALGPS2', 'BANK1', 'OSBPL10', 'MS4A1', 'CD79A'],
#    'Cluster 13': ['IGHM', 'BANK1', 'FCRL1', 'AFF3', 'PAX5'],
#    'Cluster 14': ['SLC4A10', 'PHACTR2', 'NKG7', 'KLRB1', 'RORA'],
#    'Cluster 15': ['CD74', 'CST3', 'HLA-DPB1', 'HLA-DRB1', 'HLA-DRA'],
#    'Cluster 17': ['TCF4', 'FCHSD2', 'RHEX', 'AUTS2', 'BCL11A'],
#    'Cluster 18': ['IKZF2', 'RTKN2', 'IL32', 'CASK', 'IL2RA']
#}
#
## --- 运行与保存 ---
## 请确保 Cell_marker_All.xlsx 在当前目录下
#result_summary = run_cellmarker_annotation_xlsx(my_degs, "Cell_marker_All.xlsx")
#
#if result_summary is not None:
#    print("\n--- 匹配结果预览 (Top 10) ---")
#    print(result_summary.head(10))
#    result_summary.to_csv("annotation_summary_v2.csv", index=False)
#    print("\n分析完成！结果已聚合保存至: annotation_summary_v2.csv")







import h5py
import numpy as np
import pandas as pd
import anndata as ad

#def split_multiome_h5_to_h5ad(input_file, rna_output_file, atac_output_file):
#    """
#    读取一个非标准的10x Multiome HDF5文件，并将其拆分为RNA和ATAC两个AnnData对象，
#    然后分别保存为h5ad文件。
#
#    参数:
#    input_file (str): 输入的HDF5文件路径。
#    rna_output_file (str): 输出的RNA数据h5ad文件路径。
#    atac_output_file (str): 输出的ATAC数据h5ad文件路径。
#    """
#    print(f"开始处理文件: {input_file}")
#
#    # 使用 'r' 模式（只读）打开HDF5文件
#    with h5py.File(input_file, 'r') as f:
#        # --- 1. 读取共享信息 ---
#        print("正在读取共享信息 (Barcodes, Celltypes)...")
#        # h5py读取的字符串是bytes类型，需要解码为UTF-8
#        barcodes = [b.decode('utf-8') for b in f['Barcodes'][:]]
#        celltypes = [ct.decode('utf-8') for ct in f['Celltypes'][:]]
#        
#        # 创建一个用于存储细胞注释的DataFrame
#        obs_df = pd.DataFrame(index=barcodes)
#        obs_df['cell_type'] = celltypes
#
#        # --- 2. 处理RNA数据 ---
#        print("\n正在处理RNA数据...")
#        if 'X1' in f and 'Genes' in f:
#            # 读取表达矩阵和基因名
#            rna_matrix = f['X1'][:]
#            gene_names = [g.decode('utf-8') for g in f['Genes'][:]]
#            
#            # 创建RNA的AnnData对象
#            adata_rna = ad.AnnData(
#                X=rna_matrix,
#                obs=obs_df.copy(),  # 使用共享的obs DataFrame
#                var=pd.DataFrame(index=gene_names)
#            )
#            
#            # 可选：添加一些元数据
#            adata_rna.uns['data_type'] = 'Gene Expression (RNA)'
#            
#            print(f"RNA AnnData 对象创建成功: {adata_rna}")
#            print(f"  - 细胞数: {adata_rna.n_obs}")
#            print(f"  - 基因数: {adata_rna.n_vars}")
#            
#            # 保存为h5ad文件
#            adata_rna.write_h5ad(rna_output_file)
#            print(f"RNA数据已保存至: {rna_output_file}")
#        else:
#            print("错误：未在HDF5文件中找到RNA数据 (X1, Genes)。")
#
#        # --- 3. 处理ATAC数据 ---
#        print("\n正在处理ATAC数据...")
#        if 'X2' in f and 'GeneFromPeaks' in f:
#            # 读取可及性矩阵和Peak特征名
#            atac_matrix = f['X2'][:]
#            peak_names = [p.decode('utf-8') for p in f['GeneFromPeaks'][:]]
#            
#            # 创建ATAC的AnnData对象
#            adata_atac = ad.AnnData(
#                X=atac_matrix,
#                obs=obs_df.copy(),  # 使用共享的obs DataFrame
#                var=pd.DataFrame(index=peak_names)
#            )
#            
#            # 可选：添加一些元数据
#            adata_atac.uns['data_type'] = 'Chromatin Accessibility (ATAC)'
#            
#            print(f"ATAC AnnData 对象创建成功: {adata_atac}")
#            print(f"  - 细胞数: {adata_atac.n_obs}")
#            print(f"  - Peak特征数: {adata_atac.n_vars}")
#            
#            # 保存为h5ad文件
#            adata_atac.write_h5ad(atac_output_file)
#            print(f"ATAC数据已保存至: {atac_output_file}")
#        else:
#            print("错误：未在HDF5文件中找到ATAC数据 (X2, GeneFromPeaks)。")
#            
#    print("\n处理完成！")
#
#
## --- 主程序入口 ---
#if __name__ == "__main__":
#    # 定义你的文件路径
#    input_h5_file = './input/pbmc3k2000/10XMultiomics_pbmc_3k_granulocyte_plus.h5'
#    
#    # 定义输出的h5ad文件名
#    rna_h5ad_file = './input/pbmc3k2000/pbmc3k_rna.h5ad'
#    atac_h5ad_file = './input/pbmc3k2000/pbmc3k_atac.h5ad'
#    
#    # 调用函数执行拆分
#    split_multiome_h5_to_h5ad(input_h5_file, rna_h5ad_file, atac_h5ad_file)



rna = ad.read_h5ad("./input/PBMC-10k/Pbmc10k-RNA.h5ad")
atac = ad.read_h5ad("./input/PBMC-10k/Pbmc10k-ATAC.h5ad")

print(rna)
print(atac)