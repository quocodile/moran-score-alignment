import math 
import sys
import align
import torch
import numpy as np
import paste as pst
import seaborn as sns
import matplotlib.pyplot as plt
import time
import scanpy as sc

def normalize_gene_expression(expression):
  return expression / np.linalg.norm(expression)

def apply_rbf_to_tensor(dist_matrix):
  l = 25 
  unnormalized_w = -1*(dist_matrix**2)/(2 * l**2)
  unnormalized_w = torch.exp(unnormalized_w) 
  return unnormalized_w

def get_spatial_weight_matrix(stereo_spots, reconstructed_spots, dist_matrix, top_expressed_spots):

  w_matrix = apply_rbf_to_tensor(dist_matrix)

  w_sum = torch.sum(w_matrix)

  # Element-wise normalization 
  w_matrix = w_matrix * (len(stereo_spots) / w_sum)

  return w_matrix 

def get_distance_matrix(stereo_spots, reconstructed_spots,top_expressed_spots, slice_n, reconstructed_slice_i):
  dist_matrix = torch.zeros((len(stereo_spots), len(reconstructed_spots)))
  for i, spot_stereo in enumerate(stereo_spots):
    for j, spot_reconstructed in enumerate(reconstructed_spots):
      if spot_reconstructed in top_expressed_spots:
        x1, y1 = spot_stereo
        x2, y2 = spot_reconstructed
        dij = ((x1 - x2)**2 + (y1-y2)**2)**0.5
        dist_matrix[i,j] = dij
      else:
        continue

  return dist_matrix 

def global_moran_matrix_mutliplication_version(spatial_weight_matrix, stereo_data, reconstructed_data):
  stereo_mean = np.mean(stereo_data.X)
  reconstructed_mean = np.mean(reconstructed_data.X)

  stereo_std = np.std(stereo_data.X)
  reconstructed_std = np.std(reconstructed_data.X)

  stereo_gene_expression = [exp[0] for exp in stereo_data.X] 
  reconstructed_gene_expression = [exp[0] for exp in reconstructed_data.X]

  stereo_zscore = (np.array(stereo_gene_expression) - stereo_mean) / stereo_std 
  reconstructed_zscore = (np.array(reconstructed_gene_expression) - reconstructed_mean) / reconstructed_std

  n = spatial_weight_matrix.shape[0]

  score = (stereo_zscore.T @ spatial_weight_matrix @ reconstructed_zscore) / n
  
  return score 

if __name__ == "__main__":
  h5ad_file = 'zfish_subset_stereoseq.h5ad' # data filepath 
  timing_out = open("top10highlyvariablegenes.txt", "w")
  slice_names = range(1,13)

  slices = align.load_slices(h5ad_file, [7])
  # Create dict to track stats of each gene
  test_genes = []
  stats_dict = {}
  for gene_idx, name in enumerate(slices[0].var_names):
    stats_dict[name] = {} 
    stats_dict[name]['sum'] = 0 
    stats_dict[name]['var'] = [] 
  for stereo_slice_i in range(1,13):
    print(stereo_slice_i)
    slices = align.load_slices(h5ad_file, [stereo_slice_i])
    for gene_idx, name in enumerate(slices[0].var_names):
      q = slices[0][:,name].X.toarray()
      stats_dict[name]['sum'] += np.sum(q) 
      stats_dict[name]['var'].append(np.var(q)) 
  for gene in stats_dict:
    test_genes.append((gene, stats_dict[gene]['sum'], np.mean(stats_dict[gene]['var'])))
    
  test_genes = sorted(test_genes, key=lambda x: x[2], reverse=True) 
  test_genes = [x[0] for x in test_genes[5000:]]
  #print(test_genes[:])
  #sys.exit()

  # Create a heatmap for each gene
  for gene_count, gene in enumerate(test_genes):
    # print(gene, gene_count)
    moran_heatmap = np.zeros((12, 12))
    #np.save(f"numpy_data/final/ipf_{gene}_{gene_count}_parameter_l_25_nonzero_x2_final_2001to3000genes_highly_variable_within_all_slices", moran_heatmap)
    np.save(f"numpy_data/final/CTFacTomo_{gene}_{gene_count}_parameter_l_25_nonzero_x2_final_5001toEndgenes_highly_variable_within_all_slices", moran_heatmap)

  start = time.time()
  for slice_n in slice_names:
    for reconstructed_slice_i in range(23,35):
      print(f"Aligning... stereo {slice_n} reconstructed {reconstructed_slice_i}")
      slices = align.load_slices(h5ad_file, [slice_n])
      slices[0].X = normalize_gene_expression(slices[0].X.toarray())
 
      # Loading the aligned slices
      
      # CTFacTomo aligned with Stereo-seq
      stereoseq_slice = sc.read_h5ad(f"alignments/ctfactomo/stereo_slice_{slice_n}_aligned_with_reconstructed_slice_{reconstructed_slice_i}.h5ad")
      reconstructed_slice = sc.read_h5ad(f"alignments/ctfactomo/reconstructed_slice_{reconstructed_slice_i}_aligned_with_stereo_slice_{slice_n}.h5ad")

      # IPF aligned with Stereo-seq
      #stereoseq_slice = sc.read_h5ad(f"alignments/ipf/stereo_slice_{slice_n}_aligned_with_reconstructed_slice_{reconstructed_slice_i}.h5ad")
      #reconstructed_slice = sc.read_h5ad(f"alignments/ipf/reconstructed_slice_{reconstructed_slice_i}_aligned_with_stereo_slice_{slice_n}.h5ad")

      for gene_count, gene in enumerate(test_genes):
        print("Gene:", gene, "Num:", gene_count)

        # Skip this slice if less than 5 spots of expression
        #print("Nonzero counts in stereo:", np.count_nonzero(slices[0][:,gene].X))
        if np.count_nonzero(slices[0][:,gene].X) < 5:
          continue

        #moran_heatmap = np.load(f"numpy_data/final/ipf_{gene}_{gene_count}_parameter_l_25_nonzero_x2_final_2001to3000genes_highly_variable_within_all_slices.npy")
        
        moran_heatmap = np.load(f"numpy_data/final/CTFacTomo_{gene}_{gene_count}_parameter_l_25_nonzero_x2_final_5001toEndgenes_highly_variable_within_all_slices.npy")

        stereoseq_slice_gene = stereoseq_slice[:,gene]
        reconstructed_slice_gene = reconstructed_slice[:,gene]
        stereo_spots = stereoseq_slice_gene.obsm['spatial'];
        rec_spots = reconstructed_slice_gene.obsm['spatial'];

        # Load the w matrix, either for ipf or CTFacTomo, with all spots included
        #w_matrix = np.load(f"w_matrices/ipf/all_spots/unnormalized/w_matrix_stereo_slice_{slice_n}_reconstructed_slice_{reconstructed_slice_i}.npy")
        w_matrix = np.load(f"w_matrices/ctfactomo/all_spots/unnormalized/w_matrix_stereo_slice_{slice_n}_reconstructed_slice_{reconstructed_slice_i}.npy")

        # Match the number of non-zero spots
        stereo_ratio_nonzero = np.count_nonzero(slices[0][:, gene].X) / len(slices[0][:,gene].X) * 2 
        rec_spatial = []
        rec_expression = []
        for data in reconstructed_slice_gene:
          rec_spatial.append(data.obsm['spatial'][0])
          rec_expression.append(data.X.flatten()[0])

        # Choosing which spots to consider from the reconstructed data
        #print("Ratio:", stereo_ratio_nonzero * len(rec_spatial))
        rec_sort_indices = np.argsort(rec_expression)[::-1][:int(stereo_ratio_nonzero * len(rec_spatial))]
        #rec_sort_indices = np.argsort(rec_expression)[::-1]
        top_expressed_spots = np.array(rec_spatial)[rec_sort_indices]
       
        '''
        rec_x = []
        rec_y = []
        rec_val = []
        for i, p, in enumerate(rec_spots):
          rec_x.append(p[0])
          rec_y.append(p[1])
          rec_val.append(1)
        #stereo_x = []
        #stereo_y = []
        #stereo_val = []
        #for i, p, in enumerate(stereo_spots):
        #  stereo_x.append(p[0])
        #  stereo_y.append(p[1])
        #  stereo_val.append(1)
        rec_x2 = []
        rec_y2 = []
        rec_val2 = []
        for i, p, in enumerate(rec_spots):
          if p in top_expressed_spots:
            rec_x2.append(p[0])
            rec_y2.append(p[1])
            rec_val2.append(1)
      
        plt.figure()
        plt.scatter(rec_x, rec_y, c='black', marker='o', alpha=1, label='Data Points')
        plt.scatter(rec_x2, rec_y2, c='r', marker='o', alpha=1, label='Data Points')
        plt.savefig(f"figures/rec_relevant_spots_{slice_n}_{reconstructed_slice_i}_ipf_{gene}_x3.png", dpi=300, bbox_inches='tight')  # Save as PNG with high resolution
        print(f"figures/rec_relevant_spots_{slice_n}_{reconstructed_slice_i}_ipf_{gene}_x3.png")
        '''
        # Filter entries of the w matrix to only top-expressed spots
        for i in range(w_matrix.shape[1]):
          if i not in rec_sort_indices:
            w_matrix[:,i] = 0

        # Normalize the w matrix
        w_sum = np.sum(w_matrix)
        w_matrix = w_matrix * (w_matrix.shape[0] / w_sum)
        
        global_moran_score = global_moran_matrix_mutliplication_version(w_matrix, stereoseq_slice_gene, reconstructed_slice_gene)

        print(slice_n, reconstructed_slice_i)
        print(gene, global_moran_score.item())
        moran_heatmap[slice_n-1][reconstructed_slice_i-23] = global_moran_score.item()
        
        #np.save(f"numpy_data/final/ipf_{gene}_{gene_count}_parameter_l_25_nonzero_x2_final_2001to3000genes_highly_variable_within_all_slices", moran_heatmap)
        np.save(f"numpy_data/final/CTFacTomo_{gene}_{gene_count}_parameter_l_25_nonzero_x2_final_5001toEndgenes_highly_variable_within_all_slices", moran_heatmap)
        
  end = time.time()
  print("Total Elapsed Time:", end-start)
  print("Total Elapsed Time:", end-start, file=timing_out)
  timing_out.close()
