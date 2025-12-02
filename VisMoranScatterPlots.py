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
import scipy.stats as stats


if __name__ == "__main__":
  h5ad_file = 'zfish_subset_stereoseq.h5ad' # data filepath 
  slice_names = range(1,13)

  fixed_stereo_slice = int(sys.argv[1])
  # Load in the saved genes file (sorted already)
  f = open("genes_sorted_by_variability_largest_to_smallest.txt", "r")
  lines = f.readlines()
  test_genes = []
  for gene in lines:
    test_genes.append(gene.strip())
  
  # Load the heatmaps, analyze
  ipf_data = []
  gene_name_data = []
  for gene_idx, gene in enumerate(test_genes[:2000]):
    #if gene in genes_that_pass_filter:
    heatmap = np.load(f"heatmap_matrices/ipf/ipf_{gene}_{gene_idx}_parameter_l_25_nonzero_allspots_final_2000genes_highly_variable_within_all_slices.npy")
    #max_score = np.max(heatmap)
    max_score = np.max(heatmap[fixed_stereo_slice,:])
    ipf_data.append((gene, max_score))
    gene_name_data.append(gene)
  
  ctfactomo_data = []
  for gene_idx, gene in enumerate(test_genes[:2000]):
    #if gene in genes_that_pass_filter:
    heatmap = np.load(f"heatmap_matrices/ctfactomo/CTFacTomo_{gene}_{gene_idx}_parameter_l_25_nonzero_allspots_final_2000genes_highly_variable_within_all_slices.npy")
    max_score = np.max(heatmap[fixed_stereo_slice,:])
    ctfactomo_data.append((gene, max_score))
  
  # Sort by max_score
  sorted_ipf_data = sorted(ipf_data, key=lambda x: x[1] if not math.isnan(x[1]) else float('inf'))
  f = open("analyze_ipf_all_spots_heatmaps.txt", "w")
  for r in sorted_ipf_data:
    print(r, file=f)
  f.close()
  
  sorted_ctfactomo_data = sorted(ctfactomo_data, key=lambda x: x[1] if not math.isnan(x[1]) else float('inf'))
  f = open("analyze_ctfactomo_all_spots_heatmaps.txt", "w")
  for r in sorted_ctfactomo_data:
    print(r, file=f)
  
  ctf_scores = []
  ipf_scores = []
  gene_names = []
  
  #num_ctf_better = 0
  
  for i in range(len(ctfactomo_data)):
    if math.isnan(ctfactomo_data[i][1]) or math.isnan(ipf_data[i][1]):
      continue
    else:
      ctf_scores.append(ctfactomo_data[i][1])
      ipf_scores.append(ipf_data[i][1])
      gene_names.append(gene_name_data[i])
    
  avg_ctf = np.mean(ctf_scores)
  avg_ipf = np.mean(ipf_scores)
  
  count_ctf_better = 0
  for i in range(len(ctf_scores)):
    if ctf_scores[i] > ipf_scores[i]:
      count_ctf_better += 1
  
  print(count_ctf_better)
  print(len(ctf_scores))
  
  print(avg_ctf)
  print(avg_ipf)
  
  # Get a p-value for CTFacTomo and IPF means
  diff = np.array(ctf_scores) - np.array(ipf_scores) 
  print("T-test:", stats.ttest_1samp(diff, 0))
  
  # Generate the scatter plot with gene name annotations
  plt.scatter(ctf_scores, ipf_scores, s=4)
  for i, s1  in enumerate(ctf_scores):
    for j, s2 in enumerate(ipf_scores):
      if i == j:
        plt.text(s1, s2, gene_names[i], fontsize=2, color='black')
  
  plt.plot([-0.1, 0.5], [-0.1, 0.5], 'k--', linewidth=0.5)
  plt.xlim(-0.1,0.5)
  plt.ylim(-0.1,0.5)
  plt.savefig(f"figures/stereo{fixed_stereo_slice}_ctfactomo_ipf_allspots_2000genes.png", dpi=300, bbox_inches='tight')
