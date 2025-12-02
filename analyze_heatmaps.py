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

h5ad_file = 'zfish_subset_stereoseq.h5ad' # data filepath 
slice_names = range(1,13)
fixed_stereo_slice  = 11 

'''
test_genes = []
slices = align.load_slices(h5ad_file, [7])
# Create dict to track stats of each gene
stats_dict = {}
for gene_idx, name in enumerate(slices[0].var_names):
  stats_dict[name] = {} 
  stats_dict[name]['sum'] = 0 
  stats_dict[name]['var'] = [] 

# Sorts the genes by their mean variability across all slices
for stereo_slice_i in range(1,13):
  print(stereo_slice_i)
  slices = align.load_slices(h5ad_file, [stereo_slice_i])
  for gene_idx, name in enumerate(slices[0].var_names):
    q = slices[0][:,name].X.toarray()
    stats_dict[name]['sum'] += np.sum(q) 
    stats_dict[name]['var'].append(np.var(q)) 
for gene in stats_dict:
  test_genes.append((gene, stats_dict[gene]['sum'], np.mean(stats_dict[gene]['var'])))
  
coloring = [x[1] for x in test_genes[:1000]]
test_genes = sorted(test_genes, key=lambda x: x[2], reverse=True) 
test_genes = [x[0] for x in test_genes]

genes_file = open("genes_sorted_by_variability_largest_to_smallest.txt", "w")
for g in test_genes:
  print(g, file=genes_file)
genes_file.close()

sys.exit()
'''


# Load in the saved genes file (sorted already)
f = open("genes_sorted_by_variability_largest_to_smallest.txt", "r")
lines = f.readlines()
test_genes = []
for gene in lines:
  test_genes.append(gene.strip())

'''
# Filter the gene if it doesn't have at least one slice with 10 spots expressed
genes_that_pass_filter = []
for stereo_slice_i in range(1,13):
  print(stereo_slice_i)
  slices = align.load_slices(h5ad_file, [stereo_slice_i])
  for gene_idx, name in enumerate(slices[0].var_names):
    if np.count_nonzero(slices[0][:,name].X.toarray()) > -1:
      genes_that_pass_filter.append(name)

print(len(test_genes))
genes_that_pass_filter = set(genes_that_pass_filter)
print(len(genes_that_pass_filter))

test_genes_filtered = []

for gene in test_genes:
  if gene in genes_that_pass_filter:
    test_genes_filtered.append(gene)
print(len(test_genes_filtered))
'''

# Load the heatmaps, analyze
ipf_data = []
gene_name_data = []
for gene_idx, gene in enumerate(test_genes[:2000]):
  #if gene in genes_that_pass_filter:
  heatmap = np.load(f"numpy_data/final/ipf_{gene}_{gene_idx}_parameter_l_25_nonzero_allspots_final_2000genes_highly_variable_within_all_slices.npy")
  #max_score = np.max(heatmap)
  max_score = np.max(heatmap[fixed_stereo_slice,:])
  ipf_data.append((gene, max_score))
  gene_name_data.append(gene)

ctfactomo_data = []
for gene_idx, gene in enumerate(test_genes[:2000]):
  #if gene in genes_that_pass_filter:
  heatmap = np.load(f"numpy_data/final/CTFacTomo_{gene}_{gene_idx}_parameter_l_25_nonzero_allspots_final_2000genes_highly_variable_within_all_slices.npy")
  max_score = np.max(heatmap[fixed_stereo_slice,:])
  ctfactomo_data.append((gene, max_score))

'''
for gene_idx, gene in enumerate(test_genes[2000:3000]):
  #if gene in genes_that_pass_filter:
  heatmap = np.load(f"numpy_data/final/ipf_{gene}_{gene_idx}_parameter_l_25_nonzero_allspots_final_2001to3000genes_highly_variable_within_all_slices.npy")
  max_score = np.max(heatmap)
  ipf_data.append((gene, max_score))

for gene_idx, gene in enumerate(test_genes[2000:3000]):
  #if gene in genes_that_pass_filter:
  heatmap = np.load(f"numpy_data/final/CTFacTomo_{gene}_{gene_idx}_parameter_l_25_nonzero_allspots_final_2000to3000genes_highly_variable_within_all_slices.npy")
  max_score = np.max(heatmap)
  ctfactomo_data.append((gene, max_score))

for gene_idx, gene in enumerate(test_genes[3000:4000]):
  #if gene in genes_that_pass_filter:
  heatmap = np.load(f"numpy_data/final/ipf_{gene}_{gene_idx}_parameter_l_25_nonzero_allspots_final_3001to4000genes_highly_variable_within_all_slices.npy")
  max_score = np.max(heatmap)
  ipf_data.append((gene, max_score))

for gene_idx, gene in enumerate(test_genes[3000:4000]):
  #if gene in genes_that_pass_filter:
  heatmap = np.load(f"numpy_data/final/CTFacTomo_{gene}_{gene_idx}_parameter_l_25_nonzero_allspots_final_3001to4000genes_highly_variable_within_all_slices.npy")
  max_score = np.max(heatmap)
  ctfactomo_data.append((gene, max_score))

for gene_idx, gene in enumerate(test_genes[4000:5000]):
  #if gene in genes_that_pass_filter:
  heatmap = np.load(f"numpy_data/final/ipf_{gene}_{gene_idx}_parameter_l_25_nonzero_allspots_final_4001to5000genes_highly_variable_within_all_slices.npy")
  max_score = np.max(heatmap)
  ipf_data.append((gene, max_score))

for gene_idx, gene in enumerate(test_genes[4000:5000]):
  #if gene in genes_that_pass_filter:
  heatmap = np.load(f"numpy_data/final/CTFacTomo_{gene}_{gene_idx}_parameter_l_25_nonzero_allspots_final_4001to5000genes_highly_variable_within_all_slices.npy")
  max_score = np.max(heatmap)
  ctfactomo_data.append((gene, max_score))

for gene_idx, gene in enumerate(test_genes[5000:]):
  #if gene in genes_that_pass_filter:
  heatmap = np.load(f"numpy_data/final/ipf_{gene}_{gene_idx}_parameter_l_25_nonzero_allspots_final_5001toEndgenes_highly_variable_within_all_slices.npy")
  max_score = np.max(heatmap)
  ipf_data.append((gene, max_score))

for gene_idx, gene in enumerate(test_genes[5000:]):
  #if gene in genes_that_pass_filter:
  heatmap = np.load(f"numpy_data/final/CTFacTomo_{gene}_{gene_idx}_parameter_l_25_nonzero_allspots_final_5001toEndgenes_highly_variable_within_all_slices.npy")
  max_score = np.max(heatmap)
  ctfactomo_data.append((gene, max_score))
'''

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

diff = np.array(ctf_scores) - np.array(ipf_scores) 

print("T-test:", stats.ttest_1samp(diff, 0))

plt.scatter(ctf_scores, ipf_scores, s=4, cmap='viridis')
for i, s1  in enumerate(ctf_scores):
  for j, s2 in enumerate(ipf_scores):
    if i == j:
      plt.text(s1, s2, gene_names[i], fontsize=2, color='black')
plt.plot([-0.1, 0.5], [-0.1, 0.5], 'k--', linewidth=0.5)
plt.xlim(-0.1,0.5)
plt.ylim(-0.1,0.5)
plt.savefig(f"figures/stereo{fixed_stereo_slice}_ctfactomo_ipf_allspots_2000genes.png", dpi=300, bbox_inches='tight')
