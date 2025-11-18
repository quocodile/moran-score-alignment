import numpy as np
import matplotlib.pyplot as plt
import align
import paste as pst
import colorcet as cc

print

def normalize_gene_expression(expression):
  return expression / np.linalg.norm(expression)

def load_ctfactomo_slice(stereo_slice_num, reconstructed_slice_num):
  return align.load_reconstructed_slice(f"reconstructed_zfish_shield/ctfactomo_zfish_shield_genes_matched_with_stereoseq_slice_{reconstructed_slice_num}.npy", [stereo_slice_num], reconstructed_slice_num)

def load_ipf_slice(stereo_slice_num, reconstructed_slice_num):
  return align.load_reconstructed_slice(f"reconstructed_zfish_shield/ipf_zfish_shield_genes_matched_with_stereoseq_slice_{reconstructed_slice_num}.npy", [stereo_slice_num], reconstructed_slice_num)

h5ad_file = 'zfish_subset_stereoseq.h5ad' # data filepath 

# Select the configuration 
gene = 'gsc'
reconstructed_slice_number = 32 
stereo_slice_number = 2 
#reconstruction_type = 'IPF'
reconstruction_type = 'CTFacTomo'

# Load the slices 
slices = align.load_slices(h5ad_file, [stereo_slice_number])
if reconstruction_type == 'CTFacTomo':
  # CTFacTomo
  reconstructed_slice = load_ctfactomo_slice(stereo_slice_number, reconstructed_slice_number)
elif reconstruction_type == 'IPF':
  # IPF
  reconstructed_slice = load_ipf_slice(stereo_slice_number, reconstructed_slice_number)
else:
  reconstructed_slice = None

reconstructed_slice.X = normalize_gene_expression(reconstructed_slice.X)
slices[0].X = normalize_gene_expression(slices[0].X.toarray())
slices.append(reconstructed_slice)

# Pairwise align the slices
pis = []
pi = pst.pairwise_align(slices[0], slices[1])
pis.append(pi)

# To visualize the alignment you can stack 
# the slices according to the alignment pi
new_slices = pst.stack_slices_pairwise(slices, pis)
stereoseq_slice = new_slices[0] 
reconstructed_slice = new_slices[1] 

reconstructed_slice = reconstructed_slice[:,gene]
stereoseq_slice = stereoseq_slice[:,gene]

# Visualize the reconstructed slice
rec_x = []
rec_y = []
rec_val = []
rec_nonzero_x = []
rec_nonzero_y = []
rec_zero_x = []
rec_zero_y = []
all_vals = []
for i, r, in enumerate(stereoseq_slice):
  p = r.obsm['spatial']
  rec_x.append(p[0][0])
  rec_y.append(p[0][1])
  rec_val.append(r.X.flatten()[0])
  all_vals.append(r.X.flatten()[0])
  if r.X.flatten()[0] == 0:
    rec_zero_x.append(p[0][0])
    rec_zero_y.append(p[0][1])
  else:
    rec_nonzero_x.append(p[0][0])
    rec_nonzero_y.append(p[0][1])

# Selecting the top spots to visualize
stereo_ratio_nonzero = np.count_nonzero(slices[0][:, gene].X) / len(slices[0][:,gene].X) * 2 
#rec_sort_indices = np.argsort(rec_val)[::-1][:int(stereo_ratio_nonzero * len(all_vals))]
rec_sort_indices = np.argsort(rec_val)[::-1][:cutoff_expression_idx]
#not_included_rec_sort_indices = np.argsort(rec_val)[::-1][int(stereo_ratio_nonzero * len(all_vals)):]
#rec_sort_indices = np.argsort(rec_expression)[::-1]
top_expressed_x = np.array(rec_nonzero_x)[rec_sort_indices]
top_expressed_y = np.array(rec_nonzero_y)[rec_sort_indices]

plt.figure()
plt.axis('off')
plt.scatter(rec_nonzero_x, rec_nonzero_y, c=rec_val, edgecolors='black', cmap='Greens', marker='o', linewidth=0.5, alpha=1, s=25, label='Data Points', vmin=min(all_vals))
plt.scatter(rec_zero_x, rec_zero_y, facecolors='white', edgecolors='lightgray', marker='o', linewidth=1.1, alpha=0.5, s=25, label='Data Points')
plt.scatter(top_expressed_x, top_expressed_y, c='red', edgecolors='black', marker='o', linewidth=0.5, alpha=1, s=25, label='Data Points', vmin=0)
plt.savefig(f"figures/allspots_{reconstruction_type}_slice{reconstructed_slice_number}_{gene}_aligned_with_Stereoseq_{stereo_slice_number}.png", dpi=300, bbox_inches='tight')  # Save as PNG with high resolution

# Visualize the Stereo-seq slice
plt.clf()
plt.axis('off')
stereo_x = []
stereo_y = []
stereo_val = []
stereo_zero_x = []
stereo_zero_y = []
for i, r, in enumerate(stereoseq_slice):
  p = r.obsm['spatial']
  stereo_x.append(p[0][0])
  stereo_y.append(p[0][1])
  stereo_val.append(r.X.flatten()[0])
  if r.X.flatten()[0] == 0:
    stereo_zero_x.append(p[0][0])
    stereo_zero_y.append(p[0][1])

sorted_vals = sorted(all_vals, reverse=True) 
plt.scatter(stereo_x, stereo_y, c=stereo_val, edgecolors='black', cmap='Greens', marker='o', linewidth=0.5, alpha=1, s=14, label='Data Points', vmin=0)
plt.scatter(stereo_zero_x, stereo_zero_y, facecolors='white', edgecolors='lightgray', marker='o', linewidth=1.1, alpha=1, s=15, label='Data Points')
plt.savefig(f"figures/allspots_Stereoseq_slice{stereo_slice_number}_{gene}_aligned_with_{reconstruction_type}_{reconstructed_slice_number}.png", dpi=300, bbox_inches='tight')  # Save as PNG with high resolution
