import numpy as np
import scanpy as sc
import paste as pst
import anndata as ad
from scipy.spatial import KDTree
from scipy.stats import pearsonr 
import sys

# Load Slices
# h5ad_file = 'zfish_subset_stereoseq.h5ad' # change this path to the data you wish to analyze
# mask_file = 'zebrafish_mask.npy'

def normalize_gene_expression(expression):
  return expression / np.linalg.norm(expression)

def load_slice(h5ad_file, slice_name=None):
    data = sc.read_h5ad(h5ad_file)
    if slice_name:
      cur_slice = data[data.obs["slice"] == slice_name] 
    else:
      cur_slice = data
    slice_x = cur_slice.obs["spatial_x"] 
    slice_y = cur_slice.obs["spatial_y"] 
    spatial = []
    for i in range(len(slice_x)):
      spatial.append([slice_x[i], slice_y[i]]) 
    cur_slice.obsm['spatial'] = np.array(spatial) 
    return cur_slice 

def get_min_spatial_bounds(h5ad_file, slice_names):
    data = sc.read_h5ad(h5ad_file)
    x = []
    y = []
    for slice_name in slice_names:
        slice_i = data[data.obs["slice"] == slice_name] 
        slice_i_x = slice_i.obs["spatial_x"] 
        slice_i_y = slice_i.obs["spatial_y"] 
        for i in range(len(slice_i_x)):
          x.append(slice_i_x[i])
          y.append(slice_i_y[i])
    return min(x), min(y) 

def load_reconstructed_slice(filename, slice_names, reconstructed_slice_num, flipped=False):
  
  genes = [x.strip() for x in open("match.txt", "r").readlines()]
  sorted_genes = sorted(genes)
  sort_indices = np.argsort(genes)
  np_data = np.load(filename)
  mask = np.load(mask_file)
  slice_mask = mask[:,:,reconstructed_slice_num-1]
  if flipped == True:
    print("Flipped!")
    np_data = np_data[:,:,::-1]
    slice_mask = slice_mask[:,::-1]
  slice_mask_nonzero_coords = []
  expression_within_mask = []
  for i in range(slice_mask.shape[0]):
    for j in range(slice_mask.shape[1]):
      if slice_mask[i][j] == 1:
        slice_mask_nonzero_coords.append((i, j)) 
        expression_within_mask.append(np_data[:,i,j])
  expression_within_mask = np.array(expression_within_mask)
  print(expression_within_mask[0][0], expression_within_mask[0][-1])
  print("shape", slice_mask.shape)
  print(expression_within_mask.shape)
  print(len(slice_mask_nonzero_coords))
  # anndata_input = np_data.T.reshape(-1, 5439)[:, sort_indices]
  anndata_input = expression_within_mask[:, sort_indices]
  #anndata_input = []
  #for i, val in enumerate(np_data):
  #  anndata_input.append([[i,i], val])
  anndata_input = np.array(anndata_input)
  data = ad.AnnData(X=anndata_input)
  spatial = []
  print("Slice names:", slice_names)
  min_x, min_y = get_min_spatial_bounds(h5ad_file, slice_names)
  # print()
  # print(min_x)
  # print(min_y)
  # print()
  '''
  for i in range(49):
    for j in range(50):
      spatial.append([5640 + i*25, 15315 + j*25])
      # spatial.append([min_x + i*25, min_y + j*25])
  '''
  for i in range(len(expression_within_mask)):
    spatial.append([5640 + slice_mask_nonzero_coords[i][0]*25, 15315 + slice_mask_nonzero_coords[i][1]*25])
  print("spatial", len(spatial))
  data.obsm['spatial'] = np.array(spatial)
  data.var_names = sorted_genes 

  # sc.pp.filter_genes(data, min_counts = 1e-5)
  # for row in data.X:
  #  print(np.sum(row))
  # print(data.X.shape)
  # sc.pp.filter_cells(data, min_counts = 5)
  # print(data)
  # print(data.X.shape)
  return data


if __name__ == "__main__":

  # Read h5ad files
  filename1 = sys.argv[1]
  filename2 = sys.argv[2]

  slice1 = load_slice(filename1)
  slice2 = load_slice(filename2)

  # Pairwise align the slices
  pis = []
  pi = pst.pairwise_align(slice1, slice2)
  pis.append(pi)
    
  slice1.write(f"slice1.h5ad")
  slice2.write(f"slice2.h5ad")
  



