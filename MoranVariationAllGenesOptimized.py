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

def get_spatial_weight_matrix(stereo_spots, reconstructed_spots, dist_matrix):
  w_matrix = apply_rbf_to_tensor(dist_matrix)
  #w_sum = torch.sum(w_matrix)
  # Element-wise normalization 
  #w_matrix = w_matrix * (len(stereo_spots) / w_sum)

  return w_matrix 

def get_distance_matrix(stereo_spots, reconstructed_spots):
  dist_matrix = torch.zeros((len(stereo_spots), len(reconstructed_spots)))
  for i, spot_stereo in enumerate(stereo_spots):
    for j, spot_reconstructed in enumerate(reconstructed_spots):
      x1, y1 = spot_stereo
      x2, y2 = spot_reconstructed
      dij = ((x1 - x2)**2 + (y1-y2)**2)**0.5
      dist_matrix[i,j] = dij

  return dist_matrix 

def global_moran_score(spatial_weight_matrix, stereo_data, reconstructed_data):
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

  # User inputs...
  # 1. h5ad file for slice 1
  # 2. h5ad file for slice 2

  print(sys.argv)
  filename1 = sys.argv[1]
  filename2 = sys.argv[2]

  slice1 = align.load_slice(filename1)
  slice2 = align.load_slice(filename2)

  gene = 'magi1b'
  slice1_gene = slice1[:, gene] 
  slice2_gene = slice2[:, gene] 
  slice1_spots = slice1_gene.obsm['spatial']
  slice2_spots = slice2_gene.obsm['spatial']

  # Get the distance matrix
  d = get_distance_matrix(slice1_spots, slice2_spots)
  
  # Get the unnormalized spatial weight matrix
  w = get_spatial_weight_matrix(slice1_spots, slice2_spots, d)
      
  # Match the number of non-zero spots
  stereo_ratio_nonzero = np.count_nonzero(slice1_gene.X) / len(slice1_gene.X) * 2 
  rec_spatial = []
  rec_expression = []
  for data in slice2_gene:
    rec_spatial.append(data.obsm['spatial'][0])
    rec_expression.append(data.X.flatten()[0])

  # Choosing which spots to consider from the reconstructed data
  rec_sort_indices = np.argsort(rec_expression)[::-1][:int(stereo_ratio_nonzero * len(rec_spatial))]
  #rec_sort_indices = np.argsort(rec_expression)[::-1]
  top_expressed_spots = np.array(rec_spatial)[rec_sort_indices]
       
  # Filter entries of the w matrix to only top-expressed spots
  for i in range(w.shape[1]):
    if i not in rec_sort_indices:
      w[:,i] = 0

  # Normalize the w matrix
  print(type(w))
  w_sum = torch.sum(w)
  #w_matrix = w * (w.shape[0] / w_sum)
        
  global_moran_score = global_moran_score(w.numpy(), slice1_gene, slice2_gene)
  print(global_moran_score)

