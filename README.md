# Bivariate Moran Score for Cross-Platform Tissue Slice Comparison

What's the plan?

Tool should be able to be used as an imported module.
But also as a command line tool.


What is the input format?
Probably AnnData
An AnnData object that has a 'spatial' attribute specifying  X,Y coordinates
The slices should have already been aligned by PASTE.

## Generate a Moran Score for PASTE-aligned slices

invoke the following command:

```
python3 MoranScore.py <h5ad_filename1> <h5ad_filename2> <gene_name>
```

## Visualize slices and top-contributing spots

```
python3 VisAlignedSlices.py <gene_name> <reconstructed_slice> <stereoseq_slice> <reconstruction_mode> 
```

For the slices of the gene 'gsc' that we previously generated a MoranScore for, you should get something like this.

![Image of the reconstructed CTFacTomo slice for the gene 'gsc'](figures/allspots_CTFacTomo_slice31_gsc_aligned_with_Stereoseq_2.png)

