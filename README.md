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

python3 MoranVariationAllGenesOptimized

## Visualize slices and top-contributing spots
