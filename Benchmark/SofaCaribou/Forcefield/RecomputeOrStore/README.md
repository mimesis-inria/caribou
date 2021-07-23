# Hyperelastic force field benchark (storing the F and S tensors)

Here, the following components are compared:
1. The `HyperelasticForcefieldRecomputeF` component does not store `F`
on each integration points. It will therefore compute `F`, `C` and `S` on every call to 
`addForce`, `addDForce` and `addKToMatrix`.
   
2. The `HyperelasticForcefieldStoreF` component stores `F`. It will therefore only 
recompute `C` and `S` without needing to access the position vector.
   
3. The `HyperelasticForcefieldStoreFandS` component stores all quantities (`F`, `C` and `S`).

The test scene is the usual rectangular beam cantilever bending simulation
for different discretizations of the beam.

## Results
In the following tables, the mean times (in milliseconds) are shown for various stage of 
a Newton iteration. The total time is the mean time taken to complete one Newton iteration.
The LHS and RHS times are the time taken to assemble the stiffness matrix and the force 
residual vector, respectivley. Finally, The ANA is the time taken to anaylse the stiffness
matrix pattern. This is done only one time per time step. The simulation is doing 10 "steps" 
(load increments), where each step has a maximum of 10 Newton iterations.

### Tri-linear hexahedral meshes (8 Gauss nodes per element)

```
_________________________________________________________________________________________________________________________________________________________________________________________
      Mean times per            |          Total time          |             LHS              |             RHS              |             SOL              |             ANA              
      Newton iteration          |  RecomputeF StoreF StoreF&S  |  RecomputeF StoreF StoreF&S  |  RecomputeF StoreF StoreF&S  |  RecomputeF StoreF StoreF&S  |  RecomputeF StoreF StoreF&S  
_________________________________________________________________________________________________________________________________________________________________________________________
750 elements (1116 nodes)       |    29.064   28.365  27.474   |    19.954   19.739  19.231   |    1.502    1.612   1.513    |    2.932    2.907   2.748    |    14.636   11.849  11.535   
1440 elements (2009 nodes)      |    50.619   50.794  50.667   |    35.518   35.757  35.279   |    2.433    2.610   2.746    |    5.562    5.732   5.839    |    24.161   21.899  22.177  
6000 elements (7381 nodes)      |   248.232   244.009 243.494  |   157.112   153.759 150.981  |    9.819    10.186  11.114   |    47.282   46.556  47.596   |   104.668   102.159 103.271   
16660 elements (19350 nodes)    |   834.902   818.301  839.881 |   444.668   440.269 446.005  |    27.309   28.463  30.395   |   254.845   244.138 254.750  |   303.556   293.636 297.022   
```
### Linear tetrahedral meshes (1 Gauss node per element)
```
_________________________________________________________________________________________________________________________________________________________________________________________
      Mean times per            |          Total time          |             LHS              |             RHS              |             SOL              |             ANA              
      Newton iteration          |  RecomputeF StoreF StoreF&S  |  RecomputeF StoreF StoreF&S  |  RecomputeF StoreF StoreF&S  |  RecomputeF StoreF StoreF&S  |  RecomputeF StoreF StoreF&S  
_________________________________________________________________________________________________________________________________________________________________________________________
4500 elements (1116 nodes)      |    16.323   14.874  14.990   |    10.471   9.674   9.637    |    0.687    0.688   0.758    |    1.827    1.729   1.776    |    10.533   7.696   7.746
8640 elements (2009 nodes)      |    33.644   32.596  34.052   |    20.838   20.168  20.866   |    1.478    1.534   1.763    |    5.009    5.083   5.479    |    20.142   17.733  17.506      
36000 elements (7381 nodes)     |   174.428   175.135 172.887  |    98.168   97.945   97.741  |    5.590    5.807   6.074    |    43.473   44.456  42.280   |    73.364   72.157  70.333      
99960 elements (19350 nodes)    |   639.916   632.366 635.478  |   298.308   295.349 301.426  |    14.807   14.950  16.288   |   237.195   234.478 230.474   |   215.233   210.419 209.397      
```

## Conclusion
Recomputing the `F` tensor is usually a bit slower than storing it when assembling the RHS vector.
However, that means that we cannot implement a generic method to recompute the stiffness matrix
taken has argument the current position vector. We have to rely on the `addForce` method to
store the `F` tensor. In Caribou, we wanted to provide a `update_stiffness(x)` method, hence
we choose to recompute the `F` tensor, even knowing it is a bit slower.