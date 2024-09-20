*Data Preprocessing:*
The reflectance data presented here have been processed from raw DNs to reflectance using the R package `asdreadr` (https://cran.r-project.org/web/packages/asdreader/index.html) and have undergone splice correction using the package `prospectr` (https://cran.r-project.org/web/packages/prospectr/index.html), as outlined in the manuscript.
These steps were omitted from the repository to reduce the overall size of the dataset.

`ProximalInfesatation.xlsx` variable descriptions

Variable                            Description                                                  Units

1. rep                              planting group ID                                            -
2. sampleName                       ID within rep                                                -
3. plantName                        unique ID across all reps                                    -
4. treatment                        treatment condition for pot                                  -
5. uninfestedStems                  # uninfested stems                                           count
6. infestedNotCut                   # stems with 2 or more nodes burrowed through                count
7. cutStems                         # stems cut by WSS                                           count
8. totalStems                       # stems with pot                                             count
9. totalSignificantInfestations     # (infestedNotCut - neonates) + cutStems                     count
10. neonates                        # stems with dead neonates but no further WSS damage         count
11. uninfestedYield                 total threshed grain weight of uninfested stems              grams
12. infestedYield                   total threshed grain weight of infested stems                grams
13. uninfestedNoHead                # uninfested stems with no seed head                         count
14. infestedNoHead                  # WSS infested stems with no seed head                       count
15. percentSigInfested              totalSignificantInfestations / totalStems  * 100             proportion
16. grainYieldTotal                 total yield                                                  grams
17. stemYield                       average yield per stem                                       grams
