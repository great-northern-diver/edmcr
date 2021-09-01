# edmcr
An R package for Euclidean (squared) distance matrix completion ( and determining point configurations based on the completed matrix.

## Description
Implements various general algorithms to estimate missing elements of a Euclidean (squared) distance matrix.  
Includes optimization methods based on semi-definite programming, nonparametric position, and dissimilarity parameterization formulas.
   
When the only non-missing distances are those on the minimal spanning tree, the guided random search algorithm will complete the matrix while preserving the minimal spanning tree.
   
Point configurations in specified dimensions can be determined from the completions. 

Special problems such as the sensor localization problem and reconstructing
the geometry of a molecular structure can also be solved.

Online documentation: [https://great-northern-diver.github.io/edmcr/](https://great-northern-diver.github.io/edmcr/)  

## References

- Alfakih, Khandani, and Wolkowicz (1999) "Solving Euclidean Distance Matrix Completion Problems Via Semidefinite Programming", Computational Optimization and Applications, Volume 12, pages 13–30
- Trosset (2000) "Distance Matrix Completion by Numerical Optimization", Computational Optimization and Applications, Volume 17, pages 11–22
- Krislock and Henry Wolkowicz (2010) "Explicit sensor network localization using semidefinite representations and facial reductions", SIAM Journal on Optimization, Volume 20(5), pages 2679–2708
- Fang and O'Leary (2012) "Euclidean Matrix Completion Problems", Optimization Methods and Software, Volume 27, pages 695-717, and
- Rahman and Oldford (2018) "Euclidean Distance Matrix Completion and Point Configurations from the Minimal Spanning Tree", SIAM Journal on Optimization, Volume 28, pages 528-550
