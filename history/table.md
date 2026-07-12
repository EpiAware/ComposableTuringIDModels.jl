|                                                                    | 0d15531404559d...   |
|:-------------------------------------------------------------------|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 15.3 ± 15 μs        |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 19.5 ± 0.55 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.151 ± 0.036 ms    |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0757 ± 0.0007 ms  |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.228 ± 0.021 ms    |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0852 ± 0.0012 ms  |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.187 ± 0.023 ms    |
| Model evaluation/AR latent/forward                                 | 2.15 ± 1.7 μs       |
| Model evaluation/AR latent/rand                                    | 2.88 ± 2.1 μs       |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0732 ± 0.00083 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0723 ± 0.00084 ms |
| Model evaluation/RandomWalk latent/forward                         | 0.468 ± 0.64 μs     |
| Model evaluation/RandomWalk latent/rand                            | 1.28 ± 0.81 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0791 ± 0.00092 ms |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0765 ± 0.0011 ms  |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.625 ± 0.043 s     |
| time_to_load                                                       | 4.2 ± 0.041 s       |

|                                                                    | 0d15531404559d...         |
|:-------------------------------------------------------------------|:-------------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 0.332 k allocs: 0.0522 MB |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 0.524 k allocs: 20 kB     |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 2.17 k allocs: 0.0764 MB  |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.242 k allocs: 12.3 kB   |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.68 k allocs: 0.0835 MB  |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.32 k allocs: 15.5 kB    |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 1.63 k allocs: 0.0648 MB  |
| Model evaluation/AR latent/forward                                 | 0.109 k allocs: 4.73 kB   |
| Model evaluation/AR latent/rand                                    | 0.118 k allocs: 5.77 kB   |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.35 k allocs: 15.8 kB    |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.349 k allocs: 15.1 kB   |
| Model evaluation/RandomWalk latent/forward                         | 16  allocs: 1.83 kB       |
| Model evaluation/RandomWalk latent/rand                            | 15  allocs: 2.05 kB       |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.513 k allocs: 23.5 kB   |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.509 k allocs: 22.6 kB   |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.71 M allocs: 0.212 GB   |
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   |

