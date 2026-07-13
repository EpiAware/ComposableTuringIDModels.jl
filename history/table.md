|                                                                    | 1b0040b6107b65...   |
|:-------------------------------------------------------------------|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 14.8 ± 4.7 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 19.2 ± 0.62 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.151 ± 0.036 ms    |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.076 ± 0.00067 ms  |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.228 ± 0.011 ms    |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0852 ± 0.0013 ms  |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.188 ± 0.022 ms    |
| Model evaluation/AR latent/forward                                 | 2.13 ± 1.7 μs       |
| Model evaluation/AR latent/rand                                    | 2.89 ± 2.2 μs       |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0733 ± 0.00068 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0721 ± 0.00063 ms |
| Model evaluation/RandomWalk latent/forward                         | 0.482 ± 0.63 μs     |
| Model evaluation/RandomWalk latent/rand                            | 1.13 ± 0.77 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0775 ± 0.00089 ms |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0762 ± 0.001 ms   |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.682 ± 0.11 s      |
| time_to_load                                                       | 4.3 ± 0.02 s        |

|                                                                    | 1b0040b6107b65...         |
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
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.503 k allocs: 23.4 kB   |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.499 k allocs: 22.4 kB   |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.64 M allocs: 0.203 GB   |
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   |

