|                                                                    | 7c6f0048e7382a...   |
|:-------------------------------------------------------------------|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 14.6 ± 5.2 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 19.2 ± 0.59 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.152 ± 0.036 ms    |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0761 ± 0.00065 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.227 ± 0.015 ms    |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0884 ± 0.0013 ms  |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.189 ± 0.026 ms    |
| Model evaluation/AR latent/forward                                 | 2.17 ± 1.8 μs       |
| Model evaluation/AR latent/rand                                    | 2.92 ± 2.3 μs       |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0732 ± 0.00077 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0722 ± 0.00072 ms |
| Model evaluation/RandomWalk latent/forward                         | 0.464 ± 0.64 μs     |
| Model evaluation/RandomWalk latent/rand                            | 0.625 ± 0.8 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0787 ± 0.00082 ms |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0761 ± 0.00095 ms |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.633 ± 0.063 s     |
| time_to_load                                                       | 4.12 ± 0.00072 s    |

|                                                                    | 7c6f0048e7382a...         |
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
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.6 M allocs: 0.199 GB    |
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   |

