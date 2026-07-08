|                                                                    | 544d7610610f0d...  |
|:-------------------------------------------------------------------|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 22.1 ± 10 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 19.2 ± 1.6 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.129 ± 0.026 ms   |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0723 ± 0.0027 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.227 ± 0.016 ms   |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0806 ± 0.0043 ms |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.163 ± 0.014 ms   |
| Model evaluation/AR latent/forward                                 | 1.91 ± 1.2 μs      |
| Model evaluation/AR latent/rand                                    | 2.6 ± 1.5 μs       |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.068 ± 0.0028 ms  |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0674 ± 0.0026 ms |
| Model evaluation/RandomWalk latent/forward                         | 0.463 ± 0.5 μs     |
| Model evaluation/RandomWalk latent/rand                            | 1.08 ± 0.62 μs     |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0723 ± 0.0043 ms |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0717 ± 0.0035 ms |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.576 ± 0.08 s     |
| time_to_load                                                       | 4.75 ± 0.064 s     |

|                                                                    | 544d7610610f0d...         |
|:-------------------------------------------------------------------|:-------------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 0.332 k allocs: 0.0522 MB |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 0.522 k allocs: 19.9 kB   |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 2.17 k allocs: 0.0764 MB  |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.242 k allocs: 12.3 kB   |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.677 k allocs: 0.0834 MB |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.318 k allocs: 15.4 kB   |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 1.63 k allocs: 0.0647 MB  |
| Model evaluation/AR latent/forward                                 | 0.109 k allocs: 4.73 kB   |
| Model evaluation/AR latent/rand                                    | 0.118 k allocs: 5.77 kB   |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.347 k allocs: 15.7 kB   |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.347 k allocs: 15 kB     |
| Model evaluation/RandomWalk latent/forward                         | 16  allocs: 1.83 kB       |
| Model evaluation/RandomWalk latent/rand                            | 15  allocs: 2.05 kB       |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.499 k allocs: 23.1 kB   |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.495 k allocs: 22.2 kB   |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.63 M allocs: 0.203 GB   |
| time_to_load                                                       | 0.15 k allocs: 11.7 kB    |

