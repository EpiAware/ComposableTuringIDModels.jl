|                                                                    | 83596661d1e0b0...  |
|:-------------------------------------------------------------------|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 30.7 ± 19 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 20.4 ± 0.7 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.149 ± 0.042 ms   |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0688 ± 0.0011 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.215 ± 0.028 ms   |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.078 ± 0.0013 ms  |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.172 ± 0.021 ms   |
| Model evaluation/AR latent/forward                                 | 2.55 ± 1.9 μs      |
| Model evaluation/AR latent/rand                                    | 3.33 ± 2.4 μs      |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0653 ± 0.0012 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0645 ± 0.0013 ms |
| Model evaluation/RandomWalk latent/forward                         | 0.532 ± 0.77 μs    |
| Model evaluation/RandomWalk latent/rand                            | 1.46 ± 0.95 μs     |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0689 ± 0.0012 ms |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0672 ± 0.0014 ms |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.11 ± 0.17 s      |
| time_to_load                                                       | 4.75 ± 0.19 s      |

|                                                                    | 83596661d1e0b0...         |
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
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.48 k allocs: 22.3 kB    |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.477 k allocs: 21.6 kB   |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 3.29 M allocs: 0.407 GB   |
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   |

