|                                                                    | 616713f77ba869...   |
|:-------------------------------------------------------------------|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 14.7 ± 5.1 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 18.8 ± 0.57 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.149 ± 0.035 ms    |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0745 ± 0.00068 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.222 ± 0.014 ms    |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0838 ± 0.0013 ms  |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.187 ± 0.025 ms    |
| Model evaluation/AR latent/forward                                 | 2.14 ± 1.7 μs       |
| Model evaluation/AR latent/rand                                    | 2.88 ± 2.2 μs       |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.072 ± 0.00067 ms  |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0707 ± 0.00064 ms |
| Model evaluation/RandomWalk latent/forward                         | 0.472 ± 0.62 μs     |
| Model evaluation/RandomWalk latent/rand                            | 1.08 ± 0.76 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0755 ± 0.00086 ms |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0747 ± 0.001 ms   |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.434 ± 0.19 s      |
| time_to_load                                                       | 4.16 ± 0.025 s      |

|                                                                    | 616713f77ba869...         |
|:-------------------------------------------------------------------|:-------------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 0.332 k allocs: 0.0522 MB |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 0.522 k allocs: 19.9 kB   |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 2.17 k allocs: 0.0764 MB  |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.242 k allocs: 12.3 kB   |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.68 k allocs: 0.0835 MB  |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.318 k allocs: 15.4 kB   |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 1.63 k allocs: 0.0648 MB  |
| Model evaluation/AR latent/forward                                 | 0.109 k allocs: 4.73 kB   |
| Model evaluation/AR latent/rand                                    | 0.118 k allocs: 5.77 kB   |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.35 k allocs: 15.8 kB    |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.349 k allocs: 15.1 kB   |
| Model evaluation/RandomWalk latent/forward                         | 16  allocs: 1.83 kB       |
| Model evaluation/RandomWalk latent/rand                            | 15  allocs: 2.05 kB       |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.503 k allocs: 23.4 kB   |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.499 k allocs: 22.4 kB   |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.03 M allocs: 0.128 GB   |
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   |

