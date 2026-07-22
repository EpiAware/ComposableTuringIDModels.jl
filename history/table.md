|                                                                    | be266685d53c5d...   |
|:-------------------------------------------------------------------|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 11.7 ± 17 μs        |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.31 ± 2.3 μs       |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0532 ± 0.0095 ms  |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0748 ± 0.00071 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.222 ± 0.013 ms    |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0841 ± 0.0013 ms  |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.19 ± 0.026 ms     |
| Model evaluation/AR latent/forward                                 | 0.601 ± 0.067 μs    |
| Model evaluation/AR latent/rand                                    | 0.93 ± 0.99 μs      |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.072 ± 0.00072 ms  |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.071 ± 0.00068 ms  |
| Model evaluation/RandomWalk latent/forward                         | 1.08 ± 0.64 μs      |
| Model evaluation/RandomWalk latent/rand                            | 1.26 ± 0.77 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0738 ± 0.00077 ms |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0733 ± 0.0019 ms  |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.599 ± 0.19 s      |
| time_to_load                                                       | 4.01 ± 0.035 s      |

|                                                                    | be266685d53c5d...         |
|:-------------------------------------------------------------------|:-------------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 0.056 k allocs: 0.0508 MB |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 0.04 k allocs: 4.98 kB    |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.775 k allocs: 0.0319 MB |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.242 k allocs: 12.3 kB   |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.68 k allocs: 0.0835 MB  |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.32 k allocs: 15.5 kB    |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 1.65 k allocs: 0.0654 MB  |
| Model evaluation/AR latent/forward                                 | 20  allocs: 2.41 kB       |
| Model evaluation/AR latent/rand                                    | 22  allocs: 2.83 kB       |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.35 k allocs: 15.8 kB    |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.349 k allocs: 15.1 kB   |
| Model evaluation/RandomWalk latent/forward                         | 16  allocs: 1.83 kB       |
| Model evaluation/RandomWalk latent/rand                            | 15  allocs: 2.05 kB       |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.48 k allocs: 22.3 kB    |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.477 k allocs: 21.6 kB   |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.31 M allocs: 0.163 GB   |
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   |

