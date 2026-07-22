|                                                                    | 4568cd49f3cdf8...   |
|:-------------------------------------------------------------------|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 13.6 ± 10 μs        |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 7.47 ± 1.3 μs       |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0494 ± 0.0062 ms  |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0723 ± 0.00067 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.22 ± 0.0078 ms    |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0798 ± 0.00097 ms |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.167 ± 0.014 ms    |
| Model evaluation/AR latent/forward                                 | 0.625 ± 0.12 μs     |
| Model evaluation/AR latent/rand                                    | 1.44 ± 0.72 μs      |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0707 ± 0.00083 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0695 ± 0.00069 ms |
| Model evaluation/RandomWalk latent/forward                         | 0.942 ± 0.045 μs    |
| Model evaluation/RandomWalk latent/rand                            | 1.09 ± 0.53 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0725 ± 0.00083 ms |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0715 ± 0.00092 ms |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.02 ± 0.075 s      |
| time_to_load                                                       | 4.68 ± 0.13 s       |

|                                                                    | 4568cd49f3cdf8...         |
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
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 2.58 M allocs: 0.321 GB   |
| time_to_load                                                       | 0.15 k allocs: 11.7 kB    |

