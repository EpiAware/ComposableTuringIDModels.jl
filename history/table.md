|                                                                    | 6910bff2dbe061...   |
|:-------------------------------------------------------------------|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 12.6 ± 16 μs        |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.36 ± 2.3 μs       |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0534 ± 0.0095 ms  |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.076 ± 0.00071 ms  |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.228 ± 0.013 ms    |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0857 ± 0.0013 ms  |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.191 ± 0.025 ms    |
| Model evaluation/AR latent/forward                                 | 0.625 ± 0.073 μs    |
| Model evaluation/AR latent/rand                                    | 1.45 ± 1 μs         |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0736 ± 0.00078 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0724 ± 0.00071 ms |
| Model evaluation/RandomWalk latent/forward                         | 1.09 ± 0.62 μs      |
| Model evaluation/RandomWalk latent/rand                            | 1.26 ± 0.77 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0753 ± 0.00081 ms |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0745 ± 0.0009 ms  |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.32 ± 0.073 s      |
| time_to_load                                                       | 4.25 ± 0.044 s      |

|                                                                    | 6910bff2dbe061...         |
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
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.913 M allocs: 0.114 GB  |
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   |

