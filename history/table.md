|                                                                    | 84643c0bacc202...   |
|:-------------------------------------------------------------------|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 12.6 ± 16 μs        |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.39 ± 2.3 μs       |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0537 ± 0.01 ms    |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0763 ± 0.0008 ms  |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.228 ± 0.016 ms    |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0853 ± 0.0014 ms  |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.192 ± 0.026 ms    |
| Model evaluation/AR latent/forward                                 | 0.624 ± 0.074 μs    |
| Model evaluation/AR latent/rand                                    | 1.59 ± 0.99 μs      |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0732 ± 0.00076 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0724 ± 0.0007 ms  |
| Model evaluation/RandomWalk latent/forward                         | 1.11 ± 0.65 μs      |
| Model evaluation/RandomWalk latent/rand                            | 1.3 ± 0.8 μs        |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.077 ± 0.0009 ms   |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0746 ± 0.00096 ms |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.654 ± 0.1 s       |
| time_to_load                                                       | 4.57 ± 0.019 s      |

|                                                                    | 84643c0bacc202...         |
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
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.66 M allocs: 0.207 GB   |
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   |

