|                                                                    | 7c61e740279838...   |
|:-------------------------------------------------------------------|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 14.2 ± 20 μs        |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.89 ± 2.7 μs       |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0513 ± 0.0092 ms  |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0703 ± 0.00085 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.211 ± 0.015 ms    |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0793 ± 0.0012 ms  |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.175 ± 0.022 ms    |
| Model evaluation/AR latent/forward                                 | 0.642 ± 0.091 μs    |
| Model evaluation/AR latent/rand                                    | 1.78 ± 1.2 μs       |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0672 ± 0.00098 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0665 ± 0.00096 ms |
| Model evaluation/RandomWalk latent/forward                         | 1.25 ± 0.74 μs      |
| Model evaluation/RandomWalk latent/rand                            | 1.45 ± 0.93 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0696 ± 0.0011 ms  |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0688 ± 0.0014 ms  |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.405 ± 0.21 s      |
| time_to_load                                                       | 4.28 ± 0.034 s      |

|                                                                    | 7c61e740279838...         |
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
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.05 M allocs: 0.132 GB   |
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   |

