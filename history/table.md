|                                                                    | 9c29d46270a6d9...   |
|:-------------------------------------------------------------------|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 12.8 ± 16 μs        |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.12 ± 2.1 μs       |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0536 ± 0.0098 ms  |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0751 ± 0.00069 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.223 ± 0.014 ms    |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0851 ± 0.0015 ms  |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.188 ± 0.024 ms    |
| Model evaluation/AR latent/forward                                 | 0.607 ± 0.068 μs    |
| Model evaluation/AR latent/rand                                    | 0.857 ± 0.98 μs     |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0717 ± 0.00068 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0709 ± 0.0007 ms  |
| Model evaluation/RandomWalk latent/forward                         | 1.08 ± 0.62 μs      |
| Model evaluation/RandomWalk latent/rand                            | 1.25 ± 0.77 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0741 ± 0.00094 ms |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0736 ± 0.0023 ms  |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.91 ± 0.38 s       |
| time_to_load                                                       | 4.16 ± 0.012 s      |

|                                                                    | 9c29d46270a6d9...         |
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
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 4.21 M allocs: 0.522 GB   |
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   |

