|                                                                    | 28d40868bf54fb...  |
|:-------------------------------------------------------------------|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 12.1 ± 3.7 μs      |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 7.38 ± 1.3 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0335 ± 0.0027 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0518 ± 0.0042 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.171 ± 0.022 ms   |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0587 ± 0.0054 ms |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.117 ± 0.0076 ms  |
| Model evaluation/AR latent/forward                                 | 0.454 ± 0.39 μs    |
| Model evaluation/AR latent/rand                                    | 0.706 ± 0.68 μs    |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0582 ± 0.0077 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0542 ± 0.006 ms  |
| Model evaluation/RandomWalk latent/forward                         | 0.756 ± 0.036 μs   |
| Model evaluation/RandomWalk latent/rand                            | 0.872 ± 0.48 μs    |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0567 ± 0.0094 ms |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0524 ± 0.0068 ms |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.441 ± 0.22 s     |
| time_to_load                                                       | 3.79 ± 0.028 s     |

|                                                                    | 28d40868bf54fb...         |
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
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.69 M allocs: 0.21 GB    |
| time_to_load                                                       | 0.15 k allocs: 11.7 kB    |

