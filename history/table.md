|                                                                    | c87ef18775c631...  |
|:-------------------------------------------------------------------|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 18.4 ± 8.5 μs      |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 18.7 ± 1.4 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.126 ± 0.022 ms   |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0682 ± 0.0057 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.21 ± 0.013 ms    |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0784 ± 0.0046 ms |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.155 ± 0.016 ms   |
| Model evaluation/AR latent/forward                                 | 1.89 ± 1.3 μs      |
| Model evaluation/AR latent/rand                                    | 2.61 ± 1.3 μs      |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0658 ± 0.0029 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.065 ± 0.0029 ms  |
| Model evaluation/RandomWalk latent/forward                         | 0.408 ± 0.5 μs     |
| Model evaluation/RandomWalk latent/rand                            | 1.03 ± 0.62 μs     |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0697 ± 0.0031 ms |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0663 ± 0.0029 ms |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.393 ± 0.14 s     |
| time_to_load                                                       | 4.51 ± 0.026 s     |

|                                                                    | c87ef18775c631...         |
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
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.51 k allocs: 22.7 kB    |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.507 k allocs: 22 kB     |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.05 M allocs: 0.131 GB   |
| time_to_load                                                       | 0.15 k allocs: 11.7 kB    |

