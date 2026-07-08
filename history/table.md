|                                                                    | 09262d9b0bc80a...  |
|:-------------------------------------------------------------------|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 15.8 ± 6.5 μs      |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 19.6 ± 0.88 μs     |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.153 ± 0.036 ms   |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0752 ± 0.0029 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.231 ± 0.021 ms   |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0852 ± 0.0042 ms |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.188 ± 0.025 ms   |
| Model evaluation/AR latent/forward                                 | 2.16 ± 1.8 μs      |
| Model evaluation/AR latent/rand                                    | 2.93 ± 2.1 μs      |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0723 ± 0.0029 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0712 ± 0.0026 ms |
| Model evaluation/RandomWalk latent/forward                         | 0.515 ± 0.62 μs    |
| Model evaluation/RandomWalk latent/rand                            | 1.27 ± 0.79 μs     |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0776 ± 0.0036 ms |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0754 ± 0.0031 ms |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 2.35 ± 0.05 s      |
| time_to_load                                                       | 4.55 ± 0.1 s       |

|                                                                    | 09262d9b0bc80a...         |
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
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 6.77 M allocs: 0.835 GB   |
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   |

