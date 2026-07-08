|                                                                    | 402702821a96ce...  |
|:-------------------------------------------------------------------|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 31.5 ± 21 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 20.5 ± 0.77 μs     |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.149 ± 0.039 ms   |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.071 ± 0.0071 ms  |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.219 ± 0.031 ms   |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0787 ± 0.0015 ms |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.179 ± 0.025 ms   |
| Model evaluation/AR latent/forward                                 | 2.53 ± 2.2 μs      |
| Model evaluation/AR latent/rand                                    | 3.37 ± 2.5 μs      |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0668 ± 0.0011 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0654 ± 0.001 ms  |
| Model evaluation/RandomWalk latent/forward                         | 0.549 ± 0.78 μs    |
| Model evaluation/RandomWalk latent/rand                            | 1.44 ± 0.93 μs     |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0714 ± 0.0012 ms |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.071 ± 0.0014 ms  |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.328 ± 0.36 s     |
| time_to_load                                                       | 4.56 ± 0.054 s     |

|                                                                    | 402702821a96ce...         |
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
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.627 M allocs: 0.0785 GB |
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   |

