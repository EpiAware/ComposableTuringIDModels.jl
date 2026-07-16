|                                                                    | 25949fcea99b5f...   |
|:-------------------------------------------------------------------|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 25 ± 17 μs          |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 19.5 ± 0.64 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.153 ± 0.037 ms    |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0763 ± 0.00074 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.238 ± 0.024 ms    |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0859 ± 0.0014 ms  |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.189 ± 0.024 ms    |
| Model evaluation/AR latent/forward                                 | 2.17 ± 1.7 μs       |
| Model evaluation/AR latent/rand                                    | 2.9 ± 2.1 μs        |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0736 ± 0.00077 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0723 ± 0.00068 ms |
| Model evaluation/RandomWalk latent/forward                         | 0.484 ± 0.63 μs     |
| Model evaluation/RandomWalk latent/rand                            | 1.29 ± 0.8 μs       |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.077 ± 0.00081 ms  |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0745 ± 0.001 ms   |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.222 ± 0.68 s      |
| time_to_load                                                       | 4.32 ± 0.035 s      |

|                                                                    | 25949fcea99b5f...         |
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
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.48 k allocs: 22.3 kB    |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.477 k allocs: 21.6 kB   |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.585 M allocs: 0.0739 GB |
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   |

