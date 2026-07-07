|                                                                    | b04208f6ad4af2...   |
|:-------------------------------------------------------------------|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 26 ± 16 μs          |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 19.2 ± 0.53 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.152 ± 0.038 ms    |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0763 ± 0.00073 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.238 ± 0.023 ms    |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0858 ± 0.0013 ms  |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.189 ± 0.023 ms    |
| Model evaluation/AR latent/forward                                 | 2.18 ± 1.8 μs       |
| Model evaluation/AR latent/rand                                    | 2.88 ± 2.1 μs       |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0733 ± 0.00091 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0724 ± 0.0011 ms  |
| Model evaluation/RandomWalk latent/forward                         | 0.476 ± 0.65 μs     |
| Model evaluation/RandomWalk latent/rand                            | 1.17 ± 0.78 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0765 ± 0.00099 ms |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0761 ± 0.0011 ms  |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.546 ± 0.13 s      |
| time_to_load                                                       | 4.23 ± 0.081 s      |

|                                                                    | b04208f6ad4af2...         |
|:-------------------------------------------------------------------|:-------------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 0.332 k allocs: 0.0522 MB |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 0.522 k allocs: 19.9 kB   |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 2.17 k allocs: 0.0764 MB  |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.242 k allocs: 12.3 kB   |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.677 k allocs: 0.0834 MB |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.318 k allocs: 15.4 kB   |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 1.63 k allocs: 0.0647 MB  |
| Model evaluation/AR latent/forward                                 | 0.109 k allocs: 4.73 kB   |
| Model evaluation/AR latent/rand                                    | 0.118 k allocs: 5.77 kB   |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.347 k allocs: 15.7 kB   |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.347 k allocs: 15 kB     |
| Model evaluation/RandomWalk latent/forward                         | 16  allocs: 1.83 kB       |
| Model evaluation/RandomWalk latent/rand                            | 15  allocs: 2.05 kB       |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.499 k allocs: 23.1 kB   |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.495 k allocs: 22.2 kB   |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.14 M allocs: 0.142 GB   |
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   |

