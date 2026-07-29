|                                                                    | daf499cded1e8c...  |
|:-------------------------------------------------------------------|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 11.5 ± 2.6 μs      |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 6.48 ± 0.88 μs     |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0336 ± 0.0047 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0515 ± 0.0041 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.163 ± 0.026 ms   |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0648 ± 0.0053 ms |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.115 ± 0.014 ms   |
| Model evaluation/AR latent/forward                                 | 0.435 ± 0.19 μs    |
| Model evaluation/AR latent/rand                                    | 0.678 ± 0.68 μs    |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0525 ± 0.0066 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0522 ± 0.0057 ms |
| Model evaluation/RandomWalk latent/forward                         | 0.763 ± 0.32 μs    |
| Model evaluation/RandomWalk latent/rand                            | 0.887 ± 0.48 μs    |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0536 ± 0.0085 ms |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0524 ± 0.0069 ms |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.04 ± 0.14 s      |
| time_to_load                                                       | 3.96 ± 0.069 s     |

|                                                                    | daf499cded1e8c...         |
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
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 3.57 M allocs: 0.442 GB   |
| time_to_load                                                       | 0.15 k allocs: 11.7 kB    |

