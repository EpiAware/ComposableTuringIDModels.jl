|                                                                    | 1541d9697d47e9...  |
|:-------------------------------------------------------------------|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 21 ± 11 μs         |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.84 ± 1.9 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0483 ± 0.0054 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0742 ± 0.001 ms  |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.226 ± 0.014 ms   |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0821 ± 0.0019 ms |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.165 ± 0.014 ms   |
| Model evaluation/AR latent/forward                                 | 0.566 ± 0.15 μs    |
| Model evaluation/AR latent/rand                                    | 0.966 ± 0.79 μs    |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0713 ± 0.0016 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0696 ± 0.0012 ms |
| Model evaluation/RandomWalk latent/forward                         | 0.93 ± 0.049 μs    |
| Model evaluation/RandomWalk latent/rand                            | 1.08 ± 0.6 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0753 ± 0.0016 ms |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0722 ± 0.0019 ms |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 12.4 s             |
| time_to_load                                                       | 4.84 ± 0.021 s     |

|                                                                    | 1541d9697d47e9...         |
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
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.575 k allocs: 23.8 kB   |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.572 k allocs: 23 kB     |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.0369 G allocs: 4.55 GB  |
| time_to_load                                                       | 0.15 k allocs: 11.7 kB    |

