|                                                                    | 4bb0533819c0a9...   |
|:-------------------------------------------------------------------|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 14.1 ± 8.8 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 7.3 ± 1.2 μs        |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0408 ± 0.0043 ms  |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0617 ± 0.00081 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.191 ± 0.011 ms    |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.069 ± 0.0019 ms   |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.141 ± 0.013 ms    |
| Model evaluation/AR latent/forward                                 | 0.465 ± 0.13 μs     |
| Model evaluation/AR latent/rand                                    | 0.687 ± 0.71 μs     |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0591 ± 0.00085 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0579 ± 0.00084 ms |
| Model evaluation/RandomWalk latent/forward                         | 0.811 ± 0.038 μs    |
| Model evaluation/RandomWalk latent/rand                            | 0.944 ± 0.51 μs     |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.062 ± 0.001 ms    |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.06 ± 0.0012 ms    |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.06 ± 2.4 s        |
| time_to_load                                                       | 4.11 ± 0.035 s      |

|                                                                    | 4bb0533819c0a9...         |
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
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.505 k allocs: 22.7 kB   |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.502 k allocs: 22 kB     |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 3.31 M allocs: 0.41 GB    |
| time_to_load                                                       | 0.15 k allocs: 11.7 kB    |

