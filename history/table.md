|                                                                    | 01ef7b3bea964d...  |
|:-------------------------------------------------------------------|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 17.9 ± 9.9 μs      |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.06 ± 1.4 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0477 ± 0.0062 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0723 ± 0.0039 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.219 ± 0.017 ms   |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.078 ± 0.0033 ms  |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.16 ± 0.012 ms    |
| Model evaluation/AR latent/forward                                 | 0.535 ± 0.15 μs    |
| Model evaluation/AR latent/rand                                    | 0.756 ± 0.8 μs     |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0676 ± 0.0028 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0682 ± 0.0035 ms |
| Model evaluation/RandomWalk latent/forward                         | 0.924 ± 0.11 μs    |
| Model evaluation/RandomWalk latent/rand                            | 1.05 ± 0.61 μs     |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.071 ± 0.0031 ms  |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.07 ± 0.0031 ms   |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.618 ± 0.9 s      |
| time_to_load                                                       | 4.5 ± 0.058 s      |

|                                                                    | 01ef7b3bea964d...         |
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
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.85 M allocs: 0.23 GB    |
| time_to_load                                                       | 0.15 k allocs: 11.7 kB    |

