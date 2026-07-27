|                                                                    | 3965a49cb8aa8b...   |
|:-------------------------------------------------------------------|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 14.8 ± 25 μs        |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 10.2 ± 3.2 μs       |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0591 ± 0.012 ms   |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0819 ± 0.00097 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.244 ± 0.014 ms    |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0922 ± 0.0017 ms  |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.205 ± 0.034 ms    |
| Model evaluation/AR latent/forward                                 | 0.742 ± 0.14 μs     |
| Model evaluation/AR latent/rand                                    | 2.25 ± 1.3 μs       |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0783 ± 0.0011 ms  |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.077 ± 0.0011 ms   |
| Model evaluation/RandomWalk latent/forward                         | 1.45 ± 0.9 μs       |
| Model evaluation/RandomWalk latent/rand                            | 1.68 ± 1.1 μs       |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0807 ± 0.0012 ms  |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.08 ± 0.0014 ms    |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.25 ± 1.2 s        |
| time_to_load                                                       | 4.75 ± 0.023 s      |

|                                                                    | 3965a49cb8aa8b...         |
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
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 3 M allocs: 0.372 GB      |
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   |

