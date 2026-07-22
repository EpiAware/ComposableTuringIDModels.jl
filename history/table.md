|                                                                    | c582ef676cad50...   |
|:-------------------------------------------------------------------|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 26.2 ± 17 μs        |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 7.21 ± 1 μs         |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0406 ± 0.007 ms   |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0532 ± 0.00072 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.167 ± 0.022 ms    |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0601 ± 0.001 ms   |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.137 ± 0.019 ms    |
| Model evaluation/AR latent/forward                                 | 0.534 ± 0.066 μs    |
| Model evaluation/AR latent/rand                                    | 0.687 ± 0.97 μs     |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0509 ± 0.00082 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0501 ± 0.0008 ms  |
| Model evaluation/RandomWalk latent/forward                         | 1.01 ± 0.59 μs      |
| Model evaluation/RandomWalk latent/rand                            | 1.17 ± 0.76 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0531 ± 0.00085 ms |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0522 ± 0.001 ms   |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.463 ± 0.11 s      |
| time_to_load                                                       | 3.77 ± 0.019 s      |

|                                                                    | c582ef676cad50...         |
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
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.74 M allocs: 0.217 GB   |
| time_to_load                                                       | 0.15 k allocs: 11.7 kB    |

