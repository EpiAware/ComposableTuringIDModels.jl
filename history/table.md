|                                                                    | 1eb527fdc5fe45...   |
|:-------------------------------------------------------------------|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 0.0338 ± 0.022 ms   |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.92 ± 2.6 μs       |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0516 ± 0.003 ms   |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0717 ± 0.00099 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.219 ± 0.03 ms     |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0792 ± 0.0014 ms  |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.175 ± 0.024 ms    |
| Model evaluation/AR latent/forward                                 | 0.677 ± 0.11 μs     |
| Model evaluation/AR latent/rand                                    | 1.92 ± 1.2 μs       |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0666 ± 0.0011 ms  |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0652 ± 0.0011 ms  |
| Model evaluation/RandomWalk latent/forward                         | 1.27 ± 0.69 μs      |
| Model evaluation/RandomWalk latent/rand                            | 1.46 ± 0.91 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0686 ± 0.0011 ms  |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0679 ± 0.0014 ms  |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.472 ± 5.7 s       |
| time_to_load                                                       | 4.77 ± 0.033 s      |

|                                                                    | 1eb527fdc5fe45...         |
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
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.1 M allocs: 0.138 GB    |
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   |

