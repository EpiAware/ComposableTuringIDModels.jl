|                                                                    | v0.1.0              | 70bd186f7e5436...   | v0.1.0 / 70bd186f7e5436... |
|:-------------------------------------------------------------------|:-------------------:|:-------------------:|:--------------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 13 ± 21 μs          | 12.5 ± 21 μs        | 1.04 ± 2.4                 |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.98 ± 2.6 μs       | 8.92 ± 2.7 μs       | 1.01 ± 0.42                |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0513 ± 0.0088 ms  | 0.0511 ± 0.0091 ms  | 1 ± 0.25                   |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0708 ± 0.00096 ms | 0.0694 ± 0.00087 ms | 1.02 ± 0.019               |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.22 ± 0.028 ms     | 0.206 ± 0.014 ms    | 1.07 ± 0.16                |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0799 ± 0.0018 ms  | 0.0781 ± 0.0014 ms  | 1.02 ± 0.029               |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.178 ± 0.025 ms    | 0.175 ± 0.025 ms    | 1.02 ± 0.21                |
| Model evaluation/AR latent/forward                                 | 0.659 ± 0.098 μs    | 0.658 ± 0.14 μs     | 1 ± 0.26                   |
| Model evaluation/AR latent/rand                                    | 1.88 ± 1.2 μs       | 1.84 ± 1.2 μs       | 1.02 ± 0.94                |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0678 ± 0.0011 ms  | 0.0661 ± 0.0011 ms  | 1.03 ± 0.024               |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0668 ± 0.0011 ms  | 0.0652 ± 0.001 ms   | 1.02 ± 0.023               |
| Model evaluation/RandomWalk latent/forward                         | 1.27 ± 0.62 μs      | 1.28 ± 0.72 μs      | 0.994 ± 0.74               |
| Model evaluation/RandomWalk latent/rand                            | 1.46 ± 0.95 μs      | 1.45 ± 0.94 μs      | 1.01 ± 0.93                |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.073 ± 0.0011 ms   | 0.0698 ± 0.0011 ms  | 1.05 ± 0.023               |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0698 ± 0.0014 ms  | 0.0677 ± 0.0013 ms  | 1.03 ± 0.029               |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.345 ± 0.08 s      | 0.891 ± 0.3 s       | 0.387 ± 0.16               |
| time_to_load                                                       | 4.44 ± 0.27 s       | 4.48 ± 0.087 s      | 0.99 ± 0.063               |

|                                                                    | v0.1.0                    | 70bd186f7e5436...         | v0.1.0 / 70bd186f7e5436... |
|:-------------------------------------------------------------------|:-------------------------:|:-------------------------:|:--------------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 0.056 k allocs: 0.0508 MB | 0.056 k allocs: 0.0508 MB | 1                          |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 0.04 k allocs: 4.98 kB    | 0.04 k allocs: 4.98 kB    | 1                          |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.775 k allocs: 0.0319 MB | 0.775 k allocs: 0.0319 MB | 1                          |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.242 k allocs: 12.3 kB   | 0.242 k allocs: 12.3 kB   | 1                          |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.68 k allocs: 0.0835 MB  | 0.68 k allocs: 0.0835 MB  | 1                          |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.32 k allocs: 15.5 kB    | 0.32 k allocs: 15.5 kB    | 1                          |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 1.65 k allocs: 0.0654 MB  | 1.65 k allocs: 0.0654 MB  | 1                          |
| Model evaluation/AR latent/forward                                 | 20  allocs: 2.41 kB       | 20  allocs: 2.41 kB       | 1                          |
| Model evaluation/AR latent/rand                                    | 22  allocs: 2.83 kB       | 22  allocs: 2.83 kB       | 1                          |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.35 k allocs: 15.8 kB    | 0.35 k allocs: 15.8 kB    | 1                          |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.349 k allocs: 15.1 kB   | 0.349 k allocs: 15.1 kB   | 1                          |
| Model evaluation/RandomWalk latent/forward                         | 16  allocs: 1.83 kB       | 16  allocs: 1.83 kB       | 1                          |
| Model evaluation/RandomWalk latent/rand                            | 15  allocs: 2.05 kB       | 15  allocs: 2.05 kB       | 1                          |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.545 k allocs: 23.3 kB   | 0.48 k allocs: 22.3 kB    | 1.05                       |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.542 k allocs: 22.6 kB   | 0.477 k allocs: 21.6 kB   | 1.05                       |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.912 M allocs: 0.114 GB  | 2.03 M allocs: 0.252 GB   | 0.453                      |
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   | 0.149 k allocs: 11.2 kB   | 1                          |

