|                                                                    | v0.1.1              | v0.1.0              | 25e392d4948329...   |
|:-------------------------------------------------------------------|:-------------------:|:-------------------:|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 12.4 ± 17 μs        | 12 ± 3.5 μs         | 26.8 ± 16 μs        |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 7.29 ± 0.99 μs      | 7.19 ± 0.96 μs      | 7.17 ± 1 μs         |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.04 ± 0.0014 ms    | 0.0405 ± 0.0083 ms  | 0.0399 ± 0.0013 ms  |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0542 ± 0.00068 ms | 0.0553 ± 0.00075 ms | 0.0541 ± 0.00067 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.165 ± 0.018 ms    | 0.169 ± 0.011 ms    | 0.172 ± 0.02 ms     |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.061 ± 0.00099 ms  | 0.0631 ± 0.0011 ms  | 0.0608 ± 0.00098 ms |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.137 ± 0.018 ms    | 0.14 ± 0.02 ms      | 0.134 ± 0.016 ms    |
| Model evaluation/AR latent/forward                                 | 0.529 ± 0.076 μs    | 0.517 ± 0.058 μs    | 0.508 ± 0.13 μs     |
| Model evaluation/AR latent/rand                                    | 0.655 ± 0.98 μs     | 0.662 ± 1 μs        | 0.674 ± 0.97 μs     |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0519 ± 0.00084 ms | 0.0528 ± 0.00079 ms | 0.0517 ± 0.00084 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0511 ± 0.00079 ms | 0.0521 ± 0.00079 ms | 0.0509 ± 0.00085 ms |
| Model evaluation/RandomWalk latent/forward                         | 1.02 ± 0.061 μs     | 1.02 ± 0.57 μs      | 0.999 ± 0.057 μs    |
| Model evaluation/RandomWalk latent/rand                            | 1.14 ± 0.74 μs      | 1.15 ± 0.76 μs      | 1.17 ± 0.73 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0559 ± 0.00097 ms | 0.0572 ± 0.00098 ms | 0.0557 ± 0.0009 ms  |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0537 ± 0.0011 ms  | 0.0548 ± 0.0011 ms  | 0.0538 ± 0.0011 ms  |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.748 ± 0.023 s     | 0.758 ± 0.028 s     | 0.75 ± 0.014 s      |
| time_to_load                                                       | 3.62 ± 0.0065 s     | 3.63 ± 0.013 s      | 3.66 ± 0.026 s      |

|                                                                    | v0.1.1                    | v0.1.0                    | 25e392d4948329...         |
|:-------------------------------------------------------------------|:-------------------------:|:-------------------------:|:-------------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 0.056 k allocs: 0.0508 MB | 0.056 k allocs: 0.0508 MB | 0.056 k allocs: 0.0508 MB |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 0.04 k allocs: 4.98 kB    | 0.04 k allocs: 4.98 kB    | 0.04 k allocs: 4.98 kB    |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.775 k allocs: 0.0319 MB | 0.775 k allocs: 0.0319 MB | 0.775 k allocs: 0.0319 MB |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.242 k allocs: 12.3 kB   | 0.242 k allocs: 12.3 kB   | 0.242 k allocs: 12.3 kB   |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.68 k allocs: 0.0835 MB  | 0.68 k allocs: 0.0835 MB  | 0.68 k allocs: 0.0835 MB  |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.32 k allocs: 15.5 kB    | 0.32 k allocs: 15.5 kB    | 0.32 k allocs: 15.5 kB    |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 1.65 k allocs: 0.0654 MB  | 1.65 k allocs: 0.0654 MB  | 1.65 k allocs: 0.0654 MB  |
| Model evaluation/AR latent/forward                                 | 20  allocs: 2.41 kB       | 20  allocs: 2.41 kB       | 20  allocs: 2.41 kB       |
| Model evaluation/AR latent/rand                                    | 22  allocs: 2.83 kB       | 22  allocs: 2.83 kB       | 22  allocs: 2.83 kB       |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.35 k allocs: 15.8 kB    | 0.35 k allocs: 15.8 kB    | 0.35 k allocs: 15.8 kB    |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.349 k allocs: 15.1 kB   | 0.349 k allocs: 15.1 kB   | 0.349 k allocs: 15.1 kB   |
| Model evaluation/RandomWalk latent/forward                         | 16  allocs: 1.83 kB       | 16  allocs: 1.83 kB       | 16  allocs: 1.83 kB       |
| Model evaluation/RandomWalk latent/rand                            | 15  allocs: 2.05 kB       | 15  allocs: 2.05 kB       | 15  allocs: 2.05 kB       |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.57 k allocs: 23.7 kB    | 0.57 k allocs: 23.7 kB    | 0.57 k allocs: 23.7 kB    |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.567 k allocs: 23 kB     | 0.567 k allocs: 23 kB     | 0.567 k allocs: 23 kB     |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 2.99 M allocs: 0.371 GB   | 2.99 M allocs: 0.371 GB   | 2.99 M allocs: 0.371 GB   |
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   | 0.149 k allocs: 11.2 kB   | 0.149 k allocs: 11.2 kB   |

