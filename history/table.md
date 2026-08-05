|                                                                    | v0.1.1              | v0.1.0              | 4e9c175a248d1d...   |
|:-------------------------------------------------------------------|:-------------------:|:-------------------:|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 13.3 ± 3.5 μs       | 12.7 ± 3.9 μs       | 12.9 ± 3.4 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 7.17 ± 0.83 μs      | 7.33 ± 1.1 μs       | 7.34 ± 1 μs         |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0409 ± 0.0051 ms  | 0.0407 ± 0.0046 ms  | 0.0408 ± 0.0059 ms  |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0626 ± 0.00078 ms | 0.0629 ± 0.00088 ms | 0.062 ± 0.00079 ms  |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.192 ± 0.0091 ms   | 0.192 ± 0.0098 ms   | 0.194 ± 0.0096 ms   |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0694 ± 0.0013 ms  | 0.0702 ± 0.0013 ms  | 0.0695 ± 0.0014 ms  |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.141 ± 0.014 ms    | 0.142 ± 0.014 ms    | 0.141 ± 0.014 ms    |
| Model evaluation/AR latent/forward                                 | 0.485 ± 0.14 μs     | 0.469 ± 0.12 μs     | 0.498 ± 0.14 μs     |
| Model evaluation/AR latent/rand                                    | 0.756 ± 0.72 μs     | 0.698 ± 0.72 μs     | 0.761 ± 0.71 μs     |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0596 ± 0.00084 ms | 0.0598 ± 0.00099 ms | 0.0596 ± 0.001 ms   |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0585 ± 0.00085 ms | 0.0585 ± 0.0008 ms  | 0.0588 ± 0.00088 ms |
| Model evaluation/RandomWalk latent/forward                         | 0.82 ± 0.032 μs     | 0.832 ± 0.037 μs    | 0.823 ± 0.03 μs     |
| Model evaluation/RandomWalk latent/rand                            | 0.972 ± 0.53 μs     | 0.962 ± 0.56 μs     | 0.955 ± 0.53 μs     |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0637 ± 0.00099 ms | 0.0638 ± 0.0011 ms  | 0.0638 ± 0.001 ms   |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0611 ± 0.0011 ms  | 0.0615 ± 0.0012 ms  | 0.0615 ± 0.0011 ms  |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.843 ± 0.023 s     | 0.847 ± 0.019 s     | 0.84 ± 0.022 s      |
| time_to_load                                                       | 4.99 ± 0.025 s      | 5 ± 0.016 s         | 5 ± 0.032 s         |

|                                                                    | v0.1.1                    | v0.1.0                    | 4e9c175a248d1d...         |
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
| time_to_load                                                       | 0.15 k allocs: 11.7 kB    | 0.15 k allocs: 11.7 kB    | 0.15 k allocs: 11.7 kB    |

