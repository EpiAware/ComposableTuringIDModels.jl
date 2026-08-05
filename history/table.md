|                                                                    | v0.1.1              | v0.1.0              | a74ad17466d071...   |
|:-------------------------------------------------------------------|:-------------------:|:-------------------:|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 11.1 ± 3.2 μs       | 30.5 ± 22 μs        | 14 ± 21 μs          |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.82 ± 2.6 μs       | 8.97 ± 2.7 μs       | 8.9 ± 2.6 μs        |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0508 ± 0.0094 ms  | 0.0508 ± 0.0023 ms  | 0.0499 ± 0.0089 ms  |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.069 ± 0.00078 ms  | 0.0702 ± 0.00082 ms | 0.0703 ± 0.00081 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.205 ± 0.013 ms    | 0.218 ± 0.027 ms    | 0.216 ± 0.027 ms    |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0785 ± 0.0012 ms  | 0.0791 ± 0.0011 ms  | 0.0793 ± 0.0013 ms  |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.173 ± 0.024 ms    | 0.176 ± 0.023 ms    | 0.175 ± 0.024 ms    |
| Model evaluation/AR latent/forward                                 | 0.674 ± 0.093 μs    | 0.628 ± 0.09 μs     | 0.637 ± 0.089 μs    |
| Model evaluation/AR latent/rand                                    | 1.86 ± 1.2 μs       | 0.956 ± 1.2 μs      | 1.81 ± 1.2 μs       |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0673 ± 0.00096 ms | 0.0671 ± 0.00099 ms | 0.0678 ± 0.00094 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0653 ± 0.00098 ms | 0.0661 ± 0.00093 ms | 0.0673 ± 0.00098 ms |
| Model evaluation/RandomWalk latent/forward                         | 1.25 ± 0.73 μs      | 1.28 ± 0.76 μs      | 1.28 ± 0.74 μs      |
| Model evaluation/RandomWalk latent/rand                            | 1.45 ± 0.93 μs      | 1.45 ± 0.94 μs      | 1.46 ± 0.95 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0716 ± 0.0011 ms  | 0.0722 ± 0.001 ms   | 0.0736 ± 0.001 ms   |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0685 ± 0.0012 ms  | 0.0694 ± 0.0012 ms  | 0.0704 ± 0.0012 ms  |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.966 ± 0.055 s     | 0.958 ± 0.037 s     | 0.953 ± 0.047 s     |
| time_to_load                                                       | 4.52 ± 0.0095 s     | 4.4 ± 0.03 s        | 4.42 ± 0.047 s      |

|                                                                    | v0.1.1                    | v0.1.0                    | a74ad17466d071...         |
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

