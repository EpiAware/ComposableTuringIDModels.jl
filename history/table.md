|                                                                    | v0.1.2             | v0.1.1              | v0.1.0              | 234dc4447d406a...  |
|:-------------------------------------------------------------------|:------------------:|:-------------------:|:-------------------:|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 10.9 ± 18 μs       | 10.4 ± 17 μs        | 10.4 ± 17 μs        | 10.8 ± 18 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 7.53 ± 1.7 μs      | 7.36 ± 1.6 μs       | 7.43 ± 1.4 μs       | 7.45 ± 1.3 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0478 ± 0.0062 ms | 0.0525 ± 0.0095 ms  | 0.0528 ± 0.009 ms   | 0.0471 ± 0.0083 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 5.34 ± 0.41 μs     | 0.0733 ± 0.0008 ms  | 0.0747 ± 0.00078 ms | 5.24 ± 0.42 μs     |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 16.8 ± 20 μs       | 0.226 ± 0.023 ms    | 0.23 ± 0.022 ms     | 16.6 ± 20 μs       |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 9.41 ± 0.4 μs      | 0.0819 ± 0.0011 ms  | 0.0829 ± 0.0011 ms  | 9.2 ± 0.38 μs      |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.127 ± 0.022 ms   | 0.189 ± 0.025 ms    | 0.192 ± 0.026 ms    | 0.125 ± 0.021 ms   |
| Model evaluation/AR latent/forward                                 | 0.639 ± 0.84 μs    | 0.62 ± 0.079 μs     | 0.61 ± 0.072 μs     | 0.623 ± 0.84 μs    |
| Model evaluation/AR latent/rand                                    | 1.74 ± 1 μs        | 1.29 ± 1 μs         | 0.856 ± 1 μs        | 1.76 ± 0.99 μs     |
| Model evaluation/DirectInfections+Poisson/forward                  | 2.03 ± 0.8 μs      | 0.0704 ± 0.00085 ms | 0.0713 ± 0.00072 ms | 2.05 ± 0.2 μs      |
| Model evaluation/DirectInfections+Poisson/rand                     | 1.78 ± 0.94 μs     | 0.0698 ± 0.00094 ms | 0.0705 ± 0.00068 ms | 1.83 ± 0.11 μs     |
| Model evaluation/RandomWalk latent/forward                         | 1.12 ± 0.61 μs     | 1.12 ± 0.62 μs      | 1.11 ± 0.62 μs      | 1.13 ± 0.59 μs     |
| Model evaluation/RandomWalk latent/rand                            | 1.32 ± 0.8 μs      | 1.34 ± 0.8 μs       | 1.31 ± 0.81 μs      | 1.33 ± 0.8 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 7.16 ± 2.4 μs      | 0.0756 ± 0.00091 ms | 0.0762 ± 0.00085 ms | 7.27 ± 2.4 μs      |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 4.67 ± 2.6 μs      | 0.0723 ± 0.00098 ms | 0.0732 ± 0.00093 ms | 4.68 ± 2.6 μs      |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.118 ± 0.0073 s   | 0.98 ± 0.033 s      | 0.999 ± 0.045 s     | 0.117 ± 0.0084 s   |
| time_to_load                                                       | 5.41 ± 0.095 s     | 5.19 ± 0.063 s      | 4.95 ± 0.041 s      | 5.34 ± 0.077 s     |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | 234dc4447d406a...         |
|:-------------------------------------------------------------------|:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 0.059 k allocs: 0.0512 MB | 0.056 k allocs: 0.0508 MB | 0.056 k allocs: 0.0508 MB | 0.059 k allocs: 0.0512 MB |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 0.043 k allocs: 4.64 kB   | 0.036 k allocs: 4.42 kB   | 0.036 k allocs: 4.42 kB   | 0.043 k allocs: 4.64 kB   |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.738 k allocs: 30.6 kB   | 0.775 k allocs: 0.0319 MB | 0.775 k allocs: 0.0319 MB | 0.738 k allocs: 30.6 kB   |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.034 k allocs: 4.36 kB   | 0.242 k allocs: 12.3 kB   | 0.242 k allocs: 12.3 kB   | 0.035 k allocs: 4.39 kB   |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.068 k allocs: 0.0591 MB | 0.68 k allocs: 0.0835 MB  | 0.68 k allocs: 0.0835 MB  | 0.068 k allocs: 0.0591 MB |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.038 k allocs: 4.83 kB   | 0.316 k allocs: 14.9 kB   | 0.316 k allocs: 14.9 kB   | 0.038 k allocs: 4.83 kB   |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 1.7 k allocs: 0.0643 MB   | 1.65 k allocs: 0.0654 MB  | 1.65 k allocs: 0.0654 MB  | 1.7 k allocs: 0.0643 MB   |
| Model evaluation/AR latent/forward                                 | 21  allocs: 2.44 kB       | 20  allocs: 2.41 kB       | 20  allocs: 2.41 kB       | 21  allocs: 2.44 kB       |
| Model evaluation/AR latent/rand                                    | 23  allocs: 2.86 kB       | 22  allocs: 2.83 kB       | 22  allocs: 2.83 kB       | 23  allocs: 2.86 kB       |
| Model evaluation/DirectInfections+Poisson/forward                  | 23  allocs: 2.52 kB       | 0.35 k allocs: 15.8 kB    | 0.35 k allocs: 15.8 kB    | 23  allocs: 2.52 kB       |
| Model evaluation/DirectInfections+Poisson/rand                     | 20  allocs: 2.67 kB       | 0.349 k allocs: 15.1 kB   | 0.349 k allocs: 15.1 kB   | 20  allocs: 2.67 kB       |
| Model evaluation/RandomWalk latent/forward                         | 17  allocs: 1.86 kB       | 16  allocs: 1.83 kB       | 16  allocs: 1.83 kB       | 17  allocs: 1.86 kB       |
| Model evaluation/RandomWalk latent/rand                            | 16  allocs: 2.08 kB       | 15  allocs: 2.05 kB       | 15  allocs: 2.05 kB       | 16  allocs: 2.08 kB       |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.153 k allocs: 8.28 kB   | 0.57 k allocs: 23.7 kB    | 0.57 k allocs: 23.7 kB    | 0.153 k allocs: 8.28 kB   |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.148 k allocs: 8.38 kB   | 0.567 k allocs: 23 kB     | 0.567 k allocs: 23 kB     | 0.148 k allocs: 8.38 kB   |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.479 M allocs: 0.269 GB  | 2.99 M allocs: 0.371 GB   | 2.99 M allocs: 0.371 GB   | 0.479 M allocs: 0.269 GB  |
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   | 0.149 k allocs: 11.2 kB   | 0.149 k allocs: 11.2 kB   | 0.149 k allocs: 11.2 kB   |

