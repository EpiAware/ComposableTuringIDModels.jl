|                                                                    | v0.1.2             | v0.1.1             | v0.1.0             | 81cd9e3a7a25d0... |
|:-------------------------------------------------------------------|:------------------:|:------------------:|:------------------:|:-----------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 10.9 ± 3.3 μs      | 14.8 ± 8.5 μs      | 13.1 ± 7.7 μs      | 11.8 ± 4.1 μs     |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 5.88 ± 0.83 μs     | 5.51 ± 0.64 μs     | 6.39 ± 1.1 μs      | 6.62 ± 1.2 μs     |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0342 ± 0.0044 ms | 0.0321 ± 0.0041 ms | 0.0322 ± 0.0013 ms | 29.7 ± 4.7 μs     |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 3.27 ± 0.28 μs     | 0.0546 ± 0.0043 ms | 0.0547 ± 0.0062 ms | 3.29 ± 0.5 μs     |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 20.1 ± 9.7 μs      | 0.17 ± 0.024 ms    | 0.167 ± 0.022 ms   | 18.6 ± 9.2 μs     |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 7.5 ± 0.58 μs      | 0.0605 ± 0.0052 ms | 0.0619 ± 0.007 ms  | 7.35 ± 1.3 μs     |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.0768 ± 0.012 ms  | 0.116 ± 0.016 ms   | 0.12 ± 0.0074 ms   | 0.0766 ± 0.013 ms |
| Model evaluation/AR latent/forward                                 | 0.522 ± 0.59 μs    | 0.443 ± 0.59 μs    | 0.487 ± 0.56 μs    | 0.53 ± 0.57 μs    |
| Model evaluation/AR latent/rand                                    | 1.17 ± 0.67 μs     | 1.09 ± 0.67 μs     | 0.768 ± 0.68 μs    | 1.2 ± 0.71 μs     |
| Model evaluation/DirectInfections+Poisson/forward                  | 1.44 ± 0.046 μs    | 0.0576 ± 0.0071 ms | 0.0548 ± 0.0073 ms | 1.51 ± 0.16 μs    |
| Model evaluation/DirectInfections+Poisson/rand                     | 1.18 ± 0.67 μs     | 0.0579 ± 0.0057 ms | 0.0544 ± 0.0063 ms | 1.22 ± 0.68 μs    |
| Model evaluation/RandomWalk latent/forward                         | 0.783 ± 0.071 μs   | 0.778 ± 0.059 μs   | 0.783 ± 0.047 μs   | 0.791 ± 0.058 μs  |
| Model evaluation/RandomWalk latent/rand                            | 0.904 ± 0.52 μs    | 0.891 ± 0.49 μs    | 0.913 ± 0.53 μs    | 0.926 ± 0.54 μs   |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 5.24 ± 1.1 μs      | 0.0579 ± 0.0097 ms | 0.0579 ± 0.0094 ms | 5.36 ± 1.1 μs     |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 3.27 ± 1.8 μs      | 0.0562 ± 0.011 ms  | 0.0589 ± 0.0036 ms | 3.48 ± 1.9 μs     |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.118 ± 0.029 s    | 0.747 ± 0.038 s    | 0.748 ± 0.026 s    | 0.101 ± 0.011 s   |
| time_to_load                                                       | 4.86 ± 0.021 s     | 4.81 ± 0.06 s      | 5.41 ± 0.077 s     | 5.1 ± 0.054 s     |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | 81cd9e3a7a25d0...         |
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
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   | 0.149 k allocs: 11.2 kB   | 0.15 k allocs: 11.7 kB    | 0.15 k allocs: 11.7 kB    |

