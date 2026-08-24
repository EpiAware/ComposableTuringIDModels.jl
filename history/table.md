|                                                                    | v0.1.2             | v0.1.1              | v0.1.0              | 3b6590f2c16e25... |
|:-------------------------------------------------------------------|:------------------:|:-------------------:|:-------------------:|:-----------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 12.4 ± 20 μs       | 29.9 ± 21 μs        | 12.3 ± 21 μs        | 11.8 ± 20 μs      |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 9.27 ± 1.8 μs      | 8.92 ± 2.6 μs       | 8.83 ± 2.7 μs       | 9.07 ± 2.7 μs     |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0472 ± 0.0083 ms | 0.0509 ± 0.0087 ms  | 0.0507 ± 0.0095 ms  | 0.047 ± 0.0087 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 5.89 ± 0.49 μs     | 0.0705 ± 0.00086 ms | 0.0696 ± 0.00083 ms | 6.08 ± 0.5 μs     |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 19.8 ± 25 μs       | 0.221 ± 0.028 ms    | 0.209 ± 0.017 ms    | 19.1 ± 26 μs      |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 10.4 ± 0.46 μs     | 0.0797 ± 0.0012 ms  | 0.0788 ± 0.0012 ms  | 10.3 ± 0.46 μs    |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.122 ± 0.023 ms   | 0.175 ± 0.024 ms    | 0.175 ± 0.025 ms    | 0.123 ± 0.023 ms  |
| Model evaluation/AR latent/forward                                 | 0.663 ± 1 μs       | 0.635 ± 0.1 μs      | 0.656 ± 0.098 μs    | 0.658 ± 1 μs      |
| Model evaluation/AR latent/rand                                    | 1.99 ± 1.2 μs      | 0.955 ± 1.2 μs      | 1.91 ± 1.2 μs       | 1.97 ± 1.2 μs     |
| Model evaluation/DirectInfections+Poisson/forward                  | 2.23 ± 0.92 μs     | 0.0677 ± 0.00092 ms | 0.0668 ± 0.00084 ms | 2.25 ± 0.83 μs    |
| Model evaluation/DirectInfections+Poisson/rand                     | 1.99 ± 1.2 μs      | 0.0671 ± 0.00096 ms | 0.0656 ± 0.00078 ms | 1.97 ± 1.2 μs     |
| Model evaluation/RandomWalk latent/forward                         | 1.27 ± 0.71 μs     | 1.26 ± 0.073 μs     | 1.27 ± 0.73 μs      | 1.3 ± 0.11 μs     |
| Model evaluation/RandomWalk latent/rand                            | 1.45 ± 0.95 μs     | 1.44 ± 0.92 μs      | 1.45 ± 0.94 μs      | 1.46 ± 0.95 μs    |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 7.8 ± 2.9 μs       | 0.0732 ± 0.00098 ms | 0.0721 ± 0.00096 ms | 7.84 ± 2.9 μs     |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 5.31 ± 3.1 μs      | 0.07 ± 0.0011 ms    | 0.0689 ± 0.0011 ms  | 5.33 ± 3.1 μs     |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.131 ± 0.007 s    | 0.968 ± 0.047 s     | 0.96 ± 0.059 s      | 0.134 ± 0.0067 s  |
| time_to_load                                                       | 4.86 ± 0.075 s     | 4.83 ± 0.047 s      | 4.83 ± 0.16 s       | 4.93 ± 0.046 s    |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | 3b6590f2c16e25...         |
|:-------------------------------------------------------------------|:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 0.059 k allocs: 0.0512 MB | 0.056 k allocs: 0.0508 MB | 0.056 k allocs: 0.0508 MB | 0.059 k allocs: 0.0512 MB |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 0.047 k allocs: 5.2 kB    | 0.04 k allocs: 4.98 kB    | 0.04 k allocs: 4.98 kB    | 0.047 k allocs: 5.2 kB    |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.738 k allocs: 30.6 kB   | 0.775 k allocs: 0.0319 MB | 0.775 k allocs: 0.0319 MB | 0.738 k allocs: 30.6 kB   |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.034 k allocs: 4.36 kB   | 0.242 k allocs: 12.3 kB   | 0.242 k allocs: 12.3 kB   | 0.034 k allocs: 4.36 kB   |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.068 k allocs: 0.0591 MB | 0.68 k allocs: 0.0835 MB  | 0.68 k allocs: 0.0835 MB  | 0.068 k allocs: 0.0591 MB |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.042 k allocs: 5.39 kB   | 0.32 k allocs: 15.5 kB    | 0.32 k allocs: 15.5 kB    | 0.042 k allocs: 5.39 kB   |
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

