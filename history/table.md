|                                                                    | v0.1.2             | v0.1.1              | v0.1.0              | d58c3c426c57ff...  |
|:-------------------------------------------------------------------|:------------------:|:-------------------:|:-------------------:|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 10.4 ± 17 μs       | 25.6 ± 18 μs        | 10.1 ± 17 μs        | 10.4 ± 18 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 7.56 ± 1.5 μs      | 7.39 ± 1.9 μs       | 7.4 ± 1.5 μs        | 7.51 ± 1.3 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0473 ± 0.0032 ms | 0.0518 ± 0.0033 ms  | 0.0523 ± 0.0085 ms  | 0.0472 ± 0.0028 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 5.15 ± 0.42 μs     | 0.073 ± 0.00068 ms  | 0.0732 ± 0.00076 ms | 5.39 ± 0.4 μs      |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 16.8 ± 20 μs       | 0.228 ± 0.024 ms    | 0.227 ± 0.022 ms    | 16.4 ± 20 μs       |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 9.18 ± 0.37 μs     | 0.0813 ± 0.00086 ms | 0.0815 ± 0.00097 ms | 10.1 ± 1.7 μs      |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.127 ± 0.022 ms   | 0.188 ± 0.024 ms    | 0.191 ± 0.025 ms    | 0.127 ± 0.022 ms   |
| Model evaluation/AR latent/forward                                 | 0.658 ± 0.84 μs    | 0.624 ± 0.1 μs      | 0.609 ± 0.071 μs    | 0.658 ± 0.86 μs    |
| Model evaluation/AR latent/rand                                    | 1.79 ± 1 μs        | 1.71 ± 1 μs         | 0.845 ± 1 μs        | 1.79 ± 1 μs        |
| Model evaluation/DirectInfections+Poisson/forward                  | 2.06 ± 0.82 μs     | 0.0706 ± 0.00077 ms | 0.0702 ± 0.00077 ms | 2.11 ± 0.24 μs     |
| Model evaluation/DirectInfections+Poisson/rand                     | 1.8 ± 0.97 μs      | 0.0695 ± 0.00076 ms | 0.0694 ± 0.00071 ms | 1.82 ± 0.97 μs     |
| Model evaluation/RandomWalk latent/forward                         | 1.14 ± 0.62 μs     | 1.14 ± 0.61 μs      | 1.12 ± 0.61 μs      | 1.16 ± 0.63 μs     |
| Model evaluation/RandomWalk latent/rand                            | 1.34 ± 0.83 μs     | 1.31 ± 0.79 μs      | 1.31 ± 0.81 μs      | 1.36 ± 0.83 μs     |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 7.16 ± 2.4 μs      | 0.0754 ± 0.00081 ms | 0.0751 ± 0.00085 ms | 7.2 ± 2.4 μs       |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 4.69 ± 2.6 μs      | 0.0724 ± 0.00096 ms | 0.0722 ± 0.001 ms   | 4.71 ± 2.6 μs      |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.126 ± 0.011 s    | 0.973 ± 0.031 s     | 0.977 ± 0.033 s     | 0.12 ± 0.0089 s    |
| time_to_load                                                       | 5.44 ± 0.18 s      | 5.04 ± 0.088 s      | 5.04 ± 0.071 s      | 5.53 ± 0.06 s      |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | d58c3c426c57ff...         |
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

