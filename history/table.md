|                                                                    | v0.1.2             | v0.1.1              | v0.1.0              | 0d632f2f12d430...  |
|:-------------------------------------------------------------------|:------------------:|:-------------------:|:-------------------:|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 9.65 ± 17 μs       | 9.71 ± 17 μs        | 9.88 ± 17 μs        | 9.72 ± 17 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.56 ± 2.3 μs      | 8.35 ± 2.4 μs       | 8.29 ± 2.3 μs       | 8.43 ± 2.3 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0477 ± 0.0044 ms | 0.0526 ± 0.0092 ms  | 0.0518 ± 0.0035 ms  | 0.0478 ± 0.0071 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 5.16 ± 0.4 μs      | 0.0733 ± 0.00074 ms | 0.0744 ± 0.00075 ms | 5.14 ± 0.4 μs      |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 16.4 ± 21 μs       | 0.227 ± 0.024 ms    | 0.23 ± 0.022 ms     | 16.8 ± 20 μs       |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 10.1 ± 0.37 μs     | 0.0832 ± 0.0013 ms  | 0.0842 ± 0.0011 ms  | 10.1 ± 0.39 μs     |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.126 ± 0.021 ms   | 0.189 ± 0.025 ms    | 0.19 ± 0.024 ms     | 0.126 ± 0.021 ms   |
| Model evaluation/AR latent/forward                                 | 0.604 ± 0.85 μs    | 0.586 ± 0.071 μs    | 0.613 ± 0.074 μs    | 0.619 ± 0.83 μs    |
| Model evaluation/AR latent/rand                                    | 1.72 ± 1 μs        | 1.59 ± 1 μs         | 0.878 ± 1 μs        | 1.72 ± 1 μs        |
| Model evaluation/DirectInfections+Poisson/forward                  | 2.04 ± 0.74 μs     | 0.0703 ± 0.00073 ms | 0.0714 ± 0.00082 ms | 2.02 ± 0.18 μs     |
| Model evaluation/DirectInfections+Poisson/rand                     | 1.78 ± 0.99 μs     | 0.0697 ± 0.00068 ms | 0.0706 ± 0.00078 ms | 1.75 ± 1 μs        |
| Model evaluation/RandomWalk latent/forward                         | 1.12 ± 0.62 μs     | 1.11 ± 0.64 μs      | 1.1 ± 0.64 μs       | 1.11 ± 0.6 μs      |
| Model evaluation/RandomWalk latent/rand                            | 1.29 ± 0.79 μs     | 1.3 ± 0.8 μs        | 1.31 ± 0.81 μs      | 1.3 ± 0.79 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 7.08 ± 2.3 μs      | 0.0753 ± 0.00083 ms | 0.0762 ± 0.00084 ms | 7.11 ± 2.3 μs      |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 4.54 ± 2.5 μs      | 0.0721 ± 0.00094 ms | 0.0734 ± 0.00092 ms | 4.62 ± 2.5 μs      |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.105 ± 0.0041 s   | 0.973 ± 0.034 s     | 0.986 ± 0.032 s     | 0.108 ± 0.0036 s   |
| time_to_load                                                       | 4.73 ± 0.021 s     | 4.69 ± 0.047 s      | 4.67 ± 0.042 s      | 4.71 ± 0.026 s     |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | 0d632f2f12d430...         |
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

