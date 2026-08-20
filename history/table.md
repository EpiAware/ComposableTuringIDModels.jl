|                                                                    | v0.1.2            | v0.1.1              | v0.1.0              | 1a3e4853454310... |
|:-------------------------------------------------------------------|:-----------------:|:-------------------:|:-------------------:|:-----------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 8.91 ± 10 μs      | 8.78 ± 9.5 μs       | 9.14 ± 9.9 μs       | 8.87 ± 10 μs      |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 4.9 ± 0.42 μs     | 4.85 ± 0.43 μs      | 4.9 ± 0.39 μs       | 4.96 ± 0.41 μs    |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 24 ± 5 μs         | 26.4 ± 5.7 μs       | 26.4 ± 5.4 μs       | 24.2 ± 1.5 μs     |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 2.68 ± 0.28 μs    | 0.0378 ± 0.00053 ms | 0.0376 ± 0.00047 ms | 2.7 ± 0.28 μs     |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 13.1 ± 12 μs      | 0.116 ± 0.014 ms    | 0.114 ± 0.014 ms    | 13.4 ± 12 μs      |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 6.47 ± 0.46 μs    | 0.0422 ± 0.00081 ms | 0.0422 ± 0.00073 ms | 6.48 ± 0.45 μs    |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.066 ± 0.013 ms  | 0.0935 ± 0.015 ms   | 0.0935 ± 0.015 ms   | 0.0662 ± 0.015 ms |
| Model evaluation/AR latent/forward                                 | 0.383 ± 0.57 μs   | 0.363 ± 0.035 μs    | 0.367 ± 0.036 μs    | 0.395 ± 0.59 μs   |
| Model evaluation/AR latent/rand                                    | 1.07 ± 0.73 μs    | 0.474 ± 0.71 μs     | 0.513 ± 0.72 μs     | 1.02 ± 0.73 μs    |
| Model evaluation/DirectInfections+Poisson/forward                  | 1.31 ± 0.62 μs    | 0.0358 ± 0.00051 ms | 0.036 ± 0.00074 ms  | 1.32 ± 0.61 μs    |
| Model evaluation/DirectInfections+Poisson/rand                     | 1.07 ± 0.72 μs    | 0.0352 ± 0.00058 ms | 0.0352 ± 0.00053 ms | 0.504 ± 0.72 μs   |
| Model evaluation/RandomWalk latent/forward                         | 0.315 ± 0.44 μs   | 0.301 ± 0.44 μs     | 0.302 ± 0.44 μs     | 0.328 ± 0.44 μs   |
| Model evaluation/RandomWalk latent/rand                            | 0.331 ± 0.54 μs   | 0.302 ± 0.54 μs     | 0.317 ± 0.54 μs     | 0.351 ± 0.54 μs   |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 4.19 ± 1.6 μs     | 0.0386 ± 0.00068 ms | 0.0383 ± 0.00065 ms | 4.2 ± 1.6 μs      |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 2.65 ± 1.7 μs     | 0.0371 ± 0.00072 ms | 0.0369 ± 0.00082 ms | 2.6 ± 1.6 μs      |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.0771 ± 0.0032 s | 0.528 ± 0.0083 s    | 0.527 ± 0.018 s     | 0.0765 ± 0.0037 s |
| time_to_load                                                       | 3.37 ± 0.024 s    | 2.98 ± 0.018 s      | 2.98 ± 0.012 s      | 3.4 ± 0.023 s     |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | 1a3e4853454310...         |
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

