|                                                                    | v0.1.2            | v0.1.1              | v0.1.0              | ad0815dd1e267f... |
|:-------------------------------------------------------------------|:-----------------:|:-------------------:|:-------------------:|:-----------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 10 ± 12 μs        | 9.96 ± 1.9 μs       | 8.91 ± 1.7 μs       | 9.28 ± 2.1 μs     |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 4.47 ± 0.55 μs    | 4.45 ± 0.51 μs      | 4.32 ± 0.47 μs      | 4.38 ± 0.48 μs    |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 25.2 ± 3.1 μs     | 28.6 ± 6.8 μs       | 27.3 ± 6.4 μs       | 24.7 ± 5.6 μs     |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 2.78 ± 0.3 μs     | 0.0416 ± 0.00092 ms | 0.0421 ± 0.00087 ms | 2.54 ± 0.29 μs    |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 13.3 ± 2.3 μs     | 0.129 ± 0.0063 ms   | 0.13 ± 0.006 ms     | 13.3 ± 2.2 μs     |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 5.41 ± 0.4 μs     | 0.046 ± 0.0015 ms   | 0.0465 ± 0.0011 ms  | 5.1 ± 0.36 μs     |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.0714 ± 0.017 ms | 0.0973 ± 0.017 ms   | 0.0962 ± 0.018 ms   | 0.0688 ± 0.017 ms |
| Model evaluation/AR latent/forward                                 | 0.607 ± 0.68 μs   | 0.463 ± 0.71 μs     | 0.429 ± 0.71 μs     | 0.409 ± 0.68 μs   |
| Model evaluation/AR latent/rand                                    | 0.73 ± 0.8 μs     | 0.585 ± 0.81 μs     | 0.552 ± 0.8 μs      | 0.559 ± 0.76 μs   |
| Model evaluation/DirectInfections+Poisson/forward                  | 1.44 ± 0.7 μs     | 0.0408 ± 0.0013 ms  | 0.0409 ± 0.0016 ms  | 1.39 ± 0.7 μs     |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.696 ± 0.81 μs   | 0.0403 ± 0.0015 ms  | 0.0406 ± 0.0012 ms  | 0.551 ± 0.76 μs   |
| Model evaluation/RandomWalk latent/forward                         | 0.361 ± 0.51 μs   | 0.314 ± 0.51 μs     | 0.291 ± 0.5 μs      | 0.301 ± 0.5 μs    |
| Model evaluation/RandomWalk latent/rand                            | 0.424 ± 0.6 μs    | 0.364 ± 0.65 μs     | 0.348 ± 0.63 μs     | 0.362 ± 0.61 μs   |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 4.94 ± 1.5 μs     | 0.0445 ± 0.0021 ms  | 0.0429 ± 0.0014 ms  | 4.66 ± 1.8 μs     |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 3.31 ± 1.8 μs     | 0.0436 ± 0.0021 ms  | 0.0407 ± 0.0012 ms  | 2.96 ± 1.9 μs     |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.0931 ± 0.012 s  | 0.588 ± 0.024 s     | 0.604 ± 0.0094 s    | 0.084 ± 0.009 s   |
| time_to_load                                                       | 3.41 ± 0.086 s    | 3.28 ± 0.14 s       | 3.17 ± 0.13 s       | 2.54 ± 0.16 s     |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | ad0815dd1e267f...         |
|:-------------------------------------------------------------------|:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 0.059 k allocs: 0.0512 MB | 0.056 k allocs: 0.0508 MB | 0.056 k allocs: 0.0508 MB | 0.059 k allocs: 0.0512 MB |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 0.043 k allocs: 4.64 kB   | 0.036 k allocs: 4.42 kB   | 0.036 k allocs: 4.42 kB   | 0.043 k allocs: 4.64 kB   |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.73 k allocs: 30.4 kB    | 0.764 k allocs: 0.0315 MB | 0.764 k allocs: 0.0315 MB | 0.73 k allocs: 30.4 kB    |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.036 k allocs: 4.41 kB   | 0.244 k allocs: 12.4 kB   | 0.244 k allocs: 12.4 kB   | 0.037 k allocs: 4.44 kB   |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.068 k allocs: 0.0591 MB | 0.68 k allocs: 0.0835 MB  | 0.68 k allocs: 0.0835 MB  | 0.068 k allocs: 0.0591 MB |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.038 k allocs: 4.83 kB   | 0.346 k allocs: 15.9 kB   | 0.346 k allocs: 15.9 kB   | 0.038 k allocs: 4.83 kB   |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 1.68 k allocs: 0.0637 MB  | 1.58 k allocs: 0.062 MB   | 1.58 k allocs: 0.062 MB   | 1.68 k allocs: 0.0637 MB  |
| Model evaluation/AR latent/forward                                 | 21  allocs: 2.44 kB       | 20  allocs: 2.41 kB       | 20  allocs: 2.41 kB       | 21  allocs: 2.44 kB       |
| Model evaluation/AR latent/rand                                    | 23  allocs: 2.86 kB       | 22  allocs: 2.83 kB       | 22  allocs: 2.83 kB       | 23  allocs: 2.86 kB       |
| Model evaluation/DirectInfections+Poisson/forward                  | 23  allocs: 2.52 kB       | 0.35 k allocs: 15.8 kB    | 0.35 k allocs: 15.8 kB    | 23  allocs: 2.52 kB       |
| Model evaluation/DirectInfections+Poisson/rand                     | 20  allocs: 2.67 kB       | 0.349 k allocs: 15.1 kB   | 0.349 k allocs: 15.1 kB   | 20  allocs: 2.67 kB       |
| Model evaluation/RandomWalk latent/forward                         | 17  allocs: 1.86 kB       | 16  allocs: 1.83 kB       | 16  allocs: 1.83 kB       | 17  allocs: 1.86 kB       |
| Model evaluation/RandomWalk latent/rand                            | 16  allocs: 2.08 kB       | 15  allocs: 2.05 kB       | 15  allocs: 2.05 kB       | 16  allocs: 2.08 kB       |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.153 k allocs: 8.28 kB   | 0.57 k allocs: 23.7 kB    | 0.57 k allocs: 23.7 kB    | 0.153 k allocs: 8.28 kB   |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.148 k allocs: 8.38 kB   | 0.567 k allocs: 23 kB     | 0.567 k allocs: 23 kB     | 0.148 k allocs: 8.38 kB   |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.369 M allocs: 0.252 GB  | 2.88 M allocs: 0.354 GB   | 2.88 M allocs: 0.354 GB   | 0.369 M allocs: 0.252 GB  |
| time_to_load                                                       | 0.2 k allocs: 11.8 kB     | 0.2 k allocs: 11.8 kB     | 0.2 k allocs: 11.8 kB     | 0.2 k allocs: 11.8 kB     |

