|                                                                    | v0.1.2             | v0.1.1              | v0.1.0              | 707e3309db5acb...  |
|:-------------------------------------------------------------------|:------------------:|:-------------------:|:-------------------:|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 12.8 ± 1.6 μs      | 13.3 ± 8.8 μs       | 13.2 ± 3.3 μs       | 12.5 ± 1.9 μs      |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 6.06 ± 0.82 μs     | 5.96 ± 0.84 μs      | 5.96 ± 0.81 μs      | 6.01 ± 0.82 μs     |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0346 ± 0.0073 ms | 0.0376 ± 0.0075 ms  | 0.0372 ± 0.0078 ms  | 0.0346 ± 0.0069 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 3.56 ± 0.3 μs      | 0.0538 ± 0.00096 ms | 0.0548 ± 0.00092 ms | 3.55 ± 0.3 μs      |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 17.5 ± 1.7 μs      | 0.171 ± 0.0098 ms   | 0.171 ± 0.0094 ms   | 17.5 ± 1.6 μs      |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 7.7 ± 0.79 μs      | 0.061 ± 0.0019 ms   | 0.0616 ± 0.0013 ms  | 7.4 ± 0.68 μs      |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.0926 ± 0.02 ms   | 0.129 ± 0.02 ms     | 0.13 ± 0.021 ms     | 0.0921 ± 0.019 ms  |
| Model evaluation/AR latent/forward                                 | 0.582 ± 0.88 μs    | 0.703 ± 0.86 μs     | 1.29 ± 0.86 μs      | 0.65 ± 0.87 μs     |
| Model evaluation/AR latent/rand                                    | 0.706 ± 1 μs       | 0.698 ± 0.97 μs     | 0.693 ± 0.98 μs     | 0.72 ± 0.99 μs     |
| Model evaluation/DirectInfections+Poisson/forward                  | 1.86 ± 0.17 μs     | 0.0533 ± 0.00094 ms | 0.0538 ± 0.0023 ms  | 1.83 ± 0.14 μs     |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.76 ± 1 μs        | 0.0525 ± 0.00088 ms | 0.053 ± 0.0015 ms   | 0.798 ± 0.99 μs    |
| Model evaluation/RandomWalk latent/forward                         | 1.02 ± 0.68 μs     | 1.02 ± 0.67 μs      | 1.02 ± 0.67 μs      | 1.01 ± 0.66 μs     |
| Model evaluation/RandomWalk latent/rand                            | 1.12 ± 0.78 μs     | 1.18 ± 0.79 μs      | 1.16 ± 0.78 μs      | 1.13 ± 0.78 μs     |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 6.4 ± 2.1 μs       | 0.0572 ± 0.0013 ms  | 0.0579 ± 0.0016 ms  | 6.37 ± 2.2 μs      |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 3.81 ± 2.4 μs      | 0.0551 ± 0.0011 ms  | 0.0558 ± 0.0016 ms  | 3.7 ± 2.6 μs       |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.113 ± 0.0071 s   | 0.815 ± 0.03 s      | 0.803 ± 0.026 s     | 0.122 ± 0.013 s    |
| time_to_load                                                       | 3.98 ± 0.036 s     | 3.97 ± 0.21 s       | 4.25 ± 0.15 s       | 3.41 ± 0.049 s     |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | 707e3309db5acb...         |
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

