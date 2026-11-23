|                                                                    | v0.1.2             | v0.1.1              | v0.1.0              | b0d9cee71398a8...  |
|:-------------------------------------------------------------------|:------------------:|:-------------------:|:-------------------:|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 10 ± 18 μs         | 10.2 ± 0.82 μs      | 9.96 ± 0.82 μs      | 9.71 ± 1.9 μs      |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 7.08 ± 1.8 μs      | 6.94 ± 1.8 μs       | 7.05 ± 0.94 μs      | 7.16 ± 1.9 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0446 ± 0.0086 ms | 0.0493 ± 0.01 ms    | 0.0482 ± 0.01 ms    | 0.0446 ± 0.0079 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 3.79 ± 0.41 μs     | 0.0743 ± 0.00075 ms | 0.0746 ± 0.00073 ms | 3.9 ± 0.42 μs      |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 15.8 ± 0.91 μs     | 0.229 ± 0.011 ms    | 0.228 ± 0.011 ms    | 0.0374 ± 0.024 ms  |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 9.17 ± 1.8 μs      | 0.0831 ± 0.001 ms   | 0.0831 ± 0.00097 ms | 9.06 ± 1.9 μs      |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.122 ± 0.024 ms   | 0.175 ± 0.028 ms    | 0.173 ± 0.027 ms    | 0.12 ± 0.024 ms    |
| Model evaluation/AR latent/forward                                 | 0.714 ± 0.97 μs    | 0.628 ± 0.92 μs     | 0.723 ± 0.94 μs     | 0.651 ± 0.93 μs    |
| Model evaluation/AR latent/rand                                    | 1.77 ± 1 μs        | 0.854 ± 1.1 μs      | 0.842 ± 1 μs        | 1.75 ± 1.1 μs      |
| Model evaluation/DirectInfections+Poisson/forward                  | 2.06 ± 0.56 μs     | 0.0732 ± 0.00074 ms | 0.0737 ± 0.0012 ms  | 2.04 ± 0.76 μs     |
| Model evaluation/DirectInfections+Poisson/rand                     | 1.76 ± 1.1 μs      | 0.0724 ± 0.00067 ms | 0.0731 ± 0.001 ms   | 1.73 ± 1.1 μs      |
| Model evaluation/RandomWalk latent/forward                         | 1.17 ± 0.75 μs     | 1.13 ± 0.72 μs      | 1.16 ± 0.75 μs      | 1.13 ± 0.72 μs     |
| Model evaluation/RandomWalk latent/rand                            | 1.25 ± 0.86 μs     | 1.28 ± 0.85 μs      | 1.27 ± 0.85 μs      | 1.25 ± 0.85 μs     |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 7.42 ± 2.4 μs      | 0.0786 ± 0.0017 ms  | 0.0787 ± 0.0016 ms  | 7.12 ± 2.5 μs      |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 4.81 ± 2.6 μs      | 0.0759 ± 0.0013 ms  | 0.0761 ± 0.0014 ms  | 4.82 ± 2.6 μs      |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.102 ± 0.0078 s   | 1.03 ± 0.032 s      | 1.01 ± 0.034 s      | 0.106 ± 0.0089 s   |
| time_to_load                                                       | 4.76 ± 0.077 s     | 4.84 ± 0.089 s      | 4.77 ± 0.074 s      | 3.81 ± 0.035 s     |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | b0d9cee71398a8...         |
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
| time_to_load                                                       | 0.2 k allocs: 11.8 kB     | 0.204 k allocs: 11.9 kB   | 0.2 k allocs: 11.8 kB     | 0.2 k allocs: 11.8 kB     |

