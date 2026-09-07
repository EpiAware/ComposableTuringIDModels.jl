|                                                                    | v0.1.2             | v0.1.1              | v0.1.0              | 7ce3dc1d602e86...  |
|:-------------------------------------------------------------------|:------------------:|:-------------------:|:-------------------:|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 11.8 ± 21 μs       | 14.1 ± 21 μs        | 16.1 ± 21 μs        | 12.9 ± 21 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.07 ± 2.2 μs      | 7.91 ± 2 μs         | 7.91 ± 2.2 μs       | 7.93 ± 2 μs        |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0466 ± 0.0019 ms | 0.0505 ± 0.009 ms   | 0.0509 ± 0.0082 ms  | 0.0466 ± 0.0016 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 5.9 ± 0.49 μs      | 0.0708 ± 0.00078 ms | 0.0764 ± 0.00094 ms | 5.92 ± 0.5 μs      |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 19 ± 24 μs         | 0.214 ± 0.026 ms    | 0.217 ± 0.027 ms    | 18.6 ± 23 μs       |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 9.44 ± 0.44 μs     | 0.0786 ± 0.00097 ms | 0.0788 ± 0.001 ms   | 9.45 ± 0.5 μs      |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.122 ± 0.023 ms   | 0.175 ± 0.024 ms    | 0.175 ± 0.024 ms    | 0.123 ± 0.023 ms   |
| Model evaluation/AR latent/forward                                 | 0.665 ± 1 μs       | 0.672 ± 0.13 μs     | 0.654 ± 0.16 μs     | 0.657 ± 1 μs       |
| Model evaluation/AR latent/rand                                    | 1.94 ± 1.2 μs      | 1.87 ± 1.2 μs       | 0.873 ± 1.2 μs      | 1.95 ± 1.2 μs      |
| Model evaluation/DirectInfections+Poisson/forward                  | 2.27 ± 0.97 μs     | 0.0676 ± 0.00089 ms | 0.068 ± 0.00085 ms  | 2.29 ± 0.96 μs     |
| Model evaluation/DirectInfections+Poisson/rand                     | 1.98 ± 1.1 μs      | 0.0666 ± 0.00082 ms | 0.0666 ± 0.00082 ms | 2 ± 1.2 μs         |
| Model evaluation/RandomWalk latent/forward                         | 1.28 ± 0.7 μs      | 1.27 ± 0.71 μs      | 1.25 ± 0.7 μs       | 1.3 ± 0.64 μs      |
| Model evaluation/RandomWalk latent/rand                            | 1.44 ± 0.93 μs     | 1.46 ± 0.95 μs      | 1.45 ± 0.93 μs      | 1.46 ± 0.95 μs     |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 7.68 ± 2.9 μs      | 0.0728 ± 0.00097 ms | 0.0726 ± 0.00094 ms | 7.59 ± 2.8 μs      |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 5.21 ± 3 μs        | 0.0698 ± 0.0012 ms  | 0.0702 ± 0.0012 ms  | 5.19 ± 3.1 μs      |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.128 ± 0.014 s    | 0.962 ± 0.041 s     | 0.973 ± 0.033 s     | 0.132 ± 0.0056 s   |
| time_to_load                                                       | 4.91 ± 0.11 s      | 4.91 ± 0.051 s      | 5.12 ± 0.014 s      | 4.91 ± 0.17 s      |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | 7ce3dc1d602e86...         |
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

