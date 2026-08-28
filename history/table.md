|                                                                    | v0.1.2             | v0.1.1              | v0.1.0              | 92b4c41fdc740f...  |
|:-------------------------------------------------------------------|:------------------:|:-------------------:|:-------------------:|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 10.1 ± 17 μs       | 10.2 ± 17 μs        | 9.81 ± 17 μs        | 10.1 ± 17 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 7.58 ± 1.8 μs      | 7.34 ± 1.5 μs       | 7.3 ± 1.9 μs        | 7.34 ± 1.4 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0473 ± 0.0018 ms | 0.0518 ± 0.0085 ms  | 0.0517 ± 0.0057 ms  | 0.0471 ± 0.0079 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 5.23 ± 0.41 μs     | 0.0742 ± 0.00066 ms | 0.0728 ± 0.00066 ms | 5.08 ± 0.4 μs      |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 16.8 ± 21 μs       | 0.23 ± 0.021 ms     | 0.226 ± 0.022 ms    | 16.8 ± 19 μs       |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 10.9 ± 0.42 μs     | 0.0823 ± 0.00094 ms | 0.0806 ± 0.00086 ms | 10.2 ± 1.7 μs      |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.125 ± 0.021 ms   | 0.189 ± 0.025 ms    | 0.188 ± 0.024 ms    | 0.125 ± 0.021 ms   |
| Model evaluation/AR latent/forward                                 | 0.64 ± 0.83 μs     | 0.595 ± 0.066 μs    | 0.618 ± 0.063 μs    | 0.616 ± 0.83 μs    |
| Model evaluation/AR latent/rand                                    | 1.72 ± 0.98 μs     | 0.814 ± 0.97 μs     | 1.58 ± 0.99 μs      | 1.68 ± 0.98 μs     |
| Model evaluation/DirectInfections+Poisson/forward                  | 2.04 ± 0.79 μs     | 0.0712 ± 0.00072 ms | 0.0702 ± 0.00063 ms | 2 ± 0.75 μs        |
| Model evaluation/DirectInfections+Poisson/rand                     | 1.78 ± 0.15 μs     | 0.0703 ± 0.00073 ms | 0.0691 ± 0.00066 ms | 1.73 ± 1 μs        |
| Model evaluation/RandomWalk latent/forward                         | 1.12 ± 0.59 μs     | 1.09 ± 0.62 μs      | 1.11 ± 0.64 μs      | 1.1 ± 0.6 μs       |
| Model evaluation/RandomWalk latent/rand                            | 1.3 ± 0.77 μs      | 1.28 ± 0.79 μs      | 1.28 ± 0.8 μs       | 1.28 ± 0.78 μs     |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 6.94 ± 2.3 μs      | 0.076 ± 0.00082 ms  | 0.0747 ± 0.00075 ms | 7.15 ± 2.3 μs      |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 4.52 ± 2.6 μs      | 0.073 ± 0.00092 ms  | 0.0717 ± 0.00083 ms | 4.58 ± 2.4 μs      |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.107 ± 0.011 s    | 0.985 ± 0.034 s     | 0.968 ± 0.035 s     | 0.102 ± 0.0046 s   |
| time_to_load                                                       | 5.12 ± 0.17 s      | 5.01 ± 0.36 s       | 4.94 ± 0.014 s      | 4.99 ± 0.045 s     |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | 92b4c41fdc740f...         |
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

