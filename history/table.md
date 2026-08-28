|                                                                    | v0.1.2             | v0.1.1              | v0.1.0              | 5884e17dc51af6...  |
|:-------------------------------------------------------------------|:------------------:|:-------------------:|:-------------------:|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 12.7 ± 21 μs       | 12 ± 20 μs          | 12.4 ± 22 μs        | 19.8 ± 21 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.16 ± 2.1 μs      | 7.96 ± 1.8 μs       | 7.93 ± 2.2 μs       | 8.05 ± 2.1 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0463 ± 0.0013 ms | 0.051 ± 0.01 ms     | 0.0518 ± 0.0092 ms  | 0.0467 ± 0.0012 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 5.93 ± 0.48 μs     | 0.071 ± 0.00086 ms  | 0.0706 ± 0.00085 ms | 5.93 ± 0.5 μs      |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 21 ± 23 μs         | 0.22 ± 0.029 ms     | 0.212 ± 0.022 ms    | 20.9 ± 24 μs       |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 9.45 ± 0.47 μs     | 0.0786 ± 0.001 ms   | 0.079 ± 0.00098 ms  | 9.57 ± 0.48 μs     |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.122 ± 0.023 ms   | 0.177 ± 0.024 ms    | 0.173 ± 0.023 ms    | 0.122 ± 0.023 ms   |
| Model evaluation/AR latent/forward                                 | 0.674 ± 1 μs       | 0.661 ± 0.11 μs     | 0.676 ± 0.13 μs     | 0.678 ± 1 μs       |
| Model evaluation/AR latent/rand                                    | 1.97 ± 1.2 μs      | 1.89 ± 1.3 μs       | 1.94 ± 1.2 μs       | 2 ± 1.2 μs         |
| Model evaluation/DirectInfections+Poisson/forward                  | 2.29 ± 0.94 μs     | 0.0676 ± 0.00097 ms | 0.0682 ± 0.00099 ms | 2.27 ± 0.99 μs     |
| Model evaluation/DirectInfections+Poisson/rand                     | 2.02 ± 0.13 μs     | 0.0668 ± 0.00095 ms | 0.0673 ± 0.001 ms   | 2.03 ± 0.16 μs     |
| Model evaluation/RandomWalk latent/forward                         | 1.31 ± 0.086 μs    | 1.3 ± 0.079 μs      | 1.26 ± 0.52 μs      | 1.31 ± 0.077 μs    |
| Model evaluation/RandomWalk latent/rand                            | 1.49 ± 0.95 μs     | 1.47 ± 0.96 μs      | 1.49 ± 0.94 μs      | 1.49 ± 0.95 μs     |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 7.05 ± 2.8 μs      | 0.0727 ± 0.001 ms   | 0.0733 ± 0.0011 ms  | 6.92 ± 2.8 μs      |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 5.04 ± 2.8 μs      | 0.07 ± 0.0012 ms    | 0.0707 ± 0.0013 ms  | 4.9 ± 2.9 μs       |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.142 ± 0.0092 s   | 0.97 ± 0.034 s      | 0.966 ± 0.032 s     | 0.144 ± 0.013 s    |
| time_to_load                                                       | 5.53 ± 0.036 s     | 5.46 ± 0.032 s      | 5.02 ± 0.069 s      | 5.52 ± 0.055 s     |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | 5884e17dc51af6...         |
|:-------------------------------------------------------------------|:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 0.059 k allocs: 0.0512 MB | 0.056 k allocs: 0.0508 MB | 0.056 k allocs: 0.0508 MB | 0.059 k allocs: 0.0512 MB |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 0.043 k allocs: 4.64 kB   | 0.036 k allocs: 4.42 kB   | 0.036 k allocs: 4.42 kB   | 0.043 k allocs: 4.64 kB   |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.738 k allocs: 30.6 kB   | 0.775 k allocs: 0.0319 MB | 0.775 k allocs: 0.0319 MB | 0.738 k allocs: 30.6 kB   |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.034 k allocs: 4.36 kB   | 0.242 k allocs: 12.3 kB   | 0.242 k allocs: 12.3 kB   | 0.034 k allocs: 4.36 kB   |
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

