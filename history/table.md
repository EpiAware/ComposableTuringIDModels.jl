|                                                                    | v0.1.2             | v0.1.1              | v0.1.0              | ee4f35d2c8e611...  |
|:-------------------------------------------------------------------|:------------------:|:-------------------:|:-------------------:|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 10.2 ± 17 μs       | 10 ± 5.6 μs         | 10.5 ± 17 μs        | 10.3 ± 17 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 7.44 ± 1.7 μs      | 7.4 ± 1.2 μs        | 7.45 ± 1.4 μs       | 7.54 ± 1.4 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0472 ± 0.0047 ms | 0.0532 ± 0.01 ms    | 0.0518 ± 0.0099 ms  | 0.0477 ± 0.0051 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 5.12 ± 0.41 μs     | 0.0743 ± 0.00075 ms | 0.0746 ± 0.00075 ms | 5.33 ± 0.42 μs     |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 17 ± 19 μs         | 0.224 ± 0.022 ms    | 0.223 ± 0.022 ms    | 16.3 ± 19 μs       |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 9.17 ± 0.39 μs     | 0.0827 ± 0.0011 ms  | 0.0828 ± 0.0012 ms  | 9.22 ± 0.41 μs     |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.125 ± 0.021 ms   | 0.192 ± 0.029 ms    | 0.191 ± 0.027 ms    | 0.126 ± 0.021 ms   |
| Model evaluation/AR latent/forward                                 | 0.648 ± 0.81 μs    | 0.604 ± 0.073 μs    | 0.59 ± 0.073 μs     | 0.616 ± 0.82 μs    |
| Model evaluation/AR latent/rand                                    | 1.68 ± 0.98 μs     | 0.843 ± 1 μs        | 0.809 ± 1 μs        | 1.72 ± 1 μs        |
| Model evaluation/DirectInfections+Poisson/forward                  | 2.01 ± 0.77 μs     | 0.0714 ± 0.00078 ms | 0.0719 ± 0.00097 ms | 2.02 ± 0.54 μs     |
| Model evaluation/DirectInfections+Poisson/rand                     | 1.74 ± 1 μs        | 0.0707 ± 0.00071 ms | 0.0705 ± 0.00073 ms | 1.74 ± 1 μs        |
| Model evaluation/RandomWalk latent/forward                         | 1.1 ± 0.58 μs      | 1.09 ± 0.6 μs       | 1.09 ± 0.62 μs      | 1.08 ± 0.58 μs     |
| Model evaluation/RandomWalk latent/rand                            | 1.28 ± 0.78 μs     | 1.27 ± 0.78 μs      | 1.28 ± 0.8 μs       | 1.24 ± 0.74 μs     |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 6.6 ± 2.3 μs       | 0.0765 ± 0.001 ms   | 0.0762 ± 0.00079 ms | 7.07 ± 2.3 μs      |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 4.55 ± 2.4 μs      | 0.0735 ± 0.0011 ms  | 0.0734 ± 0.00096 ms | 4.48 ± 2.4 μs      |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.113 ± 0.003 s    | 0.988 ± 0.038 s     | 1 ± 0.032 s         | 0.105 ± 0.0027 s   |
| time_to_load                                                       | 4.78 ± 0.062 s     | 4.75 ± 0.12 s       | 4.83 ± 0.096 s      | 4.74 ± 0.081 s     |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | ee4f35d2c8e611...         |
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

