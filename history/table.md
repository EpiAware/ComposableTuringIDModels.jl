|                                                                    | v0.1.2             | v0.1.1              | v0.1.0              | 2a78930a4e7ba2...  |
|:-------------------------------------------------------------------|:------------------:|:-------------------:|:-------------------:|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 9.81 ± 18 μs       | 10.2 ± 2.7 μs       | 9.61 ± 11 μs        | 9.92 ± 17 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 7.43 ± 1.4 μs      | 7.38 ± 2 μs         | 7.31 ± 1.8 μs       | 7.5 ± 1.7 μs       |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0469 ± 0.0077 ms | 0.0518 ± 0.0097 ms  | 0.0517 ± 0.009 ms   | 0.0482 ± 0.0079 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 5.29 ± 0.41 μs     | 0.0745 ± 0.00075 ms | 0.0733 ± 0.00068 ms | 5.16 ± 0.41 μs     |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 17.3 ± 21 μs       | 0.221 ± 0.021 ms    | 0.219 ± 0.014 ms    | 16.2 ± 21 μs       |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 9.16 ± 0.27 μs     | 0.083 ± 0.001 ms    | 0.0812 ± 0.00092 ms | 10.1 ± 1.7 μs      |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.126 ± 0.022 ms   | 0.19 ± 0.027 ms     | 0.188 ± 0.028 ms    | 0.126 ± 0.022 ms   |
| Model evaluation/AR latent/forward                                 | 0.615 ± 0.85 μs    | 0.601 ± 0.068 μs    | 0.602 ± 0.069 μs    | 0.611 ± 0.83 μs    |
| Model evaluation/AR latent/rand                                    | 1.73 ± 1 μs        | 0.823 ± 1 μs        | 0.813 ± 1.1 μs      | 1.71 ± 1 μs        |
| Model evaluation/DirectInfections+Poisson/forward                  | 2.03 ± 0.72 μs     | 0.0716 ± 0.00064 ms | 0.0702 ± 0.00073 ms | 2.04 ± 0.79 μs     |
| Model evaluation/DirectInfections+Poisson/rand                     | 1.77 ± 1 μs        | 0.0705 ± 0.00063 ms | 0.0693 ± 0.00068 ms | 1.76 ± 1 μs        |
| Model evaluation/RandomWalk latent/forward                         | 1.11 ± 0.64 μs     | 1.1 ± 0.65 μs       | 1.1 ± 0.65 μs       | 1.1 ± 0.65 μs      |
| Model evaluation/RandomWalk latent/rand                            | 1.31 ± 0.82 μs     | 1.28 ± 0.81 μs      | 1.31 ± 0.84 μs      | 1.29 ± 0.8 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 7.25 ± 2.4 μs      | 0.0762 ± 0.0007 ms  | 0.0748 ± 0.00079 ms | 7.19 ± 2.4 μs      |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 4.69 ± 2.6 μs      | 0.0734 ± 0.00084 ms | 0.0721 ± 0.00088 ms | 4.6 ± 2.6 μs       |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.103 ± 0.002 s    | 1 ± 0.042 s         | 0.984 ± 0.032 s     | 0.103 ± 0.0018 s   |
| time_to_load                                                       | 5.04 ± 0.012 s     | 4.65 ± 0.025 s      | 4.64 ± 0.021 s      | 5.03 ± 0.003 s     |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | 2a78930a4e7ba2...         |
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

