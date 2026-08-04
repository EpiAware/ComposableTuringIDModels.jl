|                                                                    | v0.1.1              | v0.1.0              | db982c4b493853...   |
|:-------------------------------------------------------------------|:-------------------:|:-------------------:|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 30.8 ± 22 μs        | 21.7 ± 22 μs        | 29.7 ± 21 μs        |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.82 ± 2.5 μs       | 8.73 ± 2.6 μs       | 9.03 ± 2.6 μs       |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0507 ± 0.0075 ms  | 0.0504 ± 0.0091 ms  | 0.0506 ± 0.0033 ms  |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0706 ± 0.00088 ms | 0.0716 ± 0.00084 ms | 0.0707 ± 0.00088 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.215 ± 0.029 ms    | 0.217 ± 0.026 ms    | 0.219 ± 0.027 ms    |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0793 ± 0.0012 ms  | 0.0796 ± 0.0012 ms  | 0.0795 ± 0.0012 ms  |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.175 ± 0.024 ms    | 0.176 ± 0.023 ms    | 0.175 ± 0.023 ms    |
| Model evaluation/AR latent/forward                                 | 0.64 ± 0.1 μs       | 0.635 ± 0.092 μs    | 0.642 ± 0.099 μs    |
| Model evaluation/AR latent/rand                                    | 1.83 ± 1.2 μs       | 1.87 ± 1.2 μs       | 0.881 ± 1.2 μs      |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0665 ± 0.00095 ms | 0.0678 ± 0.00094 ms | 0.0671 ± 0.00097 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0652 ± 0.00094 ms | 0.067 ± 0.00098 ms  | 0.0662 ± 0.001 ms   |
| Model evaluation/RandomWalk latent/forward                         | 1.25 ± 0.75 μs      | 1.24 ± 0.71 μs      | 1.25 ± 0.73 μs      |
| Model evaluation/RandomWalk latent/rand                            | 1.44 ± 0.93 μs      | 1.43 ± 0.91 μs      | 1.43 ± 0.92 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0713 ± 0.0011 ms  | 0.0723 ± 0.00099 ms | 0.0697 ± 0.001 ms   |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0688 ± 0.0012 ms  | 0.0707 ± 0.0012 ms  | 0.0688 ± 0.0012 ms  |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.618 ± 1.3 s       | 0.448 ± 3.9 s       | 0.654 ± 0.061 s     |
| time_to_load                                                       | 4.52 ± 0.037 s      | 4.43 ± 0.03 s       | 4.52 ± 0.056 s      |

|                                                                    | v0.1.1                    | v0.1.0                    | db982c4b493853...         |
|:-------------------------------------------------------------------|:-------------------------:|:-------------------------:|:-------------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 0.056 k allocs: 0.0508 MB | 0.056 k allocs: 0.0508 MB | 0.056 k allocs: 0.0508 MB |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 0.04 k allocs: 4.98 kB    | 0.04 k allocs: 4.98 kB    | 0.04 k allocs: 4.98 kB    |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.775 k allocs: 0.0319 MB | 0.775 k allocs: 0.0319 MB | 0.775 k allocs: 0.0319 MB |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.242 k allocs: 12.3 kB   | 0.242 k allocs: 12.3 kB   | 0.242 k allocs: 12.3 kB   |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.68 k allocs: 0.0835 MB  | 0.68 k allocs: 0.0835 MB  | 0.68 k allocs: 0.0835 MB  |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.32 k allocs: 15.5 kB    | 0.32 k allocs: 15.5 kB    | 0.32 k allocs: 15.5 kB    |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 1.65 k allocs: 0.0654 MB  | 1.65 k allocs: 0.0654 MB  | 1.65 k allocs: 0.0654 MB  |
| Model evaluation/AR latent/forward                                 | 20  allocs: 2.41 kB       | 20  allocs: 2.41 kB       | 20  allocs: 2.41 kB       |
| Model evaluation/AR latent/rand                                    | 22  allocs: 2.83 kB       | 22  allocs: 2.83 kB       | 22  allocs: 2.83 kB       |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.35 k allocs: 15.8 kB    | 0.35 k allocs: 15.8 kB    | 0.35 k allocs: 15.8 kB    |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.349 k allocs: 15.1 kB   | 0.349 k allocs: 15.1 kB   | 0.349 k allocs: 15.1 kB   |
| Model evaluation/RandomWalk latent/forward                         | 16  allocs: 1.83 kB       | 16  allocs: 1.83 kB       | 16  allocs: 1.83 kB       |
| Model evaluation/RandomWalk latent/rand                            | 15  allocs: 2.05 kB       | 15  allocs: 2.05 kB       | 15  allocs: 2.05 kB       |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.6 k allocs: 24.1 kB     | 0.525 k allocs: 23 kB     | 0.48 k allocs: 22.3 kB    |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.597 k allocs: 23.4 kB   | 0.522 k allocs: 22.3 kB   | 0.477 k allocs: 21.6 kB   |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.9 M allocs: 0.236 GB    | 1.05 M allocs: 0.131 GB   | 1.8 M allocs: 0.224 GB    |
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   | 0.149 k allocs: 11.2 kB   | 0.149 k allocs: 11.2 kB   |

