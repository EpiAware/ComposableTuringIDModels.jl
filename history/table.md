|                                                                    | v0.1.2            | v0.1.1              | v0.1.0             | cfb6177073c483...  |
|:-------------------------------------------------------------------|:-----------------:|:-------------------:|:------------------:|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 12.3 ± 20 μs      | 30.1 ± 21 μs        | 12.5 ± 19 μs       | 13.7 ± 21 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.96 ± 2.5 μs     | 8.88 ± 2.7 μs       | 8.94 ± 2.7 μs      | 9.16 ± 2.6 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.047 ± 0.0013 ms | 0.0511 ± 0.0091 ms  | 0.0507 ± 0.0019 ms | 0.0471 ± 0.0033 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 5.98 ± 0.62 μs    | 0.0728 ± 0.00088 ms | 0.0695 ± 0.0009 ms | 5.78 ± 0.46 μs     |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 20.9 ± 24 μs      | 0.217 ± 0.027 ms    | 0.206 ± 0.016 ms   | 18.9 ± 23 μs       |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 10.6 ± 0.51 μs    | 0.0793 ± 0.0015 ms  | 0.0785 ± 0.0012 ms | 10.3 ± 0.46 μs     |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.122 ± 0.022 ms  | 0.174 ± 0.023 ms    | 0.172 ± 0.021 ms   | 0.122 ± 0.022 ms   |
| Model evaluation/AR latent/forward                                 | 0.686 ± 0.99 μs   | 0.633 ± 0.087 μs    | 0.654 ± 0.11 μs    | 0.65 ± 1 μs        |
| Model evaluation/AR latent/rand                                    | 1.98 ± 1.2 μs     | 1.25 ± 1.2 μs       | 1.85 ± 1.2 μs      | 1.96 ± 1.2 μs      |
| Model evaluation/DirectInfections+Poisson/forward                  | 2.23 ± 0.93 μs    | 0.0663 ± 0.001 ms   | 0.0659 ± 0.0011 ms | 2.26 ± 0.94 μs     |
| Model evaluation/DirectInfections+Poisson/rand                     | 1.98 ± 1.2 μs     | 0.0653 ± 0.00093 ms | 0.065 ± 0.0011 ms  | 1.99 ± 0.17 μs     |
| Model evaluation/RandomWalk latent/forward                         | 1.28 ± 0.073 μs   | 1.28 ± 0.73 μs      | 1.29 ± 0.74 μs     | 1.3 ± 0.74 μs      |
| Model evaluation/RandomWalk latent/rand                            | 1.45 ± 0.93 μs    | 1.48 ± 0.97 μs      | 1.48 ± 0.96 μs     | 1.46 ± 0.95 μs     |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 7.41 ± 1.8 μs     | 0.072 ± 0.0011 ms   | 0.0714 ± 0.0011 ms | 7.81 ± 3 μs        |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 5.11 ± 2.8 μs     | 0.0683 ± 0.0012 ms  | 0.0682 ± 0.0013 ms | 5.21 ± 3.1 μs      |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.144 ± 0.01 s    | 0.949 ± 0.032 s     | 0.946 ± 0.035 s    | 0.134 ± 0.0064 s   |
| time_to_load                                                       | 4.99 ± 0.011 s    | 4.63 ± 0.05 s       | 4.73 ± 0.036 s     | 4.79 ± 0.022 s     |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | cfb6177073c483...         |
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

