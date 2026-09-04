|                                                                    | v0.1.2             | v0.1.1              | v0.1.0              | 2b9b89e2f9b4f6...  |
|:-------------------------------------------------------------------|:------------------:|:-------------------:|:-------------------:|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 13.7 ± 10 μs       | 13.2 ± 10 μs        | 14.2 ± 11 μs        | 17.7 ± 11 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 6.55 ± 1.5 μs      | 6.52 ± 1.1 μs       | 6.54 ± 1.1 μs       | 6.51 ± 1.2 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0442 ± 0.0058 ms | 0.0476 ± 0.0064 ms  | 0.0474 ± 0.0062 ms  | 0.0444 ± 0.0042 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 4.22 ± 0.32 μs     | 0.0733 ± 0.00083 ms | 0.0747 ± 0.00082 ms | 4.28 ± 0.31 μs     |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 20 ± 13 μs         | 0.223 ± 0.0081 ms   | 0.229 ± 0.015 ms    | 21 ± 14 μs         |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 9.56 ± 1.1 μs      | 0.0793 ± 0.00078 ms | 0.0807 ± 0.0008 ms  | 9.35 ± 1.2 μs      |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.113 ± 0.014 ms   | 0.165 ± 0.015 ms    | 0.167 ± 0.015 ms    | 0.112 ± 0.015 ms   |
| Model evaluation/AR latent/forward                                 | 0.725 ± 0.62 μs    | 0.634 ± 0.12 μs     | 0.623 ± 0.12 μs     | 0.648 ± 0.61 μs    |
| Model evaluation/AR latent/rand                                    | 1.48 ± 0.7 μs      | 1.42 ± 0.71 μs      | 1.4 ± 0.72 μs       | 1.47 ± 0.69 μs     |
| Model evaluation/DirectInfections+Poisson/forward                  | 1.86 ± 0.17 μs     | 0.0711 ± 0.00084 ms | 0.0724 ± 0.00076 ms | 1.9 ± 0.4 μs       |
| Model evaluation/DirectInfections+Poisson/rand                     | 1.5 ± 0.054 μs     | 0.0703 ± 0.0007 ms  | 0.0716 ± 0.00071 ms | 1.51 ± 0.06 μs     |
| Model evaluation/RandomWalk latent/forward                         | 0.961 ± 0.036 μs   | 0.948 ± 0.039 μs    | 0.953 ± 0.043 μs    | 0.95 ± 0.036 μs    |
| Model evaluation/RandomWalk latent/rand                            | 1.1 ± 0.53 μs      | 1.09 ± 0.55 μs      | 1.1 ± 0.54 μs       | 1.1 ± 0.55 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 6.85 ± 1.7 μs      | 0.0763 ± 0.00095 ms | 0.0774 ± 0.00089 ms | 6.5 ± 1.7 μs       |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 4.1 ± 1.8 μs       | 0.0732 ± 0.0011 ms  | 0.0742 ± 0.001 ms   | 4.04 ± 1.8 μs      |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.143 ± 0.015 s    | 0.977 ± 0.03 s      | 0.986 ± 0.025 s     | 0.139 ± 0.011 s    |
| time_to_load                                                       | 5.67 ± 0.0061 s    | 5.59 ± 0.038 s      | 5.58 ± 0.029 s      | 5.58 ± 0.021 s     |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | 2b9b89e2f9b4f6...         |
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
| time_to_load                                                       | 0.15 k allocs: 11.7 kB    | 0.15 k allocs: 11.7 kB    | 0.15 k allocs: 11.7 kB    | 0.15 k allocs: 11.7 kB    |

