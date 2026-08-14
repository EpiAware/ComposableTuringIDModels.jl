|                                                                    | v0.1.1              | v0.1.0              | 9eda01168efcbf...  |
|:-------------------------------------------------------------------|:-------------------:|:-------------------:|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 10.1 ± 17 μs        | 9.61 ± 3.5 μs       | 10 ± 18 μs         |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.41 ± 2.4 μs       | 8.34 ± 2.3 μs       | 8.47 ± 2.2 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0537 ± 0.0038 ms  | 0.0528 ± 0.003 ms   | 0.0475 ± 0.0075 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0776 ± 0.00075 ms | 0.0761 ± 0.00076 ms | 5.23 ± 0.41 μs     |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.238 ± 0.024 ms    | 0.225 ± 0.021 ms    | 16.6 ± 20 μs       |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.087 ± 0.0013 ms   | 0.0853 ± 0.0018 ms  | 10.1 ± 0.4 μs      |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.194 ± 0.026 ms    | 0.191 ± 0.024 ms    | 0.125 ± 0.021 ms   |
| Model evaluation/AR latent/forward                                 | 0.597 ± 0.081 μs    | 0.602 ± 0.15 μs     | 0.622 ± 0.83 μs    |
| Model evaluation/AR latent/rand                                    | 1.62 ± 1 μs         | 1.68 ± 1 μs         | 1.71 ± 1 μs        |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0743 ± 0.00084 ms | 0.0732 ± 0.0008 ms  | 2.04 ± 0.19 μs     |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0732 ± 0.00079 ms | 0.072 ± 0.00069 ms  | 1.77 ± 0.95 μs     |
| Model evaluation/RandomWalk latent/forward                         | 1.11 ± 0.63 μs      | 1.13 ± 0.62 μs      | 1.13 ± 0.6 μs      |
| Model evaluation/RandomWalk latent/rand                            | 1.29 ± 0.8 μs       | 1.3 ± 0.81 μs       | 1.3 ± 0.79 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0788 ± 0.00082 ms | 0.0776 ± 0.00081 ms | 7.24 ± 2.4 μs      |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0761 ± 0.00091 ms | 0.0747 ± 0.00088 ms | 4.64 ± 2.5 μs      |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.03 ± 0.04 s       | 1.01 ± 0.035 s      | 0.11 ± 0.011 s     |
| time_to_load                                                       | 4.66 ± 0.027 s      | 4.56 ± 0.049 s      | 4.59 ± 0.05 s      |

|                                                                    | v0.1.1                    | v0.1.0                    | 9eda01168efcbf...         |
|:-------------------------------------------------------------------|:-------------------------:|:-------------------------:|:-------------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 0.056 k allocs: 0.0508 MB | 0.056 k allocs: 0.0508 MB | 0.059 k allocs: 0.0512 MB |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 0.04 k allocs: 4.98 kB    | 0.04 k allocs: 4.98 kB    | 0.047 k allocs: 5.2 kB    |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.775 k allocs: 0.0319 MB | 0.775 k allocs: 0.0319 MB | 0.738 k allocs: 30.6 kB   |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.242 k allocs: 12.3 kB   | 0.242 k allocs: 12.3 kB   | 0.034 k allocs: 4.36 kB   |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.68 k allocs: 0.0835 MB  | 0.68 k allocs: 0.0835 MB  | 0.068 k allocs: 0.0591 MB |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.32 k allocs: 15.5 kB    | 0.32 k allocs: 15.5 kB    | 0.042 k allocs: 5.39 kB   |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 1.65 k allocs: 0.0654 MB  | 1.65 k allocs: 0.0654 MB  | 1.7 k allocs: 0.0643 MB   |
| Model evaluation/AR latent/forward                                 | 20  allocs: 2.41 kB       | 20  allocs: 2.41 kB       | 21  allocs: 2.44 kB       |
| Model evaluation/AR latent/rand                                    | 22  allocs: 2.83 kB       | 22  allocs: 2.83 kB       | 23  allocs: 2.86 kB       |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.35 k allocs: 15.8 kB    | 0.35 k allocs: 15.8 kB    | 23  allocs: 2.52 kB       |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.349 k allocs: 15.1 kB   | 0.349 k allocs: 15.1 kB   | 20  allocs: 2.67 kB       |
| Model evaluation/RandomWalk latent/forward                         | 16  allocs: 1.83 kB       | 16  allocs: 1.83 kB       | 17  allocs: 1.86 kB       |
| Model evaluation/RandomWalk latent/rand                            | 15  allocs: 2.05 kB       | 15  allocs: 2.05 kB       | 16  allocs: 2.08 kB       |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.57 k allocs: 23.7 kB    | 0.57 k allocs: 23.7 kB    | 0.153 k allocs: 8.28 kB   |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.567 k allocs: 23 kB     | 0.567 k allocs: 23 kB     | 0.148 k allocs: 8.38 kB   |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 2.99 M allocs: 0.371 GB   | 2.99 M allocs: 0.371 GB   | 0.479 M allocs: 0.269 GB  |
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   | 0.149 k allocs: 11.2 kB   | 0.149 k allocs: 11.2 kB   |

