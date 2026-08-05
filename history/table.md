|                                                                    | v0.1.1              | v0.1.0              | 378b39b0cd2d25...   |
|:-------------------------------------------------------------------|:-------------------:|:-------------------:|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 24.4 ± 18 μs        | 10 ± 17 μs          | 10.3 ± 17 μs        |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.64 ± 2.3 μs       | 8.62 ± 2.4 μs       | 8.42 ± 2.5 μs       |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0536 ± 0.0036 ms  | 0.0532 ± 0.0087 ms  | 0.0519 ± 0.0016 ms  |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.077 ± 0.00073 ms  | 0.0762 ± 0.00078 ms | 0.076 ± 0.00074 ms  |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.237 ± 0.023 ms    | 0.234 ± 0.023 ms    | 0.235 ± 0.023 ms    |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0867 ± 0.0014 ms  | 0.0857 ± 0.0013 ms  | 0.0862 ± 0.0012 ms  |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.193 ± 0.024 ms    | 0.192 ± 0.023 ms    | 0.19 ± 0.022 ms     |
| Model evaluation/AR latent/forward                                 | 0.614 ± 0.077 μs    | 0.622 ± 0.078 μs    | 0.596 ± 0.079 μs    |
| Model evaluation/AR latent/rand                                    | 0.892 ± 1 μs        | 1.61 ± 1 μs         | 1.56 ± 1 μs         |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0741 ± 0.00078 ms | 0.0727 ± 0.00075 ms | 0.0725 ± 0.00065 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.073 ± 0.00077 ms  | 0.0718 ± 0.00071 ms | 0.072 ± 0.00065 ms  |
| Model evaluation/RandomWalk latent/forward                         | 1.11 ± 0.58 μs      | 1.12 ± 0.64 μs      | 1.1 ± 0.65 μs       |
| Model evaluation/RandomWalk latent/rand                            | 1.29 ± 0.79 μs      | 1.3 ± 0.8 μs        | 1.29 ± 0.78 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0791 ± 0.00087 ms | 0.0776 ± 0.0008 ms  | 0.0777 ± 0.00076 ms |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0761 ± 0.00095 ms | 0.0747 ± 0.001 ms   | 0.0745 ± 0.00096 ms |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.04 ± 0.05 s       | 1.02 ± 0.046 s      | 1.01 ± 0.045 s      |
| time_to_load                                                       | 4.61 ± 0.087 s      | 4.4 ± 0.073 s       | 4.42 ± 0.0055 s     |

|                                                                    | v0.1.1                    | v0.1.0                    | 378b39b0cd2d25...         |
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
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.57 k allocs: 23.7 kB    | 0.57 k allocs: 23.7 kB    | 0.57 k allocs: 23.7 kB    |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.567 k allocs: 23 kB     | 0.567 k allocs: 23 kB     | 0.567 k allocs: 23 kB     |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 2.99 M allocs: 0.371 GB   | 2.99 M allocs: 0.371 GB   | 2.99 M allocs: 0.371 GB   |
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   | 0.149 k allocs: 11.2 kB   | 0.149 k allocs: 11.2 kB   |

