|                                                                    | v0.1.1              | v0.1.0              | 278810260e7711...   |
|:-------------------------------------------------------------------|:-------------------:|:-------------------:|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 12.4 ± 17 μs        | 14.1 ± 18 μs        | 10.5 ± 17 μs        |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.51 ± 2.5 μs       | 8.41 ± 2.4 μs       | 8.35 ± 2.2 μs       |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0532 ± 0.0086 ms  | 0.0529 ± 0.0042 ms  | 0.0529 ± 0.0083 ms  |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0773 ± 0.00075 ms | 0.0771 ± 0.00067 ms | 0.0758 ± 0.00067 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.243 ± 0.023 ms    | 0.237 ± 0.024 ms    | 0.235 ± 0.023 ms    |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0871 ± 0.0012 ms  | 0.0871 ± 0.0012 ms  | 0.0859 ± 0.0015 ms  |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.193 ± 0.025 ms    | 0.192 ± 0.024 ms    | 0.191 ± 0.025 ms    |
| Model evaluation/AR latent/forward                                 | 0.621 ± 0.068 μs    | 0.594 ± 0.072 μs    | 0.593 ± 0.064 μs    |
| Model evaluation/AR latent/rand                                    | 0.831 ± 1 μs        | 0.797 ± 1 μs        | 0.81 ± 1 μs         |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0743 ± 0.00075 ms | 0.0744 ± 0.00068 ms | 0.0731 ± 0.00078 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0732 ± 0.0007 ms  | 0.0731 ± 0.00063 ms | 0.0719 ± 0.00073 ms |
| Model evaluation/RandomWalk latent/forward                         | 1.13 ± 0.64 μs      | 1.1 ± 0.64 μs       | 1.1 ± 0.65 μs       |
| Model evaluation/RandomWalk latent/rand                            | 1.3 ± 0.78 μs       | 1.29 ± 0.8 μs       | 1.31 ± 0.82 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0764 ± 0.00079 ms | 0.0765 ± 0.00075 ms | 0.0758 ± 0.00076 ms |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0755 ± 0.00097 ms | 0.0753 ± 0.00089 ms | 0.0741 ± 0.0009 ms  |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.872 ± 0.36 s      | 0.688 ± 5.9 s       | 0.659 ± 0.38 s      |
| time_to_load                                                       | 4.39 ± 0.041 s      | 4.38 ± 0.038 s      | 4.43 ± 0.017 s      |

|                                                                    | v0.1.1                    | v0.1.0                    | 278810260e7711...         |
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
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.48 k allocs: 22.3 kB    | 0.48 k allocs: 22.3 kB    | 0.48 k allocs: 22.3 kB    |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.477 k allocs: 21.6 kB   | 0.477 k allocs: 21.6 kB   | 0.477 k allocs: 21.6 kB   |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.93 M allocs: 0.239 GB   | 1.88 M allocs: 0.233 GB   | 0.589 M allocs: 0.0744 GB |
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   | 0.149 k allocs: 11.2 kB   | 0.149 k allocs: 11.2 kB   |

