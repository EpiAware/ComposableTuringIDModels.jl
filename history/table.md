|                                                                    | v0.1.1              | v0.1.0              | 71e539ce009e87...   |
|:-------------------------------------------------------------------|:-------------------:|:-------------------:|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 9.55 ± 3.1 μs       | 9.78 ± 17 μs        | 9.9 ± 17 μs         |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.36 ± 2.3 μs       | 8.34 ± 2.4 μs       | 8.71 ± 2.4 μs       |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0533 ± 0.0095 ms  | 0.0527 ± 0.0092 ms  | 0.0521 ± 0.0053 ms  |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0773 ± 0.00069 ms | 0.0769 ± 0.00071 ms | 0.0759 ± 0.00073 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.229 ± 0.022 ms    | 0.231 ± 0.023 ms    | 0.233 ± 0.023 ms    |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0863 ± 0.0011 ms  | 0.0863 ± 0.0013 ms  | 0.0853 ± 0.0014 ms  |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.193 ± 0.025 ms    | 0.193 ± 0.027 ms    | 0.191 ± 0.024 ms    |
| Model evaluation/AR latent/forward                                 | 0.62 ± 0.07 μs      | 0.64 ± 0.068 μs     | 0.616 ± 0.075 μs    |
| Model evaluation/AR latent/rand                                    | 1.59 ± 0.99 μs      | 1.48 ± 1 μs         | 1.71 ± 1 μs         |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.074 ± 0.00072 ms  | 0.0742 ± 0.00074 ms | 0.0728 ± 0.0007 ms  |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0732 ± 0.0007 ms  | 0.0733 ± 0.00072 ms | 0.0718 ± 0.00066 ms |
| Model evaluation/RandomWalk latent/forward                         | 1.11 ± 0.64 μs      | 1.1 ± 0.64 μs       | 1.11 ± 0.62 μs      |
| Model evaluation/RandomWalk latent/rand                            | 1.31 ± 0.8 μs       | 1.28 ± 0.75 μs      | 1.3 ± 0.78 μs       |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0786 ± 0.00078 ms | 0.0789 ± 0.00073 ms | 0.0775 ± 0.00075 ms |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0757 ± 0.00087 ms | 0.0761 ± 0.00088 ms | 0.0746 ± 0.00087 ms |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 1.04 ± 0.029 s      | 1.03 ± 0.039 s      | 1.01 ± 0.037 s      |
| time_to_load                                                       | 4.35 ± 0.047 s      | 4.47 ± 0.028 s      | 4.47 ± 0.076 s      |

|                                                                    | v0.1.1                    | v0.1.0                    | 71e539ce009e87...         |
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

