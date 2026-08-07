|                                                                    | v0.1.1              | v0.1.0              | 011169f6594c89...   |
|:-------------------------------------------------------------------|:-------------------:|:-------------------:|:-------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 14.1 ± 22 μs        | 11.1 ± 3.4 μs       | 12 ± 22 μs          |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.81 ± 2.6 μs       | 8.94 ± 2.7 μs       | 8.87 ± 2.7 μs       |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0508 ± 0.0087 ms  | 0.0507 ± 0.0018 ms  | 0.051 ± 0.0087 ms   |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.0709 ± 0.00085 ms | 0.0705 ± 0.00085 ms | 0.0703 ± 0.00081 ms |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 0.232 ± 0.028 ms    | 0.21 ± 0.026 ms     | 0.206 ± 0.014 ms    |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 0.0798 ± 0.0013 ms  | 0.0799 ± 0.0013 ms  | 0.0786 ± 0.0012 ms  |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.176 ± 0.025 ms    | 0.174 ± 0.023 ms    | 0.174 ± 0.023 ms    |
| Model evaluation/AR latent/forward                                 | 0.633 ± 0.078 μs    | 0.643 ± 0.1 μs      | 0.639 ± 0.09 μs     |
| Model evaluation/AR latent/rand                                    | 0.883 ± 1.2 μs      | 0.865 ± 1.2 μs      | 1.23 ± 1.2 μs       |
| Model evaluation/DirectInfections+Poisson/forward                  | 0.0683 ± 0.00094 ms | 0.0678 ± 0.0012 ms  | 0.0663 ± 0.00096 ms |
| Model evaluation/DirectInfections+Poisson/rand                     | 0.0666 ± 0.00099 ms | 0.0666 ± 0.0013 ms  | 0.0653 ± 0.00096 ms |
| Model evaluation/RandomWalk latent/forward                         | 1.28 ± 0.75 μs      | 1.28 ± 0.72 μs      | 1.29 ± 0.75 μs      |
| Model evaluation/RandomWalk latent/rand                            | 1.45 ± 0.97 μs      | 1.45 ± 0.96 μs      | 1.48 ± 0.98 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 0.0733 ± 0.0011 ms  | 0.0725 ± 0.001 ms   | 0.0715 ± 0.001 ms   |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 0.0713 ± 0.0013 ms  | 0.0699 ± 0.0013 ms  | 0.0684 ± 0.0012 ms  |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.965 ± 0.037 s     | 0.973 ± 0.047 s     | 0.943 ± 0.037 s     |
| time_to_load                                                       | 4.51 ± 0.11 s       | 4.48 ± 0.031 s      | 4.49 ± 0.02 s       |

|                                                                    | v0.1.1                    | v0.1.0                    | 011169f6594c89...         |
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

