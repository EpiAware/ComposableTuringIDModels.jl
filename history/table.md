|                                                                    | v0.1.2             | v0.1.1              | v0.1.0              | 5a8d4ab9d45ce6...  |
|:-------------------------------------------------------------------|:------------------:|:-------------------:|:-------------------:|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 10.1 ± 18 μs       | 10.7 ± 18 μs        | 10.7 ± 17 μs        | 10.5 ± 17 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.41 ± 2.3 μs      | 8.45 ± 2.4 μs       | 8.32 ± 2.2 μs       | 8.54 ± 2.2 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0478 ± 0.0076 ms | 0.0531 ± 0.0039 ms  | 0.0534 ± 0.011 ms   | 0.0486 ± 0.0045 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 5.13 ± 0.41 μs     | 0.0764 ± 0.00077 ms | 0.0775 ± 0.00083 ms | 5.15 ± 0.41 μs     |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 17.2 ± 20 μs       | 0.235 ± 0.024 ms    | 0.232 ± 0.022 ms    | 17.3 ± 20 μs       |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 10.4 ± 0.45 μs     | 0.0856 ± 0.0013 ms  | 0.0868 ± 0.0016 ms  | 10.2 ± 0.46 μs     |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.125 ± 0.021 ms   | 0.192 ± 0.025 ms    | 0.194 ± 0.027 ms    | 0.125 ± 0.021 ms   |
| Model evaluation/AR latent/forward                                 | 0.624 ± 0.86 μs    | 0.632 ± 0.078 μs    | 0.621 ± 0.077 μs    | 0.648 ± 0.84 μs    |
| Model evaluation/AR latent/rand                                    | 1.74 ± 1 μs        | 0.914 ± 1 μs        | 0.859 ± 1 μs        | 1.76 ± 1 μs        |
| Model evaluation/DirectInfections+Poisson/forward                  | 2.03 ± 0.76 μs     | 0.073 ± 0.00073 ms  | 0.0742 ± 0.00079 ms | 2.05 ± 0.7 μs      |
| Model evaluation/DirectInfections+Poisson/rand                     | 1.79 ± 1 μs        | 0.0724 ± 0.00082 ms | 0.0737 ± 0.00086 ms | 1.78 ± 0.94 μs     |
| Model evaluation/RandomWalk latent/forward                         | 1.12 ± 0.64 μs     | 1.1 ± 0.31 μs       | 1.12 ± 0.58 μs      | 1.12 ± 0.57 μs     |
| Model evaluation/RandomWalk latent/rand                            | 1.32 ± 0.81 μs     | 1.31 ± 0.79 μs      | 1.31 ± 0.79 μs      | 1.3 ± 0.78 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 7.23 ± 2.4 μs      | 0.0787 ± 0.00089 ms | 0.0791 ± 0.00093 ms | 7.24 ± 2.4 μs      |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 4.65 ± 2.6 μs      | 0.075 ± 0.001 ms    | 0.0761 ± 0.001 ms   | 4.63 ± 2.5 μs      |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.112 ± 0.006 s    | 1.02 ± 0.025 s      | 1.05 ± 0.035 s      | 0.129 ± 0.0065 s   |
| time_to_load                                                       | 4.67 ± 0.027 s     | 4.8 ± 0.085 s       | 5.04 ± 0.14 s       | 5.1 ± 0.16 s       |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | 5a8d4ab9d45ce6...         |
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

