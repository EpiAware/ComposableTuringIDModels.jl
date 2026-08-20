|                                                                    | v0.1.2             | v0.1.1              | v0.1.0              | 3bbf1697cffd70...  |
|:-------------------------------------------------------------------|:------------------:|:-------------------:|:-------------------:|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 11.9 ± 20 μs       | 25.6 ± 21 μs        | 30.1 ± 21 μs        | 11.4 ± 20 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.96 ± 2.6 μs      | 8.92 ± 2.7 μs       | 9.06 ± 2.7 μs       | 8.99 ± 2.7 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0464 ± 0.0012 ms | 0.0509 ± 0.0088 ms  | 0.0511 ± 0.009 ms   | 0.0462 ± 0.0013 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 5.89 ± 0.5 μs      | 0.0696 ± 0.00086 ms | 0.0698 ± 0.0008 ms  | 6 ± 0.61 μs        |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 19.6 ± 24 μs       | 0.214 ± 0.028 ms    | 0.216 ± 0.028 ms    | 18.9 ± 24 μs       |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 10.2 ± 0.52 μs     | 0.0788 ± 0.0012 ms  | 0.0785 ± 0.0012 ms  | 10.2 ± 0.48 μs     |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.121 ± 0.023 ms   | 0.174 ± 0.024 ms    | 0.174 ± 0.025 ms    | 0.12 ± 0.021 ms    |
| Model evaluation/AR latent/forward                                 | 0.679 ± 1 μs       | 0.645 ± 0.092 μs    | 0.652 ± 0.094 μs    | 0.649 ± 1 μs       |
| Model evaluation/AR latent/rand                                    | 1.95 ± 1.2 μs      | 1.25 ± 1.2 μs       | 1.94 ± 1.2 μs       | 1.94 ± 1.2 μs      |
| Model evaluation/DirectInfections+Poisson/forward                  | 2.23 ± 0.95 μs     | 0.0664 ± 0.00091 ms | 0.0667 ± 0.00093 ms | 2.24 ± 1 μs        |
| Model evaluation/DirectInfections+Poisson/rand                     | 2 ± 0.19 μs        | 0.0657 ± 0.0009 ms  | 0.0658 ± 0.00087 ms | 2 ± 1.2 μs         |
| Model evaluation/RandomWalk latent/forward                         | 1.27 ± 0.077 μs    | 1.27 ± 0.7 μs       | 1.25 ± 0.71 μs      | 1.27 ± 0.72 μs     |
| Model evaluation/RandomWalk latent/rand                            | 1.46 ± 0.93 μs     | 1.46 ± 0.95 μs      | 1.45 ± 0.95 μs      | 1.45 ± 0.96 μs     |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 7.04 ± 2.8 μs      | 0.0717 ± 0.001 ms   | 0.0721 ± 0.001 ms   | 6.75 ± 2.8 μs      |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 4.99 ± 2.9 μs      | 0.0691 ± 0.0013 ms  | 0.0691 ± 0.0013 ms  | 5.21 ± 3.1 μs      |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.146 ± 0.012 s    | 0.954 ± 0.04 s      | 0.957 ± 0.034 s     | 0.127 ± 0.0071 s   |
| time_to_load                                                       | 5.43 ± 0.044 s     | 4.85 ± 0.095 s      | 4.75 ± 0.088 s      | 5.11 ± 0.014 s     |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | 3bbf1697cffd70...         |
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

