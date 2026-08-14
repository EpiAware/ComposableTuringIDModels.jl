|                                                                    | v0.1.2             | v0.1.1              | v0.1.0              | 9eda01168efcbf... |
|:-------------------------------------------------------------------|:------------------:|:-------------------:|:-------------------:|:-----------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 9.88 ± 18 μs       | 9.72 ± 16 μs        | 9.7 ± 17 μs         | 9.87 ± 18 μs      |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.42 ± 2.4 μs      | 8.36 ± 2.4 μs       | 8.25 ± 2.3 μs       | 8.52 ± 2.4 μs     |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0475 ± 0.0071 ms | 0.0523 ± 0.0041 ms  | 0.0522 ± 0.0085 ms  | 0.0484 ± 0.004 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 5.12 ± 0.41 μs     | 0.077 ± 0.00064 ms  | 0.0757 ± 0.00065 ms | 5.22 ± 0.41 μs    |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 16.4 ± 20 μs       | 0.238 ± 0.023 ms    | 0.233 ± 0.022 ms    | 16.6 ± 20 μs      |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 9.9 ± 0.37 μs      | 0.0859 ± 0.0011 ms  | 0.0841 ± 0.0011 ms  | 10 ± 0.35 μs      |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.124 ± 0.02 ms    | 0.189 ± 0.024 ms    | 0.19 ± 0.024 ms     | 0.127 ± 0.021 ms  |
| Model evaluation/AR latent/forward                                 | 0.631 ± 0.83 μs    | 0.608 ± 0.076 μs    | 0.602 ± 0.072 μs    | 0.624 ± 0.83 μs   |
| Model evaluation/AR latent/rand                                    | 1.73 ± 1 μs        | 1.62 ± 1 μs         | 1.65 ± 1 μs         | 1.71 ± 1 μs       |
| Model evaluation/DirectInfections+Poisson/forward                  | 2.02 ± 0.18 μs     | 0.0781 ± 0.00076 ms | 0.0725 ± 0.00068 ms | 2.02 ± 0.79 μs    |
| Model evaluation/DirectInfections+Poisson/rand                     | 1.76 ± 1 μs        | 0.0744 ± 0.00066 ms | 0.0718 ± 0.00076 ms | 1.75 ± 1 μs       |
| Model evaluation/RandomWalk latent/forward                         | 1.12 ± 0.61 μs     | 1.1 ± 0.62 μs       | 1.11 ± 0.65 μs      | 1.1 ± 0.62 μs     |
| Model evaluation/RandomWalk latent/rand                            | 1.29 ± 0.8 μs      | 1.29 ± 0.78 μs      | 1.3 ± 0.79 μs       | 1.28 ± 0.78 μs    |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 7.2 ± 2.3 μs       | 0.0801 ± 0.0008 ms  | 0.0774 ± 0.00076 ms | 7.2 ± 2.4 μs      |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 4.64 ± 2.6 μs      | 0.079 ± 0.00089 ms  | 0.0742 ± 0.00087 ms | 4.61 ± 2.5 μs     |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.112 ± 0.0072 s   | 1.03 ± 0.025 s      | 1.03 ± 0.052 s      | 0.111 ± 0.0052 s  |
| time_to_load                                                       | 4.53 ± 0.042 s     | 4.48 ± 0.0064 s     | 4.55 ± 0.07 s       | 4.44 ± 0.02 s     |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | 9eda01168efcbf...         |
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

