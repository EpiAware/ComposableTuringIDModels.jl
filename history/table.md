|                                                                    | v0.1.2             | v0.1.1              | v0.1.0              | b9bde4a305c8fa...  |
|:-------------------------------------------------------------------|:------------------:|:-------------------:|:-------------------:|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 11.7 ± 21 μs       | 29.7 ± 21 μs        | 22.1 ± 22 μs        | 12.4 ± 20 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 7.92 ± 2 μs        | 7.85 ± 2.2 μs       | 7.89 ± 1.5 μs       | 7.93 ± 1.9 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0465 ± 0.0015 ms | 0.0518 ± 0.0087 ms  | 0.0516 ± 0.0091 ms  | 0.0468 ± 0.0084 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 5.95 ± 0.48 μs     | 0.0698 ± 0.00079 ms | 0.07 ± 0.00095 ms   | 6.04 ± 0.62 μs     |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 19.3 ± 25 μs       | 0.215 ± 0.028 ms    | 0.215 ± 0.027 ms    | 18.1 ± 23 μs       |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 9.43 ± 0.39 μs     | 0.0776 ± 0.00091 ms | 0.0778 ± 0.00099 ms | 10.4 ± 2.1 μs      |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.121 ± 0.022 ms   | 0.175 ± 0.025 ms    | 0.174 ± 0.024 ms    | 0.121 ± 0.021 ms   |
| Model evaluation/AR latent/forward                                 | 0.672 ± 1 μs       | 0.637 ± 0.089 μs    | 0.635 ± 0.083 μs    | 0.653 ± 1 μs       |
| Model evaluation/AR latent/rand                                    | 1.98 ± 1.2 μs      | 0.847 ± 1.2 μs      | 0.915 ± 1.2 μs      | 1.95 ± 1.2 μs      |
| Model evaluation/DirectInfections+Poisson/forward                  | 2.25 ± 0.96 μs     | 0.0663 ± 0.00086 ms | 0.0665 ± 0.00093 ms | 2.29 ± 0.92 μs     |
| Model evaluation/DirectInfections+Poisson/rand                     | 1.98 ± 0.88 μs     | 0.0658 ± 0.00087 ms | 0.0655 ± 0.00098 ms | 2 ± 1.1 μs         |
| Model evaluation/RandomWalk latent/forward                         | 1.29 ± 0.65 μs     | 1.27 ± 0.75 μs      | 1.26 ± 0.74 μs      | 1.26 ± 0.71 μs     |
| Model evaluation/RandomWalk latent/rand                            | 1.47 ± 0.92 μs     | 1.44 ± 0.96 μs      | 1.44 ± 0.96 μs      | 1.47 ± 0.98 μs     |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 6.84 ± 2.8 μs      | 0.0714 ± 0.00092 ms | 0.0718 ± 0.001 ms   | 7.25 ± 2.7 μs      |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 5.17 ± 3 μs        | 0.0688 ± 0.001 ms   | 0.0689 ± 0.0011 ms  | 5.16 ± 3 μs        |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.134 ± 0.009 s    | 0.943 ± 0.04 s      | 0.951 ± 0.048 s     | 0.124 ± 0.0095 s   |
| time_to_load                                                       | 5.33 ± 0.04 s      | 5.26 ± 0.03 s       | 4.79 ± 0.43 s       | 5.25 ± 0.016 s     |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | b9bde4a305c8fa...         |
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
| time_to_load                                                       | 0.149 k allocs: 11.2 kB   | 0.149 k allocs: 11.2 kB   | 0.149 k allocs: 11.2 kB   | 0.149 k allocs: 11.2 kB   |

