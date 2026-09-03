|                                                                    | v0.1.2             | v0.1.1              | v0.1.0              | d92356b3816e4d...  |
|:-------------------------------------------------------------------|:------------------:|:-------------------:|:-------------------:|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 10.3 ± 17 μs       | 10.6 ± 18 μs        | 11 ± 18 μs          | 10.2 ± 17 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 7.44 ± 1.8 μs      | 7.53 ± 1.8 μs       | 7.45 ± 1.9 μs       | 7.55 ± 1.8 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0477 ± 0.0022 ms | 0.0529 ± 0.0084 ms  | 0.0526 ± 0.0035 ms  | 0.0477 ± 0.0024 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 5.36 ± 0.42 μs     | 0.0736 ± 0.00078 ms | 0.0737 ± 0.0008 ms  | 5.18 ± 0.41 μs     |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 16.6 ± 21 μs       | 0.228 ± 0.024 ms    | 0.228 ± 0.024 ms    | 17.2 ± 20 μs       |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 9.4 ± 0.38 μs      | 0.0829 ± 0.0011 ms  | 0.082 ± 0.0011 ms   | 10.3 ± 1.7 μs      |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.125 ± 0.021 ms   | 0.192 ± 0.027 ms    | 0.19 ± 0.024 ms     | 0.125 ± 0.021 ms   |
| Model evaluation/AR latent/forward                                 | 0.634 ± 0.84 μs    | 0.608 ± 0.069 μs    | 0.607 ± 0.075 μs    | 0.636 ± 0.83 μs    |
| Model evaluation/AR latent/rand                                    | 1.74 ± 1 μs        | 1.59 ± 1 μs         | 0.854 ± 1 μs        | 1.74 ± 1 μs        |
| Model evaluation/DirectInfections+Poisson/forward                  | 2.05 ± 0.78 μs     | 0.0706 ± 0.0008 ms  | 0.0703 ± 0.00075 ms | 2.04 ± 0.8 μs      |
| Model evaluation/DirectInfections+Poisson/rand                     | 1.8 ± 1 μs         | 0.07 ± 0.00074 ms   | 0.0692 ± 0.00071 ms | 1.77 ± 0.95 μs     |
| Model evaluation/RandomWalk latent/forward                         | 1.13 ± 0.61 μs     | 1.13 ± 0.61 μs      | 1.12 ± 0.63 μs      | 1.14 ± 0.6 μs      |
| Model evaluation/RandomWalk latent/rand                            | 1.33 ± 0.81 μs     | 1.31 ± 0.8 μs       | 1.29 ± 0.81 μs      | 1.32 ± 0.8 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 7 ± 2.3 μs         | 0.0753 ± 0.00081 ms | 0.0751 ± 0.00085 ms | 7.05 ± 2.4 μs      |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 4.65 ± 2.6 μs      | 0.0726 ± 0.00099 ms | 0.0722 ± 0.0009 ms  | 4.74 ± 2.7 μs      |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.126 ± 0.0061 s   | 0.997 ± 0.035 s     | 0.98 ± 0.041 s      | 0.114 ± 0.014 s    |
| time_to_load                                                       | 5.06 ± 0.021 s     | 5.05 ± 0.15 s       | 5.06 ± 0.23 s       | 5.3 ± 0.084 s      |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | d92356b3816e4d...         |
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

