|                                                                    | v0.1.2            | v0.1.1              | v0.1.0              | 71731c361e29b5...  |
|:-------------------------------------------------------------------|:-----------------:|:-------------------:|:-------------------:|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 10.9 ± 18 μs      | 10.4 ± 17 μs        | 10.3 ± 17 μs        | 11.1 ± 17 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 7.51 ± 1.7 μs     | 7.36 ± 1.9 μs       | 7.61 ± 1.8 μs       | 7.46 ± 1.7 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.047 ± 0.0014 ms | 0.052 ± 0.0028 ms   | 0.0516 ± 0.0067 ms  | 0.0472 ± 0.0028 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 5.23 ± 0.42 μs    | 0.0729 ± 0.00076 ms | 0.0734 ± 0.00085 ms | 5.07 ± 0.41 μs     |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 17.4 ± 19 μs      | 0.228 ± 0.022 ms    | 0.228 ± 0.023 ms    | 17.4 ± 20 μs       |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 9.17 ± 0.37 μs    | 0.0813 ± 0.00091 ms | 0.0818 ± 0.001 ms   | 10.2 ± 1.7 μs      |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.127 ± 0.023 ms  | 0.19 ± 0.024 ms     | 0.188 ± 0.024 ms    | 0.125 ± 0.023 ms   |
| Model evaluation/AR latent/forward                                 | 0.662 ± 0.83 μs   | 0.621 ± 0.076 μs    | 0.624 ± 0.071 μs    | 0.635 ± 0.84 μs    |
| Model evaluation/AR latent/rand                                    | 1.76 ± 1 μs       | 1.63 ± 0.98 μs      | 0.873 ± 1 μs        | 1.73 ± 0.99 μs     |
| Model evaluation/DirectInfections+Poisson/forward                  | 2.06 ± 0.82 μs    | 0.0702 ± 0.00075 ms | 0.0706 ± 0.00074 ms | 2.05 ± 0.8 μs      |
| Model evaluation/DirectInfections+Poisson/rand                     | 1.8 ± 1 μs        | 0.0693 ± 0.0007 ms  | 0.0696 ± 0.00074 ms | 1.76 ± 1 μs        |
| Model evaluation/RandomWalk latent/forward                         | 1.13 ± 0.55 μs    | 1.12 ± 0.63 μs      | 1.12 ± 0.61 μs      | 1.12 ± 0.6 μs      |
| Model evaluation/RandomWalk latent/rand                            | 1.31 ± 0.78 μs    | 1.33 ± 0.81 μs      | 1.3 ± 0.8 μs        | 1.29 ± 0.78 μs     |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 7.06 ± 2.4 μs     | 0.0751 ± 0.00087 ms | 0.0757 ± 0.00084 ms | 7.08 ± 2.4 μs      |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 4.51 ± 2.5 μs     | 0.0722 ± 0.0011 ms  | 0.0727 ± 0.0011 ms  | 4.61 ± 2.5 μs      |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.118 ± 0.0092 s  | 0.975 ± 0.034 s     | 0.979 ± 0.036 s     | 0.111 ± 0.0071 s   |
| time_to_load                                                       | 5.55 ± 0.058 s    | 5.2 ± 0.047 s       | 5.4 ± 0.033 s       | 5.4 ± 0.14 s       |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | 71731c361e29b5...         |
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

