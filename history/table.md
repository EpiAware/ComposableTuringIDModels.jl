|                                                                    | v0.1.2             | v0.1.1              | v0.1.0              | f7d6cc9f523b76...  |
|:-------------------------------------------------------------------|:------------------:|:-------------------:|:-------------------:|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 13.3 ± 21 μs       | 12.6 ± 20 μs        | 29.9 ± 21 μs        | 11.9 ± 21 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 8.11 ± 1.9 μs      | 8.05 ± 1.6 μs       | 7.99 ± 1.3 μs       | 8.04 ± 2.3 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0481 ± 0.0012 ms | 0.0507 ± 0.0092 ms  | 0.0505 ± 0.0087 ms  | 0.0466 ± 0.0012 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 5.76 ± 0.49 μs     | 0.0694 ± 0.00083 ms | 0.0692 ± 0.00081 ms | 5.73 ± 0.51 μs     |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 20.6 ± 24 μs       | 0.21 ± 0.016 ms     | 0.216 ± 0.028 ms    | 19.2 ± 24 μs       |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 9.45 ± 0.48 μs     | 0.0776 ± 0.001 ms   | 0.0776 ± 0.00098 ms | 9.56 ± 0.46 μs     |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.128 ± 0.023 ms   | 0.174 ± 0.024 ms    | 0.174 ± 0.023 ms    | 0.121 ± 0.023 ms   |
| Model evaluation/AR latent/forward                                 | 0.679 ± 1 μs       | 0.652 ± 0.099 μs    | 0.669 ± 0.16 μs     | 0.676 ± 0.98 μs    |
| Model evaluation/AR latent/rand                                    | 1.96 ± 1.2 μs      | 0.882 ± 1.2 μs      | 1.92 ± 1.2 μs       | 1.95 ± 1.2 μs      |
| Model evaluation/DirectInfections+Poisson/forward                  | 2.27 ± 0.96 μs     | 0.0667 ± 0.00097 ms | 0.0666 ± 0.00089 ms | 2.27 ± 0.97 μs     |
| Model evaluation/DirectInfections+Poisson/rand                     | 1.98 ± 1.2 μs      | 0.066 ± 0.0011 ms   | 0.0655 ± 0.00083 ms | 2.03 ± 0.14 μs     |
| Model evaluation/RandomWalk latent/forward                         | 1.3 ± 0.078 μs     | 1.29 ± 0.46 μs      | 1.26 ± 0.74 μs      | 1.27 ± 0.087 μs    |
| Model evaluation/RandomWalk latent/rand                            | 1.47 ± 0.93 μs     | 1.46 ± 0.91 μs      | 1.46 ± 0.93 μs      | 1.49 ± 0.92 μs     |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 7.49 ± 2.7 μs      | 0.0718 ± 0.001 ms   | 0.0716 ± 0.00098 ms | 6.97 ± 2.8 μs      |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 5 ± 2.9 μs         | 0.0692 ± 0.0013 ms  | 0.0687 ± 0.0012 ms  | 4.96 ± 2.9 μs      |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.139 ± 0.014 s    | 0.952 ± 0.037 s     | 0.944 ± 0.044 s     | 0.138 ± 0.0099 s   |
| time_to_load                                                       | 5.45 ± 0.11 s      | 4.99 ± 0.24 s       | 5.25 ± 0.4 s        | 5.32 ± 0.1 s       |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | f7d6cc9f523b76...         |
|:-------------------------------------------------------------------|:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 0.059 k allocs: 0.0512 MB | 0.056 k allocs: 0.0508 MB | 0.056 k allocs: 0.0508 MB | 0.059 k allocs: 0.0512 MB |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 0.043 k allocs: 4.64 kB   | 0.036 k allocs: 4.42 kB   | 0.036 k allocs: 4.42 kB   | 0.043 k allocs: 4.64 kB   |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.738 k allocs: 30.6 kB   | 0.775 k allocs: 0.0319 MB | 0.775 k allocs: 0.0319 MB | 0.738 k allocs: 30.6 kB   |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 0.034 k allocs: 4.36 kB   | 0.242 k allocs: 12.3 kB   | 0.242 k allocs: 12.3 kB   | 0.034 k allocs: 4.36 kB   |
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

