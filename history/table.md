|                                                                    | v0.1.2             | v0.1.1              | v0.1.0              | 377cd8b7592759...  |
|:-------------------------------------------------------------------|:------------------:|:-------------------:|:-------------------:|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 11.9 ± 21 μs       | 31.5 ± 22 μs        | 21.5 ± 22 μs        | 30.1 ± 21 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 9.23 ± 2.7 μs      | 9.13 ± 2.6 μs       | 8.84 ± 2.6 μs       | 9.12 ± 2.7 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0468 ± 0.0012 ms | 0.0508 ± 0.009 ms   | 0.0515 ± 0.0095 ms  | 0.0474 ± 0.0013 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 6.05 ± 0.49 μs     | 0.0712 ± 0.00095 ms | 0.0693 ± 0.00085 ms | 5.81 ± 0.5 μs      |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 19.7 ± 26 μs       | 0.221 ± 0.028 ms    | 0.214 ± 0.028 ms    | 25.5 ± 25 μs       |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 10.6 ± 0.49 μs     | 0.0799 ± 0.0013 ms  | 0.0789 ± 0.0013 ms  | 10.3 ± 0.48 μs     |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.121 ± 0.023 ms   | 0.176 ± 0.024 ms    | 0.175 ± 0.024 ms    | 0.122 ± 0.022 ms   |
| Model evaluation/AR latent/forward                                 | 0.689 ± 1 μs       | 0.647 ± 0.14 μs     | 0.648 ± 0.11 μs     | 0.674 ± 1 μs       |
| Model evaluation/AR latent/rand                                    | 1.98 ± 1.2 μs      | 1.87 ± 1.2 μs       | 1.94 ± 1.2 μs       | 1.98 ± 1.2 μs      |
| Model evaluation/DirectInfections+Poisson/forward                  | 2.29 ± 0.96 μs     | 0.0676 ± 0.001 ms   | 0.0662 ± 0.0011 ms  | 2.27 ± 0.95 μs     |
| Model evaluation/DirectInfections+Poisson/rand                     | 2.01 ± 0.19 μs     | 0.0667 ± 0.001 ms   | 0.0653 ± 0.00099 ms | 2 ± 1.2 μs         |
| Model evaluation/RandomWalk latent/forward                         | 1.31 ± 0.24 μs     | 1.28 ± 0.75 μs      | 1.29 ± 0.74 μs      | 1.3 ± 0.093 μs     |
| Model evaluation/RandomWalk latent/rand                            | 1.48 ± 0.95 μs     | 1.48 ± 0.95 μs      | 1.46 ± 0.95 μs      | 1.49 ± 0.96 μs     |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 6.89 ± 2.7 μs      | 0.0733 ± 0.001 ms   | 0.0714 ± 0.001 ms   | 7.28 ± 2.8 μs      |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 5.17 ± 3 μs        | 0.07 ± 0.0012 ms    | 0.0685 ± 0.0012 ms  | 5.1 ± 3 μs         |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.138 ± 0.0084 s   | 0.973 ± 0.036 s     | 0.944 ± 0.035 s     | 0.143 ± 0.013 s    |
| time_to_load                                                       | 4.8 ± 0.16 s       | 4.92 ± 0.09 s       | 4.79 ± 0.044 s      | 4.78 ± 0.095 s     |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | 377cd8b7592759...         |
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

