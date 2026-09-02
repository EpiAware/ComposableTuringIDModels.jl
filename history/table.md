|                                                                    | v0.1.2            | v0.1.1              | v0.1.0              | 5423cd88248dfe...  |
|:-------------------------------------------------------------------|:-----------------:|:-------------------:|:-------------------:|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 10.1 ± 17 μs      | 10.2 ± 18 μs        | 25.3 ± 18 μs        | 10.2 ± 17 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 7.59 ± 1.9 μs     | 7.61 ± 1.7 μs       | 7.37 ± 1.6 μs       | 7.56 ± 1.7 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.048 ± 0.0018 ms | 0.0523 ± 0.0085 ms  | 0.0523 ± 0.0066 ms  | 0.0474 ± 0.0061 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 5.29 ± 0.41 μs    | 0.0748 ± 0.00077 ms | 0.0734 ± 0.00077 ms | 5.4 ± 0.41 μs      |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 16.7 ± 20 μs      | 0.23 ± 0.022 ms     | 0.227 ± 0.023 ms    | 17.2 ± 20 μs       |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 9.35 ± 0.43 μs    | 0.0826 ± 0.00092 ms | 0.0822 ± 0.001 ms   | 9.16 ± 0.37 μs     |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.127 ± 0.021 ms  | 0.189 ± 0.025 ms    | 0.19 ± 0.024 ms     | 0.127 ± 0.022 ms   |
| Model evaluation/AR latent/forward                                 | 0.625 ± 0.85 μs   | 0.606 ± 0.067 μs    | 0.621 ± 0.07 μs     | 0.608 ± 0.84 μs    |
| Model evaluation/AR latent/rand                                    | 1.74 ± 1 μs       | 1.63 ± 1 μs         | 0.876 ± 1 μs        | 1.75 ± 1 μs        |
| Model evaluation/DirectInfections+Poisson/forward                  | 2.06 ± 0.77 μs    | 0.0715 ± 0.00069 ms | 0.0701 ± 0.0007 ms  | 2.04 ± 0.82 μs     |
| Model evaluation/DirectInfections+Poisson/rand                     | 1.78 ± 1 μs       | 0.0708 ± 0.00073 ms | 0.0694 ± 0.00069 ms | 1.79 ± 1 μs        |
| Model evaluation/RandomWalk latent/forward                         | 1.14 ± 0.63 μs    | 1.12 ± 0.62 μs      | 1.12 ± 0.63 μs      | 1.13 ± 0.61 μs     |
| Model evaluation/RandomWalk latent/rand                            | 1.32 ± 0.8 μs     | 1.3 ± 0.8 μs        | 1.31 ± 0.8 μs       | 1.32 ± 0.8 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 6.9 ± 2.3 μs      | 0.0767 ± 0.00081 ms | 0.0748 ± 0.00079 ms | 7.15 ± 2.3 μs      |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 4.59 ± 2.5 μs     | 0.0735 ± 0.00086 ms | 0.0721 ± 0.00094 ms | 4.62 ± 2.5 μs      |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.105 ± 0.0032 s  | 1 ± 0.036 s         | 0.972 ± 0.033 s     | 0.105 ± 0.0069 s   |
| time_to_load                                                       | 5.14 ± 0.0074 s   | 5.09 ± 0.4 s        | 4.77 ± 0.0094 s     | 5.15 ± 0.013 s     |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | 5423cd88248dfe...         |
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

