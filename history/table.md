|                                                                    | v0.1.2            | v0.1.1              | v0.1.0              | 076865ba7f9383... |
|:-------------------------------------------------------------------|:-----------------:|:-------------------:|:-------------------:|:-----------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 15.1 ± 11 μs      | 13.8 ± 10 μs        | 16.1 ± 11 μs        | 14.7 ± 11 μs      |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 7.5 ± 0.97 μs     | 7.36 ± 0.95 μs      | 7.48 ± 1.3 μs       | 7.57 ± 0.94 μs    |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.044 ± 0.0053 ms | 0.0479 ± 0.0059 ms  | 0.0475 ± 0.0065 ms  | 0.044 ± 0.0037 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 4.02 ± 0.32 μs    | 0.0739 ± 0.00067 ms | 0.074 ± 0.00064 ms  | 4 ± 0.32 μs       |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 21.3 ± 13 μs      | 0.227 ± 0.01 ms     | 0.228 ± 0.015 ms    | 20.6 ± 13 μs      |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 9.12 ± 0.41 μs    | 0.0814 ± 0.00091 ms | 0.081 ± 0.00099 ms  | 9.44 ± 0.37 μs    |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.11 ± 0.013 ms   | 0.167 ± 0.014 ms    | 0.166 ± 0.015 ms    | 0.11 ± 0.015 ms   |
| Model evaluation/AR latent/forward                                 | 0.682 ± 0.59 μs   | 0.621 ± 0.15 μs     | 0.623 ± 0.13 μs     | 0.684 ± 0.6 μs    |
| Model evaluation/AR latent/rand                                    | 1.47 ± 0.69 μs    | 1.54 ± 0.1 μs       | 1.44 ± 0.7 μs       | 1.48 ± 0.55 μs    |
| Model evaluation/DirectInfections+Poisson/forward                  | 1.85 ± 0.23 μs    | 0.072 ± 0.00075 ms  | 0.0721 ± 0.00075 ms | 1.83 ± 0.4 μs     |
| Model evaluation/DirectInfections+Poisson/rand                     | 1.49 ± 0.074 μs   | 0.0715 ± 0.00078 ms | 0.0715 ± 0.00076 ms | 1.49 ± 0.055 μs   |
| Model evaluation/RandomWalk latent/forward                         | 0.95 ± 0.034 μs   | 0.937 ± 0.041 μs    | 0.941 ± 0.042 μs    | 0.956 ± 0.038 μs  |
| Model evaluation/RandomWalk latent/rand                            | 1.09 ± 0.54 μs    | 1.09 ± 0.54 μs      | 1.09 ± 0.56 μs      | 1.09 ± 0.53 μs    |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 6.81 ± 1.7 μs     | 0.0773 ± 0.00084 ms | 0.0772 ± 0.00084 ms | 6.45 ± 1.6 μs     |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 4.07 ± 1.8 μs     | 0.0738 ± 0.00087 ms | 0.0742 ± 0.0009 ms  | 4.02 ± 1.7 μs     |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.149 ± 0.013 s   | 0.981 ± 0.022 s     | 0.979 ± 0.029 s     | 0.147 ± 0.012 s   |
| time_to_load                                                       | 5.71 ± 0.045 s    | 5.64 ± 0.045 s      | 5.64 ± 0.06 s       | 5.56 ± 0.013 s    |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | 076865ba7f9383...         |
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
| time_to_load                                                       | 0.15 k allocs: 11.7 kB    | 0.149 k allocs: 11.2 kB   | 0.15 k allocs: 11.7 kB    | 0.149 k allocs: 11.2 kB   |

