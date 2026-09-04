|                                                                    | v0.1.2            | v0.1.1              | v0.1.0              | 658bf8fa28e84b... |
|:-------------------------------------------------------------------|:-----------------:|:-------------------:|:-------------------:|:-----------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 10.4 ± 1.6 μs     | 10.4 ± 1.8 μs       | 9.87 ± 1.9 μs       | 9.9 ± 10 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 4.42 ± 0.41 μs    | 4.53 ± 0.48 μs      | 4.25 ± 0.43 μs      | 4.42 ± 0.41 μs    |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 24.7 ± 2.1 μs     | 27.4 ± 5.2 μs       | 27.5 ± 3.7 μs       | 24.9 ± 4.9 μs     |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 2.7 ± 0.27 μs     | 0.0386 ± 0.00082 ms | 0.0378 ± 0.0013 ms  | 2.81 ± 0.29 μs    |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 13.8 ± 12 μs      | 0.124 ± 0.015 ms    | 0.115 ± 0.012 ms    | 13.7 ± 11 μs      |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 6.34 ± 0.45 μs    | 0.0425 ± 0.00081 ms | 0.0413 ± 0.00085 ms | 5.85 ± 0.41 μs    |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.0686 ± 0.016 ms | 0.0965 ± 0.017 ms   | 0.0936 ± 0.013 ms   | 0.0681 ± 0.016 ms |
| Model evaluation/AR latent/forward                                 | 0.42 ± 0.59 μs    | 0.377 ± 0.046 μs    | 0.39 ± 0.043 μs     | 0.383 ± 0.57 μs   |
| Model evaluation/AR latent/rand                                    | 1.14 ± 0.68 μs    | 1.09 ± 0.69 μs      | 0.635 ± 0.68 μs     | 1.11 ± 0.69 μs    |
| Model evaluation/DirectInfections+Poisson/forward                  | 1.38 ± 0.091 μs   | 0.0363 ± 0.00076 ms | 0.0357 ± 0.0011 ms  | 1.32 ± 0.6 μs     |
| Model evaluation/DirectInfections+Poisson/rand                     | 1.1 ± 0.68 μs     | 0.0361 ± 0.0015 ms  | 0.0351 ± 0.00098 ms | 1.05 ± 0.69 μs    |
| Model evaluation/RandomWalk latent/forward                         | 0.326 ± 0.43 μs   | 0.709 ± 0.44 μs     | 0.314 ± 0.43 μs     | 0.324 ± 0.44 μs   |
| Model evaluation/RandomWalk latent/rand                            | 0.379 ± 0.52 μs   | 0.394 ± 0.53 μs     | 0.338 ± 0.53 μs     | 0.367 ± 0.53 μs   |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 4.15 ± 1.5 μs     | 0.039 ± 0.0011 ms   | 0.0383 ± 0.001 ms   | 4.15 ± 1.6 μs     |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 2.62 ± 1.7 μs     | 0.0376 ± 0.001 ms   | 0.0368 ± 0.00099 ms | 2.65 ± 1.6 μs     |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.0879 ± 0.008 s  | 0.545 ± 0.021 s     | 0.524 ± 0.012 s     | 0.0817 ± 0.0042 s |
| time_to_load                                                       | 4.32 ± 0.16 s     | 4.01 ± 0.079 s      | 3.92 ± 0.039 s      | 4.01 ± 0.052 s    |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | 658bf8fa28e84b...         |
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
| time_to_load                                                       | 0.15 k allocs: 11.7 kB    | 0.15 k allocs: 11.7 kB    | 0.15 k allocs: 11.7 kB    | 0.149 k allocs: 11.2 kB   |

