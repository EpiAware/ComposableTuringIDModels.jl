|                                                                    | v0.1.2             | v0.1.1              | v0.1.0              | ad9da788950112...  |
|:-------------------------------------------------------------------|:------------------:|:-------------------:|:-------------------:|:------------------:|
| AD gradients/AR latent logjoint/ForwardDiff                        | 10.4 ± 17 μs       | 25.5 ± 18 μs        | 10.1 ± 18 μs        | 10.2 ± 17 μs       |
| AD gradients/AR latent logjoint/Mooncake reverse                   | 7.54 ± 1.9 μs      | 7.55 ± 1.9 μs       | 7.36 ± 2 μs         | 7.56 ± 1.5 μs      |
| AD gradients/AR latent logjoint/ReverseDiff (tape)                 | 0.0482 ± 0.0076 ms | 0.0532 ± 0.0018 ms  | 0.0522 ± 0.0038 ms  | 0.0478 ± 0.0082 ms |
| AD gradients/DirectInfections+Poisson posterior/Enzyme reverse     | 5.19 ± 0.4 μs      | 0.0751 ± 0.0008 ms  | 0.0746 ± 0.00071 ms | 5.14 ± 0.41 μs     |
| AD gradients/DirectInfections+Poisson posterior/ForwardDiff        | 17.2 ± 20 μs       | 0.232 ± 0.023 ms    | 0.222 ± 0.017 ms    | 16.7 ± 20 μs       |
| AD gradients/DirectInfections+Poisson posterior/Mooncake reverse   | 9.44 ± 0.42 μs     | 0.0827 ± 0.00088 ms | 0.0829 ± 0.00086 ms | 10.2 ± 1.6 μs      |
| AD gradients/DirectInfections+Poisson posterior/ReverseDiff (tape) | 0.127 ± 0.022 ms   | 0.19 ± 0.024 ms     | 0.192 ± 0.023 ms    | 0.125 ± 0.021 ms   |
| Model evaluation/AR latent/forward                                 | 0.641 ± 0.84 μs    | 0.606 ± 0.086 μs    | 0.605 ± 0.071 μs    | 0.645 ± 0.84 μs    |
| Model evaluation/AR latent/rand                                    | 1.75 ± 1.1 μs      | 0.837 ± 1 μs        | 1.58 ± 1 μs         | 1.77 ± 1 μs        |
| Model evaluation/DirectInfections+Poisson/forward                  | 2.06 ± 0.8 μs      | 0.0718 ± 0.0008 ms  | 0.0716 ± 0.00069 ms | 2.05 ± 0.44 μs     |
| Model evaluation/DirectInfections+Poisson/rand                     | 1.78 ± 0.97 μs     | 0.071 ± 0.00075 ms  | 0.0707 ± 0.00065 ms | 1.78 ± 1 μs        |
| Model evaluation/RandomWalk latent/forward                         | 1.13 ± 0.062 μs    | 1.13 ± 0.07 μs      | 1.12 ± 0.66 μs      | 1.12 ± 0.61 μs     |
| Model evaluation/RandomWalk latent/rand                            | 1.31 ± 0.8 μs      | 1.31 ± 0.81 μs      | 1.3 ± 0.79 μs       | 1.31 ± 0.8 μs      |
| Model evaluation/Renewal+NegativeBinomial/forward                  | 7.28 ± 2.4 μs      | 0.0763 ± 0.00083 ms | 0.0762 ± 0.0008 ms  | 7.18 ± 2.3 μs      |
| Model evaluation/Renewal+NegativeBinomial/rand                     | 4.69 ± 2.5 μs      | 0.0738 ± 0.0011 ms  | 0.0734 ± 0.00086 ms | 4.67 ± 2.6 μs      |
| Sampling/NUTS (DirectInfections+Poisson, 50 draws)                 | 0.111 ± 0.0069 s   | 1 ± 0.043 s         | 0.998 ± 0.063 s     | 0.106 ± 0.0033 s   |
| time_to_load                                                       | 5.43 ± 0.16 s      | 4.9 ± 0.086 s       | 4.76 ± 0.055 s      | 5.12 ± 0.03 s      |

|                                                                    | v0.1.2                    | v0.1.1                    | v0.1.0                    | ad9da788950112...         |
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

