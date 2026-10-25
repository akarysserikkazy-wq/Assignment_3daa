- Python: 3.12.3 (CPython)
- OS: Linux-6.18.44-fc-v153-x86_64-with-glibc2.39
- CPU: Intel(R) Xeon(R) Processor @ 2.30GHz (1 logical cores visible)
- Seed: 20260101 (graph for size n uses random.Random(20260101 + n)); repeats per method: 5

| n | m | MSF cost | components | Kruskal runs (ms) | Kruskal median (ms) | Prim runs (ms) | Prim median (ms) |
|---:|---:|---:|---:|---|---:|---|---:|
| 1000 | 4000 | -682182526 | 1 | 1.8, 1.7, 1.7, 1.7, 1.7 | 1.7 | 2.8, 2.8, 2.7, 2.6, 2.8 | 2.8 |
| 10000 | 40000 | -7019931301 | 1 | 29.2, 37.4, 42.8, 36.8, 37.0 | 37.0 | 86.8, 119.1, 85.8, 85.7, 84.6 | 85.8 |
| 100000 | 250000 | -53176181339 | 1 | 482.2, 417.0, 493.7, 347.0, 342.5 | 417.0 | 1149.9, 981.5, 1113.4, 769.8, 731.7 | 981.5 |
