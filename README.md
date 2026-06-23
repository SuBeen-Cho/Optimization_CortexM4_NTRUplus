# NTRU+ KEM on ARM Cortex-M4: Optimized Implementation and Benchmarks

Optimized implementations of NTRU+KEM (768/864/1152) on ARM Cortex-M4, with cycle and stack benchmarks measured under pqm4.

Source code will be added later; this repository currently provides the benchmark results only. The implementation folders (clean, m4opt, m4opt_mem) under each KEM are placeholders.

## Abstract

NTRU+ is an NTT-friendly lattice-based KEM. This repository reports the performance of optimized NTRU+KEM (768/864/1152) implementations on ARM Cortex-M4. The number-theoretic transform (NTT) used in polynomial multiplication is written in assembly to improve speed (m4opt), and a memory-optimized variant (m4opt_mem) additionally reuses and compacts the working buffers of key generation, encapsulation, and decapsulation to lower the runtime stack usage. Compared with the reference C implementation (clean), speed improves by about 1.2x to 1.9x depending on the operation, and the memory-optimized variant reduces stack usage by 74% to 84%. All measurements were taken with the pqm4 framework on an STM32L4R5 board.

## Variants

- clean: reference C implementation (NTRU+ / KpqClean ver2), no optimization
- m4opt: assembly optimization (asm NTT)
- m4opt_mem: assembly + memory optimization (asm NTT + stack-memory optimization)

## Benchmark Setup

Measurements were taken with pqm4 (commit a24bb4b), a PQC benchmarking framework for ARM Cortex-M4. The target board is the NUCLEO-L4R5ZI (STM32L4R5, Cortex-M4F); the compiler is arm-none-eabi-gcc 15.2 with -O3. Flashing and serial communication use OpenOCD 0.12 and the ST-Link virtual COM port.

randombytes uses the STM32 on-chip hardware TRNG provided by pqm4, and the hash functions (SHAKE, SHA-2) use pqm4's shared assembly implementation across all variants. The measured differences between variants therefore come only from the NTT and the memory layout.

### Cycle measurement
Each operation (key generation, encapsulation, decapsulation) runs on the board through pqm4's speed.bin and is measured with the Cortex-M4 cycle counter. Each operation is run 10,000 times and the median is reported. Encapsulation and decapsulation are constant time and show almost no variance; key generation includes a rejection step (invertibility check), so the median is used.

### Stack measurement
Stack usage is measured with pqm4's stack.bin. The stack region is filled with a known pattern, and after an operation runs, the bytes left untouched are subtracted to obtain the actual usage in bytes. Buffers allocated outside the procedures (keys, ciphertexts) are not included.

## Speed (cycles, median of 10,000 runs)

Values in parentheses are the speedup over clean.

### NTRU+KEM768
| Operation | clean | m4opt | m4opt_mem |
|---|---:|---:|---:|
| KeyGen | 662,103 | 346,287 (1.91x) | 361,977 (1.83x) |
| Encaps | 556,926 | 344,607 (1.62x) | 362,019 (1.54x) |
| Decaps | 709,579 | 410,895 (1.73x) | 544,747 (1.30x) |

### NTRU+KEM864
| Operation | clean | m4opt | m4opt_mem |
|---|---:|---:|---:|
| KeyGen | 767,950 | 446,897 (1.72x) | 485,023 (1.58x) |
| Encaps | 625,560 | 447,068 (1.40x) | 481,172 (1.30x) |
| Decaps | 835,673 | 535,811 (1.56x) | 708,331 (1.18x) |

### NTRU+KEM1152
| Operation | clean | m4opt | m4opt_mem |
|---|---:|---:|---:|
| KeyGen | 978,949 | 578,205 (1.69x) | 604,823 (1.62x) |
| Encaps | 831,399 | 585,890 (1.42x) | 611,965 (1.36x) |
| Decaps | 1,080,465 | 694,140 (1.56x) | 897,616 (1.20x) |

## Stack (bytes)

Values in parentheses are the reduction over clean. m4opt keeps the reference memory layout; m4opt_mem is the memory-optimized variant.

### NTRU+KEM768
| Operation | clean | m4opt | m4opt_mem |
|---|---:|---:|---:|
| KeyGen | 9,440 | 9,440 | 2,032 (-78%) |
| Encaps | 9,216 | 9,216 | 2,320 (-75%) |
| Decaps | 15,880 | 15,880 | 2,736 (-83%) |

### NTRU+KEM864
| Operation | clean | m4opt | m4opt_mem |
|---|---:|---:|---:|
| KeyGen | 10,568 | 10,568 | 2,232 (-79%) |
| Encaps | 10,312 | 10,312 | 2,728 (-74%) |
| Decaps | 17,816 | 17,816 | 3,200 (-82%) |

### NTRU+KEM1152
| Operation | clean | m4opt | m4opt_mem |
|---|---:|---:|---:|
| KeyGen | 13,952 | 13,952 | 2,840 (-80%) |
| Encaps | 13,584 | 13,584 | 3,232 (-76%) |
| Decaps | 23,608 | 23,608 | 3,840 (-84%) |

## License

MIT. See the LICENSE file.
