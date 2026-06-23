# NTRU+ KEM on ARM Cortex-M4 — Optimized Implementation & Benchmarks

Optimized implementations of **NTRU+KEM (768 / 864 / 1152)** on the ARM Cortex-M4,
with cycle and stack benchmarks measured using the **[pqm4](https://github.com/mupq/pqm4)** framework.

> ⚠️ **Source code will be added later.** This repository currently publishes the
> **benchmark results**; the implementation sources (`clean`, `m4opt`, `m4opt_mem` per KEM)
> are placeholders and will be populated in a future update.

---

## Variants

| Folder | Description |
|---|---|
| `clean` | Reference C implementation (NTRU+ / KpqClean ver2), no optimization |
| `m4opt` | **Assembly-optimized** — asm NTT (no memory optimization) |
| `m4opt_mem` | **Assembly + Memory-optimized** — asm NTT + stack-memory optimization |

---

## Benchmark Setup

| Item | Value |
|---|---|
| Framework | **pqm4** (commit `a24bb4b`) |
| Board | **NUCLEO-L4R5ZI** (STM32L4R5, ARM Cortex-M4F) |
| Toolchain | arm-none-eabi-gcc 15.2, OpenOCD 0.12 |
| Build | `-O3` (`-o speed`) |
| randombytes | On-chip **hardware TRNG** (pqm4 default) |
| Hash (SHAKE/SHA-2) | pqm4 shared **assembly Keccak / SHA-2** — common to all variants |

> All variants share pqm4's assembly Keccak/SHA-2 and the hardware TRNG, so the comparison
> isolates the **NTT + memory-optimization** contributions. This is the same measurement
> basis as other KpqC Cortex-M4 results (e.g., via pqm4), enabling direct comparison.

### Cycle Measurement
Cycles are measured with the Cortex-M4 cycle counter via pqm4's `speed.bin`,
**`-i 10000`** (10,000 executions). Encapsulation/Decapsulation are constant-time
(median = min); KeyGen uses rejection sampling, so the **median over 10,000 runs** is reported.

### Stack Usage Measurement
Stack consumption is measured with pqm4's `stack.bin` (watermark/painting), in bytes,
excluding memory allocated outside the procedures (keys, ciphertexts).

---

## Performance — Speed (cycles, median of 10,000)

Parenthesis = speedup vs `clean`.

### NTRU+KEM768
| Operation | clean | m4opt (ASM) | m4opt_mem (ASM+MEM) |
|---|---:|---:|---:|
| KeyGen | 662,103 | 346,287 (1.91×) | 361,977 (1.83×) |
| Encaps | 556,926 | 344,607 (1.62×) | 362,019 (1.54×) |
| Decaps | 709,579 | 410,895 (1.73×) | 544,747 (1.30×) |

### NTRU+KEM864
| Operation | clean | m4opt (ASM) | m4opt_mem (ASM+MEM) |
|---|---:|---:|---:|
| KeyGen | 767,950 | 446,897 (1.72×) | 485,023 (1.58×) |
| Encaps | 625,560 | 447,068 (1.40×) | 481,172 (1.30×) |
| Decaps | 835,673 | 535,811 (1.56×) | 708,331 (1.18×) |

### NTRU+KEM1152
| Operation | clean | m4opt (ASM) | m4opt_mem (ASM+MEM) |
|---|---:|---:|---:|
| KeyGen | 978,949 | 578,205 (1.69×) | 604,823 (1.62×) |
| Encaps | 831,399 | 585,890 (1.42×) | 611,965 (1.36×) |
| Decaps | 1,080,465 | 694,140 (1.56×) | 897,616 (1.20×) |

---

## Performance — Stack (bytes)

Parenthesis = reduction vs `clean`. (`m4opt` keeps the reference memory layout; `m4opt_mem` is the memory-optimized variant.)

### NTRU+KEM768
| Operation | clean | m4opt | m4opt_mem |
|---|---:|---:|---:|
| KeyGen | 9,440 | 9,440 | 2,032 (−78%) |
| Encaps | 9,216 | 9,216 | 2,320 (−75%) |
| Decaps | 15,880 | 15,880 | 2,736 (−83%) |

### NTRU+KEM864
| Operation | clean | m4opt | m4opt_mem |
|---|---:|---:|---:|
| KeyGen | 10,568 | 10,568 | 2,232 (−79%) |
| Encaps | 10,312 | 10,312 | 2,728 (−74%) |
| Decaps | 17,816 | 17,816 | 3,200 (−82%) |

### NTRU+KEM1152
| Operation | clean | m4opt | m4opt_mem |
|---|---:|---:|---:|
| KeyGen | 13,952 | 13,952 | 2,840 (−80%) |
| Encaps | 13,584 | 13,584 | 3,232 (−76%) |
| Decaps | 23,608 | 23,608 | 3,840 (−84%) |

---

## License
See [LICENSE](LICENSE) (MIT).
