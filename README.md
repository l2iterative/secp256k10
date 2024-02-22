## Application-specific optimization for secp256k1 ECDSA key recovery in RISC Zero

<img src="https://github.com/l2iterative/secp256k10/blob/main/title.png?raw=true" align="right" width="300">

It is necessary to give some background why we work on secp256k1 ECDSA key recovery for RISC Zero all of a sudden.

- **Succinct released [SP1](https://github.com/succinctlabs/sp1).** SP1 is another RISC-V zkVM and is a direct 
competitor to RISC Zero. On the one hand, we are glad because it validates our fund's vision about RISC-V alt-VM for Ethereum;
Paradigm's opinions matter a lot in the industry given that it is a sophisticated, very technically strong Ethereum developer team.
On the other hand, since RISC Zero is our portfolio company, we do feel a responsibility to "protect RISC Zero" especially when 
Succinct reports better performance in ECDSA signature verification, thanks to customized gates. The observation is that 
Succinct's approach is able to shrink the number of cycles by a factor of 3. We want RISC Zero to achieve the same speedup.
- **Demands from [Zeth](https://github.com/risc0/zeth), [Eclipse](https://github.com/Eclipse-Laboratories-Inc/zk-bpf), and others.** Optimism has been working with RISC Zero to bring ZK (fraud) proofs to OP Stack. 
Zeth is production-ready, but it is desirous to reduce the proof generation cost further. The main bottlenecks are now 
in ECDSA signature verification and Keccak. We believe that our technique can be generalized to ed25519, which would be 
useful for [Eclipse](https://github.com/Eclipse-Laboratories-Inc/zk-bpf), another portfolio company of ours, which works on bringing Solana VM to Ethereum ecosystem through ZK.
- **Study of the sandiness of application-specific proofs.** It is academically interesting to study application-specific 
proof systems, particularly, how to design them. SP1 creates more specialized application-specific optimization toward 
secp256k1 ECDSA signature verification. However, we feel that it can be similarly fruitful if we provide optimization for 
something more generic, such as the [big integer syscalls](https://dev.risczero.com/api/zkvm/acceleration) in RISC Zero that 
works for any 256-bit numbers. As a result of the implementation in this repository, we submitted a few requests ([#1](https://github.com/risc0/risc0/issues/1432), [#2](https://github.com/risc0/risc0/issues/1443)) 
for features to RISC Zero for even more efficient signature verification through some loosely application-specific primitives.

This repository provides an implementation that improves over existing implementation (see [here](https://github.com/risc0/risc0/tree/main/examples/ecdsa)), which is a patched k256 Rust crate, 
from RISC Zero team, that has been used in Zeth, for a factor of about 3 in terms of compute cycles. 
- k256 takes 925,349 to perform a linear combination of two points
- our implementation takes 300,098 to perform a linear combination of two points where one of the points is fixed and known

The benchmark focuses on the specific step within ECDSA key recovery that reconstructs the public key, which we find it representative.

### Background: existing approach used in the patched k256 Rust crate
The k256 Rust crate patched by RISC Zero implements ECDSA public key recovery with the following highlights.
- It uses the [GLV](https://www.iacr.org/archive/crypto2001/21390189.pdf) endomorphism.
- Modular multiplications are done with RISC Zero's bigint syscalls.
- For each base point, it (pre-)computes a lookup table of size 8.
- The decomposed scalar is sliced into 32 windows each of signed 4-bit.
- During the GLV endomorphism, the doubling is merged together.

Our implementation also uses the GLV endomorphism (but with a more transferable version of the implementation), inherits 
the same lookup table and bigint syscall, and performs slicing and doubling in a similar fashion, but with differences here 
and there. To summarize, our differences are as follows.
- For the fixed group generator point G, a larger lookup table is precomputed, which eliminates the need of doubling when computing 
scalar multiplication involving the generator point.
- For the dynamic group generator point R, a joint lookup table of size 15 is used, and the slicing becomes unsigned 4-bit.
- Apply the original GLV endomorphism algorithm with a difference that we ceil the number instead of round the number when 
decomposing the scalar. Doing so has a cost, as the resultant vector can be larger, but we handle it in our implementation.

The improvement comes from the following techniques that I will now go over.

### Technique: adders with few RISC-V cycles

We reuse an adder implementation of ours that uses 32-bit limbs that seems to result in very few cycles in RISC-V. 
The idea is to use u64 `wrapping_add` to add three numbers, and to use u32 `overflowing_add` to add two numbers. Then, 
we leave it to the Rust compiler to figure out how to prepare assembly for them. We are not aware of a more efficient 
implementation, and would be more than happy if anyone has a better solution.

```rust
#[inline(always)]
pub fn add32_and_overflow(a: u32, b: u32, carry: u32) -> (u32, u32) {
    let v = (a as u64).wrapping_add(b as u64).wrapping_add(carry as u64);
    ((v >> 32) as u32, (v & 0xffffffff) as u32)
}

#[inline(always)]
pub fn carry32_and_overflow(a: u32, carry: u32) -> (u32, u32) {
    let (v, carry) = a.overflowing_add(carry);
    (carry as u32, v)
}

#[inline]
pub fn add<const I: usize, const J: usize>(accm: &mut [u32; I], new: &[u32; J]) -> u32 {
    let mut carry = 0;
    (carry, accm[0]) = add32_and_overflow(accm[0], new[0], carry);
    for i in 1..J {
        (carry, accm[i]) = add32_and_overflow(accm[i], new[i], carry);
    }
    for i in J..I {
        (carry, accm[i]) = carry32_and_overflow(accm[i], carry);
    }
    carry
}
```

We do, however, seek to request RISC Zero to add syscalls for modular additions and subtractions, which we believe is the 
bottleneck right now for big integer arithmetics, and it might be a low-hanging fruit in Zirgen, and the current bigint syscall does leave future space of development, 
such as supporting more OP codes for this syscall. Before that happens, we will stick to this approach.

### Technique: hinted arithmetics

Scalar multiplication for plain ECDSA is generally done on projective coordinates to reduce the necessity of modular inversions, 
and this is among one of the reasons why affine coordinates are not used in the middle of the computation.

The situation is, however, very different when you are allowed to provide hints. 

Originally, in the code, we need to compute a few things:
- **Modular inverse: r_inv * r = 1.** We need to compute the inverse of `r` which is the x coordinate of the point R
involved in the computation.
- **Legendre symbol and square root.** We need to compute the Legendre symbol and square-root in order 
to recover the y coordinate.
- **Slope.** For each point doubling and addition, we need to compute the slope 
which necessitates a modular division. 

But with hints, these become different. 
- **Modular inverse: r_inv * r = 1.** We present `r_inv` as a hint and verify this hint in the guest.
- **Legendre symbol and square root.** We provide a square root of an imaginary square root, which itself serves 
as a proof of Legendre symbol, and verify that this is indeed a square root.
- **Slope.** We present the slope and verify that the slope is correct (without doing inversions) in the guest.

To implement the hints, we let the host supply the hints as input. The guest will be reading the hints in a streaming manner.

### Technique: fast modular quotient

When computing the GLV endomorphism, we need to perform a modular division where we want the quotient instead of the remainder. 
The current bigint syscall does not offer the quotient.

We identify a trick as follows to compute `a * b // n`.

- `r = a * b % n`
- `r' = a * b % (n-1)`
- if `r > r'`, then `r' = r + q + 1 - n`
- if `r <= r'`, then `r' = r + q`

Computing `r` and `r'` can be done with the bigint syscall. Since the application is security-sensitive, we find 
it necessary to prove its correctness. We prove in Coq the correctness of this trick. See [here](coq/quotient.v). The main theorem is as follows:
```coq
Theorem MAIN: forall n q q' r r' a b: nat,
  nonnegative_nat(a) /\ nonnegative_nat(b) /\ positive_nat(n) /\ nonnegative_nat(q) /\ nonnegative_nat(q') 
  /\ nonnegative_nat(r) /\ nonnegative_nat(r') /\ a * b = n * q + r
  /\ a * b = (n - 1) * q' + r' /\ a < n /\ b < n /\ r < n /\ r' < n - 1 /\ n > 1 
  -> (r > r' -> r' = r + q + 1 - n) /\ (r <= r' -> r' = r + q).
```

### License

Since the `k256` crate is dual-licensed under MIT or Apache 2.0, we adopt the same license model.
