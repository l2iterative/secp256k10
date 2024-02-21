## Application-specific optimization for secp256k1 ECDSA key recovery in RISC Zero

<img src="https://github.com/l2iterative/secp256k10/blob/main/title.png?raw=true" align="right" width="300">

It is necessary to give some background why we work on secp256k1 ECDSA key recovery for RISC Zero all of a sudden.

- **Succinct/Paradigm released [SP1](https://github.com/succinctlabs/sp1).** SP1 is another RISC-V zkVM and is a direct 
competitor to RISC Zero. On the one hand, we are glad because it validates our fund's vision about RISC-V alt-VM for Ethereum;
Paradigm's opinions matter a lot in the industry given that it is a sophisticated, very technically strong Ethereum developer team.
On the other hand, since RISC Zero is our portfolio company, we do feel a responsibility to "protect RISC Zero" especially when 
Succinct reports better performance in ECDSA signature verification, thanks to customized gates. The observation is that 
Succinct's approach is able to shrink the number of constraints by a factor of 3. We want RISC Zero to achieve the same speedup.
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
from RISC Zero team, that has been used in Zeth, for a factor of about 3 in terms of compute cycles. The benchmark is about the number of cycles
involved in the step within ECDSA key recovery that reconstructs the public key, which we find it representative.

| **k256** | **this** |
|----------|----------|
| 925,349  | 300,098  |

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

### Technique: adders with few RISC-V cycles

### Technique: hinted arithmetics

### Technique: fast modular quotient