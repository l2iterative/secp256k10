## Application-specific optimization for secp256k1 ECDSA key recovery in RISC Zero

It is necessary to give some background why we work on secp256k1 ECDSA key recovery for RISC Zero.

- **Succinct/Paradigm released [SP1](https://github.com/succinctlabs/sp1).** SP1 is another RISC-V zkVM and is a direct 
competitor to RISC Zero. On the one hand, we are glad because it validates our fund's vision about RISC-V alt-VM for Ethereum;
Paradigm's opinions matter a lot in the industry given that it is a sophisticated, very technically strong Ethereum developer team.
On the other hand, since RISC Zero is our portfolio company, we do feel a responsibility to "protect RISC Zero" especially when 
Succinct reports better performance in ECDSA signature verification, thanks to customized gates. The observation is that 
Succinct's approach is able to shrink the number of constraints by a factor of 3. We want RISC Zero to achieve the same speedup.
- **Demands from Zeth and others.** 
- **Sandiness of application-specific optimization.**

### Approach

There are a large number of techniques used in this implementation. But the general idea follows the k256 Rust crate's 
optimization strategy.

#### The k256 Rust crate's approach
- GLV endomorphism
- lookup table of size 8
- the decomposed scalar is sliced into 32 windows each of signed 4-bit
- merged doubling between the GLV endomorphism