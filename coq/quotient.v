Require Import Arith.
Require Import Lia.


Definition nonnegative_nat (n : nat) := n >= 0.
Definition positive_nat (n : nat) := n >= 1.

Theorem QPRIME_GE_Q: forall n q q' r r': nat,
  positive_nat(n) /\ nonnegative_nat(q) /\ nonnegative_nat(q') /\ nonnegative_nat(r) /\ nonnegative_nat(r') /\ n * q + r = (n - 1) * q' + r' /\ n >= q + 2 /\ r < n /\ r' < n - 1 -> q' >= q.
Proof.
  intros n q q' r r' H.
  destruct (lt_eq_lt_dec q' q) as [[Hlt | Heq] | Hgt].
  exfalso.
  assert (q' <= q - 1) by (lia).
  destruct (lt_eq_lt_dec q' 0) as [[Hqneg | Hqzero] | Hqpos].
  exfalso.
  easy.
  assert (q = 1) by (lia).
  assert (Hn: n < n) by (lia).
  apply Nat.lt_irrefl in Hn.
  contradiction.
  assert ((n - 1) * q' <= (n - 1) * (q - 1)) by (apply Nat.mul_le_mono_nonneg_l; lia).
  assert (Hr: r' < r') by (lia).
  apply Nat.lt_irrefl in Hr.
  contradiction.
  rewrite Heq.
  auto.
  apply le_S_n.
  apply le_S.
  assumption.
Qed.

Theorem QPRIME_IS_AT_MOST_Q_PLUS_1: forall n q q' r r': nat,
  positive_nat(n) /\ nonnegative_nat(q) /\ nonnegative_nat(q') /\ nonnegative_nat(r) /\ nonnegative_nat(r') /\ n * q + r = (n - 1) * q' + r' /\ n >= q + 2 /\ r < n /\ r' < n - 1 -> q' <= q + 1.
Proof.
  intros n q q' r r' H.
  pose (gap := q' - q).
  assert (q' >= q) by (apply QPRIME_GE_Q in H; lia).
  destruct (lt_eq_lt_dec 1 gap) as [[Hlt | Heq] | Hgt].
  exfalso.
  assert (q' >= q + 2) by (lia).
  assert ((n - 1) * q' >= (n - 1) * (q + 2)) by (apply Nat.mul_le_mono_nonneg_l; lia).
  assert (Hr: r > r) by (lia).
  apply Nat.lt_irrefl in Hr.
  contradiction.
  lia.
  lia.
Qed.

Theorem Q_MUST_BE_EITHER_Q_OR_Q_PLUS_1: forall q q': nat, q' >= q /\ q' <= q + 1 -> q' = q \/ q' = q + 1.
Proof.
  lia.
Qed.

Theorem MAIN: forall n q q' r r' a b: nat,
  nonnegative_nat(a) /\ nonnegative_nat(b) /\ positive_nat(n) /\ nonnegative_nat(q) /\ nonnegative_nat(q') /\ nonnegative_nat(r) /\ nonnegative_nat(r') /\ a * b = n * q + r
  /\ a * b = (n - 1) * q' + r' /\ a < n /\ b < n /\ r < n /\ r' < n - 1 /\ n > 1 -> (r > r' -> r' = r + q + 1 - n) /\ (r <= r' -> r' = r + q).
Proof.
  intros n q q' r r' a b H.
  pose (product := a * b).
  destruct (lt_eq_lt_dec product 0) as [[Hlt | Heq] | Hgt].
  exfalso.
  lia.
  assert (q = 0 /\ r = 0 /\ q' = 0 /\ r' = 0 /\ r = r') by (lia).
  split.
  - lia.
  - lia.
  -
  assert (a > 0 /\ b > 0) by (lia).
  assert (b * a < b * n) by (apply Nat.mul_lt_mono_pos_l; lia).
  assert (n * b < n * n) by (apply Nat.mul_lt_mono_pos_l; lia).
  assert (n * q < n * n -> q < n) by (apply Nat.mul_lt_mono_pos_l; lia).
  assert (n >= q + 1) by (lia).
  pose (q_plus_one := q + 1).
  destruct (lt_eq_lt_dec n q_plus_one) as [[Hlt2 | Heq2] | Hgt2].
  -- exfalso. lia.
  --
  assert (a > 0 /\ a <= n - 1 /\ b > 0 /\ b <= n - 1) by (lia).
  assert (a * b <= (n - 1) * b) by (apply Nat.mul_le_mono_pos_r; lia).
  assert ((n - 1) * b <= (n - 1) * (n - 1)) by (apply Nat.mul_le_mono_pos_l; lia).
  assert (n-1 <= 0) by (lia).
  exfalso.
  lia.
  --
  assert (n >= q + 2) by (lia).
  assert (n * q + r = (n - 1) * q' + r') by (lia).
  assert (q' <= q + 1).
  apply QPRIME_IS_AT_MOST_Q_PLUS_1 with (n:=n) (r:=r) (r':=r').
  easy.
  assert (H8:q' = q \/ q' = q + 1).
  apply Q_MUST_BE_EITHER_Q_OR_Q_PLUS_1.
  split.
  apply QPRIME_GE_Q  with (n:=n) (r:=r) (r':=r').
  easy.
  easy.
  destruct H8 as [H8a | H8b].
  assert (r' = r + q /\ r <= r') by (lia).
  split.
  lia.
  easy.
  assert (r - n = r' - q - 1) by (lia).
  lia.
Qed.