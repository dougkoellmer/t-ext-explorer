# Current Algorithms 1–13

Algorithms 1–9 and the full Ext normal coordinates remain unchanged. New
Algorithms 10–13 follow main.tex (2026-10-12), as in the corrected Mathematica
package. They work over the specified exact field and with formal Frobenius.

| Algorithm | Public entry point | Result |
| --- | --- | --- |
| 1–6 | CompatibleBasis routines | Tagged jump basis and compatible order |
| 7 | FJReductionData / ReductionData | Validated original presentations and jump data |
| 8 | FJReduce / FJReduceReport | Full normal representative and inner witness |
| 9 | FJExtension | Full action Pi in the original monomial order |
| 10 ConstantData | FJConstantData | Column-vectorized D; rank-based h; B,Z,U,P0 |
| 11 SequenceCoordinates | FJSequenceCoordinates | Ext0/quotient bases and constant K,R,P |
| 12 SequenceActions | FJSequenceActions(Pi,data) | Pi0 and general Q_t |
| 13 ExactSequence | FJExactSequence / FJReport | Assembly, ranks, Ore identities, exactness, Anderson checks |

D = NPhi^T tensor I_e - I_d tensor NPsi. The implementation obtains these
same columns by applying H -> H NPhi - NPsi H to raw matrix units ordered
by source column, then target row. Independent original columns of D form B.
Independent rows of B determine which standard vectors complete B to [B Z].
The inverse has rows U and P0, exactly as in Algorithm 10.

The full monomials retain their original row/column/degree order. Algorithm 11
uses one scan to locate constants and positive monomials. It inserts B at the
column-vectorized constant positions and the positive identity at the existing
positive positions. Thus the returned K,R,P refer to precisely the same full
coordinates as Pi. ExtTau0Basis consists of image-D matrix descriptors first,
then positive monomial descriptors with one-based tags.

Algorithm 12 computes R C_a K^(a) coefficient by coefficient, generating each
successive K twist by one Frobenius application. Q_t is theta I_h plus the
induced target nilpotent on coker D. The quotient need not be scalar.

Algorithm 13 reuses the full action and successive twisted inclusion products.
It checks RK=I, PK=0, rank K=M-h, rank P=h, Pi K=K Pi0, P Pi=Q_t P, and the
e-th powers of the differential nilpotent parts. Additional certificates are
P Section=I and K R+Section P=I. R and Section are additive coordinate maps,
not asserted t-module splittings. Empty blocks retain their explicit shapes.

Ext0 class reduction uses full reduction and P x=0, then reads R x. The full
normal representative can have a removable constant coefficient. The separate
ZeroConstantRepresentative and ZeroConstantInnerWitness certify its removal.
Formal coefficient families are supported as well as ordinary representatives.

The partial basis shifted construction is retained as FJPartial* and
ReductionData.partial_* methods, with its old applicability checks. It appears
only as the final partial reduction example. FJGeneral* retains its older
behavior for compatibility. GaPart now returns rank-based h; epsilon and
PseudocodeGaCount are diagnostic. basis(zero_part=True) retains the partial basis
shifted monomial family; use FJSequenceCoordinates for the current Ext0 basis.

## Descriptive report adapter (0.4.5)

FJExampleCompute constructs data and the true sequence. FJExampleReport numbers
bases and documents each matrix and checked identity. FJPartialBasisComparison
computes the Ore maps W=R T_partial and WInverse=L_partial K when applicable.
FJJ2Example uses Psi_t upper-right 1+tau and keeps the general result when
partial basis is unavailable. BEGIN/END comments identify the actual true Ext0
calculation blocks; raw algorithm data remain available.
