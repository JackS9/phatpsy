# A PROJECTED HAMILTONIAN APPROACH TO POLYATOMIC SYSTEMS

**By**  
**JACK A. SMITH**

A DISSERTATION PRESENTED TO THE GRADUATE COUNCIL OF  
THE UNIVERSITY OF FLORIDA  
IN PARTIAL FULFILLMENT OF THE REQUIREMENTS FOR THE  
DEGREE OF DOCTOR OF PHILOSOPHY

UNIVERSITY OF FLORIDA  
1978

---

*To Malissa*

---

## ACKNOWLEDGEMENTS

I would like to acknowledge the guidance and assistance given by Professor Yngve Öhrn in the course of this work. I also wish to take this opportunity to thank all the members of the Quantum Theory Project for such a stimulating environment in which to work and study. I owe special thanks to Professor Per-Olov Löwdin for providing me the many opportunities to attend summer schools, winter schools, symposia and the like in such places as Sweden, Norway, Belgium, Sanibel Island, Palm Coast and last, but not least, Gainesville.

I also acknowledge the financial support of the Northeast Regional Data Center, the Air Force Office of Scientific Research (AFOSR - 74-2656C), and the National Science Foundation (NSF - CHE74-03948, SHI76-23036).

---

## TABLE OF CONTENTS

- ACKNOWLEDGEMENTS ... iii
- LIST OF TABLES ... v
- LIST OF FIGURES ... vi
- ABSTRACT ... vii

### CHAPTERS

I. INTRODUCTION ... 1

II. THE EFFECTIVE ONE-ELECTRON OPERATOR AND THE PARTITIONING OF ITS EIGENSPACE ... 5

III. THE ENERGY-DEPENDENT PSEUDOPOTENTIAL ... 19

IV. A MODEL APPROXIMATION TO THE EFFECTIVE ONE-ELECTRON INTERATOMIC COULOMB POTENTIAL ... 28

V. AN APPROXIMATE NONLOCAL INTERATOMIC EXCHANGE POTENTIAL ... 37

VI. A MODEL HAMILTONIAN AND AN APPROXIMATE ELECTRON PROPAGATOR ... 39

VII. THE GRAND-CANONICAL HARTREE-FOCK SCHEME ... 45

VIII. A SELF-CONSISTENT CHARGE AND CONFIGURATION PROCEDURE ... 50

IX. APPLICATIONS ... 58
- A. Nitric Oxide Molecule ... 58
- B. Nickel(100) Surface ... 78
- C. Nitric Oxide on Nickel(100) Surface ... 88

X. DISCUSSION ... 107

BIBLIOGRAPHY ... 111

BIOGRAPHICAL SKETCH ... 115

---

## ABSTRACT

A method of treating polyatomic systems, finite or extended, is presented which fully exploits their "atom-in-molecule" nature. Within an independent-particle model a partitioning technique is applied to the projection of the full polyatomic space into many atomic subspaces. The subspaces are then each coupled to one another approximately through second-order in overlap in a piecewise self-consistent fashion. The inherent localized picture allows for an approximate form of the interatomic interactions without resorting to neglect of any differential overlap or use of any empirical parameters. Molecular orbitals and energies may be obtained from an approximate electron propagator which is based on a model Hamiltonian built up from atomic one-electron Hamiltonians in diagonal form. Symmetric orthonormalization of these orbitals gives a density matrix which can be used as a guide for charge redistribution within the system. A generalization of the Hartree-Fock approximation based on a statistical ensemble is employed which permits the use of fractional occupation numbers in the atomic configurations. Application of the method here includes the nitric oxide molecule and its chemisorption on a nickel(100) surface.

---

## CHAPTER I: INTRODUCTION

In the quantum mechanical treatment of polyatomic systems, whether they be finite or extended, one feature consistently emerges, and that is they appear to be built up from "atoms". The reasonable success of the cellular (Slater, 1934), valence-bond (VB) (Gerratt, 1974), atom-in-molecule (AIM) (Moffitt, 1951), and similar building-block (Adams, 1971; Gilbert, 1972) approaches is indebted to this feature. However, as one progresses from atoms to molecules serious complications are encountered in the numerical and analytic techniques which have been so successfully applied to atoms. These complications are largely due to the loss of spherical symmetry and one-center expansions. For numerical techniques, the threat of multidimensional integrations and quadratures has led to cellular approximations to the potential (Slater, 1934) and local approximations to the exchange interactions (Slater, 1971). On the other hand, analytic methods are faced with the tremendous number of complicated multicenter two-electron integrals, and one is usually forced to sacrifice both quality and quantity in the selection of basis sets (Dunning and Hay, 1977). The exclusion of core electrons by means of effective potentials (Kahn et al., 1976) has become a popular means of reducing the number of integrals. A whole gamut of semiempirical methods have also been devised and constantly revised to help overcome these difficulties.

One wonders, though, whether or not the ability to even describe the "atoms" within a system is lost in all these sacrifices. This brings us to the theme of the present work. We wish to maintain as much as possible a complete and accurate description of the "atoms" in a system at the controlled expense of an approximate, but sufficiently rigorous, nonempirical treatment of the interatomic interactions. Some corollaries to this theme will be to keep intact the molecular Hamiltonian while partitioning the discrete eigenspace (the finite space spanned by the discrete eigenfunctions) of an effective one-electron operator; to employ analytic techniques enabling proper treatment of the exchange interactions; to use suitable basis sets, such as Slater-type orbitals (STO's) of at least "double-zeta" (DZ) quality with polarization functions (DZP) when needed; to correctly compute and retain all one-center integrals and all two-center one-electron integrals; to reduce all multicenter two-electron integrals to two-center one-electron integrals by exploiting the localized picture, but without resorting to empiricism or neglect of any differential overlap; and to establish atomic valency as a better-defined and workable concept.

---

## CHAPTER II: THE EFFECTIVE ONE-ELECTRON OPERATOR AND THE PARTITIONING OF ITS EIGENSPACE

In order to illustrate the partitioning of a polyatomic system into many atomic subspaces it is sufficient, at least for the moment, to consider a diatomic molecule, A-B. The average value of the nonrelativistic many-electron Hamiltonian, with respect to a proper choice of density operator, for the diatomic molecule can be written in second quantization as:

$$\langle H \rangle = \sum_i h_i^{(A)} \langle a_i^+ a_i \rangle + \frac{1}{2} \sum_{ij} (J_{ij}^{AA} - K_{ij}^{AA}) \langle a_i^+ a_j^+ a_j a_i \rangle$$

$$+ \sum_i h_i^{(B)} \langle b_i^+ b_i \rangle + \frac{1}{2} \sum_{ij} (J_{ij}^{BB} - K_{ij}^{BB}) \langle b_i^+ b_j^+ b_j b_i \rangle$$

$$+ \frac{1}{2} \sum_i \sum_j (J_{ij}^{AB} - K_{ij}^{AB}) \langle a_i^+ b_j^+ b_j a_i \rangle$$

$$+ \frac{1}{2} \sum_j \sum_i (J_{ji}^{BA} - K_{ji}^{BA}) \langle b_j^+ a_i^+ a_i b_j \rangle \quad \text{(II-1)}$$

where the $a$'s and $b$'s are field operators corresponding to spin-orbitals assigned to atoms A and B, respectively, and the summations have been restricted accordingly.

Such a total energy functional leads to the effective one-electron operator:

$$F = -\frac{1}{2}\nabla^2 - Z^A r_{1A}^{-1} + J^A - K^A - Z^B r_{1B}^{-1} + J^B - K^B \quad \text{(II-2)}$$

where the Coulomb and exchange operators associated with center A are defined by:

$$J^A \phi(1) = \sum_i q_i^{(A)} \langle \phi_i^A | r_{12}^{-1} | \phi_i^A \rangle \phi(1) \quad \text{(II-3)}$$

$$K^A \phi(1) = \sum_i q_i^{(A)} \langle \phi_i^A | r_{12}^{-1} | \phi \rangle \phi_i^A(1) \quad \text{(II-4)}$$

where the Dirac brackets here denote integration over the coordinates of electron 2 and summation over its spin. The $q$'s are the spin-orbital occupation numbers. The spin-orbitals are now restricted to be eigensolutions to a pseudoeigenvalue problem:

$$F \phi(1) = \varepsilon \phi(1) \quad \text{(II-5)}$$

At this point we have done nothing but invoke the Hartree-Fock approximation with the only twist being the localized restriction on the spin-orbitals. We shall use such localized equivalent orbitals in our treatment and restrict our spin-orbitals to be centered on a single atomic site.

To effect the partitioning, let us look at the diatomic molecule A-B at very large interatomic separation where the eigenfunctions of F approach the atomic canonical Hartree-Fock spin-orbitals associated with the two eigenvalue problems:

$$F^A \phi_i^A[0] = \{-\frac{1}{2}\nabla^2 - Z^A r_{1A}^{-1} + J^A[0] - K^A[0]\} \phi_i^A[0] \quad \text{(II-6)}$$

$$F^B \phi_i^B[0] = \{-\frac{1}{2}\nabla^2 - Z^B r_{1B}^{-1} + J^B[0] - K^B[0]\} \phi_i^B[0] \quad \text{(II-7)}$$

where the [0] denotes this separated-atom case.

Now as we allow the two atoms to approach each other and interact, let us focus our attention on a specific atom, say A. Suppose that instead of solving (II-2) for all of its eigensolutions, we seek only those solutions which we associate with atom A:

$$F \phi_i^A[1] = \{-\frac{1}{2}\nabla^2 - Z^A r_{1A}^{-1} + J^A[1] - K^A[1] - Z^B r_{1B}^{-1} + J^B[0] - K^B[0]\} \phi_i^A[1] \quad \text{(II-8)}$$

In order for the functions on both centers A and B to be simultaneous eigenfunctions of F, they must be noninteracting through F, that is:

$$\langle \phi_i^A | F | \phi_j^B \rangle = 0 \quad \text{(II-9)}$$

must be satisfied. It is sufficient, though, that they be orthogonal:

$$\langle \phi_i^A[1] | \phi_j^B[0] \rangle = 0, \quad \forall i,j \quad \text{(II-10)}$$

We obviously could have just as easily chosen to concentrate on atom B and ended up with the eigenvalue problem:

$$F \phi_i^B[1] = \{-\frac{1}{2}\nabla^2 - Z^B r_{1B}^{-1} + J^B[1] - K^B[1] - Z^A r_{1A}^{-1} + J^A[0] - K^A[0]\} \phi_i^B[1] \quad \text{(II-11)}$$

subject to the constraints:

$$\langle \phi_i^A[0] | \phi_j^B[1] \rangle = 0, \quad \forall i,j \quad \text{(II-12)}$$

Getting back now to the problem of handling the additional constraints on the pseudoeigenvalue problem, we intend to treat this problem with the projection operator technique (Löwdin, 1961). Let us start by considering a general unconstrained function whose projection is the desired function:

$$\phi_i^A = \tilde{\phi}_i^A - \sum_j \langle \phi_j^B | \tilde{\phi}_i^A \rangle \phi_j^B = (1 - P_B) \tilde{\phi}_i^A \quad \text{(II-15)}$$

where:

$$P_B = \sum_j |\phi_j^B \rangle \langle \phi_j^B| \quad \text{(II-16)}$$

is the projection operator which projects out the "B-part" of the functions on A.

The eigenvalue problem associated with atom A in terms of these unconstrained functions would be:

$$F(1 - P_B) \tilde{\phi}_i^A = \varepsilon_i^A (1 - P_B) \tilde{\phi}_i^A \quad \text{(II-18)}$$

and the secular problem would then become:

$$|\langle \tilde{\phi}_i^A | (1 - P_B)(F - \varepsilon)(1 - P_B) | \tilde{\phi}_i^A \rangle| = 0 \quad \text{(II-19)}$$

Thus, we have transferred the restriction on the eigenfunctions to a modification of the operator. We therefore seek instead the unconstrained solutions to a projected Hamiltonian.

If we rewrite equations (II-19) and (II-20) with the unprojected operator F separated out, we get:

$$|\langle \tilde{\phi}_i^A | F + V_A - \varepsilon | \tilde{\phi}_i^A \rangle| = 0 \quad \text{(II-21)}$$

$$|\langle \tilde{\phi}_j^B | F + V_B - \varepsilon | \tilde{\phi}_j^B \rangle| = 0 \quad \text{(II-22)}$$

where the pseudopotentials are given by:

$$V_A = -(P_B F - \varepsilon P_B)(2 - P_B) - [F, P_B] \quad \text{(II-23)}$$

$$V_B = -(P_A F - \varepsilon P_A)(2 - P_A) - [F, P_A] \quad \text{(II-24)}$$

---

## CHAPTER III: THE ENERGY-DEPENDENT PSEUDOPOTENTIAL

In the previous chapter, we have essentially defined an effective atom-in-molecule one-electron Hamiltonian. Its associated eigenvalue problem has the form:

$$Q^+ F Q \phi = \varepsilon \phi \quad \text{(III-1)}$$

where F is the polyatomic analog of (II-2) and Q is a product of projection operators:

$$Q = (1 - P_B)(1 - P_{B'})(1 - P_{B''})... \quad \text{(III-2)}$$

If we expand Q and retain only terms up to second-order in differential overlap (second-order in P), and then invoke the Mulliken approximation for the terms involving interatomic differential overlap:

$$|\phi_i^B \rangle \langle \phi_j^{B'}| = \frac{1}{2} \langle \phi_i^B | \phi_j^{B'} \rangle (|\phi_i^B \rangle \langle \phi_i^B| + |\phi_j^{B'} \rangle \langle \phi_j^{B'}|) \quad \text{(III-4)}$$

we get the general form for Q:

$Q = 1 - P \quad \text{(III-5)}$

where:

$P = \sum_i^{(B)} (1 - \sum_{j \neq i}^{(B')} \langle \phi_i^B | \phi_j^{B'} \rangle^2)^{1/2} |\phi_i^B \rangle \langle \phi_i^B| \quad \text{(III-6)}$

The secular problem which we then wish to solve is of the same form as (II-19):

$|\langle \phi_i | (1-P)(F-\varepsilon)(1-P) | \phi_i \rangle| = 0 \quad \text{(III-9)}$

or in terms of a pseudopotential:

$|\langle \phi_i | F + V - \varepsilon | \phi_i \rangle| = 0 \quad \text{(III-10)}$

where, analogous to (II-23):

$V = -(PF - \varepsilon P)(2 - P) - [F,P] \quad \text{(III-11)}$

In order to simplify the form of V, we need to examine FP and PP in some detail. Let us first consider the former:

$FP = F \sum_i^{(B)} \gamma_i^B |\phi_i^B \rangle \langle \phi_i^B| \quad \text{(III-12)}$

We must recall that the eigenfunctions in (III-12) are not eigenfunctions of F but rather:

$(F + V^B) \phi_i^B = \varepsilon_i^B \phi_i^B \quad \text{(III-13)}$

with the corresponding pseudopotential associated with the environment of atom B. Rearranging (III-13) gives:

$F\phi_i^B = \varepsilon_i^B \phi_i^B - V^B \phi_i^B \quad \text{(III-14)}$

which through first-order in perturbation theory becomes:

$F\phi_i^B = \tilde{\varepsilon}_i^B \phi_i^B \quad \text{(III-15)}$

where:

$\tilde{\varepsilon}_i^B = \varepsilon_i^B - \langle \phi_i^B | V^B | \phi_i^B \rangle \quad \text{(III-16)}$

Using this now in (III-12) gives:

$FP = \sum_i^{(B)} \tilde{\varepsilon}_i^B \gamma_i^B |\phi_i^B \rangle \langle \phi_i^B| \quad \text{(III-17)}$

and since F is a self-adjoint operator, (III-15) also grants:

$PF = FP \quad \text{(III-18)}$

causing the commutator in (III-11) to vanish, leaving us with:

$V = -(F - \varepsilon)(2P - PP) \quad \text{(III-19)}$

Looking now at PP we have that:

$PP = \sum_i^{(B)} \sum_j^{(B')} \gamma_i^B \gamma_j^{B'} |\phi_i^B \rangle \langle \phi_i^B | \phi_j^{B'} \rangle \langle \phi_j^{B'}| \quad \text{(III-20)}$

If we again invoke Mulliken's approximation (III-4) then:

$PP = \sum_i^{(B)} \tilde{\gamma}_i^B \gamma_i^B |\phi_i^B \rangle \langle \phi_i^B| \quad \text{(III-21)}$

where:

$\tilde{\gamma}_i^B = \gamma_i^B + \sum_{j \neq i}^{(B')} \gamma_j^{B'} \langle \phi_i^B | \phi_j^{B'} \rangle^2 \quad \text{(III-22)}$

Putting this into (III-19) we obtain:

$V = \sum_i^{(B)} (\varepsilon - \tilde{\varepsilon}_i^B)(2 - \tilde{\gamma}_i^B)\gamma_i^B |\phi_i^B \rangle \langle \phi_i^B| \quad \text{(III-23)}$

which in the limit of an orthonormal "projection" manifold reduces to the Phillips-Kleinman (1959) pseudopotential:

$V = \sum_i^{(B)} (\varepsilon - \varepsilon_i^B) |\phi_i^B \rangle \langle \phi_i^B| \quad \text{(III-24)}$

One should note the explicit dependence of V on the energy eigenvalue ε. This energy dependence is the only remaining complication in writing down the final matrix equations for solving the secular problem in (III-10).

For the sake of clarity and simplicity, we consider an orthonormal basis centered on the atom (A) with which we are currently concerned, such that:

$\phi_i = \sum_{k=1}^M u_k C_{ki} \quad \text{(III-25)}$

or in matrix notation:

$\phi = uC \quad \text{(III-26)}$

We then consider a composite basis for all the other eigenfunctions centered on their respective atomic centers:

$\phi_i^B = \sum_k^{(B)} u_k^B C_{ki}^B \quad \text{(III-27)}$

If we now define the matrices:

$S_{kl} = \langle u_k | u_l \rangle \quad \text{(III-29)}$

$E_{kl} = \sum_i^{(B)} \varepsilon_i^B (2 - \tilde{\gamma}_i^B) \gamma_i^B C_{ki}^B C_{li}^{B*} \quad \text{(III-30)}$

$B_{kl} = \sum_i^{(B)} \tilde{\varepsilon}_i^B (2 - \tilde{\gamma}_i^B) \gamma_i^B C_{ki}^B C_{li}^{B*} \quad \text{(III-31)}$

then:

$V_{kl} = \langle u_k | V | u_l \rangle = \varepsilon(SES^+)_{kl} - (SBS^+)_{kl} \quad \text{(III-32)}$

Thus, aside from the energy dependence of V, the secular problem (III-10) reduces to finding the unitary matrix C such that:

$C^+(F + V)C = \varepsilon \quad \text{(III-33)}$

where:

$F_{kl} = \langle u_k | F | u_l \rangle \quad \text{(III-34)}$

and ε is now a diagonal matrix.

As long as this energy dependence is there, we cannot seek simultaneous eigenfunctions of the same operator. We shall now derive an approximate form for the pseudopotential which is explicitly energy independent and Hermitean. First, let us rewrite (III-33) in terms of a modified Fock matrix:

$C^+ \bar{F} C = \varepsilon \quad \text{(III-35)}$

where:

$\bar{F} = F + \varepsilon P - B \quad \text{(III-36)}$

$P = SES^+ \quad \text{(III-37)}$

$B = SBS^+ \quad \text{(III-38)}$

Now putting (III-36) back into (III-35) we have:

$C^+ FC + C^+ PC\varepsilon - C^+ BC = \varepsilon \quad \text{(III-40)}$

which upon substitution becomes:

$C^+(F + PF - B)C = \varepsilon \quad \text{(III-42)}$

So comparing this with (III-35) we have that:

$\bar{F} = F + PF - B \quad \text{(III-43)}$

which upon rearranging gives:

$\bar{F} = (1 - P)^{-1}(F - B) \quad \text{(III-44)}$

If we now expand the inverse, we generate the series:

$\bar{F} = F - B + PF - PB + PPF - ... \quad \text{(III-45)}$

from which we gather that:

$\varepsilon P = PF = BF - B + BPF - ... \quad \text{(III-46)}$

which is similar to a perturbation expansion. Since the pseudopotential matrices are already second-order in overlap, this series should be rapidly convergent. In fact, we shall take only the first two terms of (III-46) in our approximation. The non-Hermiticity is also apparent in this expansion form, and so we choose the approximate Hermitean form:

$\varepsilon P = \frac{1}{2}(PF - BP + FP - BP) \quad \text{(III-47)}$

Our modified Fock matrix then becomes:

$\bar{F} = F - B + \frac{1}{2}(PF - BP + FP - BP) \quad \text{(III-48)}$

---

## CHAPTER IV: A MODEL APPROXIMATION TO THE EFFECTIVE ONE-ELECTRON INTERATOMIC COULOMB POTENTIAL

If one can substantially reduce the number of multicenter two-electron integrals encountered in conventional analytic ab initio methods, the computational savings would be overwhelming. In the last chapter, the multicenter integrals were all reduced to at most two-center integrals by restricting the eigenfunctions to one-center expansions. The actual number of integrals, however, has not been changed. We shall in this chapter exploit the localized picture and substantially reduce the number of two-center Coulomb integrals.

The integrals which we wish to approximate are of the form:

$\langle \phi_i | J | \phi_i \rangle \quad \text{(IV-1)}$

which is the Coulombic interaction of the i-th spin-orbital with all the other spin-orbitals in the polyatomic system except those on the same atomic center (A). The effective one-electron operator is given more explicitly by:

$J = \sum_i^{(B)} q_i^B \langle \phi_i^B | r_{12}^{-1} | \phi_i^B \rangle \quad \text{(IV-2)}$

Let us restrict ourselves to the case of a diatomic molecule for the moment, where J is centered on a single atom (B). Then for large interatomic separation J approaches the approximate form:

$J = N^B r_{1B}^{-1} \quad \text{(IV-4)}$

where the total electronic charge is reduced to a point charge on the distant atomic center. At lesser interatomic separation, where the charge density can no longer be considered a point charge, the interelectronic distance must undergo an effective dilation to account for the more diffuse charge. The following model potential function, which we now choose, has these desired properties:

$J = N^B \chi^B(\vec{r}_{1B}) r_{1B}^{-1} \quad \text{(IV-6)}$

where χ is a screening function with the asymptotic behavior:

$\chi(0) = 1 \quad \text{(IV-7)}$

$\chi(\infty) = 0 \quad \text{(IV-8)}$

Such a screening function is inherent to the Thomas-Fermi model (Thomas, 1928; Fermi, 1928) of the atom. A significant feature of the Thomas-Fermi screening function is that it is universal for all neutral atoms with respect to the dimensionless variable:

$x = r_{1B}[(3\pi/8)Z^B]^{-1/3} \quad \text{(IV-9)}$

We make the choice:

$\chi^B(\vec{r}_{1B}) = \sum_t A_t^B \chi_t(n,l,m,a;\vec{r}_{1B}) \quad \text{(IV-10)}$

where:

$\chi_t(n,l,m,a;\vec{r}_{1B}) = r_{1B}^{-1} \exp(-ax_{1B}) Y_{lm}(\theta_{1B}, \phi_{1B}) \quad \text{(IV-11)}$

and where the exponential factors are, in principle, universal and could be used for all atoms. Any angular distortions are taken into account with the real normalized spherical harmonics. The linear coefficients are to be determined by fitting the approximate form of J in (IV-6) to its true form in (IV-5) along with the constraint in (IV-7).

The linear coefficients are computationally quite simple to handle since they can be taken outside the integrals in which they appear. The coefficients can be determined by substituting (IV-6) for (IV-5) in the one-center Coulomb integrals occurring on atom B and matching their values:

$\langle \phi_i^B | J^B | \phi_i^B \rangle = \langle \phi_i^B | \tilde{J}^B | \phi_i^B \rangle \quad \text{(IV-12)}$

In more detail we have:

$N^B \langle \phi_i^B | r_{1B}^{-1} | \phi_i^B \rangle - \sum_t A_t^B \langle \phi_i^B | \chi_t(\vec{r}_{1B}) r_{1B}^{-1} | \phi_i^B \rangle = \sum_j^{(B)} q_j^B \langle \phi_i^B \phi_j^B | r_{12}^{-1} | \phi_j^B \phi_i^B \rangle \quad \text{(IV-13)}$

where the Coulomb integrals on the righthand side of the equation are one-center integrals which we properly compute and retain. We have as many equations like (IV-13) as there are eigenfunctions on that atom.

---

## CHAPTER V: AN APPROXIMATE NONLOCAL INTERATOMIC EXCHANGE POTENTIAL

In the last chapter, we substantially reduced the number of two-center two-electron Coulomb integrals by exploiting the localized picture. We shall now attempt a similar reduction in the analogous exchange integrals.

The two-center exchange integrals are of the form:

$\langle \phi_i | K | \phi_i \rangle \quad \text{(V-1)}$

where K is the effective nonlocal one-electron operator for the exchange interaction with all the other spin-orbitals (of same spin) in the polyatomic system except those associated with this center (A). The operator K is defined such that:

$K\phi_i(1) = \sum_j^{(B)} q_j^B \langle \phi_j^B | r_{12}^{-1} | \phi_i \rangle \phi_j^B(1) \quad \text{(V-2)}$

Let us consider the asymptotic behaviour when the average interelectronic distance approaches the internuclear distance:

$\langle \phi_i | K | \phi_i \rangle = r_{AB}^{-1} \sum_j^{(B)} q_j^B \langle \phi_i | \phi_j^B \rangle \langle \phi_j^B | \phi_i \rangle \quad \text{(V-3)}$

$K = r_{AB}^{-1} \sum_j^{(B)} q_j^B |\phi_j^B \rangle \langle \phi_j^B| \quad \text{(V-4)}$

Since K is nonlocal, we can no longer use any of its "local" properties to help determine an intermediate form for its approximation. Our endeavor thus far has been to retain a theory which is valid through at least second-order in diatomic overlap, and so it is felt that retaining K in its asymptotic form (V-4), which is second-order in overlap, would not be inconsistent with any other approximations made thus far.

---

## CHAPTER VI: A MODEL HAMILTONIAN AND AN APPROXIMATE ELECTRON PROPAGATOR

According to the Heisenberg equation of motion for the electron propagator (Linderberg and Öhrn, 1973), in the energy representation, we have in matrix form that:

$G_{ij}(E) = \langle\langle a_i ; a_j^+ \rangle\rangle_E = E^{-1}[\langle [a_i, a_j^+]_+ \rangle + \langle\langle [a_i, H] ; a_j^+ \rangle\rangle_E] \quad \text{(VI-1)}$

where:

$a_i = \langle \phi_i | \psi \rangle \quad \text{and} \quad a_i^+ = \langle \psi^+ | \phi_i \rangle \quad \text{(VI-2)}$

are the annihilation and creation operators, respectively, with the anticommutation relations:

$[a_i, a_j]_+ = [a_i^+, a_j^+]_+ = 0 \quad \text{and} \quad [a_i, a_j^+]_+ = \delta_{ij} \quad \text{(VI-3)}$

The many-electron Hamiltonian, in second quantized form, is:

$H = \sum_{ij} h_{ij} a_i^+ a_j + \frac{1}{2} \sum_{ijkl} (kl|ij) a_i^+ a_k^+ a_l a_j \quad \text{(VI-4)}$

where:

$h_{ij} = \langle \phi_i | h | \phi_j \rangle \quad \text{(VI-5)}$

$(kl|ij) = \langle \phi_k \phi_l | r_{12}^{-1} | \phi_i \phi_j \rangle \quad \text{(VI-6)}$

Equation (VI-1) can be iterated to yield:

$G^{-1}(E) = ES^{-1} + FS^{-1}E^{-2} + ... \quad \text{(VI-7)}$

where:

$F_{ij} = \langle\langle [a_i, H], a_j^+ \rangle\rangle \quad \text{(VI-8)}$

is the first moment and so on. Substituting the Hamiltonian (VI-4) into (VI-8), the first moment becomes:

$F_{ij} = h_{ij} + \sum_{kl} [(ij|kl) - (il|kj)] \langle a_k^+ a_l \rangle \quad \text{(VI-9)}$

which has the same form as the effective one-electron Hamiltonian as originally presented by Fock (1932) and will henceforth be called the Fock matrix.

Suppose now that instead of using the correct many-electron Hamiltonian, one substitutes into equation (VI-8) an approximate one-electron model Hamiltonian of the form:

$H = \sum_k \varepsilon_k a_k^+ a_k = \sum_k \varepsilon_k n_k \quad \text{(VI-10)}$

where the n's are just the occupation number operators and the ε's are real, negative energies of the noninteracting electrons. The Fock matrix then takes the approximate form:

$F_{ij} = \langle\langle [a_i, H], a_j^+ \rangle\rangle = \sum_k S_{ik} \varepsilon_k S_{kj} \quad \text{(VI-11)}$

or in matrix form:

$F = S\varepsilon S \quad \text{(VI-12)}$

where ε is a diagonal matrix. Equation (VI-7) in matrix form becomes:

$G(E) = ES^{-1} - E^{-2}S\varepsilon S + ... = (ES^{-1} - \varepsilon)^{-1} \quad \text{(VI-13)}$

where one can see that the poles of the propagator will occur at the eigenvalues of the Fock matrix. The corresponding eigenvalue problem:

$FC = SC\varepsilon \quad \text{(VI-15)}$

can be reduced to the simpler one:

$F'C' = C'\varepsilon \quad \text{(VI-16)}$

where:

$C'^+ C' = C'C'^+ = 1 \quad \text{(VI-17)}$

$F' = X^{-1}SX^{-1} \quad \text{(VI-18)}$

with X being the diagonal matrix with elements:

$X_{kk} = (-\varepsilon_k)^{1/2} \quad \text{(VI-19)}$

This can be considered as the diagonalization of a Fock matrix in a basis which is energy-weighted Löwdin orthogonalized, and this is the reason for the term "Energy-Weighted Maximum Overlap" (EWMO) used by Linderberg et al. (1976) to describe this method. The molecular orbital coefficients in the original basis are given by:

$C_{ki} = C'_{ki}(\varepsilon_i/\varepsilon_k)^{1/2} \quad \text{(VI-20)}$

---

## CHAPTER VII: THE GRAND-CANONICAL HARTREE-FOCK SCHEME

The notion of an atom in a molecule led Slater (1970) to consider an energy functional obtained as an average over multiplets arising from a given configuration. His extension of this idea to configurations with fractional occupations, known as the Hyper-Hartree-Fock method, however, has been subject to some criticism because of the appearance of off-diagonal Lagrangian multipliers. In this same spirit we wish to introduce an average based on the statistical mechanical concept of an ensemble (Abdulnur et al., 1972).

A system of noninteracting electrons described by the Hamiltonian:

$H = \sum_i \varepsilon_i n_i \quad \text{(VII-1)}$

can be described by the Grand Canonical partition function related to the density operator:

$\rho = \prod_i (1 - n_i + z_i n_i)/(1 + z_i) \quad \text{(VII-2)}$

where:

$z_i = \exp[-(\varepsilon_i - \mu)/T] \quad \text{(VII-3)}$

The parameters μ and T are the thermodynamic chemical potential and absolute temperature. The average value of an operator A for such a system is then given in terms of its trace with respect to the density operator:

$\langle A \rangle = \text{Tr}\{A\rho\} \quad \text{(VII-4)}$

In particular, the average value of the number operator is:

$\langle n_i \rangle = z_i/(1 + z_i) \quad \text{(VII-5)}$

Using this in expression (VII-2) gives a unique way of defining the density operator from a given set of occupation numbers:

$\rho = \prod_i [1 - \langle n_i \rangle + (2\langle n_i \rangle - 1)n_i] \quad \text{(VII-6)}$

If one considers an electron-conserving (Canonical) average value of the many-electron Hamiltonian in the orthonormal basis which diagonalizes the density matrix, one has:

$\langle H \rangle = \sum_i h_{ii} \langle n_i \rangle + \frac{1}{2} \sum_{ij} (J_{ij} - K_{ij}) \langle n_i n_j \rangle \quad \text{(VII-7)}$

If, however, the average is performed with respect to the Grand-Canonical density operator then:

$\langle n_i n_j \rangle = \langle n_i \rangle \langle n_j \rangle \quad \text{(VII-8)}$

and equation (VII-7) can be written:

$\langle H \rangle = \sum_i \varepsilon_i \langle n_i \rangle - \frac{1}{2} \sum_{ij} V_{ij} \langle n_i \rangle \langle n_j \rangle \quad \text{(VII-9)}$

where the effective one-electron energy is given by:

$\varepsilon_i = h_{ii} + \sum_j V_{ij} \langle n_j \rangle \quad \text{(VII-10)}$

and the effective interaction energy by:

$V_{ij} = J_{ij} - K_{ij} \quad \text{(VII-11)}$

One recognizes that these energies can be interpreted as the first and second partial derivatives of the total energy functional:

$\partial \langle H \rangle / \partial \langle n_i \rangle = \varepsilon_i \quad \text{(VII-12)}$

$\partial^2 \langle H \rangle / \partial \langle n_i \rangle \partial \langle n_j \rangle = \partial \varepsilon_i / \partial \langle n_j \rangle = V_{ij} \quad \text{(VII-13)}$

Such interpretations lead to the finite-difference approximations to ionization energies:

$\Delta E^{(i)} = \varepsilon_i \Delta \langle n_i \rangle \quad \text{(VII-14)}$

and excitation energies:

$\Delta E^{(i \to j)} = \varepsilon_j \Delta \langle n_j \rangle + \varepsilon_i \Delta \langle n_i \rangle + V_{ij} \Delta \langle n_i \rangle \Delta \langle n_j \rangle \quad \text{(VII-15)}$

for which Koopmans' (1935) theorem is a special case.

---

## CHAPTER VIII: A SELF-CONSISTENT CHARGE AND CONFIGURATION PROCEDURE

Up to now we have said very little about the manner in which the spin-orbital occupations are assigned. This is somewhat of a sensitive subject, here and in many other "building-block" approaches. We would like to make an analogy between our method and the single-configuration of nonorthogonal orbitals method. This method amounts to a condensation of many configurations to one, built up from distorted (hybridized) atomic orbitals.

The concept of "atomic valency" arises when molecules are allowed to separate into their constituent atoms. When the most general linear combination of spin couplings is formed, the resulting separated atomic states can be regarded as well-defined "valence states". In our treatment we use spin-orbitals, which from the beginning puts space and spin on similar footing. A single configuration of spin-orbitals will, in general, not correspond to a proper spin-coupled state, but the atomic "valence states" remain rather well-defined, since each valence electron is accommodated in a distinct atomic orbital with a distinct spin.

The procedure is to start out with a configuration which corresponds to a covalent structure that one generates from valence shell electron-pair repulsion theory (VSEPR) as described in any general chemistry text. The spins are assigned such that the electrons on each atom are coupled to give maximum resultant and such that the resultant spin on each atom alternates in sign with respect to each of its bonded neighbors.

Consider the example of a carbon monoxide molecule. It has the valid structure:

$\text{:C≡O:} \quad \text{(VIII-1)}$

corresponding to the configuration:

$(C1s)(C1s')(N1s)(N1s')(C2s)(C2s')(N2s)(N2s')$
$(C2p_z)(N2p_z')(C2p_x)(N2p_x')(C2p_y)(C2p_y')(N2p_y') \quad \text{(VIII-2)}$

where the primes denote opposite spin. The main point here is that this structure separates into the atomic states represented by the configurations:

$(C1s)(C1s')(C2s)(C2s')(C2p_x)(C2p_x')(C2p_y)(C2p_z) \quad \text{(VIII-3)}$

$(N1s)(N1s')(N2s)(N2s')(N2p_z')(N2p_y') \quad \text{(VIII-4)}$

each satisfying Hund's rule.

As mentioned earlier there are other valid (ionic) structures one can write down for carbon monoxide besides the covalent one. By employing fractional occupations in our configuration we can go smoothly from one structure to another. In general, the molecular state itself might best be described by such a configuration with fractional occupations which goes "adiabatically" into the separated atomic states with integral occupations.

---

## CHAPTER IX: APPLICATIONS

### A. Nitric Oxide Molecule

Despite the important role of nitric oxide and its positive ion in the upper atmosphere and its ecologically undesirable presence in the exhaust emissions of our ever-so-popular automobiles down here on earth, there has been relatively little theoretical electronic structure work reported for this first-row diatomic molecule.

For nitric oxide in its ground (X²Π) state, we will use the configuration:

$(O1s)(O1s')(N1s)(N1s')(O2s)(O2s')(N2s)(N2s')$
$(O2p_z)(N2p_z')(O2p_x)(N2p_x')(O2p_y)(O2p_y')(N2p_y') \quad \text{(IX-1)}$

in terms of "perfect-paired" atomic spin-orbitals where the primes denote the minority spin. This is unitarily equivalent to a molecular configuration:

$(1\sigma)(1\sigma')(2\sigma)(2\sigma')(3\sigma)(3\sigma')(4\sigma)(4\sigma')$
$(5\sigma)(5\sigma')(1\pi_x)(1\pi_x')(1\pi_y)(1\pi_y')(2\pi_y') \quad \text{(IX-2)}$

The basis sets used for nitrogen and oxygen are STO's of double-zeta quality plus d-type polarization functions, for a total of 15 functions on each atom.

### Table 1: Basis sets for nitrogen and oxygen

| Type | Nitrogen | Oxygen |
|------|----------|--------|
| 1s   | 6.66     | 7.66   |
| 2s   | 1.95     | 2.25   |
| 2s'  | 3.78     | 4.55   |
| 2p   | 2.00     | 2.25   |
| 2p'  | 3.75     | 4.55   |
| 3d   | 3.00     | 3.50   |

### Table 2: Model potential exponents for nitrogen and oxygen

| Type | Nitrogen | Oxygen |
|------|----------|--------|
| s₁   | 0.652    | 0.715  |
| s₂   | 1.303    | 1.431  |
| s₃   | 2.607    | 2.861  |
| s₄   | 5.213    | 5.722  |
| p    | 1.303    | 1.431  |

### Table 3: Calculated and experimental ionization potentials for nitric oxide (X²Π)

| Ion State | MO    | Expt⁴ | Expt^b | Expt^c | This Calc | Calc^d | Calc^e | Calc^f |
|-----------|-------|-------|--------|--------|-----------|--------|--------|--------|
| ¹Σ⁺       | 2π    | 9.26  | 9.25   | 9.3    | 9.9       | 10.3   | 9.3    | 9.3    |
| ³Π        | 5σ    | 15.7  | --     | 15.5   | 17.2      | 13.2   | 15.8   | 15.9   |
| ¹Π        | 5σ    | 16.6  | 16.5   | 16.5   | 17.5      | 14.0   | 16.8   | 16.8   |
| ³Σ⁻       | 1π    | 16.9  | 16.9   | 16.9   | 18.1      | 16.9   | 17.5   | 17.5   |
| ¹Σ⁻       | 1π    | 18.3  | 18.2   | 18.3   | 18.4      | 17.8   | 18.8   | 18.8   |
| ³Σ⁺       | 4σ    | 20.1  | --     | 20.1   | 21.9      | 19.0   | 20.3   | 20.4   |
| ¹Σ⁺       | 4σ    | 21.7  | 21.9   | 21.9   | 21.4      | 20.2   | 21.9   | 21.9   |

⁴Turner et al. (1970)  
^b Collin and Natalis (1968)  
^c Thulstrup and Öhrn (1972)  
^d Green (1973), Restricted HF  
^e Lefebvre-Brion and Moser (1966), CI  
^f Thulstrup and Öhrn (1972), CI

### Table 4: Calculated and experimental spectroscopic constants for nitric oxide (X²Π)

| Constant      | Expt⁴  | This Calc | Calc^b |
|---------------|--------|-----------|--------|
| rₑ (Å)        | 1.15   | 1.10      | 1.25   |
| ωₑ (cm⁻¹)     | 1904   | 1993      | 1731   |
| ωₑxₑ (cm⁻¹)   | 14     | 20        | 24     |
| Bₑ (cm⁻¹)     | 1.70   | 1.86      | 1.45   |
| αₑ (cm⁻¹)     | 0.018  | 0.032     | 0.023  |
| Dₑ (eV)       | 6.5    | 6.0       | 3.8    |

⁴Herzberg (1950)  
^b Thulstrup and Öhrn (1972)

### B. Nickel(100) Surface

Much has been said in these energy-conscious times about the importance of the study of surfaces and of chemisorption on them in the field of catalysis. On the theoretical side, various approaches stem from both solid-state band theory and theories generally applied to molecular complexes; however, we are still far from possessing a reliable computational method for a satisfactory description of the chemisorption bond.

The approach we take here is, of course, a localized one with respect to the role of the surface in the chemisorption bond. In order to avoid the effects of cluster boundaries, we shall compose a finite atom-cluster model built up from atoms which retain their "bulk", "surface" and "adsorbate" identities.

### Table 5: Basis set for nickel

| Type | Exponent |
|------|----------|
| 1s   | 27.70    |
| 2s   | 20.50    |
| 2p   | 17.90    |
| 3s   | 9.86     |
| 3s'  | 5.52     |
| 3p   | 8.65     |
| 3p'  | 4.76     |
| 3d   | 5.35     |
| 3d'  | 2.14     |
| 4s   | 1.61     |
| 4s'  | 0.63     |
| 4p   | 1.00     |

### Table 6: Model potential exponents for nickel

| Type | Exponent |
|------|----------|
| s₁   | 1.14     |
| s₂   | 2.29     |
| s₃   | 4.57     |
| s₄   | 9.14     |

### C. Nitric Oxide on Nickel(100) Surface

The catalytic reduction of nitric oxide (by hydrogen) over noble metal surfaces is a reasonably well-known process. It is, in fact, the underlying process of the catalytic converters installed in newer emission-controlled automobiles.

Recent X-ray (XPS) and ultraviolet (UPS) studies of nitric oxide and nitrogen dioxide interactions with nickel, and infrared (IR) data for nitric oxide on nickel have posed a couple of interesting questions. There is reasonable evidence to suggest that nitrogen dioxide dissociates on nickel even at very low temperatures (80K) leaving nitric oxide and atomic oxygen adsorbed on the surface.

### Table 7: Calculated and experimental core(1s) binding energies for nitric oxide on nickel(100) surface

| (eV)     | N1s    | O1s    |
|----------|--------|--------|
| **Ni-N-O** |        |        |
| Calc     | 415.5  | 548.0  |
| Expt⁴    | 400    | 531    |
| **Ni-O-N** |        |        |
| Calc     | 417.4  | 548.9  |
| Expt^b   | 403    | 532    |

⁴Brundle (1976)  
^b Attributed to different state of nitric oxide

### Table 8: Calculated spectroscopic constants for nitric oxide on nickel(100) surface

| Constant      | Calc  | Expt⁴  |
|---------------|-------|--------|
| rₑ (Å)        | 1.24  | --     |
| ωₑ (cm⁻¹)     | 1621  | 1590   |
| ωₑxₑ (cm⁻¹)   | 68    | --     |
| Dₑ (eV)       | 1.0   | --     |

⁴Blyholder and Allen (1965), IR spectrum

---

## CHAPTER X: DISCUSSION

The method we have tried to present here is still in a youthful stage of its development, having been conceived, given vitality in the ever-amendable form of computer software, and ready to withstand the test of time. Its maturity will most certainly come only with growth and inevitable change, for nearly every assumption or approximation that has been made in its early stages can be improved upon to some extent if necessary.

It would be worthwhile, perhaps, to review the assumptions and approximations which are to be placed under such scrutiny. The most fundamental approximation, outside the independent-particle model, that must be examined carefully for its justification is the use of the Mulliken approximation at various stages. It was the use of this approximate form of the differential overlap (accurate to second-order in diatomic overlap) which allows for the replacement of the constrained eigenfunctions by their renormalized unconstrained counterparts. And this is what leads to the single-centered expansions.

The next approximation which comes to mind is the Hermitean truncated (through second-order in overlap) expansion of the energy-dependent part of the pseudopotential. The integral approximations are perhaps the easiest to criticize, but these approximations are rather straightforward and can be improved if the need is indicated.

The model Hamiltonian used to obtain molecular orbitals is perhaps the most disputable aspect of our treatment. The problem is not so much with the molecular orbitals or energies themselves but rather with the subsequent population analysis which one would like to employ in the self-consistent charge and configuration procedure.

To reiterate, all of these points demand further testing before much can be said in their defense or otherwise. Let us end this report of a new and promising approach to the study of polyatomic systems, as though it really were in its youth and eagerly awaiting its chance to prove itself in a world which is equally anxious for the arrival of such a method.

---

## BIBLIOGRAPHY

Abdulnur, S.F., Linderberg, J., Öhrn, Y., and Thulstrup, P.W. (1972). Phys. Rev. A 6, 889.

Adams, W.H. (1965). J. Chem. Phys. 42, 4030.

Adams, W.H. (1967). Phys. Rev. 156, 109.

Adams, W.H. (1971). Chem. Phys. Letters 11, 441.

Anderson, A.B., and Hoffmann, R. (1974). J. Chem. Phys. 61, 4545.

Bardo, R.D., and Ruedenberg, K. (1973). J. Chem. Phys. 59, 5956.

Batra, I.P., and Brundle, C.R. (1976). Surf. Sci. 57, 12.

Blyholder, G., and Allen, M.C. (1965). J. Phys. Chem. 69, 3998.

Brundle, C.R. (1976). J. Vac. Sci. Technol. 13, 301.

Bullett, D.W., and Cohen, M.L. (1977a). J. Phys. C 10, 2083.

Bullett, D.W., and Cohen, M.L. (1977b). J. Phys. C 10, 2101.

Collin, J.E., and Natalis, P. (1968). Chem. Phys. Letters 2, 194.

Connolly, J.W.D., Siegbahn, H., Gelius, U., and Nordling, C. (1973). J. Chem. Phys. 58, 4265.

Coulson, C.A., and Fischer, I. (1949). Phil. Mag. 40, 386.

Csavinsky, P. (1968). Phys. Rev. 166, 53.

Csavinsky, P. (1973). Phys. Rev. A 8, 1688.

Das, G., and Wahl, A.C. (1976). J. Chem. Phys. 64, 4672.

Dunning, T.H., and Hay, P. Jeffrey, Jr. (1977). in "Modern Theoretical Chemistry" Vol. 3 (H.F. Schaefer, ed.), Chap. 1. Plenum Press, New York.

Fermi, E. (1928). Z. Physik 48, 73.

Fock, V. (1932). Z. Physik 75, 622.

Gerratt, J. (1974). in "Theoretical Chemistry" Vol. 1 (R.N. Dixon, ed.). The Chemical Society, London.

Gerratt, J. (1971). Adv. At. Mol. Phys. 7, 141.

Gilbert, T.L. (1972). Phys. Rev. A6, 580.

Green, S. (1973). Chem. Phys. Letters 23, 115.

Green, S. (1972). Chem. Phys. Letters 13, 552.

Gunnarsson, O., Harris, J., and Jones, R.O. (1977). J. Chem. Phys. 67, 3970.

Heitler, W. (1934). Marx Handb. d. Radiologie II, 485.

Hellmann, H. (1935). J. Chem. Phys. 3, 61.

Herzberg, G. (1950). "Spectra of Diatomic Molecules" Van Nostrand, New York.

Hodges, L., Watson, R.E., and Ehrenreich, H. (1972). Phys. Rev. B5, 3953.

Hollister, C., and Sinanoglu, O. (1966). J. Am. Chem. Soc. 88, 13.

Hulbert, H.M., and Hirschfelder, J.O. (1941). J. Chem. 9, 61.

Hurley, A.C., Lennard-Jones, J., and Pople, J.A. (1953). Proc. Roy. Soc. London A220, 446.

Huzinaga, S., McWilliams, D., and Cantu, A.A. (1973). Adv. Quant. Chem. 7, 187.

Jørgensen, P., and Öhrn, Y. (1973). Phys. Rev. A8, 112.

Kahn, L.R., Baybutt, P., and Truhlar, D.G. (1976). J. Chem. Phys. 65, 3826.

Kobylinski, T.P., and Taylor, B.W. (1974). J. Catalysis 33, 376.

Koopmans, T.A. (1935). Physica 1, 104.

Kouba, J., and Öhrn, Y. (1971). Int. J. Quant. Chem. 5, 539.

Landman, U., and Adams, D.L. (1974). J. Vac. Sci. Technol. 11, 191.

Lefebvre-Brion, H., and Moser, C.M. (1966). J. Chem. Phys. 44, 2951.

Linderberg, J., Öhrn, Y., and Thulstrup, P.W. (1976). in "Quantum Science" Plenum Press, New York.

Linderberg, J., and Öhrn, Y. (1973). "Propagators in Quantum Chemistry" Academic Press, New York.

Löwdin, P.O. (1961). J. Math. Phys. 2, 969.

Löwdin, P.O. (1964). J. Mol. Spec. 14, 112.

Löwdin, P.O. (1970). Adv. Quant. Chem. 5, 185.

McIntosh, H.V. (1975). "Plot 75" I.N.E.N., Mexico.

Moffitt, W. (1951). Proc. Roy. Soc. London A210, 247.

Morrison, R.T., and Boyd, R.N. (1966). "Organic Chemistry" (2nd ed.) Allyn and Bacon, Boston.

Mulliken, R.S. (1949). J. Chim. Phys. 46, 497.

Mulliken, R.S. (1955). J. Chem. Phys. 23, 1833, 1841, 2338, 2343.

O'Malley, T.F. (1971). Adv. At. Mol. Phys. 7, 223.

Pauncz, R. (1967). "Alternate Molecular Orbital Method" Saunders, Philadelphia.

Phillips, J., and Kleinman, L. (1959). Phys. Rev. 116, 287.

Roetti, C., and Clementi, E. (1974). Atomic Data and Nuclear Data Tables 14, 177.

Siegbahn, K., Nordling, C., Johansson, G., Hedman, J., Heden, P.F., Hamrin, K., Gelius, U., Bergmark, T., Werme, L.O., Manne, R., and Baer, Y. (1969). "ESCA Applied to Free Molecules" North-Holland, Amsterdam.

Slater, J.C. (1934). Phys. Rev. 45, 794.

Slater, J.C., Mann, J.B., Wilson, T.M., and Wood, J.H. (1969). Phys. Rev. 184, 672.

Slater, J.C. (1970). Int. J. Quant. Chem. 3S, 727.

Slater, J.C., and Wood, J.H. (1971). Int. J. Quant. Chem. 4S, 3.

Slater, J.C. (1971). Int. J. Quant. Chem. 5, 403.

Thomas, L.H. (1928). Proc. Camb. Philos. Soc. 23, 542.

Thulstrup, P.W., and Öhrn, Y. (1972). J. Chem. Phys. 57, 3716.

Turner, D.W., Baker, C., Baker, A.D., and Brundle, C.R. (1970). "Molecular Photoelectron Spectroscopy" Wiley-Interscience, London.

Weeks, J.D., Hazi, A., and Rice, S.A. (1969). Adv. Chem. Phys. 16, 283.

Wolfsberg, M., and Helmholtz, L. (1952). J. Chem. Phys. 20, 837.

---

## BIOGRAPHICAL SKETCH

Jack A. Smith was born October 16, 1949, in Indianapolis, Indiana. He graduated from Arlington High School at Indianapolis in June, 1967. In June, 1971, he received the degree of Bachelor of Science with majors in chemistry and mathematics from Indiana University - Purdue University at Indianapolis. In September, 1971, he enrolled in the Graduate School at the University of Florida. In August, 1972, he joined the Quantum Theory Project at the University of Florida. In December, 1973, he received the degree of Master of Science in chemistry. Since that time he has pursued the degree of Doctor of Philosophy at the University of Florida and is presently still a member of the Quantum Theory Project. On July 3, 1976, he married the former Miss (Susan) Malissa Brockett.