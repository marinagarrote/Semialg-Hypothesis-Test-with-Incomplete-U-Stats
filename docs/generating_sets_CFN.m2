-- =============================================================================
--                   Generating sets for the CFN ideal 
-- =============================================================================
-- This script computes different generating sets for the ideal associated with the 
-- Cavender-Farris-Neyman (CFN) two-state model of evolution on a 4-leaf unrooted 
-- tree with topology ((1,2),(3,4)).
--
-- The script describes the way of obtaining the different generating sets:
--   1. Partially/Completely Distinguishing Minimal (PDM & CDM): 
--      As a minimal basis of the kernel of a parametrization map.
--   2. Completely Distinguishing Determinantal (CDD):
--      Through determinantal equations from a Fourier transform.
--   3. Partially/Completely Distinguishing Rank (PDR & CDR): 
--      from 3x3 minors of a certain flattening of the distribution p.
--
-- It also shows the equivalence between these different generating sets systems.
-- =============================================================================



-- -----------------------------------------------------------------------------
-- I. SETUP: DEFINE POLYNOMIAL RINGS
-- -----------------------------------------------------------------------------

-- R: The ring of parameters. There are 5 parameters (θi), one for each edge
--    of the unrooted 4-leaf tree T₁₂₋₃₄. These parameters relate to edge lengths.
R = QQ[θ1,θ2,θ3,θ4,θ5];

-- Sp: The ring of symmetrized pattern probabilities. For 4 binary leaves, there are
--     2^4 = 16 possible patterns (0000, 0001, ..., 1111). Due to the symmetry
--     of the model, p(ijkl) = p(1-i,1-j,1-k,1-l), so we only need 2^(4-1) = 8
--     coordinates. `pxxxx` represents p(0000)+p(1111), `pxxxy` represents
--     p(0001)+p(1110), and so on.
Sp = QQ[pxxxx, pxxxy, pxxyx, pxxyy, pxyxx, pxyxy, pxyyx, pxyyy];



-- -----------------------------------------------------------------------------
-- II. PARAMETRIZATION: Define the Model of Evolution
-- -----------------------------------------------------------------------------

-- Define the transition probabilities for the CFN model on an edge.
β = θ -> (1-θ)/2; -- Probability of state change (e.g., 0 -> 1)
α = θ -> (1+θ)/2; -- Probability of state conservation (e.g., 0 -> 0)

-- Define the parametrization for the symmetrized pattern probabilities (p̄)
-- for the CFN model on the tree T = T₁₂₋₃₄.
Pxxxx = α(θ1)*α(θ2)*α(θ3)*α(θ4)*α(θ5) + α(θ1)*α(θ2)*β(θ3)*β(θ4)*β(θ5) + α(θ3)*α(θ4)*β(θ1)*β(θ2)*β(θ5) + α(θ5)*β(θ1)*β(θ2)*β(θ3)*β(θ4); -- p0000 + p1111
Pxxxy = α(θ1)*α(θ2)*α(θ3)*α(θ5)*β(θ4) + α(θ1)*α(θ2)*α(θ4)*β(θ3)*β(θ5) + α(θ3)*β(θ1)*β(θ2)*β(θ4)*β(θ5) + α(θ4)*α(θ5)*β(θ1)*β(θ2)*β(θ3); -- p0001 + p1110
Pxxyx = α(θ1)*α(θ2)*α(θ3)*β(θ4)*β(θ5) + α(θ1)*α(θ2)*α(θ4)*α(θ5)*β(θ3) + α(θ3)*α(θ5)*β(θ1)*β(θ2)*β(θ4) + α(θ4)*β(θ1)*β(θ2)*β(θ3)*β(θ5); -- p0010 + p1101
Pxxyy = α(θ1)*α(θ2)*α(θ3)*α(θ4)*β(θ5) + α(θ1)*α(θ2)*α(θ5)*β(θ3)*β(θ4) + α(θ3)*α(θ4)*α(θ5)*β(θ1)*β(θ2) + β(θ1)*β(θ2)*β(θ3)*β(θ4)*β(θ5); -- p0011 + p1100
Pxyxx = α(θ1)*α(θ3)*α(θ4)*α(θ5)*β(θ2) + α(θ1)*β(θ2)*β(θ3)*β(θ4)*β(θ5) + α(θ2)*α(θ3)*α(θ4)*β(θ1)*β(θ5) + α(θ2)*α(θ5)*β(θ1)*β(θ3)*β(θ4); -- p0100 + p1011
Pxyxy = α(θ1)*α(θ3)*α(θ5)*β(θ2)*β(θ4) + α(θ1)*α(θ4)*β(θ2)*β(θ3)*β(θ5) + α(θ2)*α(θ3)*β(θ1)*β(θ4)*β(θ5) + α(θ2)*α(θ4)*α(θ5)*β(θ1)*β(θ3); -- p0101 + p1010
Pxyyx = α(θ1)*α(θ3)*β(θ2)*β(θ4)*β(θ5) + α(θ1)*α(θ4)*α(θ5)*β(θ2)*β(θ3) + α(θ2)*α(θ3)*α(θ5)*β(θ1)*β(θ4) + α(θ2)*α(θ4)*β(θ1)*β(θ3)*β(θ5); -- p0110 + p1001
Pxyyy = α(θ1)*α(θ3)*α(θ4)*β(θ2)*β(θ5) + α(θ1)*α(θ5)*β(θ2)*β(θ3)*β(θ4) + α(θ2)*α(θ3)*α(θ4)*α(θ5)*β(θ1) + α(θ2)*β(θ1)*β(θ3)*β(θ4)*β(θ5); -- p0111 + p1000

-- A list containing all the parametrization polynomials.
P = {Pxxxx, Pxxxy, Pxxyx, Pxxyy, Pxyxx, Pxyxy, Pxyyx, Pxyyy};



-- -----------------------------------------------------------------------------
-- III. MINIMAL GENERATORS DESCRIPTION (PDM & CDM)
-- -----------------------------------------------------------------------------

-- This section finds the implicit algebraic relations (invariants) among the
-- probability distributions `p...` by eliminating the parameters `θi`.

f_map = map(R, Sp, P);
I = kernel f_map;
M = mingens I;
-- netList flatten entries M -- Uncomment to view the minimal generators

-- The first generator is the linear relation stating that probabilities sum to 1.
f_sum = M_0_0 -- pxxxx + pxxxy + ... + pxyyy - 1 == 0

-- The other two generators, `h0` and `h1`, are the non-linear invariants.
h0 = M_1_0
h1 = M_2_0

-- PDM: The ideal of non-linear invariants in the standard probability coordinates.
PDM = ideal{h0, h1}

-- CDM: An equivalent ideal obtained by linear combinations.
CDM = ideal{h0 + h1, h0 - h1}


-- -----------------------------------------------------------------------------
-- IV. DETERMINENTAL DESCRIPTION (CDD)
-- -----------------------------------------------------------------------------

-- This section defines the invariants using determinants after a Fourier transform.

-- Define the Fourier (Hadamard) Coordinates `q`.
qxxxx = pxxxx + pxxxy + pxxyx + pxxyy + pxyxx + pxyxy + pxyyx + pxyyy;
qxxyy = pxxxx - pxxxy - pxxyx + pxxyy + pxyxx - pxyxy - pxyyx + pxyyy;
qxyxy = pxxxx - pxxxy + pxxyx - pxxyy - pxyxx + pxyxy - pxyyx + pxyyy;
qxyyx = pxxxx + pxxxy - pxxyx - pxxyy - pxyxx - pxyxy + pxyyx + pxyyy;
qyxxy = pxxxx - pxxxy + pxxyx - pxxyy + pxyxx - pxyxy + pxyyx - pxyyy;
qyxyx = pxxxx + pxxxy - pxxyx - pxxyy + pxyxx + pxyxy - pxyyx - pxyyy;
qyyxx = pxxxx + pxxxy + pxxyx + pxxyy - pxyxx - pxyxy - pxyyx - pxyyy;
qyyyy = pxxxx - pxxxy - pxxyx + pxxyy - pxyxx + pxyxy + pxyyx - pxyyy;

-- Construct 2x2 matrices that are rank 1 on the model variety.
M1 = matrix{{qxxxx, qxxyy},
            {qyyxx, qyyyy}}

M2 = matrix{{qxyxy, qyxxy},
            {qxyyx, qyxyx}}

-- The determinants of M1 and M2 vanish on the variety and are the invariants.
F1 = det(M1)
F2 = det(M2)

-- CDD: The ideal generated by the determinantal invariants.
CDD = ideal{F1, F2}

-- --- Optional: Verification of the Determinantal Equations ---
-- Sq = QQ[qxxxx, qxxyy, qxyxy, qxyyx, qyxxy, qyxyx, qyyxx, qyyyy];
-- Qxxxx = Pxxxx + Pxxxy + Pxxyx + Pxxyy + Pxyxx + Pxyxy + Pxyyx + Pxyyy
-- Qxxyy = Pxxxx - pxxxy - pxxyx + Pxxyy + Pxyxx - pxyxy - pxyyx + Pxyyy
-- Qxyxy = Pxxxx - pxxxy + pxxyx - Pxxyy - Pxyxx + pxyxy - pxyyx + Pxyyy
-- Qxyyx = Pxxxx + pxxxy - pxxyx - Pxxyy - Pxyxx - pxyxy + pxyyx + Pxyyy
-- Qyxxy = Pxxxx - pxxxy + pxxyx - Pxxyy + Pxyxx - pxyxy + pxyyx - Pxyyy
-- Qyxyx = Pxxxx + pxxxy - pxxyx - Pxxyy + Pxyxx + pxyxy - pxyyx - Pxyyy
-- Qyyxx = Pxxxx + pxxxy + pxxyx + pxxyy - pxyxx - pxyxy - pxyyx - pxyyy
-- Qyyyy = Pxxxx - pxxxy - pxxyx + pxxyy - pxyxx + pxyxy + pxyyx - pxyyy
-- Q = {Qxxxx, Qxxyy, Qxyxy, Qxyyx, Qyxxy, Qyxyx, Qyyxx, Qyyyy};
-- g_map = map(R, Sq, Q);
-- netList entries gens kernel g_map



-- -----------------------------------------------------------------------------
-- V. RANK DESCRIPTION (PDR & CDR)
-- -----------------------------------------------------------------------------

-- This section defines the invariants using rank constraints on a flattening matrix.

-- The matrix `Flat1234` is a "flattening" that splits leaves (1,2) against (3,4).
-- The model requires that rank(Flat1234) <= 2.
Flat1234 = matrix{{pxxxx, pxxxy, pxxyx, pxxyy},
                  {pxyxx, pxyxy, pxyyx, pxyyy},
                  {pxyyy, pxyyx, pxyxy, pxyxx},
                  {pxxyy, pxxyx, pxxxy, pxxxx}};

-- The ideal of 3x3 minors of the flattening matrix.
Iminors = minors(3, Flat1234);
-- netList primaryDecomposition I_minors -- The ideal of minors is not prime.

-- Two specific linear combinations of the minors generate the relevant prime component.
g1 = (primaryDecomposition Iminors)_0_1 -- g1 = pxxyx*pxyxx - pxxyy*pxyxy - pxxxx*pxyyx + pxxxy*pxyyy;
g2 = (primaryDecomposition Iminors)_0_1 -- g2 = pxxxy*pxyxx - pxxxx*pxyxy - pxxyy*pxyyx + pxxyx*pxyyy; 

-- PDR: The ideal generated by the two rank-based invariants.
PDR = ideal{g1, g2}

-- CDR: An equivalent ideal in a different basis.
CDR = ideal{g1 + g2, g1 - g2}


-- -----------------------------------------------------------------------------
-- VI. VERIFICATION OF EQUIVALENCE
-- -----------------------------------------------------------------------------
-- The following lines assert that the ideals of non-linear invariants become equal
-- when restricted to the affine plane where probabilities sum to 1.

linear_constraint = ideal f_sum

-- Check if the different generator bases are equivalent on the model space.
(linear_constraint + PDM) == (linear_constraint + CDM) -- Expected: true
(linear_constraint + PDR) == (linear_constraint + CDR) -- Expected: true

-- Check if the main descriptions (Minimal, Determinantal, Rank) are equivalent.
(linear_constraint + PDM) == (linear_constraint + CDD) -- Expected: true
(linear_constraint + PDM) == (linear_constraint + PDR) -- Expected: true
(linear_constraint + PDR) == (linear_constraint + CDD) -- Expected: true
