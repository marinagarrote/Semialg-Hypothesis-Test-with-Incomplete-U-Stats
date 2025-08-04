-- This script computes the singular locus of the variety associated to the
-- CFN model on a binary 4 leaf tree with three topology 12|34.
-- =============================================================================

-- -----------------------------------------------------------------------------
-- I. SETUP: DEFINE POLYNOMIAL RINGS
-- -----------------------------------------------------------------------------

-- Sp: The ring of symmetrized pattern probabilities.
Sp = QQ[pxxxx, pxxxy, pxxyx, pxxyy, pxyxx, pxyxy, pxyyx, pxyyy];

-- R: The ring of parameters (θ), one for each edge of the unrooted 4-leaf tree.
R = QQ[θ1,θ2,θ3,θ4,θ5];

-- -----------------------------------------------------------------------------
-- II. PARAMETRIZATION: Define the parametrization of the model
-- -----------------------------------------------------------------------------

-- Transition probabilities for the CFN model on an edge.
α = θ -> (1+θ)/2; -- Probability of state conservation (e.g., 0 -> 0)
β = θ -> (1-θ)/2; -- Probability of state change (e.g., 0 -> 1)

-- A list of polynomials defining the parametrization from R to Sp.
-- for the CFN model on the tree T = T₁₂₋₃₄.
Pxxxx = α(θ1)*α(θ2)*α(θ3)*α(θ4)*α(θ5) + α(θ1)*α(θ2)*β(θ3)*β(θ4)*β(θ5) + β(θ1)*β(θ2)*α(θ3)*α(θ4)*β(θ5) + β(θ1)*β(θ2)*β(θ3)*β(θ4)*α(θ5); -- p0000 + p1111
Pxxxy = α(θ1)*α(θ2)*α(θ3)*β(θ4)*α(θ5) + α(θ1)*α(θ2)*β(θ3)*α(θ4)*β(θ5) + β(θ1)*β(θ2)*α(θ3)*β(θ4)*β(θ5) + β(θ1)*β(θ2)*β(θ3)*α(θ4)*α(θ5); -- p0001 + p1110
Pxxyx = α(θ1)*α(θ2)*α(θ3)*β(θ4)*β(θ5) + α(θ1)*α(θ2)*β(θ3)*α(θ4)*α(θ5) + β(θ1)*β(θ2)*α(θ3)*β(θ4)*α(θ5) + β(θ1)*β(θ2)*β(θ3)*α(θ4)*β(θ5); -- p0010 + p1101
Pxxyy = α(θ1)*α(θ2)*α(θ3)*α(θ4)*β(θ5) + α(θ1)*α(θ2)*β(θ3)*β(θ4)*α(θ5) + β(θ1)*β(θ2)*α(θ3)*α(θ4)*α(θ5) + β(θ1)*β(θ2)*β(θ3)*β(θ4)*β(θ5); -- p0011 + p1100
Pxyxx = α(θ1)*β(θ2)*α(θ3)*α(θ4)*α(θ5) + α(θ1)*β(θ2)*β(θ3)*β(θ4)*β(θ5) + β(θ1)*α(θ2)*α(θ3)*α(θ4)*β(θ5) + β(θ1)*α(θ2)*β(θ3)*β(θ4)*α(θ5); -- p0100 + p1011
Pxyxy = α(θ1)*β(θ2)*α(θ3)*β(θ4)*α(θ5) + α(θ1)*β(θ2)*β(θ3)*α(θ4)*β(θ5) + β(θ1)*α(θ2)*α(θ3)*β(θ4)*β(θ5) + β(θ1)*α(θ2)*β(θ3)*α(θ4)*α(θ5); -- p0101 + p1010
Pxyyx = α(θ1)*β(θ2)*α(θ3)*β(θ4)*β(θ5) + α(θ1)*β(θ2)*β(θ3)*α(θ4)*α(θ5) + β(θ1)*α(θ2)*α(θ3)*β(θ4)*α(θ5) + β(θ1)*α(θ2)*β(θ3)*α(θ4)*β(θ5); -- p0110 + p1001
Pxyyy = α(θ1)*β(θ2)*α(θ3)*α(θ4)*β(θ5) + α(θ1)*β(θ2)*β(θ3)*β(θ4)*α(θ5) + β(θ1)*α(θ2)*α(θ3)*α(θ4)*α(θ5) + β(θ1)*α(θ2)*β(θ3)*β(θ4)*β(θ5); -- p0111 + p1000

-- A list containing all the parametrization polynomials.
P = {Pxxxx, Pxxxy, Pxxyx, Pxxyy, Pxyxx, Pxyxy, Pxyyx, Pxyyy};


-- -----------------------------------------------------------------------------
-- III. IMPLICIT IDEAL & SINGULAR LOCUS IN PROBABILITY SPACE
-- -----------------------------------------------------------------------------

-- The parametrization map from the parameter ring R to the probability ring Sp.
f_map = map(R, Sp, P);

-- The implicit ideal I of the model, which defines the variety in Sp.
I = kernel f_map;

-- The Jacobian matrix of the ideal I. 
J = jacobian(I);

-- The singular locus is where the rank of J drops below the codimension of the variety.
-- rank(Jacobian_I) < codimension V(I)
MinorsIdeal = minors(codim I, J);

-- The ideal (M + I) gives the singular locus.
primaryDecomposition (M + I)

--------------------------------------------------------------------------------
-- IV. SINGULAR LOCUS IN PARAMETER SPACE
-- -----------------------------------------------------------------------------

-- The combined ring for implicitization and parameter-space analysis.
RS = QQ[flatten entries vars(R) | flatten entries vars(Sp)];

-- Define the parametrization ideal in the combined ring RS.
Ip = ideal({
    sub(Pxxxx, RS) - pxxxx, 
    sub(Pxxxy, RS) - pxxxy, 
    sub(Pxxyx, RS) - pxxyx, 
    sub(Pxxyy, RS) - pxxyy, 
    sub(Pxyxx, RS) - pxyxx, 
    sub(Pxyxy, RS) - pxyxy, 
    sub(Pxyyx, RS) - pxyyx, 
    sub(Pxyyy, RS) - pxyyy
});

-- The elimination ideal, E, represents the implicit equations of the singular locus in the parameter space.
E = eliminate({pxxxx, pxxxy, pxxyx, pxxyy, pxyxx, pxyxy, pxyyx, pxyyy}, sub(MinorsIdeal, RS) + Ip);

-- The primary decomposition of E gives the components of the singular locus in the parameter space.
primaryDecomposition E
